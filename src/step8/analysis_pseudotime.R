# ###############################################################
# This script aims to use Monocle method to order cells
# with pseudotime
# ###############################################################

# This script is inpired by the Monocle tutorial:
# http://cole-trapnell-lab.github.io/monocle-release/docs/

## 
# Define function to get a good color palette for cell counts
##
build_color_palette <- function( colors, local_emb, gradient.range.quantile=0.95){
  
  if(all(sign(colors)>=0)) {
    gradientPalette <- colorRampPalette(c('gray80','red'), space = "Lab")(1024)
  } else {
    gradientPalette <- colorRampPalette(c("blue", "grey70", "red"), space = "Lab")(1024)
  }
  if(all(sign(colors)>=0)) {
    zlim <- as.numeric(quantile(colors,p=c(1-gradient.range.quantile,gradient.range.quantile)))
    if(diff(zlim)==0) {
      zlim <- as.numeric(range(colors))
    }
  } else {
    zlim <- c(-1,1)*as.numeric(quantile(abs(colors),p=gradient.range.quantile))
    if(diff(zlim)==0) {
      zlim <- c(-1,1)*as.numeric(max(abs(colors)))
    }
  }
  # restrict the values
  colors[colors<zlim[1]] <- zlim[1]; colors[colors>zlim[2]] <- zlim[2];
  
  colors <- (colors-zlim[1])/(zlim[2]-zlim[1])
  cols <- gradientPalette[colors[match(rownames(local_emb),names(colors))]*(length(gradientPalette)-1)+1]
  
  return( cols)
}


## @knitr pseudotime_analysis
set.seed( 123)

cat("<H4>Computing pseudo-time representation</H4>")
# Order cells by pseudotime using markers genes from previous analysis
filtered_cds_pseudotime_supervised <- setOrderingFilter( filtered_cds_for_pseudotime, signature_genes_for_pseudotime)
plot_ordering_genes( filtered_cds_pseudotime_supervised) 
filtered_cds_pseudotime_supervised <- reduceDimension( filtered_cds_pseudotime_supervised, max_components = 2,
                            method = 'DDRTree')
filtered_cds_pseudotime_supervised <- orderCells( filtered_cds_pseudotime_supervised)


#
# Plot Ccr7 on t-SNE map because it was not possible on RNA-velocity analysis (Ccr7 is not expresses
# by enough cells to be kept)
#
cat("<H5>Plotting Ccr7 gene expression on t-SNE</H5>")

emb_Ccr7 = data.frame( emb)
row.names( emb_Ccr7) = gsub( "x", "", gsub( "10635173:", "", rownames( emb_Ccr7)))
expressions =  t( exprs( filtered_cds_pseudotime_supervised))
ccr_cells = intersect( rownames( expressions), row.names( emb_Ccr7))
emb_Ccr7$Ccr7 = rep( NA, nrow( emb_Ccr7))
emb_Ccr7[ ccr_cells, "Ccr7"] = expressions[ ccr_cells, "Ccr7"]
emb_Ccr7$Ccr7[ which( emb_Ccr7$Ccr7 == 0)] = NA
plot_tsne_Ccr7 = ggplot( data = NULL) +
  geom_point( data =  emb_Ccr7[  which( is.na( emb_Ccr7$Ccr7)), ], aes( x = X1, y = X2, color = Ccr7)) +
  geom_point( data =  emb_Ccr7[  which( !is.na( emb_Ccr7$Ccr7)), ], aes( x = X1, y = X2, color = Ccr7)) +
  theme(legend.position = "none") +
  scale_colour_gradient2( low = "blue", mid = "lightblue", high = "red",
                          midpoint =  mean( emb_Ccr7$Ccr7,  na.rm = TRUE), space = "Lab",
                          na.value = "grey", guide = "colourbar")
print( plot_tsne_Ccr7)


#
# Plot pseudotime tree by Pseudotime
#
cat("<H5>Plotting pseudotime tree by Pseudotime</H5>")
plot_cell_trajectory( filtered_cds_pseudotime_supervised, color = "Pseudotime")

#
# Plot pseudotime tree by CellType
#
cat("<H5>Plotting pseudotime tree by Cell Type</H5>")
plot_cell_trajectory( filtered_cds_pseudotime_supervised, color = "CellType") + facet_wrap(~CellType) + theme(legend.position="none")

#
# Plot pseudotime tree by State
#
cat("<H5>Plotting pseudotime tree by State</H5>")
state_colors = c( "firebrick2", "firebrick", "dodgerblue1", "dodgerblue4", "chartreuse3")
plot_cell_trajectory( filtered_cds_pseudotime_supervised, color = "State") +
  scale_color_manual(breaks = c("1", "2", "3", "4", "5"), values=state_colors) +
  theme(legend.position = "none")

#
# Show table of distribution of cell type among pseudotime states
#
cat("<H5>Table of distribution of cell type among pseudotime states</H5>")
cell_state_df = dcast( data.frame( table( pData( filtered_cds_pseudotime_supervised)[ , c( "CellType", "State")]))
                       , CellType ~ State, value.var="Freq")
row.names( cell_state_df) = cell_state_df$CellType
cell_state_df = cell_state_df[ , -which( names( cell_state_df) == "CellType")]
cell_state_df = cell_state_df[ -which( row.names( cell_state_df) == "Unknown"),]
datatable(  cell_state_df , caption = "Distribution of cell type among pseudotime states")

# Keep trace of original pseudotime
pData( filtered_cds_pseudotime_supervised)[ , "original.pseudotime"] = pData( filtered_cds_pseudotime_supervised)[ , "Pseudotime"]


#
# Compute some usefull information on pseudotime per state in order prepare change of pseudotime
# origin
#
cell_on_branch_list = list()
pseudotime_max = list()
pseudotime_min = list()
for( current_state in sort( unique( pData( filtered_cds_pseudotime_supervised)[ , "State"]))) {
  cell_on_branch_list[[ current_state]] = rownames( pData( filtered_cds_pseudotime_supervised))[ which( pData( filtered_cds_pseudotime_supervised)[ , "State"] == current_state)]
  pseudotime_max[[ current_state]] = max( pData( filtered_cds_pseudotime_supervised)[ cell_on_branch_list[[ current_state]], "original.pseudotime"])
  pseudotime_min[[ current_state]] = min( pData( filtered_cds_pseudotime_supervised)[ cell_on_branch_list[[ current_state]], "original.pseudotime"])
}

#
# Change pseudotime origin. Zero will be set to the end of state 3.
#  Pseudotime at the other branches will be set accordingly to respect the length of branches
#  provided by the original pseudotime
#
 
 # For state 3, put the max to 0 and the min to positive value
 current_state = 3
 cells_of_state = cell_on_branch_list[[ current_state]]
 pData( filtered_cds_pseudotime_supervised)[ cells_of_state, "Pseudotime"] = abs( pData( filtered_cds_pseudotime_supervised)[ cells_of_state, "original.pseudotime"] - pseudotime_max[[ current_state]])
 pseudotime_at_node = max( pData( filtered_cds_pseudotime_supervised)[ cells_of_state, "Pseudotime"])
 
 # For state 4, add the pseudotime at node to the difference to the min of the state 4
 current_state = 4
 cells_of_state = cell_on_branch_list[[ current_state]]
 pData( filtered_cds_pseudotime_supervised)[ cells_of_state, "Pseudotime"] = pseudotime_at_node + (pData( filtered_cds_pseudotime_supervised)[ cells_of_state, "original.pseudotime"] - pseudotime_min[[ current_state]])
 
 # For state 2, add the pseudotime at node to the difference to the max of the state 2 (since original pseudotime
 #  at this state is at the opposite direction)
 current_state = 2
 cells_of_state = cell_on_branch_list[[ current_state]]
 pData( filtered_cds_pseudotime_supervised)[ cells_of_state, "Pseudotime"] = pseudotime_at_node + ( pseudotime_max[[ current_state]] - pData( filtered_cds_pseudotime_supervised)[ cells_of_state, "original.pseudotime"])
 pseudotime_at_node = max( pData( filtered_cds_pseudotime_supervised)[ cells_of_state, "Pseudotime"])
 
 # For state 1, add the pseudotime at node to the difference to the min of the state 4
 current_state = 1
 cells_of_state = cell_on_branch_list[[ current_state]]
 pData( filtered_cds_pseudotime_supervised)[ cells_of_state, "Pseudotime"] = pseudotime_at_node + (pseudotime_max[[ current_state]] - pData( filtered_cds_pseudotime_supervised)[ cells_of_state, "original.pseudotime"])
 
 # For state 5, add the pseudotime at node to the difference to the min of the state 4
 current_state = 5
 cells_of_state = cell_on_branch_list[[ current_state]]
 pData( filtered_cds_pseudotime_supervised)[ cells_of_state, "Pseudotime"] = pseudotime_at_node + (pData( filtered_cds_pseudotime_supervised)[ cells_of_state, "original.pseudotime"] - pseudotime_min[[ current_state]])
 

cat("<H5>Plotting pseudotime tree by Centered Pseudotime</H5>")
plot_cell_trajectory( filtered_cds_pseudotime_supervised, color = "Pseudotime") +
  scale_colour_gradient2( low = "blue", mid = "lightblue", high = "red", 
                          midpoint = mean( filtered_cds_pseudotime_supervised$Pseudotime,  na.rm = TRUE), space = "Lab",
                          na.value = "grey", guide = "colourbar")

# Plot the cells in the t-SNE embedding (coming from previous analysis) with colors
# corresponding to the pseudotime and pseudotime states
emb_extention = data.frame( emb)
names( emb_extention) = c( "X", "Y")
row.names( emb_extention) = unlist( gsub( "x", "", gsub( "10635173:", "", row.names( emb_extention))))
emb_extention$State = rep( NA, nrow(emb_extention))
emb_extention$Pseudotime = rep( NA, nrow(emb_extention))
common_cells = intersect( row.names( pData( filtered_cds_pseudotime_supervised)), row.names( emb_extention))
emb_extention[ common_cells, "State"] = pData( filtered_cds_pseudotime_supervised)[ common_cells, "State"]
emb_extention$State = as.factor( emb_extention$State)
emb_extention[ common_cells, "Pseudotime"] = pData( filtered_cds_pseudotime_supervised)[ common_cells, "Pseudotime"]

# -- plot the t-SNE embedding with centered pseudotime
cat("<H5>Plotting pseudotime (centered) onto t-SNE embedding</H5>")
plot_tsne_pseudotime = ggplot( data = NULL) + 
  geom_point( data = emb_extention, aes( x = X, y = Y, color = Pseudotime)) + 
  theme(legend.position = "none") +
  scale_colour_gradient2( low = "blue", mid = "lightblue", high = "red", 
                          midpoint =  mean( filtered_cds_pseudotime_supervised$Pseudotime,  na.rm = TRUE), space = "Lab",
                          na.value = "grey", guide = "colourbar")

print( plot_tsne_pseudotime)

# -- plot the t-SNE embedding with pseudotime states
cat("<H5>Plotting all states onto t-SNE embedding</H5>")

# -- -- plot the states all in one
plot_tsne_state = ggplot( data = NULL) + 
  geom_point( data = emb_extention, aes( x = X, y = Y, color = State)) +
  scale_color_manual(breaks = c("1", "2", "3", "4", "5", "6", "7"), values = state_colors) +
  theme(legend.position = "none")
print( plot_tsne_state)

cat("<H5>Plotting states onto t-SNE embedding collored by state identity </H5>")

# -- -- plot the states on t-SNE mapin facet_wrap mode with color linked to state identity
plot_tsne_state_wrap = ggplot( data = NULL) + 
  geom_point( data = emb_extention, aes( x = X, y = Y, color = State)) +
  scale_color_manual(breaks = c("1", "2", "3", "4", "5", "6", "7"), values = state_colors) +
  theme(legend.position = "none") + facet_wrap(~State)
print( plot_tsne_state_wrap)

# -- -- plot the states on t-SNE map one by one with color linked to state identity
for( current_state in sort( unique( emb_extention$State))){
  current_df = emb_extention
  current_df$State = as.character( current_df$State)
  current_df$State[ which( current_df$State != current_state)] = "other"
  current_df$State = factor( current_df$State, levels = c( current_state, "other"))
  plot_tsne_one_state = ggplot( data = NULL) + 
    geom_point( data = current_df[ which( current_df$State == "other"),], aes( x = X, y = Y, color = State)) +
    geom_point( data = current_df[ which( current_df$State != "other"),], aes( x = X, y = Y, color = State)) +
    scale_color_manual(breaks = c( current_state, "other"), values = c( state_colors[ as.numeric( current_state)], "grey" )) +
    theme(legend.position = "none") + 
    ggtitle( paste( "Pseudotime state", current_state))
  print( plot_tsne_one_state)
}


cat("<H5>Plotting states onto t-SNE embedding collored by pseudotime </H5>")

# -- -- plot the states on t-SNE map in facet wrap mode with colors linked to pseudotime
plot_tsne_state_wrap = ggplot( data = NULL) + 
  geom_point( data = emb_extention, aes( x = X, y = Y, color = Pseudotime)) +
  scale_colour_gradient2( low = "blue", mid = "lightblue", high = "red", 
                          midpoint =  mean( filtered_cds_pseudotime_supervised$Pseudotime,  na.rm = TRUE), space = "Lab",
                          na.value = "grey", guide = "colourbar") +
  theme(legend.position = "none") + facet_wrap(~State)
print( plot_tsne_state_wrap)

# -- -- plot the states on t-SNE map one by one with colors linked to pseudotime
for( current_state in sort( unique( emb_extention$State))){
  current_df = emb_extention
  
  current_df$State = as.character( current_df$State)
  current_df$State[ which( current_df$State != current_state)] = "other"
  current_df$State = factor( current_df$State, levels = c( current_state, "other"))
  
  current_df$pseudotime = pData( filtered_cds_pseudotime_supervised)[ row.names( emb_extention), "Pseudotime"]
  
  plot_tsne_one_state = ggplot( data = NULL) + 
    geom_point( data = current_df[ which( current_df$State == "other"),], aes( x = X, y = Y, fill = "lightgrey"), alpha = 0.02) +
    geom_point( data = current_df[ which( current_df$State != "other"),], aes( x = X, y = Y, color = pseudotime)) +
    scale_colour_gradient2( low = "blue", mid = "lightblue", high = "red", 
                            midpoint = mean( filtered_cds_pseudotime_supervised$Pseudotime,  na.rm = TRUE), space = "Lab",
                            na.value = "grey", guide = "colourbar") +
    
    theme(legend.position = "none") + 
    ggtitle( paste( "Pseudotime state", current_state))
  print( plot_tsne_one_state)
}

cat("<H5>Plotting basins of emission and attractions onto t-SNE embedding</H5>")

# -- plot the t-SNE embedding with both clusters and basins of emission and attractions
basins = data.frame( xmin = c( -16*20/43,   -16*20/43,  -68*20/43,    -86*20/43,    38*20/43,  -110*20/43),
                     xmax = c( 17*20/43,    0*20/43,    -52*20/43,    -68*20/43,    61*20/43, -92*20/43),
                     ymin = c( -27*20/27,   6*20/27,    0*20/27,      -22*20/27,    6*20/27,  -29*20/27),
                     ymax = c( -10*20/27,   22*20/27,   12*20/27,     -13*20/27,    22*20/27, 0*20/27),
                     type = c( "emission",  "emission", "attraction", "attraction", "attraction", "attraction"))
basins$name = c( "BE1", "BE2", "BA1", "BA2", "BA3", "BA4")
row.names( basins) = basins$name

ggplot( data = NULL) + 
  geom_point( data = emb_extention, aes( x = X, y = Y, color = State)) +
  geom_rect( data = basins, aes( xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill = type), color="black", alpha=0.5) +
  geom_text( data = basins, aes( x=xmax, y=ymin, label=name))


cat("<H5>Analysis of pseudotime state dispersion over basins</H5>")

# -- look at the dispersion of cell pseudotime states in the basins of emission and attractions
bassin_dispersion_df = data.frame( stringsAsFactors = FALSE)
for( selected_bassin in c( "BE1", "BE2", "BA1", "BA2", "BA3", "BA4")){
  selected_bassin_cells = row.names( emb_extention)[ which( emb_extention$X >= basins[ selected_bassin, "xmin"] & emb_extention$X <= basins[ selected_bassin, "xmax"] &
       emb_extention$Y >= basins[ selected_bassin, "ymin"] & emb_extention$Y <= basins[ selected_bassin, "ymax"] )]

  current_dispersion_df = data.frame( t( data.frame( table( emb_extention[ selected_bassin_cells, "State"]), stringsAsFactors = FALSE)),stringsAsFactors = FALSE)
  names( current_dispersion_df) = c( "State1", "State2", "State3", "State4", "State5")
  row.names( current_dispersion_df) = c( "State" , selected_bassin)
  bassin_dispersion_df = rbind( bassin_dispersion_df, current_dispersion_df[2,])
}

bassin_dispersion_df <- sapply( bassin_dispersion_df, as.numeric )
row.names( bassin_dispersion_df) = c( "BE1", "BE2", "BA1", "BA2", "BA3", "BA4")
d3heatmap( bassin_dispersion_df, colors = "Blues", Colv=FALSE, Rowv = FALSE, scale = "none" )

#
# Compute Branched Expression Analysis Modeling (BEAM)
#

# -- Compute BEAM itself
cat("<H4>Computing Branched Expression Analysis Modeling</H4>")
BEAM_res <- BEAM( filtered_cds_pseudotime_supervised, branch_point = 1, cores = 1)
BEAM_res <- BEAM_res[ order(BEAM_res$qval),]
BEAM_res <- BEAM_res[, c("gene_short_name", "pval", "qval")]

# -- Plot BEAM for top 50 marker genes
cat("<H5>Plot branching analysis of top 50 marker genes</H5>")
ph_res1 = plot_multiple_branches_heatmap( filtered_cds_pseudotime_supervised[ row.names(BEAM_res)[1:50],],
                                branches = c( 1, 4, 5),
                                branches_name = c( "TP-State1", "DP-state4", "TN-State5"),
                                num_clusters = 3,
                                show_rownames = TRUE,
                                return_heatmap = TRUE)

pdf( file= file.path( OUTPUT_DIR, "BEAM_Top50Markers.pdf"))
print( ph_res1)
dev.off()

# -- Plot BEAM for selected marker genes
cat("<H5>Plot branching analysis of selected markers genes</H5>")
gene_set = sort( unique( c("Emb" ,"F11r", "Plet1", "Tspan10", "Treml4", "Ccr2", "Clec4d", "Trem3", "Il12b",
                           "Ccr7", "Ccl22", "Cenpa", "Ccnb2", "Mcm6", "Ccne1", "Ccne2", "Birc5", "S100a6",
                           "Pdia6", "Ranbp1", "Actb", "Gm2a", "Itm2b", "Mif", "Cd63", "Fth1",
                           "Apoc4", "Clec4a2", "Susd2", "Clec4e", "Cd9", "Nedd4", "Clec4a4",
                           "Tnnt1", "Fpr2", "Ccnd3", "Colq", "Frmd4b", "Itga9", "Wfdc17", "Igf1", "Ptgs1",
                           "Apoc4", "Tlr12", "Il1b", "Ceacam19", "Ifitm1", "Klra17", "Cd24a"
)))

ph_res2 = plot_multiple_branches_heatmap( filtered_cds_pseudotime_supervised[ gene_set,],
                                          branches = c( 1, 4, 5),
                                          branches_name = c( "TP-State1", "DP-state4", "TN-State5"),
                                          num_clusters = 3,
                                          show_rownames = TRUE,
                                          return_heatmap = TRUE)

pdf( file= file.path( OUTPUT_DIR, "BEAM_SelectedMarkers.pdf"))
print( ph_res2)
dev.off()

# -- Compute DEG and plot BEAM for best DEG
diff_test_res <- differentialGeneTest(filtered_cds_pseudotime_supervised[ signature_genes_for_pseudotime,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")

diff_test_res = diff_test_res[ order( diff_test_res$qval, decreasing = FALSE), ]

sig_gene_names <- row.names(subset(diff_test_res, qval < 0.00001))
selected_sig_gene_names = sig_gene_names[ 1:min( 50, length( sig_gene_names))]
ph_res3 = plot_multiple_branches_heatmap( filtered_cds_pseudotime_supervised[ selected_sig_gene_names,],
                                          branches = c( 1, 4, 5),
                                          branches_name = c( "TP-State1", "DP-state4", "TN-State5"),
                                          num_clusters = 3,
                                          show_rownames = TRUE,
                                          return_heatmap = TRUE)

pdf( file= file.path( OUTPUT_DIR, "BEAM_TOP50DEG.pdf"))
print( ph_res3)
dev.off()


