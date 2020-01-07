# ##################################################################
# This script aims to plot a high quality image of the
# RNA velocity result
# ##################################################################

## @knitr improve_rnavelocity_plot

## 
# Define function to get a good color palette for cell counts
##
build_color_palette <- function( colors, emb, gradient.range.quantile=0.95){
  
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
  cols <- gradientPalette[colors[match(rownames(emb),names(colors))]*(length(gradientPalette)-1)+1]
  
  return( cols)
}

##############
# MAIN
##############

NCELLS = 300
ARROW_SCALE = 20
GRIDN = 30
ARROW_LWD = 1

cat("<HR>")
par( mfrow = c(1,1))
show.velocity.on.embedding.cor(emb,
                               rvel.cd,
                               scale='sqrt',
                               cell.colors=ac(cell.colors,alpha=0.5),
                               cex=0.8,
                               show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5,
                               fixed.arrow.length = FALSE,
                               do.par=F,
                               cell.border.alpha = 0.1,
                               arrow.lwd=ARROW_LWD,
                               n= NCELLS,
                               arrow.scale=ARROW_SCALE,
                               grid.n=GRIDN
)


# Visualize a fit for maker genes (we reuse rvel.cd to save on calcualtions here):
cat("<HR>")
rgb.palette <- colorRampPalette(c("grey", "red"), space = "rgb")(1000)

gene_set = sort( c( "Emb" ,"F11r", "Plet1", "Tspan10", "Treml4", "Clec4d", "Trem3", "Il12b",
                    "Ccr7", "Ccl22", "Cenpa", "Ccnb2", "Mcm6", "Ccne1", "Birc5", "S100a6",
                    "Actb", "Gm2a", "Itm2b", "Mif", "Fth1", "Apoc4", "Clec4a2", "Susd2", 
                    "Cd9", "Clec4a4", "Tnnt1", "Ccnd3", "Wfdc17", "Igf1", 
                    "Apoc4", "Tlr12", "Il1b", "Ceacam19", "Ifitm1", "Klra17", "Cd24a"
                    

))

for( gene_name in gene_set){
  
  cat("<HR><H5>GENE", gene_name, "</H5>")
  if( gene_name %in% colnames( r_filtered$counts)){
    par( mfrow = c(1,1))
    # Plot the expression of genes without velocity grid
    r_filtered$plotEmbedding( type='PCA',embeddingType='tSNE',colors=r_filtered$counts[, gene_name], main= gene_name)
    
    # Plot the expression of genes with velocity grid
    cols = build_color_palette( r_filtered$counts[ , gene_name], emb)
    names( cols) = row.names( r_filtered$counts)
    show.velocity.on.embedding.cor(emb, rvel.cd, scale='sqrt', cell.colors= ac(cols, alpha=0.5),
                                   cex=0.8, show.grid.flow=TRUE,
                                   min.grid.cell.mass=0.5,
                                   fixed.arrow.length = FALSE, 
                                   do.par=F, cell.border.alpha = 0.1,
                                   arrow.lwd=ARROW_LWD,
                                   n=NCELLS, 
                                   arrow.scale=ARROW_SCALE, 
                                   grid.n=GRIDN)
    
  }
}

#
# Select basins of emission and attraction and look at markers genes of
# the different basins
#

# -- get the t-SNE embedding and add cluster colors to it
extended_embedding = data.frame( emb)
names( extended_embedding) = c( "X", "Y")
extended_embedding$color = cell.colors[ row.names( extended_embedding)]

# -- define the corrdinates of the basins of emission and attraction
basins = data.frame( xmin = c( -16*20/43,   -16*20/43,  -68*20/43,    -86*20/43,    38*20/43,  -110*20/43),
                      xmax = c( 17*20/43,    0*20/43,    -52*20/43,    -68*20/43,    61*20/43, -92*20/43),
                      ymin = c( -27*20/27,   6*20/27,    0*20/27,      -22*20/27,    6*20/27,  -29*20/27),
                      ymax = c( -10*20/27,   22*20/27,   12*20/27,     -13*20/27,    22*20/27, 0*20/27),
                      type = c( "emission",  "emission", "attraction", "attraction", "attraction", "attraction"))
basins$name = c( "BE1", "BE2", "BA1", "BA2", "BA3", "BA4")
row.names( basins) = basins$name

# -- plot the t-SNE embedding with both clusters and basins
ggplot( data = NULL) + 
  geom_point( data = extended_embedding, aes( x = X, y = Y, color = color)) +
  geom_rect( data = basins, aes( xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill = type), color="black", alpha=0.5) +
  geom_text( data = basins, aes( x=xmax, y=ymin, label=name))

# -- extract the cells associated to each basin
cells_in_basins = data.frame()
for( basin_name in row.names( basins)){
  cat("<BR>Basin:", basin_name)
  cell_indexes = which( extended_embedding$X >= basins[ basin_name, "xmin"] &
                        extended_embedding$X <= basins[ basin_name, "xmax"] &
                        extended_embedding$Y >= basins[ basin_name, "ymin"] &
                        extended_embedding$Y <= basins[ basin_name, "ymax"])
  cat("<BR>|-- Number of cells in basin=", length( cell_indexes))
  current_df = extended_embedding[ cell_indexes,]
  current_df$basin = rep( basin_name, nrow( current_df))
  cells_in_basins = rbind( cells_in_basins, current_df)
}

# -- write to file the list of cells in basins
write.table( cells_in_basins, file = file.path( OUTPUT_DIR, "cells_in_basins.tsv"), sep="\t",
             col.names = NA, row.names=TRUE, quote= FALSE)

# -- get the expression counts of the cells in the basins
cells_in_basins_counts_table = r_filtered$counts[ row.names( cells_in_basins),]

# -- create a Seurat object with the expression counts, adding the cell identity to its own basin name
seurat_count = CreateSeuratObject( t(cells_in_basins_counts_table), min.cells = 10, min.genes = 200, project="LysoDC")
for( current_basin in unique( cells_in_basins$basin)){
  cells_ids = row.names( cells_in_basins)[ which( cells_in_basins$basin == current_basin)]
  seurat_count = SetIdent( seurat_count, cells.use = cells_ids, ident.use = current_basin)  
}

# -- compute the marker genes between basins
all_markers = data.frame()
for( current_basin_1 in unique( cells_in_basins$basin)){
  for( current_basin_2 in unique( cells_in_basins$basin)){
    if( current_basin_1 != current_basin_2){
      current_markers = FindMarkers( seurat_count, ident.1 = "BE1", ident.2 = NULL, test.use = "bimod")
      current_markers$test = rep( paste( current_basin_1, "vs", current_basin_2, sep="_"), nrow( current_markers))
      all_markers = rbind( all_markers, current_markers)
    }
  }
}

# -- show the datatable of the results
datatable( all_markers, caption="List of DE genes between basins")

# -- write to file the result of markers genes between basins
write.csv( all_markers, file = file.path( OUTPUT_DIR, "basins_markers.csv"), row.names = TRUE, quote = FALSE)

