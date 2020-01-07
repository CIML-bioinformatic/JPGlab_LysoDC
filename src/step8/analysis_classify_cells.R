# ###############################################################
# This script aims to use Monocle method to classify cells
# with or without using predefined markers.
# ###############################################################

# This script is inpired by the Monocle tutorial:
# http://cole-trapnell-lab.github.io/monocle-release/docs/

# ###############################
# Classify cells by type
# ###############################

## @knitr classify_cells_by_type
set.seed( 123)

# Define the markers genes used to classify the cells
APOC4_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Apoc4"))
CCNB2_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Ccnb2"))
CCNE1_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Ccne1"))
CD24A_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Cd24a"))
CENPA_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Cenpa"))
CLEC4D_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Clec4d"))
EMB_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Emb"))
F11R_id <- row.names(subset(fData( filtered_cds), gene_short_name == "F11r"))
MCM6_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Mcm6"))
PLET1_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Plet1"))
CD9_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Cd9"))
TLR12_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Tlr12"))
TREM3_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Trem3"))
TSPAN10_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Tspan10"))

# Build the cell type hierarchy
cth_FILTERED<- newCellTypeHierarchy()

cat("<BR>SP= Ccne1-Mcm6-Cenpa-Ccnb2-Emb+Trem3+F11r-Cd24a-Plet1-Tspan10-Apoc4-")
cth_FILTERED <- addCellType( cth_FILTERED, "SP", 
                             classify_func = function(x) { 
                               x[CCNE1_id,] < 1 & x[ MCM6_id,] < 1 & x[ CENPA_id,] < 1 & x[ CCNB2_id,] <1 &
                                 x[ EMB_id,] >= 1 & x[ TREM3_id,] >= 1 & x[F11R_id,] < 1 & x[ CD24A_id,] < 1 &
                                 x[ PLET1_id,] < 1 & x[ TSPAN10_id,] <1 & x[APOC4_id,] < 1 })  

cat("<BR>DP= Emb+F11r+Plet1-Tspan10-")
cth_FILTERED <- addCellType( cth_FILTERED, cell_type_name="DP", 
                               classify_func = function(x) { 
                                 x[ EMB_id,] >= 1 & x[F11R_id,] >= 1 &
                                 x[ PLET1_id,] < 1 & x[ TSPAN10_id,] <1})

cat("<BR>TP= Emb+ F11r+ Plet1+ Tspan10-")
cth_FILTERED <- addCellType( cth_FILTERED, cell_type_name="TP",
                             classify_func = function(x) {
                               x[EMB_id,] >= 1 & x[ F11R_id,] >= 1 &
                               x[ PLET1_id,] >= 1 & x[ TSPAN10_id,] <1})

cat("<BR>TN1= Ccne1-Mcm6-Cenpa-Ccnb2-Emb-Trem3-F11r-Cd24a-Plet1-Clec4d-Tspan10+Apoc4+")
cth_FILTERED <- addCellType( cth_FILTERED, cell_type_name="TN",
                             classify_func = function(x) { 
                               x[CCNE1_id,] < 1 & x[ MCM6_id,] < 1 & x[ CENPA_id,] < 1 & x[ CCNB2_id,] <1 &
                               x[ EMB_id,] < 1 & x[ TREM3_id,] < 1 & x[F11R_id,] < 1 & x[ CD24A_id,] < 1 &
                               x[ PLET1_id,] < 1 & x[ CLEC4D_id,] < 1 & x[ TSPAN10_id,] >= 1 & x[APOC4_id,] >= 1 })


SELECTED_HIERARCHY = cth_FILTERED

# Apply the cell type hierarchy to classify cells
filtered_cds_cell_type <- classifyCells( filtered_cds, SELECTED_HIERARCHY)

# Plot the proportion of classified cells
cat("<HR>")
pie <- ggplot(pData( filtered_cds_cell_type), aes(x = factor(1), fill = factor( CellType))) + geom_bar( width = 1) +
      coord_polar( theta = "y") + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
print( pie)

cat("<HR>")
cat( "<BR>Total number of cell to classify:", length( pData( filtered_cds_cell_type)$CellType))
cat( "<BR>Number of cells classified as SP:", length( which( pData( filtered_cds_cell_type)$CellType == "SP")))
cat( "<BR>Number of cells classified as DP:", length( which( pData( filtered_cds_cell_type)$CellType == "DP")))
cat( "<BR>Number of cells classified as TP:", length( which( pData( filtered_cds_cell_type)$CellType == "TP")))
cat( "<BR>Number of cells classified as TN:", length( which( pData( filtered_cds_cell_type)$CellType == "TN")))
cat( "<BR>Number of cells not classified:", length( which( pData( filtered_cds_cell_type)$CellType == "Unknown")))



# ##########################################  
# Classify cells using markers
# ##########################################

## @knitr classify_cells_with_markers
set.seed( 123)

# Select the genes that co-vary with the markers genes used to classify cells
marker_diff <- markerDiffTable( filtered_cds_cell_type[ EXPRESSED_GENES,],
                                SELECTED_HIERARCHY,
                                # residualModelFormulaStr = "~Media + num_genes_expressed",
                                cores = 1)

# Look at the most specific gene for each cell type
cell_type_candidate_clustering_genes <- row.names( subset( marker_diff, qval < 0.01))
marker_spec <- calculateMarkerSpecificity( filtered_cds_cell_type[ cell_type_candidate_clustering_genes,], SELECTED_HIERARCHY)
print( datatable( selectTopMarkers( marker_spec, 10)))

# Keep the 500 most specific genes to cluster the cells
semisup_clustering_genes <- unique( selectTopMarkers( marker_spec, 200)$gene_id)
filtered_cds_cell_type <- setOrderingFilter( filtered_cds_cell_type, semisup_clustering_genes)
plot_ordering_genes( filtered_cds_cell_type)
plot_pc_variance_explained( filtered_cds_cell_type, return_all = F)

# Build the T-SNE map
set.seed( 123)
filtered_cds_cell_type <- reduceDimension( filtered_cds_cell_type, max_components = 2, num_dim = 6,
                        norm_method = 'log',
                        reduction_method = 'tSNE',
                        # residualModelFormulaStr = "~Media + num_genes_expressed",
                        verbose = T)

# Find the cell clusters
filtered_cds_cell_type <- clusterCells( filtered_cds_cell_type, method = "densityPeak", num_clusters = 5)

# Plot the data by cell types with one color per cell type
plot_cell_clusters( filtered_cds_cell_type, 1, 2, color = "CellType") +
  scale_colour_manual(values = c("violet", "blue", "green", "red", "grey")) +
  facet_wrap(~CellType)

# Plot the data by cell types with one color per cluster
plot_cell_clusters( filtered_cds_cell_type, 1, 2, color = "Cluster") +
  facet_wrap(~CellType)

# Plot the expression of markers in data
plot_cell_clusters( filtered_cds_cell_type, 1, 2, cell_size=1, color = "CellType", markers = c("Emb", "F11r", "Plet1", "Tspan10")) + 
                       scale_colour_gradient( low = "grey", high = "blue")

# Plot the expression of markers in data per cell type
plot_cell_clusters( filtered_cds_cell_type, 1, 2, cell_size=1, color = "CellType", markers = c( "Emb")) + 
    scale_colour_gradient( low = "grey", high = "blue") + facet_wrap(~CellType) +
    ggtitle( paste( "Dispersion of", "Emb", "in identified populations"))
plot_cell_clusters( filtered_cds_cell_type, 1, 2, cell_size=1, color = "CellType", markers = c( "F11r")) + 
  scale_colour_gradient( low = "grey", high = "blue") + facet_wrap(~CellType) +
  ggtitle( paste( "Dispersion of", "F11r", "in identified populations"))
plot_cell_clusters( filtered_cds_cell_type, 1, 2, cell_size=1, color = "CellType", markers = c( "Plet1")) + 
  scale_colour_gradient( low = "grey", high = "blue") + facet_wrap(~CellType) +
  ggtitle( paste( "Dispersion of", "Plet1", "in identified populations"))
plot_cell_clusters( filtered_cds_cell_type, 1, 2, cell_size=1, color = "CellType", markers = c( "Tspan10")) + 
  scale_colour_gradient( low = "grey", high = "blue") + facet_wrap(~CellType) +
  ggtitle( paste( "Dispersion of", "Tspan10", "in identified populations"))

