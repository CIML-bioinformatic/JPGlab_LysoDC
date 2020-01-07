# ################################################
# This script aims to read and filter data
# ################################################

# READ THE DATA
# -----------------------

## @knitr load_data

PATH_10x_DATA = file.path( RAW_DATA_DIR, "10635173", "outs", "filtered_gene_bc_matrices", "mm10")
ORIGINAL_SC10X.data <- Read10X(PATH_10x_DATA)
ORIGINAL_SC10X <-CreateSeuratObject(raw.data = ORIGINAL_SC10X.data, min.cells = 3, min.genes = 200, project = "10X_DATA")

cat("<br>Analyzed data are located on:", PATH_10x_DATA)

# FILTER THE DATA
# -----------------------

## @knitr filtering_data

ALL_GENES = rownames(ORIGINAL_SC10X@data)

# Identify the mitocondrial genes and compute their percentage in each cells
mito.genes <- grep(pattern = "^mt-|^MT-", x = rownames(x = ORIGINAL_SC10X@data), value = TRUE)
percent.mito <- colSums(as.matrix(ORIGINAL_SC10X@raw.data[mito.genes, ]))/colSums(as.matrix(ORIGINAL_SC10X@raw.data))

# Add the mitocondrial gene percentage as meta information in the Seurat object 
ORIGINAL_SC10X <- AddMetaData(object = ORIGINAL_SC10X, metadata = percent.mito, col.name = "percent.mito")

# Plot the scatterplot with number of UMI and percentage of mito genes
vln_plot = VlnPlot(object = ORIGINAL_SC10X, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
print(vln_plot)

# Locate the cells with low number of UMI
nUMI.drop <- isOutlier_MAD_function(ORIGINAL_SC10X@meta.data$nUMI, nmads=3, type="both", log=TRUE, title="Nb of UMI")
nUMI.drop = (ORIGINAL_SC10X@meta.data$nUMI < 10^3.3)
cat("<br>Number of cells removed from Nb of UMI filter:", sum( nUMI.drop))

# Locate the cells with low number of expressed genes
# nGene.drop <- isOutlier_MAD_function(ORIGINAL_SC10X@meta.data$nGene, nmads=4, type="both", log=TRUE, title="Nb of genes")
nGene.drop = (ORIGINAL_SC10X@meta.data$nGene > 5000)
cat("<br>Number of cells removed from Nb of genes filter:", sum( nGene.drop))

# Locate the cells with high percentage of mitocondrial genes
# mito.drop <- isOutlier_MAD_function(ORIGINAL_SC10X@meta.data$percent.mito, nmads=4, type="both", log=TRUE, title="Percent Mito Genes")
mito.drop = (ORIGINAL_SC10X@meta.data$percent.mito > 0.4)
cat("<br>Number of cells removed from Percentage of mito genes filter:", sum( mito.drop))

# Identify the cells to exclude as union of cells with low nb UMI, low nb of expressed genes or high percentage of mitocondrial genes
ORIGINAL_SC10X@meta.data$outlier = nUMI.drop | nGene.drop | mito.drop
cat("<br>Removed cells after filters:", length( which( ORIGINAL_SC10X@meta.data$outlier == TRUE)))
cat("<br>Remaining cells after filters:", length( which( ORIGINAL_SC10X@meta.data$outlier == FALSE)))

# Export the excluded cells to file
write.table( data.frame( cellid = colnames(ORIGINAL_SC10X@data)[nUMI.drop]), file= file.path( OUTPUT_DIR, "excluded_cells_LowUMINb.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")
write.table( data.frame( cellid = colnames(ORIGINAL_SC10X@data)[nGene.drop]), file= file.path( OUTPUT_DIR, "excluded_cells_LowGeneNb.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")
write.table( data.frame( cellid = colnames(ORIGINAL_SC10X@data)[mito.drop]), file= file.path( OUTPUT_DIR, "excluded_cells_HighMitoGenePerc.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")

# Plot the dispersion of excluded and non-excluded cells
gg_plot = ggplot(ORIGINAL_SC10X@meta.data, aes(nGene,percent.mito, color=outlier))+geom_point()+ggtitle("Outlier cells (we will remove them)")
print(gg_plot)

# Filter the excluded cells in the Seurat object
ORIGINAL_SC10X <- FilterCells(object = ORIGINAL_SC10X, subset.names = c(), cells.use = colnames(ORIGINAL_SC10X@data)[!ORIGINAL_SC10X@meta.data$outlier])

# NORMALIZE THE DATA
# -----------------------

## @knitr normalizing_data

NORMALIZED_FILTERED_SC10X <- NormalizeData(object = ORIGINAL_SC10X, normalization.method = "LogNormalize", scale.factor = 10000, display.progress = FALSE)
NORMALIZED_FILTERED_SC10X <- ScaleData(object = NORMALIZED_FILTERED_SC10X,
                                       do.center=TRUE, do.scale=FALSE,
                                       vars.to.regress = c( "nUMI"),
                                       display.progress = FALSE)

rm( "ORIGINAL_SC10X")
gc()
