# ################################################
# This script aims to read and filter data
# ################################################

# This script is inpired by the Monocle tutorial:
# http://cole-trapnell-lab.github.io/monocle-release/docs/

# READ THE DATA
# -----------------------

## @knitr load_data
input_file = file.path( INPUT_DIR, "filtered_raw_expression_matrix.csv")
filtered_expression_df = read.table( file = input_file, header = TRUE, sep="\t")
cat("<BR>Analyzed data are located on:", input_file)
cat("<BR>Number of cells:", ncol( filtered_expression_df))
cat("<BR>Number of genes:", nrow( filtered_expression_df))

# Import the pagoda object with clusters and t_SNE embedding
r_filtered = readRDS( file.path( INPUT_DIR, "r_filtered.rds"))
# Take t-SNE embedding from Pagoda2 pre-processing
emb <- r_filtered$embeddings$PCA$tSNE
# Load the clusters
cluster_label <- r_filtered$clusters$PCA[[1]]

# BUILD THE DATA
# ................

## @knitr build_data
cat("<H5>Building the required data for Monocle</H5>")

feature_df = AnnotatedDataFrame( data.frame( gene_short_name = row.names( filtered_expression_df)))
row.names( feature_df) = row.names( filtered_expression_df)

# Build the CDS (Cell data set) required by monocle
filtered_cds <- newCellDataSet( as.matrix( filtered_expression_df), 
                                featureData = feature_df,
                                expressionFamily=negbinomial.size())

# Compute size factors and dispersion
filtered_cds = estimateSizeFactors( filtered_cds)
filtered_cds = estimateDispersions( filtered_cds)

# FILTER CELLS AND GENES
# ......................

## @knitr filter_data
cat("<H5>Filtering the low quality data</H5>")

# Select the genes with at least 10 cells expressing it at a minimum of 0.1
filtered_cds <- detectGenes( filtered_cds, min_expr = 0.1)
EXPRESSED_GENES <- row.names( subset( fData( filtered_cds), num_cells_expressed >= 10))
print( datatable( fData( filtered_cds)[ EXPRESSED_GENES, ]))

# Look at the total number of mRNA that display the limits of mean +/- 2sd
pData( filtered_cds)$Total_UMI <- Matrix::colSums( Biobase::exprs( filtered_cds))
upper_bound <- 10^(mean(log10(pData( filtered_cds)$Total_UMI)) +
                     2*sd(log10(pData( filtered_cds)$Total_UMI)))
lower_bound <- 10^(mean(log10(pData( filtered_cds)$Total_UMI)) -
                     2*sd(log10(pData( filtered_cds)$Total_UMI)))

print( qplot( Total_UMI, data = pData( filtered_cds), geom = "density") +
         geom_vline(xintercept = lower_bound) +
         geom_vline(xintercept = upper_bound)
)

# Look if the distribution of the expression of expressed genes is lognormal
# -- Log-transform each value in the expression matrix.
log_exprs <- log( Biobase::exprs( filtered_cds)[ EXPRESSED_GENES,])

# -- Standardize each gene, so that they are all on the same scale,
# -- Then melt the data with plyr so we can plot it easily
melted_dens_df <- melt( Matrix::t( scale( Matrix::t( log_exprs))))

# -- Plot the distribution of the standardized gene expression values.
print( qplot(value, geom = "density", data = melted_dens_df) +
         stat_function(fun = dnorm, size = 0.5, color = 'red') +
         xlab("Standardized log( expr)") +
         ylab("Density")
)

# PLOT THE VALUES OF SOME GENES ALONG CLUSTERS
# ............................................

cat("<H5>Plot values of some genes across clusters</H5>")

# Get the expression values
# -- To use raw values
#counts_table_with_clusters = data.frame( t( exprs( filtered_cds)))
# -- To use normalized values(scaled with computed size factor)
counts_table_with_clusters = data.frame( t(exprs( filtered_cds)) /  pData( filtered_cds)[, 'Size_Factor'])

# Get the cell clusters with right cell names
names( cluster_label) = gsub( "x", "", gsub( "10635173:", "", names( cluster_label)))
# Get cell associated color
cell_colors <- pagoda2:::fac2col( cluster_label)

# Get the cells in common between clusters and expression
selected_cells = intersect( names( cluster_label), row.names( counts_table_with_clusters))

# Build the dataframe with merged information
counts_table_with_clusters = counts_table_with_clusters[ selected_cells, ]
counts_table_with_clusters$cluster = cluster_label[ row.names( counts_table_with_clusters)]
counts_table_with_clusters$cell.color = adjustcolor( cell_colors[ row.names( counts_table_with_clusters)], alpha.f = 1)

# Kepp only cluster 2, 3, 4 and 5
counts_table_with_clusters = counts_table_with_clusters[ which( counts_table_with_clusters$cluster %in% c( 2,3,4,5)), ]

# Reorder the cluster to have progression 1->2->3->4->6->5
counts_table_with_clusters$cluster = factor( counts_table_with_clusters$cluster, levels= c( 4, 3, 2, 5))

# Get the colors of clusters
cluster_colors = unique( counts_table_with_clusters$cell.color)[ order( unique( counts_table_with_clusters$cluster), decreasing = FALSE)]

# Define the list of gene to be plotted
gene_set = sort( c( "Emb" ,"F11r", "Plet1", "Tspan10", "Treml4", "Clec4d", "Trem3", "Il12b",
                    "Ccr7", "Ccl22", "Cenpa", "Ccnb2", "Mcm6", "Ccne1", "Birc5", "S100a6",
                    "Actb", "Gm2a", "Itm2b", "Mif", "Fth1", "Apoc4", "Clec4a2", "Susd2", 
                    "Cd9", "Clec4a4", "Tnnt1", "Ccnd3", "Wfdc17", "Igf1", 
                    "Apoc4", "Tlr12", "Il1b", "Ceacam19", "Ifitm1", "Klra17", "Cd24a"
))

# Plot the expression of selected genes along clusters
for( gene_name in gene_set){
  current_plot = ggplot( data = counts_table_with_clusters, aes_string( x="cluster", y=gene_name)) + 
                      geom_violin( aes( col=cluster)) +
                      geom_jitter( aes( col=cluster)) +
                      scale_fill_manual( values = cluster_colors) +
                      ggtitle( paste( gene_name,"expression along clusters")) +
                      theme(legend.position="none") + 
                      scale_x_discrete(labels=c("4" = "Cluster 3", "3" = "Cluster4", "2" = "Cluster 5", "5" = "Cluster 6"))
                      
  print( current_plot)
}
