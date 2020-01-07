# ########################################################
# This script aims to analyse the heterogeneity of cells
#  of cells using dimension reduction techniques
# ########################################################

## @knitr heterogeneity_analysis_by_pca
set.seed( 123)
cat("<H4>Cluster identification and analysis</H4>")
sc10x <- RunPCA(object = sc10x, pc.genes = selected_variable_genes, do.print = FALSE, pcs.print = 1:5, genes.print = 5)

# Compute the score of the cells according to various gene signatures.
#for( signature in names( SIGNATURES_LIST)){
#  sc10x <- AddModuleScore( object = sc10x, genes.list = list( SIGNATURES_LIST[[ signature]]),
#                           ctrl.size = length( SIGNATURES_LIST[[ signature]]), enrich.name = signature)
#}

# Compute the clusters
sc10x <- FindClusters(object = sc10x, reduction.type = "pca", dims.use = 1:10,
                      resolution = 0.45, print.output = 0, save.SNN = TRUE, force.recalc=TRUE,
                      temp.file.location = "/tmp/")

# Build a dataframe with merged information to be able to use ggplot		  
sc10x_merged_df = as.data.frame(t(as.matrix(sc10x@data)))
sc10x_merged_df$UniqueCellID = rownames(sc10x_merged_df)
sc10x_merged_df$nUMI = sc10x@meta.data$nUMI
sc10x_merged_df$nGene = sc10x@meta.data$nGene
sc10x_merged_df$outCells = sc10x@meta.data$outlier
sc10x_merged_df$SeuratClusters = as.vector(sc10x@ident)

#pca_result = pca_hist_function(sc10x_merged_df[,selected_variable_genes], PlotHist=FALSE,
#                               VectorToColorPCA=sc10x_merged_df$SignatureScore,
#                               NameOfVectorToColor="SignatureScore", dim_all=TRUE,
#                               title="data filtered PCA")

# Show the number of cells in each clusters
cluster_count_df = data.frame( table( sc10x_merged_df$SeuratClusters))
names( cluster_count_df) = c( "ClusterNumber", "CellCount")
datatable( cluster_count_df, caption="Distribution of cells in clusters")

# Compute the PCA and plot the PCA representation in PC1/PC2/PC3
cat("\n<H4>Ditribution of clusters on PCA map and PCA loadings</H4>")

pca_result_nUMI = pca_hist_function(sc10x_merged_df[,selected_variable_genes], PlotHist=FALSE,
                                    VectorToColorPCA=sc10x_merged_df$nUMI, NameOfVectorToColor="nUMI", dim_all=TRUE,
                                    title="data filtered PCA")

pca_result_nGene = pca_hist_function(sc10x_merged_df[,selected_variable_genes], PlotHist=FALSE,
                                     VectorToColorPCA=sc10x_merged_df$nGene, NameOfVectorToColor="nGene", dim_all=TRUE,
                                     title="data filtered PCA")

pca_result = pca_hist_function(sc10x_merged_df[,selected_variable_genes], PlotHist=FALSE,
                               VectorToColorPCA=sc10x_merged_df$SeuratClusters, NameOfVectorToColor="SeuratClusters", dim_all=TRUE,
                               title="data filtered PCA")

# Plot the loadings of PCA variables (genes) in PC1/PC2/PC3
gg_plots = pca_gene_information_function( pca_result$pca_info, dim_all=TRUE, threshold_norm=0.5)

# Compute the t-SNE mapping
tsne_coord_sc10x_merged_df = run_tsne_function(sc10x_merged_df[, selected_variable_genes], perplexity_value=30)

# Show the cluster on the t-SNE embedding
cat("\n<H4>Ditribution of clusters on t-SNE map</H4>")
plot_tsne_function(sc10x_merged_df[,selected_variable_genes], 
                   tsne_coord_sc10x_merged_df,
                   VectorToColor = sc10x_merged_df$SeuratClusters,
                   NameOfVectorToColor = "SeuratClusters",
                   title="data filtered t-SNE")
		   
# Show the expression of monitored genes on the t-SNE embedding and on a splitted histogram
cat("\n<H4>Ditribution of expression of monitored genes on t-SNE map</H4>")
max_expression = max( sc10x_merged_df[,selected_variable_genes])
for( gene_name in MONITORED_GENES){
  # Plot the t-SNE map with gradient color with the expression of selected gene
  tsne_gene_plot = plot_tsne_function(sc10x_merged_df[, selected_variable_genes], 
                     tsne_coord_sc10x_merged_df,
                     VectorToColor = sc10x_merged_df[, gene_name],
                     NameOfVectorToColor = gene_name,
                     title="Data filtered t-SNE")
  print( tsne_gene_plot)
  
  # Plot the distribution of selected gene expression along clusters
  gene_expression = sc10x_merged_df[ which( sc10x_merged_df[ , gene_name] > 0), c( gene_name, "SeuratClusters")]
  names( gene_expression) = c( "gene_exp", "SeuratClusters")
  print( ggplot( data = gene_expression) + 
           geom_histogram( aes( x=gene_exp, fill= SeuratClusters)) + facet_grid( SeuratClusters ~ .) +
           xlim( 0, max_expression) +
           xlab( paste( gene_name,"expression level")) + 
           theme( legend.position = "None"))
}

gc()

# Plot the heatmap of top markers genes
cat("\n<H4>Heatmap of top markers genes through clusters</H4>")
sc10x.markers <- FindAllMarkers(object = sc10x, only.pos = FALSE, min.pct = 0.25,
                                thresh.use = 0.25, print.bar = FALSE)

top_markers_df = sc10x.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

DoHeatmap(object = sc10x, genes.use = top_markers_df$gene, slim.col.label = TRUE, remove.key = TRUE, cex.row = 5)

# Plot the heatmap of monitored genes
cat("\n<H4>Heatmap of monitored genes through clusters</H4>")
DoHeatmap(object = sc10x, genes.use = MONITORED_GENES, slim.col.label = TRUE, remove.key = TRUE, cex.row = 5)


