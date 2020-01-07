# ###############################################################
# This script aims to remove the cells that were identified
# previously as undesired in the analysis
# ###############################################################

# This script has been inspired by the R Velocyto tutorial page:
# http://pklab.med.harvard.edu/velocyto/notebooks/R/DG1.nb.html

## @knitr filter_data

# Get the cells to filter that were identified in primary analysis
contamination_cells = read.table( file= file.path( INPUT_DIR, "excluded_cells_contamination.txt"), header=TRUE, stringsAsFactors = FALSE)
highmito_cells = read.table( file= file.path( INPUT_DIR, "excluded_cells_HighMitoGenePerc.txt"), header=TRUE, stringsAsFactors = FALSE)
lowgenenb_cells = read.table( file= file.path( INPUT_DIR, "excluded_cells_LowGeneNb.txt"), header=TRUE, stringsAsFactors = FALSE)
lowuminb_cells = read.table( file= file.path( INPUT_DIR, "excluded_cells_LowUMINb.txt"), header=TRUE, stringsAsFactors = FALSE)

cat("<BR>Cell to remove due to high percentage of mitochondrial genes:", nrow( highmito_cells))
cat("<BR>Cell to remove due to low number of genes:", nrow( lowgenenb_cells))
cat("<BR>Cell to remove due to low number of UMI:", nrow( lowuminb_cells))
cat("<BR>Cell to remove due to contamination:", nrow( contamination_cells))

cells_id_to_excludes = unique( c( contamination_cells$cellid, highmito_cells$cellid, lowgenenb_cells$cellid, lowuminb_cells$cellid))
cells_name_to_exclude = paste0( "10635173:", cells_id_to_excludes, "x")

cat("<BR>Total number of cells excluded using those 4 criterias:", length( cells_id_to_excludes))

# Get the cells that are in undesired states (identified thank to previous clustering and mapping of gene expression)

# -- Load the files with cluster assignation from Velcyto analysis
cell_cluster_mapping_df = read.table( file= file.path( INPUT_DIR, "cell_cluster_mapping.tsv"), header=TRUE, stringsAsFactors = FALSE)

# -- Show the distribution of cells along all Velocyto clusters
cells_of_cluster_df = data.frame( table( cell_cluster_mapping_df[ , "cluster.number"]))
names( cells_of_cluster_df) = c( "ClusterID", "Number of cells")
datatable( cells_of_cluster_df, rownames = FALSE, caption="Cells in Velocyto clusters")

# -- Show the distribution of cells along all clusters to remove
cells_of_cluster_to_remove_df = data.frame( table( cell_cluster_mapping_df[ which( cell_cluster_mapping_df$cluster.number %in% CLUSTER_TO_REMOVE), "cluster.number"]))
names( cells_of_cluster_to_remove_df) = c( "ClusterID", "Number of cells")
datatable( cells_of_cluster_to_remove_df, rownames = FALSE, caption="Cells to remove by Velocyto clusters")

# -- get the cells to remove
clustered_cells_to_remove = cell_cluster_mapping_df$cell.id[ which( cell_cluster_mapping_df$cluster.number %in% CLUSTER_TO_REMOVE)]
cat("<BR>Total number of cells excluded using clustering criterias:", length( clustered_cells_to_remove))

# Groups all the undesired cells names
cells_name_to_exclude = unique( append( cells_name_to_exclude, clustered_cells_to_remove))
cat("<BR>Total number of cells excluded using all criterias:", length( cells_name_to_exclude))

cat("<HR>")
# Load the Velocyto data
input_file = file.path( VELOCYTO_INPUT_DIR,"10635173.loom")
ldat <- read.loom.matrices( input_file)
cat("<BR>Analyzed data are located on:", input_file)

# Get expression matrix
emat <- ldat$spliced

# Remove the undesired cells
filtered_emat = emat[ , -which( colnames( emat) %in% cells_name_to_exclude)]

# This is where one would do some filtering
filtered_emat <- filtered_emat[ , colSums(filtered_emat)>=1e3]

# PAGODA2 PRE-PROCESSING
# .......................

## @knitr build_filtered_data
cat("<BR>Building the required data for velocity analysis")

# Create the Pagoda2 object
rownames( filtered_emat) <- make.unique( rownames( filtered_emat))
r_filtered <- Pagoda2$new(filtered_emat,modelType='plain',trim=10,log.scale=T)
# Adjust the variances
r_filtered$adjustVariance(plot=T,do.par=T,gam.k=10)

# Run basic analysis steps to generate cell embedding and clustering, visualize:

# - Run the PCA
set.seed( 123)
cat("<HR>")
r_filtered$calculatePcaReduction(nPcs=100,n.odgenes=3e3,maxit=300)
cat("<HR>")
# - Cluster the cells
set.seed( 8457)
r_filtered$makeKnnGraph(k=50,type='PCA',center=T,distance='cosine');
set.seed( 5847)
r_filtered$getKnnClusters(method=walktrap.community, type='PCA', name='multilevel')
set.seed( 7984)
r_filtered$getEmbedding(type='PCA',embeddingType='tSNE',perplexity=10,verbose=FALSE)

# Compute differentially expressed genes
r_filtered$getDifferentialGenes(type='PCA', verbose=T, clusterType='multilevel')
#de <- r_filtered$diffgenes$PCA[[1]][['4']];
de = data.frame()
for( cluster_name in names( r_filtered$diffgenes[[ "PCA"]]$multilevel)){
  current_df = head( r_filtered$diffgenes[[ "PCA"]]$multilevel[[ cluster_name]], 10)
  current_df$cluster = rep( paste0( "cluster", cluster_name), nrow( current_df))
  de = rbind( de, current_df)
}
r_filtered$plotGeneHeatmap(genes = rownames( de), groups = r_filtered$clusters$PCA[[1]])
datatable( de, caption = "Differentially expressed genes by cluster")

# Plot embedding, labeling clusters
par( mfrow = c(1,1))
cat("<HR>")
r_filtered$plotEmbedding(type='PCA',embeddingType='tSNE',show.legend=F,mark.clusters=T,min.group.size=10,shuffle.colors=F,mark.cluster.cex=1,alpha=0.3,main='cell clusters')
cat("<HR>")

# Plot embedding, labeling marker genes expression
par(mfrow=c(2,2))
for( gene_name in c( "Emb", "F11r", "Plet1", "Tspan10", "Treml4", "Ccr2", "Clec4d", "Trem3", "Il12b", 
                     "Ccr7", "Ccl22", "Cenpa", "Ccnb2", "Mcm6", "Ccne1", "Ccne2",
                     "Birc5", "S100a6", "Pdia6", "Ranbp1", "Actb", "Gm2a", "Itm2b", "Mif", "Cd63", "Fth1", "Cd24a",
                     "Clec4a1", "Clec4a2", "Clec4a4", "Clec4e", "Clec4n", "Clec5a")){
  if( gene_name %in% colnames( r_filtered$counts)){
    r_filtered$plotEmbedding(type='PCA',embeddingType='tSNE',colors=r_filtered$counts[, gene_name], main= gene_name)
  }else{
    cat("<BR>The gene", gene_name, "is not present in filtered matrix<BR>")
  }
}
par(mfrow=c(1,1))
cat("<HR>")


