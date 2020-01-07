# ################################################
# This script aims to read and filter data
# ################################################

# This script has been inspired by the R Velocyto tutorial page:
# http://pklab.med.harvard.edu/velocyto/notebooks/R/DG1.nb.html

# READ THE DATA
# -----------------------

## @knitr load_data
input_file = file.path( VELOCYTO_INPUT_DIR,"10635173.loom")
ldat <- read.loom.matrices( input_file)
cat("<BR>Analyzed data are located on:", input_file)

# Get expression matrix
emat <- ldat$spliced

# Filter some very low values
emat <- emat[,colSums(emat)>=1e3]

# Pagoda2 pre-processing
# ------------------

## @knitr build_data
cat("<BR>Building the required data for velocity analysis")

# Create the Pagoda2 object
rownames( emat) <- make.unique( rownames( emat))
r <- Pagoda2$new(emat,modelType='plain',trim=10,log.scale=T)
# Adjust the variances
r$adjustVariance(plot=T,do.par=T,gam.k=10)

# Run basic analysis steps to generate cell embedding and clustering, visualize:

# - Run the PCA
set.seed( 123)
cat("<HR>")
r$calculatePcaReduction(nPcs=100,n.odgenes=3e3,maxit=300)
cat("<HR>")
# - Cluster the cells
r$makeKnnGraph(k=40,type='PCA',center=T,distance='cosine');
r$getKnnClusters(method=multilevel.community, type='PCA', name='multilevel')
r$getEmbedding(type='PCA',embeddingType='tSNE',perplexity=50,verbose=FALSE)

# - Look at the distribution of cells along cluster
cluster_df = data.frame( table( r$clusters[[ "PCA"]]$multilevel))
names( cluster_df) = c( "Cluster ID", "Cell Number")
datatable( cluster_df, caption="Distribution of cells along clusters", rownames=FALSE)

# Compute differentially expressed genes
r$getDifferentialGenes(type='PCA', verbose=T, clusterType='multilevel')
de = data.frame()
for( cluster_name in names( r$diffgenes[[ "PCA"]]$multilevel)){
  current_df = head( r$diffgenes[[ "PCA"]]$multilevel[[ cluster_name]], 10)
  current_df$cluster = rep( paste0( "cluster", cluster_name), nrow( current_df))
  de = rbind( de, current_df)
}
cat("<HR>")
r$plotGeneHeatmap(genes = rownames( de), groups = r$clusters$PCA[[1]])
cat("<HR>")
datatable( de, caption = "Differentially expressed genes by cluster")
cat("<HR>")

# Plot embedding, labeling clusters
par( mfrow = c(1,1))
cat("<HR>")
r$plotEmbedding(type='PCA',embeddingType='tSNE',show.legend=F,mark.clusters=T,min.group.size=10,shuffle.colors=F,mark.cluster.cex=1,alpha=0.3,main='cell clusters')
cat("<HR>")

# Plot embedding, labeling marker genes expression
par(mfrow=c(2,2))
cat("<HR>")
r$plotEmbedding(type='PCA',embeddingType='tSNE',colors=r$counts[,"Emb"], main='Emb')
r$plotEmbedding(type='PCA',embeddingType='tSNE',colors=r$counts[,"F11r"], main='F11r')
r$plotEmbedding(type='PCA',embeddingType='tSNE',colors=r$counts[,"Plet1"], main='Plet1')
r$plotEmbedding(type='PCA',embeddingType='tSNE',colors=r$counts[,"Tspan10"], main='Tspan10')
cat("<HR>")
r$plotEmbedding(type='PCA',embeddingType='tSNE',colors=r$counts[,"Treml4"], main='Treml4')
r$plotEmbedding(type='PCA',embeddingType='tSNE',colors=r$counts[,"Ccr2"], main='Ccr2')
r$plotEmbedding(type='PCA',embeddingType='tSNE',colors=r$counts[,"Clec4d"], main='Clec4d')
r$plotEmbedding(type='PCA',embeddingType='tSNE',colors=r$counts[,"Trem3"], main='Trem3')
cat("<HR>")
r$plotEmbedding(type='PCA',embeddingType='tSNE',colors=r$counts[,"Il12b"], main='Il12b')
r$plotEmbedding(type='PCA',embeddingType='tSNE',colors=r$counts[,"Ccr7"], main='Ccr7')
r$plotEmbedding(type='PCA',embeddingType='tSNE',colors=r$counts[,"Ccl22"], main='Ccl22')
par(mfrow=c(1,1))
cat("<HR>")

