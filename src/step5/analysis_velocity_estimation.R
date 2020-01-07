# ####################################################
# This script aim to use the Velocyto R pipeline
# to analyse the RNA velocity on LysoDC data
# ####################################################

# This script has been inspired by the R Velocyto tutorial page:
# http://pklab.med.harvard.edu/velocyto/notebooks/R/DG1.nb.html


## @knitr estimate_velocity

# Get the raw data for spliced and unspliced RNA
emat <- ldat$spliced
nmat <- ldat$unspliced

# Restrict to cells that passed pagoda2 filter
emat <- emat[,rownames(r$counts)]
nmat <- nmat[,rownames(r$counts)]

# Take cluster labels from Pagoda2 pre-processing
cluster.label <- r$clusters$PCA[[1]]
cell.colors <- pagoda2:::fac2col(cluster.label)

# Take t-SNE embedding from Pagoda2 pre-processing
emb <- r$embeddings$PCA$tSNE

# In addition to clustering and the t-SNE embedding, from the pagoda2 processing 
# we will also take a cell-cell distance, which will be better than the default 
# whole-transcriptome correlation distance that velocyto.R would normally use.
cell.dist <- as.dist(1-armaCor(t(r$reductions$PCA)))

# Filter genes based on the minimum average expresion magnitude 
# (in at least one of the clusters), output total number of resulting valid genes
emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.5)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)
cat("<BR>Number of validated genes for analysis =",length(intersect(rownames(emat),rownames(nmat))), "<BR>")

# Estimate RNA velocity (using gene-relative model with k=20 cell kNN pooling and 
# using top/bottom 2% quantiles for gamma fit)
fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=20,cell.dist=cell.dist,fit.quantile=fit.quantile)

# Visualize velocity on the t-SNE embedding, using velocity vector fields
cat("<HR>")
par( mfrow = c(1,1))
show.velocity.on.embedding.cor(emb,
                               rvel.cd,
                               n=300,
                               scale='sqrt',
                               cell.colors=ac(cell.colors,alpha=0.5),
                               cex=0.8,
                               arrow.scale=5,
                               show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5,
                               grid.n=40,
                               arrow.lwd=1,
                               do.par=F,
                               cell.border.alpha = 0.1)

# Visualize a fit for maker genes (we reuse rvel.cd to save on calcualtions here):
cat("<HR>")
par( mfrow = c(2,2))
excluded_genes = vector()
for( gene_name in c("Emb", "F11r", "Plet1", "Tspan10", "Treml4", "Ccr2", "Clec4d", "Trem3", "Il12b", 
                    "Ccr7", "Ccl22", "Cenpa", "Ccnb2", "Mcm6", "Ccne1", "Ccne2",
                    "Birc5", "S100a6", "Pdia6", "Ranbp1", "Actb", "Gm2a", "Itm2b", "Mif", "Cd63", "Fth1")){
  if( gene_name %in% rownames( emat) && gene_name %in% rownames( nmat)){
	cat("<HR>")
    gene.relative.velocity.estimates(emat, nmat, deltaT=1, kCells = 20, kGenes=1,
                                     fit.quantile=fit.quantile, cell.emb=emb, cell.colors=cell.colors,
                                     cell.dist=cell.dist, old.fit=rvel.cd, do.par=FALSE,
                                     show.gene= gene_name)
  }else{
    cat("<HR><BR>The gene", gene_name, "is not present in the matrices : can't display graphs.")
  }
}
par( mfrow = c(1,1))
