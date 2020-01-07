# ####################################################
# This script aim to use the Velocyto R pipeline
# to analyse the RNA velocity on LysoDC data
# ####################################################

## @knitr estimate_filtered_velocity

# Get the raw data for spliced and unspliced RNA
emat <- ldat$spliced
nmat <- ldat$unspliced

# Restrict to cells that passed pagoda2 filter
emat <- emat[,rownames(r_filtered$counts)]
nmat <- nmat[,rownames(r_filtered$counts)]

# Filtered out the undesired cells
emat_columns_to_exclude = which( colnames( emat) %in% cells_name_to_exclude)
if( length( emat_columns_to_exclude) > 0){ 
  filtered_emat = emat[ , -emat_columns_to_exclude]
}else{
  filtered_emat = emat
}
nmat_columns_to_exclude = which( colnames( nmat) %in% cells_name_to_exclude)
if( length( nmat_columns_to_exclude) > 0){ 
  filtered_nmat = nmat[ , -nmat_columns_to_exclude]
}else{
  filtered_nmat = nmat
}

# Take cluster labels from Pagoda2 pre-processing
cluster.label <- r_filtered$clusters$PCA[[1]]
cell.colors <- pagoda2:::fac2col(cluster.label)

# Take t-SNE embedding from Pagoda2 pre-processing
emb <- r_filtered$embeddings$PCA$tSNE

# In addition to clustering and the t-SNE embedding, from the pagoda2 processing 
# we will also take a cell-cell distance, which will be better than the default 
# whole-transcriptome correlation distance that velocyto.R would normally use.
cell.dist <- as.dist(1-armaCor(t(r_filtered$reductions$PCA)))

# Filter genes based on the minimum average expresion magnitude 
# (in at least one of the clusters), output total number of resulting valid genes
filtered_emat <- filter.genes.by.cluster.expression(filtered_emat,cluster.label,min.max.cluster.average = 0.05)
filtered_nmat <- filter.genes.by.cluster.expression(filtered_nmat,cluster.label,min.max.cluster.average = 0.05)
cat("<BR>Number of validated genes for analysis =",length(intersect(rownames(filtered_emat),rownames(filtered_nmat))), "<BR>")

# Estimate RNA velocity (using gene-relative model with k=20 cell kNN pooling and 
# using top/bottom 2% quantiles for gamma fit)
fit.quantile <- 0.02
set.seed( 3549)
rvel.cd <- gene.relative.velocity.estimates(filtered_emat,
                                            filtered_nmat,
                                            deltaT=1,
                                            kCells=20,
                                            cell.dist=cell.dist,
                                            fit.quantile=fit.quantile,
                                            cell.emb=emb)

# Visualize velocity on the t-SNE embedding, using velocity vector fields
cat("<HR>")
par( mfrow = c(1,1))
cat("<H5>n=300, arrow.scale=5, grid.n=40</H5>")
show.velocity.on.embedding.cor(emb,
                               rvel.cd,
                               scale='sqrt',
                               cell.colors=ac(cell.colors,alpha=0.5),
                               cex=0.8,
                               show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5,
                               fixed.arrow.length = TRUE,
                               arrow.lwd=1,
                               do.par=F,
                               cell.border.alpha = 0.1,
                               n=300, # Test
                               arrow.scale=5, # Test
                               grid.n=40 # Test
                               )

cat("<H5>n=100, arrow.scale=5, grid.n=40</H5>")
show.velocity.on.embedding.cor(emb,
                               rvel.cd,
                               scale='sqrt',
                               cell.colors=ac(cell.colors,alpha=0.5),
                               cex=0.8,
                               show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5,
                               fixed.arrow.length = TRUE,
                               arrow.lwd=1,
                               do.par=F,
                               cell.border.alpha = 0.1,
                               n=100, # Test
                               arrow.scale=5, # Test
                               grid.n=40 # Test
)
cat("<H5>n=300, arrow.scale=10, grid.n=40</H5>")
show.velocity.on.embedding.cor(emb,
                               rvel.cd,
                               scale='sqrt',
                               cell.colors=ac(cell.colors,alpha=0.5),
                               cex=0.8,
                               show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5,
                               fixed.arrow.length = TRUE,
                               arrow.lwd=1,
                               do.par=F,
                               cell.border.alpha = 0.1,
                               n=300, # Test
                               arrow.scale=10, # Test
                               grid.n=40 # Test
)

cat("<H5>n=300, arrow.scale=5, grid.n=20</H5>")
show.velocity.on.embedding.cor(emb,
                               rvel.cd,
                               scale='sqrt',
                               cell.colors=ac(cell.colors,alpha=0.5),
                               cex=0.8,
                               show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5,
                               fixed.arrow.length = TRUE,
                               arrow.lwd=1,
                               do.par=F,
                               cell.border.alpha = 0.1,
                               n=300, # Test
                               arrow.scale=5, # Test
                               grid.n=20 # Test
)

# Visualize a fit for maker genes (we reuse rvel.cd to save on calcualtions here):
cat("<HR>")
rgb.palette <- colorRampPalette(c("grey", "red"), space = "rgb")(1000)
excluded_genes = vector()
for( gene_name in c("Emb", "F11r", "Plet1", "Tspan10", "Treml4", "Ccr2", "Clec4d", "Trem3", "Il12b", 
                    "Ccr7", "Ccl22", "Cenpa", "Ccnb2", "Mcm6", "Ccne1", "Ccne2",
                    "Birc5", "S100a6", "Pdia6", "Ranbp1", "Actb", "Gm2a", "Itm2b", "Mif", "Cd63", "Fth1",
                    "Clec4a1", "Clec4a2", "Clec4a4", "Clec4e", "Clec4n", "Clec5a")){

  cat("<HR><H5>GENE", gene_name, "</H5>")
  if( gene_name %in% rownames(filtered_emat) && gene_name %in% rownames(filtered_nmat)){
    par( mfrow = c(2,2))
    gene.relative.velocity.estimates(filtered_emat, filtered_nmat, deltaT=1, kCells = 20, kGenes=1,
                                 fit.quantile=fit.quantile, cell.emb=emb, cell.colors=cell.colors,
                                 cell.dist=cell.dist, old.fit=rvel.cd, do.par=FALSE,
                                 show.gene= gene_name)
  }else{    
    cat("<HR><BR>The gene", gene_name, "is not present in velocity matrices<HR>")
    excluded_genes = append( excluded_genes, gene_name)
  }
  
  if( gene_name %in% colnames( r_filtered$counts)){
    par( mfrow = c(1,1))
     cell.colors.2 = rgb.palette[ 1+floor( 999 * r_filtered$counts[ , gene_name] / max( r_filtered$counts[ , gene_name]))]
     names( cell.colors.2) = names( r_filtered$counts[ ,gene_name])
     show.velocity.on.embedding.cor(emb,
                                    rvel.cd,
                                    n=300,
                                    scale='sqrt',
                                    cell.colors=ac(cell.colors.2,alpha=0.5),
                                    cex=0.8,
                                    arrow.scale=5,
                                    show.grid.flow=TRUE,
                                    min.grid.cell.mass=0.5,
                                    grid.n=40,
                                    arrow.lwd=1,
                                    do.par=F,
                                    cell.border.alpha = 0.1)

  }
}
cat("<HR>Note: the following genes were not displayed because they were filtered out from velocity analysis:")
cat("<BR>", paste( excluded_genes, collapse="<BR>"))
cat("<HR>")
par( mfrow = c(1,1))
