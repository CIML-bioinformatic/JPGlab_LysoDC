# ###############################################################
# This script aims to load R object previously saved 
# in senary analysis
# ###############################################################

## @knitr reload_data

# Import the pagoda object with clusters and t_SNE embedding
r_filtered = readRDS( file.path( INPUT_DIR, "r_filtered.rds"))
filtered_emat = readRDS( file.path( INPUT_DIR, "filtered_emat.rds"))
filtered_nmat = readRDS( file.path( INPUT_DIR, "filtered_nmat.rds"))
rvel.cd = readRDS( file.path( INPUT_DIR, "rvel.cd.rds"))

# Take cluster labels from Pagoda2 pre-processing
cluster.label <- r_filtered$clusters$PCA[[1]]
cell.colors <- pagoda2:::fac2col(cluster.label)

# Take t-SNE embedding from Pagoda2 pre-processing
emb <- r_filtered$embeddings$PCA$tSNE

# Take the count table from Pagoda object
counts_table = r_filtered$counts

# Show the number of cells in each cluster
datatable( data.frame( table( cluster.label)), rownames=FALSE, colnames= c( "Cluster ID", "Nb cells"), 
           caption="Number of cells in clusters")

