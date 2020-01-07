# ##################################################
# The aim of this script is to vinf the most
# variable genes among the cell population
# ##################################################


## FIND THE MOST VARIABLE GENES USING THE SEURAT PACKAGE METHOD

## @knitr find_variable_genes_seurat_method

cat("\n<H4>Using Seurat package method</H4>")

sc10x <- FindVariableGenes(object = sc10x, mean.function = ExpMean, dispersion.function = LogVMR,
                           x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, display.progress = FALSE, do.plot=TRUE)
SEURAT_VARIABLE_GENES = sc10x@var.genes

cat("<br>Number of variable genes (per average bin) found:", length( SEURAT_VARIABLE_GENES))

# Determine the order of top most variable genes by normalizing ((x-mean)/sd) the gene mean and
# dispersion and maximazing their sum. Then intersect this ordered list with the most variables
# genes reported by Seurat
mean_gene.dispersion.scaled = mean( sc10x@hvg.info$gene.dispersion.scaled)
sd_gene.dispersion.scaled = sd( sc10x@hvg.info$gene.dispersion.scaled)
mean_gene.mean = mean( sc10x@hvg.info$gene.mean)
sd_gene.mean = sd( sc10x@hvg.info$gene.mean)
top_var_gene_order = order( (sc10x@hvg.info$gene.dispersion.scaled-mean_gene.dispersion.scaled)/sd_gene.dispersion.scaled + (sc10x@hvg.info$gene.mean - mean_gene.mean)/sd_gene.mean, decreasing = TRUE)
top_var_gene_order_selected = intersect( top_var_gene_order, which( names( sc10x@hvg.info$gene.dispersion.scaled) %in% SEURAT_VARIABLE_GENES))[1:VAR_GENES_MAX_NUMBER]
top_var_gene_df = data_frame( 
  name = names( sc10x@hvg.info$gene.dispersion.scaled)[ top_var_gene_order_selected],
  "gene mean" = sc10x@hvg.info$gene.mean[ top_var_gene_order_selected],
  "gene dispersion per bin"= sc10x@hvg.info$gene.dispersion.scaled[ top_var_gene_order_selected]
)

datatable( top_var_gene_df, caption = paste( "Top", VAR_GENES_MAX_NUMBER, "most variable genes (per average bin)"))

selected_variable_genes = SEURAT_VARIABLE_GENES

gc()

# ## FIND THE MOST VARIABLE GENES USING THE OLD MARIONI LAB METHOD
# http://pklab.med.harvard.edu/scw2014/subpop_tutorial.html
# ############################################################

## @knitr find_variable_genes_marioni_method

library(statmod)

cat("\n<H4>Using Marioni lab method</H4>")

# I used a normalized matrix where each UMIcount is devided by the total nbr of UMI per cell (not the seurat normalization which use log)
# this matrix should already be filtered
raw_allData_df = as.data.frame(t(as.matrix(sc10x@raw.data)))
all_genes_v = rownames(sc10x@data)
norm_allData_df = as.data.frame(raw_allData_df[,all_genes_v])/rowSums(raw_allData_df[,all_genes_v]) 

PercDropout = apply( as.matrix(norm_allData_df[,all_genes_v]),2,function(x) sum(x==0)/length(x) )
MeanAll = apply( as.matrix(norm_allData_df[,all_genes_v]),2,function(x) log( mean( x)) )
CV2All = apply( as.matrix(norm_allData_df[,all_genes_v]),2,function(x) log( (sd( x) / mean(x))^2 ))
MeanPostives = apply( as.matrix(norm_allData_df[,all_genes_v]),2,function(x) log( mean( x[ x!=0]) ))
CV2Postives = apply( as.matrix(norm_allData_df[,all_genes_v]),2,function(x) log( (sd( x[ x!=0]) / mean(x[ x!=0]))^2 ))
stats_df = data.frame( PercDropout=PercDropout, MeanAll=MeanAll, CV2All=CV2All, MeanPostives=MeanPostives, CV2Postives=CV2Postives )

# ed is a cell-gene matrix (cells as columns and gene as rows)
ed=t(norm_allData_df[,all_genes_v])

means <- rowMeans(ed)
vars <- apply(ed,1,var)
cv2 <- vars/means^2
par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9)
cat("<b>Plot of log(CV^2) versus log(mean) used to estimate the variables genes</b>")
smoothScatter(log(means),log(cv2), xlab = "Log of Mean on all cells", ylab = "log of CV^2 on all cells")

# Now fit a regression line based on the controls:
minMeanForFit <- unname( quantile( means[ which( cv2 >0.3) ], .95 ) )
useForFit <- means >= minMeanForFit # & spikeins
fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/means[useForFit] ),cv2[useForFit] )
a0 <- unname( fit$coefficients["a0"] )
a1 <- unname( fit$coefficients["a1tilde"])
# Now add the fit and the 95% confidence interval to our plot:
# repeat previous plot
par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9); smoothScatter(log(means),log(cv2));
xg <- exp(seq( min(log(means[means>0])), max(log(means)), length.out=1000 ))
vfit <- a1/xg + a0
# add fit line
lines( log(xg), log(vfit), col="black", lwd=3 )
df <- sum(useForFit) - 1
# add confidence interval
lines(log(xg),log(vfit * qchisq(0.975,df)/df),lty=2,col="black")
lines(log(xg),log(vfit * qchisq(0.025,df)/df),lty=2,col="black")

# Rank genes by the significance of deviation from the fit
afit <- a1/means+a0
varFitRatio <- vars/(afit*means^2)
varorder <- order(varFitRatio,decreasing=T)
oed <- ed[varorder,]

# We can also evaluate statistical significance of the deviation
pval <- pchisq(varFitRatio*df,df=df,lower.tail=F)
adj.pval <- p.adjust(pval,"fdr")
sigVariedGenes <- adj.pval<1e-3;
datatable( data.frame( table(sigVariedGenes)), colnames = c( "Variable gene", "Number"))

# repeat previous plot with selected genes shown
cat("<b>Plot of log(CV^2) versus log(mean) with in red selected genes (adj. p-val < 1e-3)</b>")
par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9)
smoothScatter(log(means),log(cv2))
lines( log(xg), log(vfit), col="black", lwd=3 )
lines(log(xg),log(vfit * qchisq(0.975,df)/df),lty=2,col="black")
lines(log(xg),log(vfit * qchisq(0.025,df)/df),lty=2,col="black");
# # add top 100 genes
# points(log(means[varorder[1:100]]),log(cv2[varorder[1:100]]),col=3)
points(log(means),log(cv2),col=1)
points(log(means[adj.pval<1e-3]),log(cv2[adj.pval<1e-3]),col=2)

# print the table in HTML widget
MARIONI_VARIABLE_GENES = names(sigVariedGenes)[sigVariedGenes]
datatable( stats_df[ MARIONI_VARIABLE_GENES,], caption = paste( "Selected most variable genes"))

gc()


