# ############################################################
# All the functions which will be used in a scRNA-Seq Analysis
# ############################################################

## @knitr load_functions

library( ade4)
library( Rtsne)
library( data.table)
library( ggrepel)


TWO_COLORS_VECTOR=c("#FBE9E7", "#FF3D00")
TWO_COLORS_GRADIANT=c("grey80",colorRampPalette(TWO_COLORS_VECTOR)(n = 400))
color_NA_geom_point="grey50"

col_heatmap = c("grey",colorRampPalette(c("white",  "orange"))(n = 299))

# #####################################################
# #####################################################
# #####################################################
# beautiful colors (same as in ggplot2)
# #####################################################
# #####################################################
# #####################################################
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 200)[1:n]
}

# #####################################################
# #####################################################
# #####################################################
# pca function which generate only the histogram
# with the percentage of information per PC
# #####################################################
# #####################################################
# #####################################################
pca_hist_function <- function( raw_data_for_pca, PlotHist = FALSE, VectorToColorPCA = NULL, NameOfVectorToColor = "", 
                               dim_all=FALSE, dim_x="PC1", dim_y="PC2", bool_facet_wrap=FALSE, VectorFacetWrap=NA, title="",bool_scalePCA=FALSE, bool_centerPCA=TRUE){
  limits_allFigs=c(0,27)
  
  if(length(VectorFacetWrap)<=1){
    VectorFacetWrap = VectorToColorPCA
  }
  
  # special case if we want to color the pca by gene expression (continuous variable)
  if(NameOfVectorToColor %in% c(colnames(raw_data_for_pca),"Kmt2d_ex_4_5") ) {
    VectorToColorPCA[which(VectorToColorPCA==0)]=NA
  }
  
  # raw_data_for_pca must have only the Gene Expression and no additional columns
  
  # Execute a Centered pca with ade4pca_function
  # ############################################
  result__pca <- dudi.pca( raw_data_for_pca, scannf = FALSE, nf = ncol(raw_data_for_pca), scale=bool_scalePCA, center=bool_centerPCA)
  
  # Generate the histogram ggplot with the percentage of information for each principal composants
  # ##############################################################################################
  percentageInfo_byPC_df = data.frame(PC=seq(1,length(result__pca$eig),1), PercentInfo=100*result__pca$eig/sum( result__pca$eig))
  
  if(PlotHist==TRUE){
    percentage_information_histogram <- ggplot(percentageInfo_byPC_df, aes(PC, PercentInfo)) +
      geom_bar(stat="identity", color="black", fill="gray90") +
      xlab("Principal Components") +
      ylab("Percentage of variants") +
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits=c(0,max(percentageInfo_byPC_df$PercentInfo)+2))
    print(percentage_information_histogram)
  }
  
  # #####################################
  # Information for the Sweave file
  # #####################################
  if(PlotHist==TRUE){
    cat("\n****************************************************")
    cat("\n Centered PCA:")
    cat("\n****************************************************")  
  }
  
  # Look at the most informative genes along the axes
  threshold_norm = 0.75
  threshold_percent = 0.05 # 0.05
  names_best_on_axis_all_axes = vector()
  for( cs_index in c(1,2,3)){ # c(1,2,3)
    # Compute on the axis, the best variables, i.e. the variable which have the max coordinates on that axis
    current_vector = result__pca$c1[[cs_index]]
    number_best = floor( length( current_vector) * threshold_percent)
    decreasing_order = order( abs(current_vector), decreasing = TRUE)
    names_best_on_axis = row.names( result__pca$c1)[decreasing_order][1: number_best]
    names_best_on_axis_all_axes = unique( append( names_best_on_axis_all_axes, names_best_on_axis))
    
    if(PlotHist==TRUE){
      cat("\n\nOn CS", cs_index,",", threshold_percent*100, "% of best variables are\n")
      print( names_best_on_axis)
      
      # Look at the variables with the maximal cumulative information on each axes
      cat("\nOn CS", cs_index, ", direct and cummulative information carried by best variables:\n")
      perc_info_vector=c() 
      genes_names=c() 
      perc_info_sum = 0
      for( gene_index in seq(1, length( current_vector))){
        perc_info = result__pca$c1[[ cs_index]][decreasing_order][gene_index]^2
        perc_info_sum = perc_info_sum + perc_info
        cat("\n", row.names( result__pca$c1)[decreasing_order][gene_index], "\t:\t", signif((perc_info)*100, 3), "%\t:\t", signif( perc_info_sum*100 , 3), " %")
        if( perc_info_sum >= threshold_norm){
          break 
        }
      }
    }
  }
  
  
  # ##########
  # plot PCA
  # ##########
  if(!is.null(VectorToColorPCA)){
    raw_data_for_pca$DataToColorPCA = VectorToColorPCA
    # take the Principal Composants Coordinates
    # ########################################
    raw_data_for_pca$PC1=result__pca$li$Axis1
    raw_data_for_pca$PC2=result__pca$li$Axis2
    raw_data_for_pca$PC3=result__pca$li$Axis3
    raw_data_for_pca$PC4=result__pca$li$Axis4
    raw_data_for_pca$PC5=result__pca$li$Axis5
    raw_data_for_pca$PC6=result__pca$li$Axis6
    pca_percentageInfoByAxis = as.integer(100*result__pca$eig/sum( result__pca$eig)+0.5)
    names(pca_percentageInfoByAxis)=paste0("PC",1:length(pca_percentageInfoByAxis))
    
    if(dim_all==FALSE){
      if(title!=""){
        title=paste(title,"\n")
      }
      PCX=dim_x
      PCY=dim_y
      # we change the scale of the plot's axis to have a perfect plot 
      distPCX= max(raw_data_for_pca[,PCX])-min(raw_data_for_pca[,PCX])
      distPCY= max(raw_data_for_pca[,PCY])-min(raw_data_for_pca[,PCY])
      addIntervalPCX=0
      addIntervalPCY=0
      if(distPCX>distPCY){
        addIntervalPCY=(distPCX-distPCY)/2
      } else {
        addIntervalPCX=(distPCY-distPCX)/2
      }
      if(bool_facet_wrap == FALSE){
        if(NameOfVectorToColor %in% c(colnames(raw_data_for_pca),"Kmt2d_ex_4_5") ) {
          plot_pca <- ggplot(raw_data_for_pca, aes_string(PCX, PCY, color="DataToColorPCA")) +
            geom_point() +
            coord_fixed() +
            scale_colour_gradient(low =TWO_COLORS_VECTOR[1], high = TWO_COLORS_VECTOR[2], na.value = color_NA_geom_point, name="Et") +
            xlim(min(raw_data_for_pca[,PCX])-addIntervalPCX, max(raw_data_for_pca[,PCX])+addIntervalPCX) +
            ylim(min(raw_data_for_pca[,PCY])-addIntervalPCY, max(raw_data_for_pca[,PCY])+addIntervalPCY) +
            xlab(paste0(PCX," (",pca_percentageInfoByAxis[PCX],"%)")) +
            ylab(paste0(PCY," (",pca_percentageInfoByAxis[PCY],"%)")) +
            ggtitle(paste(title,"PCA colored by",NameOfVectorToColor))
          print(plot_pca)
        }else{
          plot_pca <- ggplot(raw_data_for_pca, aes_string(PCX, PCY, color="DataToColorPCA")) +
            geom_point() +
            coord_fixed() +
            xlim(min(raw_data_for_pca[,PCX])-addIntervalPCX, max(raw_data_for_pca[,PCX])+addIntervalPCX) +
            ylim(min(raw_data_for_pca[,PCY])-addIntervalPCY, max(raw_data_for_pca[,PCY])+addIntervalPCY) +
            xlab(paste0(PCX," (",pca_percentageInfoByAxis[PCX],"%)")) +
            ylab(paste0(PCY," (",pca_percentageInfoByAxis[PCY],"%)")) +
            ggtitle(title,paste("PCA colored by",NameOfVectorToColor))
          print(plot_pca)
        }
      }else{
        # we split the plot into multiple plots
        raw_data_for_pca$DataToFacetWrap = VectorFacetWrap
        raw_data_for_pca$DataToFacetWrap2 =raw_data_for_pca$DataToFacetWrap
        if(NameOfVectorToColor %in% c(colnames(raw_data_for_pca),"Kmt2d_ex_4_5") ) {
          plot_pca <- ggplot(raw_data_for_pca, aes_string(PCX, PCY, color="DataToColorPCA")) +
            geom_point(data=raw_data_for_pca[,-match("DataToFacetWrap" , colnames(raw_data_for_pca))],aes_string(x=PCX,y=PCY,group="DataToFacetWrap2"),colour="grey")+
            geom_point() +
            coord_fixed() +
            scale_colour_gradient(low =TWO_COLORS_VECTOR[1], high = TWO_COLORS_VECTOR[2], na.value = color_NA_geom_point, name="Et") +
            xlim(min(raw_data_for_pca[,PCX])-addIntervalPCX, max(raw_data_for_pca[,PCX])+addIntervalPCX) +
            ylim(min(raw_data_for_pca[,PCY])-addIntervalPCY, max(raw_data_for_pca[,PCY])+addIntervalPCY) +
            xlab(paste0(PCX," (",pca_percentageInfoByAxis[PCX],"%)")) +
            ylab(paste0(PCY," (",pca_percentageInfoByAxis[PCY],"%)")) +
            ggtitle(paste(title,"PCA splited by",NameOfVectorToColor)) + 
            facet_wrap(~DataToFacetWrap)
          print(plot_pca)
        }else{
          plot_pca <- ggplot(raw_data_for_pca, aes_string(PCX, PCY, color="DataToColorPCA")) +
            geom_point(data=raw_data_for_pca[,-match("DataToFacetWrap" , colnames(raw_data_for_pca))],aes_string(x=PCX,y=PCY,group="DataToFacetWrap2"),colour="grey")+
            geom_point() +
            coord_fixed() +
            xlim(min(raw_data_for_pca[,PCX])-addIntervalPCX, max(raw_data_for_pca[,PCX])+addIntervalPCX) +
            ylim(min(raw_data_for_pca[,PCY])-addIntervalPCY, max(raw_data_for_pca[,PCY])+addIntervalPCY) +
            xlab(paste0(PCX," (",pca_percentageInfoByAxis[PCX],"%)")) +
            ylab(paste0(PCY," (",pca_percentageInfoByAxis[PCY],"%)")) +
            ggtitle(paste(title,"PCA splited by",NameOfVectorToColor)) + 
            facet_wrap(~DataToFacetWrap)
          print(plot_pca)
        }
      }
      
    } else {
      ggplot_list = list()
      PC_v = c("PC1","PC2","PC3")
      index=0
      for(i in 1:(length(PC_v)-1)){
        for(j in (i+1):length(PC_v)){
          index=index+1
          gg_title=""
          if(index==1){
            gg_title=title
          }
          PCX=PC_v[i]
          PCY=PC_v[j]
          # we change the scale of the plot's axis to have a perfect plot 
          distPCX= max(raw_data_for_pca[,PCX])-min(raw_data_for_pca[,PCX])
          distPCY= max(raw_data_for_pca[,PCY])-min(raw_data_for_pca[,PCY])
          addIntervalPCX=0
          addIntervalPCY=0
          if(distPCX>distPCY){
            addIntervalPCY=(distPCX-distPCY)/2
          } else {
            addIntervalPCX=(distPCY-distPCX)/2
          }
          if(NameOfVectorToColor %in% c(colnames(raw_data_for_pca),"Kmt2d_ex_4_5") ) {
            plot_pca <- ggplot(raw_data_for_pca, aes_string(PCX, PCY, color="DataToColorPCA")) +
              geom_point(size=0.6) +
              coord_fixed() +
              scale_colour_gradient(low =TWO_COLORS_VECTOR[1], high = TWO_COLORS_VECTOR[2], na.value = color_NA_geom_point, name="Et") +
              xlim(min(raw_data_for_pca[,PCX])-addIntervalPCX, max(raw_data_for_pca[,PCX])+addIntervalPCX) +
              ylim(min(raw_data_for_pca[,PCY])-addIntervalPCY, max(raw_data_for_pca[,PCY])+addIntervalPCY) +
              xlab(paste0(PCX," (",pca_percentageInfoByAxis[PCX],"%)")) +
              ylab(paste0(PCY," (",pca_percentageInfoByAxis[PCY],"%)")) +
              ggtitle(gg_title) +
              guides(colour=guide_legend(title=NameOfVectorToColor))
            print(plot_pca)
            ggplot_list[[index]]=plot_pca
          }else{
            plot_pca <- ggplot(raw_data_for_pca, aes_string(PCX, PCY, color="DataToColorPCA")) +
              geom_point(size=0.6) +
              coord_fixed() +
              xlim(min(raw_data_for_pca[,PCX])-addIntervalPCX, max(raw_data_for_pca[,PCX])+addIntervalPCX) +
              ylim(min(raw_data_for_pca[,PCY])-addIntervalPCY, max(raw_data_for_pca[,PCY])+addIntervalPCY) +
              xlab(paste0(PCX," (",pca_percentageInfoByAxis[PCX],"%)")) +
              ylab(paste0(PCY," (",pca_percentageInfoByAxis[PCY],"%)")) +
              ggtitle(gg_title) +
              guides(colour=guide_legend(title=NameOfVectorToColor))
            print(plot_pca)
            ggplot_list[[index]]=plot_pca
          }
        }
      }
      
    }
    
  }
  
  return(list(pca_info = result__pca, genes_information = names_best_on_axis_all_axes))
}

# #####################################################
# #####################################################
# #####################################################
# PCA genes information distribution
# #####################################################
# #####################################################
# #####################################################
top_genes_fromPCA_function <-function(pca_info, dim_min=1, dim_max=3, threshold_norm = 0.75){
  pca=pca_info
  variable_genes_withThreshold_v = c()
  #threshold_norm = 0.75
  # take the top gene which drive the PC1 and the PC2 from the pca report
  # #########################################################
  for(i in dim_min:dim_max){
    genes_vector = c()
    
    current_vector = pca$c1[[i]]
    decreasing_order = order( abs(current_vector), decreasing = TRUE)
    for( gene_index in seq(1, length( current_vector))){
      current_norm = sqrt (sum( (pca$c1[[ i]][decreasing_order][1: gene_index])^2))
      genes_vector=c(genes_vector,row.names( pca$c1)[decreasing_order][gene_index])
      if( current_norm >= threshold_norm){
        break
      }
    }
    variable_genes_withThreshold_v = c(variable_genes_withThreshold_v, unique(genes_vector) )
    
  }
  variable_genes_withThreshold_v = unique(variable_genes_withThreshold_v)
  return(variable_genes_withThreshold_v)
}

# #####################################################
# #####################################################
# #####################################################
# PCA genes information distribution
# #####################################################
# #####################################################
# #####################################################
pca_gene_information_function <-function(pca_info, dim_all=FALSE, dim_x=1, dim_y=2, threshold_norm = 0.75, title=""){
  pca=pca_info
  ggplot_list = list()
  index=0
  threshold_norm_top_genes = 0.6
  
  variable_genes_withThreshold_v = c()
  #threshold_norm = 0.75
  # take the top gene which drive the PC1 and the PC2 from the pca report
  # #########################################################
  for(PC_i in 1:2){
    for(PC_j in (PC_i+1):3){
      index=index+1
      gg_title=""
      if(index==1){
        gg_title=title
      }
      genes_vector = c()
      for(i in c(PC_i,PC_j)){
        current_vector = pca$c1[[i]]
        decreasing_order = order( abs(current_vector), decreasing = TRUE)
        for( gene_index in seq(1, length( current_vector))){
          current_norm = sqrt (sum( (pca$c1[[ i]][decreasing_order][1: gene_index])^2))
          genes_vector=c(genes_vector,row.names( pca$c1)[decreasing_order][gene_index])
          if( current_norm >= threshold_norm){
            break
          }
        }
      }
      genes_to_study = unique(genes_vector)
      
      pca_percentageInfoByAxis = as.integer(100*pca$eig/sum( pca$eig)+0.5)
      
      # Create the article figure AR
      # #############################
      pca_to_use_df = pca$co[genes_to_study,]
      pca_to_use_df$Labels = rownames(pca_to_use_df)
      pca_genes_plot <- ggplot(pca_to_use_df, aes_string(paste0("Comp",PC_i),paste0("Comp",PC_j), label="Labels")) +
        geom_segment(aes_string(x=0, y=0, xend=paste0("Comp",PC_i), yend=paste0("Comp",PC_j)), arrow = arrow(length = unit(0.01, "npc")), color="grey" ) +
        #geom_text(size=3.5)
        #theme_bw(base_size = 1) + #theme(legend.position = "bottom") +
        #theme( plot.title = element_text( face="bold", size=1)) +
        geom_text_repel(segment.color = NA, size=2, fontface="italic")+#,
        #box.padding = unit(0.35, "lines"),
        #point.padding = unit(0.25, "lines")) +
        xlab(paste0("PC",PC_i," (",pca_percentageInfoByAxis[PC_i],"%)")) +
        ylab(paste0("PC",PC_j," (",pca_percentageInfoByAxis[PC_j],"%)")) +
        coord_fixed()+
        ggtitle(gg_title)
      print(pca_genes_plot)
      ggplot_list[[index]]=pca_genes_plot
      
      variable_genes_withThreshold_v = c(variable_genes_withThreshold_v, rownames(pca_to_use_df))
    }
  }
  
  if(FALSE){
    # Look at the genes with the maximal cumulative information on each axes
    threshold_norm = threshold_norm_top_genes
    genes_maxCumulativeInfo_vector = c()
    for(cs_index in c(1,2)){
      cat("\n\n Genes information on Principal Composant ",cs_index," :")
      cat("\n--------------------------------------------------")
      current_vector=pca$c1[[cs_index]]
      decreasing_order = order( abs(current_vector), decreasing = TRUE)
      perc_info_vector=c() 
      genes_names=c() 
      perc_info_sum = 0
      for( gene_index in seq(1, length( current_vector))){
        perc_info = pca$c1[[ cs_index]][decreasing_order][gene_index]^2
        perc_info_sum = perc_info_sum + perc_info
        cat("\n", row.names( pca$c1)[decreasing_order][gene_index], "\t:\t", signif((perc_info)*100, 3), "%\t:\t", signif( perc_info_sum*100 , 3), " %")
        if(! row.names( pca$c1)[decreasing_order][gene_index] %in% genes_maxCumulativeInfo_vector){
          genes_maxCumulativeInfo_vector = c(genes_maxCumulativeInfo_vector, row.names( pca$c1)[decreasing_order][gene_index])
        }
        if( perc_info_sum >= threshold_norm){
          break
        }
      }
    }
  }
  # do.call("grid.arrange", c(ggplot_list, nrow=3))
  # for( index in names( ggplot_list)){
  #   print( ggplot_list[[ index]])
  # }
  
  variable_genes_withThreshold_v = unique(variable_genes_withThreshold_v)
  #return(ggplot_list)
  return(variable_genes_withThreshold_v)
}

# #####################################################
# #####################################################
# #####################################################
# Function that execute a t-SNE (bh-SNE)
# #####################################################
# #####################################################
# #####################################################
tsne_function <- function( raw_data_for_tsne, VectorToColor = NULL, NameOfVectorToColor = "" ,perplexity_value=30, bool_facet_wrap=FALSE, VectorFacetWrap=NA,title=""){
  limits_allFigs=c(0,27)
  
  if(length(VectorFacetWrap)<=1){
    VectorFacetWrap = VectorToColor
  }
  
  # special case if we want to color the pca by gene expression (continuous variable)
  if(NameOfVectorToColor %in% colnames(raw_data_for_tsne) ) {
    VectorToColor[which(VectorToColor==0)]=NA
  }
  
  #raw_data_for_tsne=raw_data_all_df[,GENES_TO_ANALYSE]
  if(min(rowSums(raw_data_for_tsne))==0){
    VectorToColor = VectorToColor[-which(rowSums(raw_data_for_tsne)==0)]
    raw_data_for_tsne = raw_data_for_tsne[-which(rowSums(raw_data_for_tsne)==0),]
  }
  
  
  # bh-SNE: t-SNE using Barnes-Hut implementation
  # #############################################
  set.seed(42) # Sets seed for reproducibility
  tsne_out <- Rtsne(as.matrix(raw_data_for_tsne),initial_dims = 20, perplexity=perplexity_value) #initial_dims = 50
  
  # take the tSNE Coordinates
  # #########################
  raw_data_for_tsne$tsne_x = tsne_out$Y[,1]
  raw_data_for_tsne$tsne_y = tsne_out$Y[,2]
  raw_data_for_tsne$DataToColor = VectorToColor
  
  if(bool_facet_wrap==FALSE){
    if(NameOfVectorToColor %in% colnames(raw_data_for_tsne) ) {
      tSNE_plot <- ggplot(raw_data_for_tsne, aes(tsne_x, tsne_y, color=DataToColor)) +
        geom_point() +
        coord_fixed() +
        xlab("t-SNE 1") +
        ylab("t-SNE 2") +
        theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank()) +
        scale_colour_gradient(low =TWO_COLORS_VECTOR[1], high = TWO_COLORS_VECTOR[2], na.value = color_NA_geom_point, name="Et") +
        ggtitle(paste(title,"\n t-SNE colored by",NameOfVectorToColor))
      print(tSNE_plot)
        
    } else{
      tSNE_plot <- ggplot(raw_data_for_tsne, aes(tsne_x, tsne_y, color=DataToColor)) +
        geom_point() +
        coord_fixed() +
        xlab("t-SNE 1") +
        ylab("t-SNE 2") +
        theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank()) +
        ggtitle(paste(title,"\n t-SNE colored by",NameOfVectorToColor))
      print(tSNE_plot)
    }
  }
  else{
    raw_data_for_tsne$DataToFacetWrap = VectorFacetWrap
    raw_data_for_tsne$DataToFacetWrap2 =raw_data_for_tsne$DataToFacetWrap
    if(NameOfVectorToColor %in% colnames(raw_data_for_tsne) ) {
      tSNE_plot <- ggplot(raw_data_for_tsne, aes(tsne_x, tsne_y, color=DataToColor)) +
        geom_point(data=raw_data_for_tsne[,-match("DataToFacetWrap" , colnames(raw_data_for_tsne))],aes(x=tsne_x,y=tsne_y,group=DataToFacetWrap2),colour="grey")+
        geom_point() +
        coord_fixed() +
        xlab("t-SNE 1") +
        ylab("t-SNE 2") +
        theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank()) +
        scale_colour_gradient(low =TWO_COLORS_VECTOR[1], high = TWO_COLORS_VECTOR[2], na.value = color_NA_geom_point, name="Et") +
        facet_wrap(~DataToFacetWrap) +
        ggtitle(paste(title,"\n t-SNE colored by",NameOfVectorToColor))
      print(tSNE_plot)
    } else{
      tSNE_plot <- ggplot(raw_data_for_tsne, aes(tsne_x, tsne_y, color=DataToColor)) +
        geom_point(data=raw_data_for_tsne[,-match("DataToFacetWrap" , colnames(raw_data_for_tsne))],aes(x=tsne_x,y=tsne_y,group=DataToFacetWrap2),colour="grey")+
        geom_point() +
        coord_fixed() +
        xlab("t-SNE 1") +
        ylab("t-SNE 2") +
        theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank()) +
        facet_wrap(~DataToFacetWrap) +
        ggtitle(paste(title,"\n t-SNE colored by",NameOfVectorToColor))
      print(tSNE_plot)
    }
    
  }
}

# #####################################################
# #####################################################
# #####################################################
# Function that run t-SNE (bh-SNE)
# #####################################################
# #####################################################
# #####################################################
run_tsne_function <- function( raw_data_for_tsne ,perplexity_value=30){
  
  #raw_data_for_tsne=raw_data_all_df[,GENES_TO_ANALYSE]
  if(min(rowSums(raw_data_for_tsne))==0){
    raw_data_for_tsne = raw_data_for_tsne[-which(rowSums(raw_data_for_tsne)==0),]
  }
  
  
  # bh-SNE: t-SNE using Barnes-Hut implementation
  # #############################################
  set.seed(42) # Sets seed for reproducibility
  tsne_out <- Rtsne(as.matrix(raw_data_for_tsne),initial_dims = 20, perplexity=perplexity_value) #initial_dims = 50
  
  # take the tSNE Coordinates
  # #########################
  raw_data_for_tsne$tsne_x = tsne_out$Y[,1]
  raw_data_for_tsne$tsne_y = tsne_out$Y[,2]
  return(raw_data_for_tsne[,c("tsne_x","tsne_y")])
}
# #####################################################
# #####################################################
# #####################################################
# Function that execute a t-SNE (bh-SNE)
# #####################################################
# #####################################################
# #####################################################
plot_tsne_function <- function( raw_data_for_tsne, tsne_coord_df,VectorToColor = NULL, NameOfVectorToColor = "" , VectorFacetWrap=NA,title=""){
  limits_allFigs=c(0,27)
  
  #if(length(VectorFacetWrap)<=1){
  #  VectorFacetWrap = VectorToColor
  #}
  
  # special case if we want to color the pca by gene expression (continuous variable)
  if(NameOfVectorToColor %in% colnames(raw_data_for_tsne) ) {
    VectorToColor[which(VectorToColor==0)]=NA
  }
  
  #raw_data_for_tsne=raw_data_all_df[,GENES_TO_ANALYSE]
  if(min(rowSums(raw_data_for_tsne))==0){
    VectorToColor = VectorToColor[-which(rowSums(raw_data_for_tsne)==0)]
    raw_data_for_tsne = raw_data_for_tsne[-which(rowSums(raw_data_for_tsne)==0),]
  }
  
  
  # take the tSNE Coordinates
  # #########################
  raw_data_for_tsne$tsne_x = tsne_coord_df$tsne_x
  raw_data_for_tsne$tsne_y = tsne_coord_df$tsne_y
  raw_data_for_tsne$DataToColor = VectorToColor
  
  if(length(VectorFacetWrap)<=1){
    if(NameOfVectorToColor %in% colnames(raw_data_for_tsne) ) {
      tSNE_plot <- ggplot(raw_data_for_tsne, aes(tsne_x, tsne_y, color=DataToColor)) +
        geom_point() +
        coord_fixed() +
        xlab("t-SNE 1") +
        ylab("t-SNE 2") +
        scale_colour_gradient(low =TWO_COLORS_VECTOR[1], high = TWO_COLORS_VECTOR[2], na.value = color_NA_geom_point, name="Et") +
        ggtitle(paste(title,"\n t-SNE colored by",NameOfVectorToColor))
      return(tSNE_plot)
      
    } else{
      tSNE_plot <- ggplot(raw_data_for_tsne, aes(tsne_x, tsne_y, color=DataToColor)) +
        geom_point() +
        coord_fixed() +
        xlab("t-SNE 1") +
        ylab("t-SNE 2") +
        ggtitle(paste(title,"\n t-SNE colored by",NameOfVectorToColor))
      return(tSNE_plot)
    }
  }
  else{
    raw_data_for_tsne$DataToFacetWrap = VectorFacetWrap
    raw_data_for_tsne$DataToFacetWrap2 =raw_data_for_tsne$DataToFacetWrap
    if(NameOfVectorToColor %in% colnames(raw_data_for_tsne) ) {
      tSNE_plot <- ggplot(raw_data_for_tsne, aes(tsne_x, tsne_y, color=DataToColor)) +
        geom_point(data=raw_data_for_tsne[,-match("DataToFacetWrap" , colnames(raw_data_for_tsne))],aes(x=tsne_x,y=tsne_y,group=DataToFacetWrap2),colour="grey")+
        geom_point() +
        coord_fixed() +
        xlab("t-SNE 1") +
        ylab("t-SNE 2") +
        scale_colour_gradient(low =TWO_COLORS_VECTOR[1], high = TWO_COLORS_VECTOR[2], na.value = color_NA_geom_point, name="Et") +
        facet_wrap(~DataToFacetWrap) +
        ggtitle(paste(title,"\n t-SNE colored by",NameOfVectorToColor))
      return(tSNE_plot)
    } else{
      tSNE_plot <- ggplot(raw_data_for_tsne, aes(tsne_x, tsne_y, color=DataToColor)) +
        geom_point(data=raw_data_for_tsne[,-match("DataToFacetWrap" , colnames(raw_data_for_tsne))],aes(x=tsne_x,y=tsne_y,group=DataToFacetWrap2),colour="grey")+
        geom_point() +
        coord_fixed() +
        xlab("t-SNE 1") +
        ylab("t-SNE 2") +
        facet_wrap(~DataToFacetWrap) +
        ggtitle(paste(title,"\n t-SNE colored by",NameOfVectorToColor))
      return(tSNE_plot)
    }
    
  }
}


# #####################################################
# #####################################################
# #####################################################
# volcano plot - combine LRT expression genes (modified)
# #####################################################
# #####################################################
# #####################################################
volcano_plot_LRT_function <- function(raw_data_for_volcano_plot, VectorToSplit, FigureLetterName=NULL){
  
  variable_to_return=list()
  
  cell_attribute_toSplit = "ParameterToSplit"
  genes_vector = colnames(raw_data_for_volcano_plot)
  
  # Build the dataframe with 1 gene per line to build the FluidigmAssay object from data
  # -- Build a basic data frame with the correct number of lines
  one_gene_per_line_df = data.frame( count=seq( 1,ncol(raw_data_for_volcano_plot)*nrow( raw_data_for_volcano_plot)), stringsAsFactors=FALSE)
  
  # -- Add the vector_to_split and UniqueCellID
  one_gene_per_line_df = cbind( one_gene_per_line_df, rep( VectorToSplit, each= ncol(raw_data_for_volcano_plot)))
  one_gene_per_line_df = cbind( one_gene_per_line_df, rep( rownames(raw_data_for_volcano_plot), each= ncol(raw_data_for_volcano_plot)))
  
  # -- Add the duplication of the genes names
  one_gene_per_line_df$Gene = rep( genes_vector, nrow( raw_data_for_volcano_plot))
  
  # -- Add the value of expressions
  gene_et = vector()
  for( row_index in seq( 1, nrow( raw_data_for_volcano_plot))){
    gene_et = append( gene_et, unlist( raw_data_for_volcano_plot[ row_index, ]))
  }
  one_gene_per_line_df$Et = (gene_et) #exp(gene_et) for normalized scRNAseq data
  
  # -- Set the name of the dataframe
  names( one_gene_per_line_df) = c("count",cell_attribute_toSplit,"UniqueCellID","Gene","Et")
  
  # Build the FluidigmAssay object from MAST library
  # using the header of dataframe
  vbeta.fa = FluidigmAssay(one_gene_per_line_df, idvars="UniqueCellID", 
                           primerid="Gene",
                           measurement="Et",
                           #ncells= "CellNumber",
                           geneid="Gene", 
                           cellvars=c("UniqueCellID"),
                           phenovars=c(cell_attribute_toSplit),
                           id= raw_data_filename)
  
  
  cell_attribute_toSplit_vector = sort(unique(VectorToSplit))
  
  for( index1 in 1:(length( cell_attribute_toSplit_vector)-1)){
    current_ref_state =  cell_attribute_toSplit_vector[ index1]
    
    lrt_test <- LRT(vbeta.fa, cell_attribute_toSplit, referent = current_ref_state)
    
    for( index2 in (index1+1):length( cell_attribute_toSplit_vector)){
      current_alt_state =  cell_attribute_toSplit_vector[ index2]
      # Get the separated two groups of cells
      vbeta.1cell.nozero.nobad.ref = subset(vbeta.fa, eval(parse(text=cell_attribute_toSplit)) == current_ref_state)
      vbeta.1cell.nozero.nobad.alt = subset(vbeta.fa, eval(parse(text=cell_attribute_toSplit)) == current_alt_state)
      
      
      # Compute the difference of mean expression on the two groups and build the mean and pvalue set for the volcano plot
      mean_diff_set = vector()

	  lrt_pval_set = vector()
      for( gene_name in genes_vector){
        mean_ref = mean( exprs( vbeta.1cell.nozero.nobad.ref)[ ,gene_name])
        mean_alt = mean( exprs( vbeta.1cell.nozero.nobad.alt)[ ,gene_name])
        pval = lrt_test[ which( lrt_test[[cell_attribute_toSplit]] == current_alt_state & lrt_test$primerid == gene_name),"p.value"]
        mean_diff_set = append( mean_diff_set, mean_ref - mean_alt)
        lrt_pval_set = append( lrt_pval_set, pval)
      }
      lrt_pval_set_adjustMultiComparison = p.adjust(lrt_pval_set, method="BH")
      lrt_pval_set_adjustMultiComparison = -log10(lrt_pval_set_adjustMultiComparison)
      lrt_pval_set = -log10(lrt_pval_set)
      mean_thresholds = quantile(mean_diff_set, c(.05, .95))
      mean_thresholds_2 = c( mean( mean_diff_set) - sd( mean_diff_set), mean( mean_diff_set) + sd( mean_diff_set))
      
      
      # Compute the list of genes that are the most extreme in volcano plot (upper left and upper right)    
      mean_index_list = which( mean_diff_set <= mean_thresholds[1] | mean_diff_set >= mean_thresholds[2])
      mean_index_list_2 = which( mean_diff_set <= mean_thresholds_2[1] | mean_diff_set >= mean_thresholds_2[2])
      pval_index_list = which( lrt_pval_set_adjustMultiComparison >= -log10(0.05))
      final_index_list = intersect( mean_index_list, pval_index_list)
      final_index_list_2 = intersect( mean_index_list_2, pval_index_list)
      
      
      # -- THE BRAND NEW & TOTALY AWESOME VOLCANO PLOT --
      # #################################################
      volcano_plot_data_frame=data.frame(row.names=genes_vector, MeanDiffSet=mean_diff_set, LRTpVal=lrt_pval_set_adjustMultiComparison) #lrt_pval_set
      
      volcano_plot_data_frame$Significant <- ifelse(rownames(volcano_plot_data_frame) %in% rownames(volcano_plot_data_frame)[final_index_list], "< 0.05 Perc or > 0.95 Perc", "Not Sig")
      volcano_plot_data_frame$Significant2 <- ifelse(rownames(volcano_plot_data_frame) %in% rownames(volcano_plot_data_frame)[final_index_list_2], "abs( Z-score) >= 1", "Not Sig")
      
      if(TRUE){ #!is.null(FigureLetterName)
        volcano_plot_ggplot <- ggplot(volcano_plot_data_frame, aes(x = MeanDiffSet, y = LRTpVal)) +
          xlab(paste( "Diff of exprs means between", current_ref_state, "and", current_alt_state)) +
          ylab("-log10 p.value")+ 
          ggtitle(paste(current_alt_state,"(overexprs left) vs",current_ref_state,"(overexprs right)"))+
          geom_hline(yintercept=-log10( 0.05), colour="grey50", alpha=0.4) +
          geom_vline(xintercept = mean_thresholds_2[1], colour="firebrick2", alpha=0.4) +
          geom_vline(xintercept = mean_thresholds_2[2], colour="firebrick2", alpha=0.4) +
          scale_size_manual(values = c(3,1.5), guide="none") +
          geom_point(aes(color = Significant2),size=3) +
          scale_color_manual(values = c("firebrick2", "gray70"), name="Significant genes") +
          theme( plot.title = element_text( face="bold", size=12),
                 aspect.ratio=1,
                 legend.position="none") +
          geom_text_repel(
            data = subset(volcano_plot_data_frame[which(volcano_plot_data_frame$MeanDiffSet>0),], Significant2=="abs( Z-score) >= 1"),
            aes(label = rownames(subset(volcano_plot_data_frame[which(volcano_plot_data_frame$MeanDiffSet>0),], Significant2=="abs( Z-score) >= 1"))),
            size = 3, fontface="italic", nudge_x=3,nudge_y=1,segment.color = '#cccccc',
            box.padding = unit(0.25, "lines"),
            point.padding = unit(0.25, "lines")) +
          geom_text_repel(
            data = subset(volcano_plot_data_frame[which(volcano_plot_data_frame$MeanDiffSet<0),], Significant2=="abs( Z-score) >= 1"),
            aes(label = rownames(subset(volcano_plot_data_frame[which(volcano_plot_data_frame$MeanDiffSet<0),], Significant2=="abs( Z-score) >= 1"))),
            size = 3, fontface="italic", nudge_x=-2, nudge_y=1,segment.color = '#cccccc',
            box.padding = unit(0.25, "lines"),
            point.padding = unit(0.25, "lines")) +
        scale_x_continuous(limits=c(-max(abs(volcano_plot_data_frame$MeanDiffSet)),max(abs(volcano_plot_data_frame$MeanDiffSet)) )) +
          scale_y_continuous(expand=c(0,0),limits=c(0,max(abs(volcano_plot_data_frame$LRTpVal))+0.5))
        print(volcano_plot_ggplot)
      }
      
      if(FALSE){
        cat("\nMost significant genes in both mean difference and LRT p-value (", current_ref_state, "vs", current_alt_state, "):")
        cat( "\n\n lesser than 5 percentile or greater than 95 percentile:\n")
        print( genes_vector[final_index_list], quote=FALSE, row.names=FALSE)
        cat( "\n\n abs( Z-score) >= 1 :\n")
        print( genes_vector[final_index_list_2], quote=FALSE, row.names=FALSE)
      }
      
      variable_to_return[[paste(current_ref_state, "vs", current_alt_state)]]=list()
      # put the genes which are more expressed in current_ref_state in one vector
      variable_to_return[[paste(current_ref_state, "vs", current_alt_state)]][[current_ref_state]]=genes_vector[final_index_list_2[which(mean_diff_set[final_index_list_2]>0)]]
      # put the genes which are more expressed in current_alt_state in one vector  
      variable_to_return[[paste(current_ref_state, "vs", current_alt_state)]][[current_alt_state]]=genes_vector[final_index_list_2[which(mean_diff_set[final_index_list_2]<0)]]
      
      
    }
  }
  
  for(i in 1:length(variable_to_return)){
    cat("\n\nMost significant genes in both mean difference and LRT p-value (", names(variable_to_return)[i], "):\n")
    for(j in 1:length(variable_to_return[[i]])){
      cat("Genes overexpressed", names(variable_to_return[[i]])[j], ":\n")
      print(variable_to_return[[i]][[j]])
    }
  }
  
  all_categories_v = sort(unique(VectorToSplit))
  if(length(all_categories_v)>2){
    cat("\n\n*****************************************************************************")
    cat("\n*****************************************************************************")
    cat("\n*****************************************************************************")
    cat("\n\nList of the genes which are overexpressed or underexpressed by CATEGORY vs ALL OTHER CATEGORIES:\n\n")
    variable_to_return
    diff_expr_list=list()
    
    for(i in 1:length(all_categories_v)){
      up_exprs_genes_v = c()
      down_exprs_genes_v = c()
      for(j in 1:length(variable_to_return)){
        if(all_categories_v[i]==names(variable_to_return[[j]])[1] || all_categories_v[i]==names(variable_to_return[[j]])[2]){
          if(length(up_exprs_genes_v)==0){
            up_exprs_genes_v = variable_to_return[[j]][[all_categories_v[i] ]]
          } else{
            up_exprs_genes_v = intersect(up_exprs_genes_v,variable_to_return[[j]][[all_categories_v[i] ]])
          }
          if(length(down_exprs_genes_v)==0){
            down_exprs_genes_v = variable_to_return[[j]][[which(names(variable_to_return[[j]])!=all_categories_v[i]) ]]
          } else{
            down_exprs_genes_v = intersect(down_exprs_genes_v,variable_to_return[[j]][[which(names(variable_to_return[[j]])!=all_categories_v[i]) ]])
          }
        }
      }
      
      diff_expr_list[[all_categories_v[i] ]] = list(UP=up_exprs_genes_v,DOWN=down_exprs_genes_v)
      
    }
    #print(diff_expr_list)
    for(i in 1:length(diff_expr_list)){
      cat("\n\n",names(diff_expr_list)[i],":")
      cat("\n**************************")
      for(j in 1:length(diff_expr_list[[i]])){
        cat("\n",names(diff_expr_list[[i]])[j],":\n")
        print(diff_expr_list[[i]][[j]])
      }
    }
  }
  
  return(variable_to_return)
}

# ######################################################
# ######################################################
# ######################################################
# volcano plot - combine LRT expression genes (modified)
# ######################################################
# ######################################################
# ######################################################
volcano_plot_LRT_function_v2 <- function(raw_data_for_volcano_plot, VectorToSplit, FigureLetterName=NULL){
  
  variable_to_return=list()
  
  cell_attribute_toSplit = "ParameterToSplit"
  genes_vector = colnames(raw_data_for_volcano_plot)
  
  # Build the dataframe with 1 gene per line to build the FluidigmAssay object from data
  # -- Build a basic data frame with the correct number of lines
  one_gene_per_line_df = data.frame( count=seq( 1,ncol(raw_data_for_volcano_plot)*nrow( raw_data_for_volcano_plot)), stringsAsFactors=FALSE)
  
  # -- Add the vector_to_split and UniqueCellID
  one_gene_per_line_df = cbind( one_gene_per_line_df, rep( VectorToSplit, each= ncol(raw_data_for_volcano_plot)))
  one_gene_per_line_df = cbind( one_gene_per_line_df, rep( rownames(raw_data_for_volcano_plot), each= ncol(raw_data_for_volcano_plot)))
  
  # -- Add the duplication of the genes names
  one_gene_per_line_df$Gene = rep( genes_vector, nrow( raw_data_for_volcano_plot))
  
  # -- Add the value of expressions
  gene_et = vector()
  for( row_index in seq( 1, nrow( raw_data_for_volcano_plot))){
    gene_et = append( gene_et, unlist( raw_data_for_volcano_plot[ row_index, ]))
  }
  one_gene_per_line_df$Et = gene_et
  
  # -- Set the name of the dataframe
  names( one_gene_per_line_df) = c("count",cell_attribute_toSplit,"UniqueCellID","Gene","Et")
  
  # Build the FluidigmAssay object from MAST library
  # using the header of dataframe
  vbeta.fa = FluidigmAssay(one_gene_per_line_df, idvars="UniqueCellID", 
                           primerid="Gene",
                           measurement="Et",
                           #ncells= "CellNumber",
                           geneid="Gene", 
                           cellvars=c("UniqueCellID"),
                           phenovars=c(cell_attribute_toSplit),
                           id= raw_data_filename)
  
  
  cell_attribute_toSplit_vector = sort(unique(VectorToSplit))
  
  for( index1 in 1:(length( cell_attribute_toSplit_vector)-1)){
    current_ref_state =  cell_attribute_toSplit_vector[ index1]
    
    lrt_test <- LRT(vbeta.fa, cell_attribute_toSplit, referent = current_ref_state)
    
    for( index2 in (index1+1):length( cell_attribute_toSplit_vector)){
      current_alt_state =  cell_attribute_toSplit_vector[ index2]
      # Get the separated two groups of cells
      vbeta.1cell.nozero.nobad.ref = subset(vbeta.fa, eval(parse(text=cell_attribute_toSplit)) == current_ref_state)
      vbeta.1cell.nozero.nobad.alt = subset(vbeta.fa, eval(parse(text=cell_attribute_toSplit)) == current_alt_state)
      
      
      # Compute the difference of mean expression on the two groups and build the mean and pvalue set for the volcano plot
      mean_diff_set_noZero = vector()
      mean_diff_set = vector()
      perc_diff_set = vector()

	  lrt_pval_set = vector()
      for( gene_name in genes_vector){
        mean_ref_noZero = mean( exprs( vbeta.1cell.nozero.nobad.ref)[ ,gene_name][exprs( vbeta.1cell.nozero.nobad.ref)[ ,gene_name]!=0])
        mean_alt_noZero = mean( exprs( vbeta.1cell.nozero.nobad.alt)[ ,gene_name][exprs( vbeta.1cell.nozero.nobad.alt)[ ,gene_name]!=0])
        if(is.nan(mean_ref)){mean_ref=0}
        if(is.nan(mean_alt)){mean_alt=0}
        perc_ref = sum(exprs( vbeta.1cell.nozero.nobad.ref)[ ,gene_name]==0)/length(exprs( vbeta.1cell.nozero.nobad.ref)[ ,gene_name])
        perc_alt = sum(exprs( vbeta.1cell.nozero.nobad.alt)[ ,gene_name]==0)/length(exprs( vbeta.1cell.nozero.nobad.alt)[ ,gene_name])
        mean_ref = mean( exprs( vbeta.1cell.nozero.nobad.ref)[ ,gene_name])
        mean_alt = mean( exprs( vbeta.1cell.nozero.nobad.alt)[ ,gene_name])
        pval = lrt_test[ which( lrt_test[[cell_attribute_toSplit]] == current_alt_state & lrt_test$primerid == gene_name),"p.value"]
        mean_diff_set_noZero = append( mean_diff_set_noZero, mean_ref_noZero - mean_alt_noZero)
        mean_diff_set = append( mean_diff_set, mean_ref - mean_alt)
        perc_diff_set = append( perc_diff_set, perc_ref - perc_alt)
        lrt_pval_set = append( lrt_pval_set, -log10(pval))
      }
      mean_thresholds_noZero = quantile(mean_diff_set_noZero, c(.05, .95))
      mean_thresholds = quantile(mean_diff_set, c(.05, .95))
      mean_thresholds_2 = c( mean( mean_diff_set) - sd( mean_diff_set), mean( mean_diff_set) + sd( mean_diff_set))
      
      
      # Compute the list of genes that are the most extreme in volcano plot (upper left and upper right)    
      mean_index_list = which( mean_diff_set <= mean_thresholds[1] | mean_diff_set >= mean_thresholds[2])
      mean_index_list_2 = which( mean_diff_set <= mean_thresholds_2[1] | mean_diff_set >= mean_thresholds_2[2])
      pval_index_list = which( lrt_pval_set >= -log10(0.05))
      final_index_list = intersect( mean_index_list, pval_index_list)
      final_index_list_2 = intersect( mean_index_list_2, pval_index_list)
      
      final_index_list_basedOnlyOnPval = pval_index_list
      
      # -- THE BRAND NEW & TOTALY AWESOME VOLCANO PLOT --
      # #################################################
      volcano_plot_data_frame=data.frame(row.names=genes_vector, MeanDiffSet=mean_diff_set_noZero, PercDiffSet=perc_diff_set, LRTpVal=lrt_pval_set)
      
      volcano_plot_data_frame$Significant <- ifelse(rownames(volcano_plot_data_frame) %in% rownames(volcano_plot_data_frame)[final_index_list], "< 0.05 Perc or > 0.95 Perc", "Not Sig")
      volcano_plot_data_frame$Significant2 <- ifelse(rownames(volcano_plot_data_frame) %in% rownames(volcano_plot_data_frame)[final_index_list_2], "abs( Z-score) >= 1", "Not Sig")
      
      volcano_plot_data_frame$SignificantOnlyPval <- ifelse(rownames(volcano_plot_data_frame) %in% rownames(volcano_plot_data_frame)[final_index_list_basedOnlyOnPval], "pval>0.05","Not Sig")
      
      if(TRUE){ 
        volcano_plot_ggplot <- ggplot(volcano_plot_data_frame, aes(x = MeanDiffSet, y = PercDiffSet)) +
          xlab(paste( "Diff of expr means (of values>0)")) +
          ylab("Difference of percentage of zero values")+ 
          ggtitle(paste(current_alt_state,"(overexpressed left) vs",current_ref_state,"(overexpressed right)"))+
          geom_hline(yintercept=0, colour="grey50", alpha=0.4) +
          geom_vline(xintercept=0, colour="grey50", alpha=0.4) +
          scale_size_manual(values = c(3,1.5), guide="none") +
          geom_point(aes(color = Significant2),size=3) +
          scale_color_manual(values = c("firebrick2", "gray70"), name="Significant genes") +
          theme( plot.title = element_text( face="bold", size=12),
                 aspect.ratio=1,
                 legend.position="none") +
          geom_text_repel(
            data = subset(volcano_plot_data_frame[which(volcano_plot_data_frame$MeanDiffSet>0),], Significant2=="abs( Z-score) >= 1"),
            aes(label = rownames(subset(volcano_plot_data_frame[which(volcano_plot_data_frame$MeanDiffSet>0),], Significant2=="abs( Z-score) >= 1"))),
            size = 3, fontface="italic", nudge_x=3,nudge_y=1,segment.color = '#cccccc',
            box.padding = unit(0.25, "lines"),
            point.padding = unit(0.25, "lines")) +
          geom_text_repel(
            data = subset(volcano_plot_data_frame[which(volcano_plot_data_frame$MeanDiffSet<0),], Significant2=="abs( Z-score) >= 1"),
            aes(label = rownames(subset(volcano_plot_data_frame[which(volcano_plot_data_frame$MeanDiffSet<0),], Significant2=="abs( Z-score) >= 1"))),
            size = 3, fontface="italic", nudge_x=-2, nudge_y=1,segment.color = '#cccccc',
            box.padding = unit(0.25, "lines"),
            point.padding = unit(0.25, "lines"))
        print(volcano_plot_ggplot)
        
        volcano_plot_ggplot <- ggplot(volcano_plot_data_frame, aes(x = MeanDiffSet, y = PercDiffSet)) +
          xlab(paste( "Difference of expression means (of values>0)")) +
          ylab("Difference of percentage of zero values")+ 
          ggtitle(paste(current_alt_state,"(overexpressed left) vs",current_ref_state,"(overexpressed right)"))+
          geom_hline(yintercept=0, colour="grey50", alpha=0.4) +
          geom_vline(xintercept=0, colour="grey50", alpha=0.4) +
          scale_size_manual(values = c(3,1.5), guide="none") +
          geom_point(aes(color = SignificantOnlyPval),size=3) +
          scale_color_manual(values = c("gray70", "firebrick2"), name="Significant genes") +
          theme( plot.title = element_text( face="bold", size=12),
                 aspect.ratio=1,
                 legend.position="none") +
          geom_text_repel(
            data = subset(volcano_plot_data_frame[which(volcano_plot_data_frame$MeanDiffSet>0),], SignificantOnlyPval=="pval>0.05"),
            aes(label = rownames(subset(volcano_plot_data_frame[which(volcano_plot_data_frame$MeanDiffSet>0),], SignificantOnlyPval=="pval>0.05"))),
            size = 3, fontface="italic", nudge_x=3,nudge_y=1,segment.color = '#cccccc',
            box.padding = unit(0.25, "lines"),
            point.padding = unit(0.25, "lines")) +
          geom_text_repel(
            data = subset(volcano_plot_data_frame[which(volcano_plot_data_frame$MeanDiffSet<0),], SignificantOnlyPval=="pval>0.05"),
            aes(label = rownames(subset(volcano_plot_data_frame[which(volcano_plot_data_frame$MeanDiffSet<0),], SignificantOnlyPval=="pval>0.05"))),
            size = 3, fontface="italic", nudge_x=-2, nudge_y=1,segment.color = '#cccccc',
            box.padding = unit(0.25, "lines"),
            point.padding = unit(0.25, "lines"))
        print(volcano_plot_ggplot)
        
        
        volcano_plot_ggplot <- ggplot(volcano_plot_data_frame, aes(x = MeanDiffSet, y = PercDiffSet, color=LRTpVal)) +
          xlab(paste( "Difference of expression means (of values>0)")) +
          ylab("Difference of percentage of zero values")+ 
          ggtitle(paste(current_alt_state,"(overexpressed left) vs",current_ref_state,"(overexpressed right)"))+
          geom_hline(yintercept=0, colour="grey50", alpha=0.4) +
          geom_vline(xintercept=0, colour="grey50", alpha=0.4) +
          scale_size_manual(values = c(3,1.5), guide="none") +
          geom_point(size=3) +
          scale_colour_gradientn(colours = rainbow(7))
        print(volcano_plot_ggplot)
        
        
        volcano_plot_ggplot <- ggplot(volcano_plot_data_frame, aes(x = PercDiffSet, y = LRTpVal, color=LRTpVal)) +
          xlab(paste( "Difference of expression means (of values>0)")) +
          ylab("Difference of percentage of zero values")+ 
          ggtitle(paste(current_alt_state,"(overexpressed left) vs",current_ref_state,"(overexpressed right)"))+
          geom_hline(yintercept=0, colour="grey50", alpha=0.4) +
          geom_vline(xintercept=0, colour="grey50", alpha=0.4) +
          scale_size_manual(values = c(3,1.5), guide="none") +
          geom_point(size=3) +
          scale_colour_gradientn(colours = rainbow(7))
        print(volcano_plot_ggplot)
      }
      
      if(FALSE){
        cat("\nMost significant genes in both mean difference and LRT p-value (", current_ref_state, "vs", current_alt_state, "):")
        cat( "\n\n lesser than 5 percentile or greater than 95 percentile:\n")
        print( genes_vector[final_index_list], quote=FALSE, row.names=FALSE)
        cat( "\n\n abs( Z-score) >= 1 :\n")
        print( genes_vector[final_index_list_2], quote=FALSE, row.names=FALSE)
      }
      
      variable_to_return[[paste(current_ref_state, "vs", current_alt_state)]]=list()
      # put the genes which are more expressed in current_ref_state in one vector
      variable_to_return[[paste(current_ref_state, "vs", current_alt_state)]][[current_ref_state]]=genes_vector[final_index_list_2[which(mean_diff_set[final_index_list_2]>0)]]
      # put the genes which are more expressed in current_alt_state in one vector  
      variable_to_return[[paste(current_ref_state, "vs", current_alt_state)]][[current_alt_state]]=genes_vector[final_index_list_2[which(mean_diff_set[final_index_list_2]<0)]]
      
      
    }
  }
  
  for(i in 1:length(variable_to_return)){
    cat("\n\nMost significant genes in both mean difference and LRT p-value (", names(variable_to_return)[i], "):\n")
    for(j in 1:length(variable_to_return[[i]])){
      cat("Genes overexpressed", names(variable_to_return[[i]])[j], ":\n")
      print(variable_to_return[[i]][[j]])
    }
  }
  
  return(variable_to_return)
}

# #####################################################
# #####################################################
# #####################################################
# Biclustering backSPIN
# #####################################################
# #####################################################
# #####################################################
biclustering_function <- function (raw_data_for_biclustering, side_bar_colored_variable="Plate", genes_for_clustering=ALL_FILTERED_GENES, fig_title=""){
  
  colors_to_use=gg_color_hue(length(unique(raw_data_for_biclustering[,side_bar_colored_variable])))
  
  
  title_ggplot=paste("Biclustering:",fig_title)
  
  # -- data for the study
  raw_data_transposed=t(raw_data_for_biclustering[, genes_for_clustering])
  # - colnames based on Pheno
  colnames(raw_data_transposed)=raw_data_for_biclustering$UniqueCellID[1:ncol(raw_data_transposed)]
  
  
  # we remove the genes with no variation, otherwise we can not do the hierarchical clustering on genes
  genes_no_variation=c()
  for (i in 1:nrow(raw_data_transposed)){
    sum_row=sum(raw_data_transposed[i,])
    if(sum_row==0){
      genes_no_variation=c(genes_no_variation,rownames(raw_data_transposed)[i])
    }
  }
  raw_data_transposed=raw_data_transposed[which(!rownames(raw_data_transposed) %in% genes_no_variation),]
  
  print(paste("We remove the genes with no variation: ",paste(genes_no_variation,collapse=" "),sep=""))
  
  cells_no_variation=c()
  index_cells_no_variation=c()
  for (i in 1:ncol(raw_data_transposed)){
    sum_col=sum(raw_data_transposed[,i])
    if(sum_col==0){
      cells_no_variation=c(cells_no_variation,colnames(raw_data_transposed)[i])
      index_cells_no_variation=c(index_cells_no_variation,i)
    }
  }
  raw_data_transposed=raw_data_transposed[,which(!colnames(raw_data_transposed) %in% cells_no_variation)]
  
  print(paste("We remove the cells with no variation: ",paste(cells_no_variation,collapse=" "),sep=""))
  
  # #####################
  # backSPIN biclustering
  # https://github.com/linnarsson-lab/BackSPIN
  # http://hms-dbmi.github.io/scw/analysis-of-heterogeneity-and-subpopulations.html
  # #####################
  
  # set the backSpin algorithm dir
  backSPIN_algorithm_python="~/workspace/Logiciels/BackSPIN-1.0/backSPIN.py"
  
  data_dir="~/workspace/"
  # Delete CEF files if they already exist
  file_backSPIN=paste(data_dir,"datafile_backSPIN.cef",sep="")
  if (file.exists(file_backSPIN)){
    file.remove(file_backSPIN)
  }
  file_backSPIN_clustered=paste(data_dir,"datafile_backSPIN.clustered.cef",sep="")
  if (file.exists(file_backSPIN_clustered)){
    file.remove(file_backSPIN_clustered)
  }
  
  # creation of the CEF file (the file used by backSPIN)
  # ####################################################
  write(file=file_backSPIN, paste("CEF\t0\t1\t1", nrow(raw_data_transposed), ncol(raw_data_transposed), "0", sep="\t"))
  write(file=file_backSPIN, paste("\tcell", paste(colnames(raw_data_transposed), collapse="\t"), sep='\t'), append=T)
  x <- cbind(rep("",nrow(raw_data_transposed)), raw_data_transposed)
  write(file=file_backSPIN, paste(c("gene"), sep="\t"), append=T)
  write.table(file=file_backSPIN, x, append=T, col.names=F, row.names=T, quote=F, sep="\t")
  
  # call the backSPIN algorithm in python which create a new CEF file where cells and genes are biclustered
  # #######################################################################################################
  system(paste(backSPIN_algorithm_python," -i ", file_backSPIN, " -o ", file_backSPIN_clustered," -v -d 3",sep=""))
  
  # retrieve the data (to create the heatmap) generated by the backSPIN algorithm 
  # #############################################################################
  backSPIN_clustered_data_frame <- read.delim(file_backSPIN_clustered,stringsAsFactors=F,sep="\t",skip=1,header=F)
  index_first_group_gene=match(0,backSPIN_clustered_data_frame[-c(1,2,3,4,5,6),3])
  index_second_group_gene=match(1,backSPIN_clustered_data_frame[-c(1,2,3,4,5,6),3])
  index_first_group_cell=match(0,backSPIN_clustered_data_frame[3,-c(1,2,3,4,5,6)])
  index_second_group_cell=match(1,backSPIN_clustered_data_frame[3,-c(1,2,3,4,5,6)])
  
  names_row=backSPIN_clustered_data_frame[-c(1,2,3,4,5,6),1] # lines specific to our type of data (check the raw CEF file in the Terminal to understand: "< datafile_backSPIN.clustered.cef cef view")
  names_col=backSPIN_clustered_data_frame[1,-c(1,2,3,4,5,6)]
  backSPIN_clustered_data_frame <- backSPIN_clustered_data_frame[-c(1:6),-c(1:6)]
  backSPIN_clustered_data_frame=matrix(as.numeric(unlist(backSPIN_clustered_data_frame)),nrow=nrow(backSPIN_clustered_data_frame))
  rownames(backSPIN_clustered_data_frame)=names_row
  colnames(backSPIN_clustered_data_frame)=names_col
  
  
  # #######################################
  # --- PLOT USING ggplot2: + geom_raster()
  # #######################################
  
  # reorganization of the data for ggplot
  
  # Build the dataframe with 1 gene per line to build the FluidigmAssay object from data
  # -- Build a basic data frame with the correct number of lines
  one_gene_per_line_df = data.frame( count=seq( 1,nrow(backSPIN_clustered_data_frame)*ncol(backSPIN_clustered_data_frame) ), stringsAsFactors=FALSE)
  
  # -- Add the duplication of the genes names
  one_gene_per_line_df$Genes = structure( rep(rev(1:nrow(backSPIN_clustered_data_frame)),ncol(backSPIN_clustered_data_frame)), .Label=rev(rownames(backSPIN_clustered_data_frame)), class="factor" )
  
  # -- Add the duplication of the genes names
  cell_names_repetition_vector=c()
  for(i in 1:ncol(backSPIN_clustered_data_frame)){
    cell_names_repetition_vector=c(cell_names_repetition_vector,rep(colnames(backSPIN_clustered_data_frame)[i],times=nrow(backSPIN_clustered_data_frame)))
  }
  one_gene_per_line_df$UniqueCellID = structure( rep(1:ncol(backSPIN_clustered_data_frame),each=nrow(backSPIN_clustered_data_frame) ), .Label=colnames(backSPIN_clustered_data_frame), class="factor")
  
  
  # -- Add the value of expressions
  gene_et = vector()
  for( i in 1:ncol(backSPIN_clustered_data_frame)){
    gene_et = append( gene_et, unlist( backSPIN_clustered_data_frame[ , i ]))
  }
  maxEt=max(gene_et)
  for(i in 1:length(gene_et)){
    if(gene_et[i]==0){
      gene_et[i]=NA
    }
  }
  
  # split the continuous data to make discrete data
  gene_et_cut = cut(gene_et, breaks= seq(0,maxEt+1, by=(maxEt+1)/6),labels=1:6)
  
  one_gene_per_line_df$Et = gene_et_cut
  
  #data frame for the side bar colored by the phenotype
  side_bar_df = data.frame(count=1:ncol(backSPIN_clustered_data_frame))
  side_bar_df$UniqueCellID =  structure(1:ncol(backSPIN_clustered_data_frame),  .Label=colnames(backSPIN_clustered_data_frame), class="factor")
  side_bar_df$Y = rep(-1,times=ncol(backSPIN_clustered_data_frame))
  side_bar_df = merge(side_bar_df, raw_data_for_biclustering[,c("UniqueCellID",side_bar_colored_variable)],by="UniqueCellID", all=FALSE)
  
  biclust_ggplot<-ggplot() +
    geom_raster(data=one_gene_per_line_df,aes(x=UniqueCellID,y=Genes,fill=Et)) +
    geom_raster(data=side_bar_df,aes_string(x="UniqueCellID",y="Y",fill=side_bar_colored_variable)) + #facet_wrap(c(side_bar_colored_variable)) +
    scale_fill_manual(values=c(scales::seq_gradient_pal("#132B43","#56B1F7","Lab")(seq(0,1,length.out=6)), "gray18",gg_color_hue(length(unique(as.vector(t(side_bar_df[side_bar_colored_variable]))))-1) ),na.value="grey50",breaks=c("1","2","3","4","5","6",sort(unique(as.vector(t(side_bar_df[side_bar_colored_variable]))))), labels=c("0 to 4.3","4.4 to 8.7","8.8 to 13.0","13.1 to 17.3","17.4 to :21.6","21.7 to 26", sort(unique(as.vector(t(side_bar_df[side_bar_colored_variable]))))) )+
    geom_rect(size=0.5,fill=NA,colour="black",aes(xmin=index_first_group_cell-0.5,xmax=(index_second_group_cell-1)+0.5, ymin=(nrow(backSPIN_clustered_data_frame)-(index_second_group_gene-1)+1)-0.5,ymax=(nrow(backSPIN_clustered_data_frame)-index_first_group_gene+1)+0.5)) +
    geom_rect(size=0.5,fill=NA,colour="black",aes(xmin=index_second_group_cell-0.5,xmax=ncol(backSPIN_clustered_data_frame)+0.5, ymin=1-0.5,ymax=(nrow(backSPIN_clustered_data_frame)-index_second_group_gene+1)+0.5)) +
    theme(axis.text.x  = element_blank(), axis.text.y  = element_text(size=7)) + #element_text(angle=90, size=5)
    ggtitle(title_ggplot)
  
  print(biclust_ggplot)
  
  
}

# #####################################################
# #####################################################
# #####################################################
# Pseudotime
# #####################################################
# #####################################################
# #####################################################
pseudotimeDDRTree_function <- function(raw_data_for_pseudotimeDDRTree, genes_to_study_v, cell_attribute_colorBy=c("GroupPheno"), cell_attribute_facetwrapBy=NA ,bool_facet_wrap=FALSE, DDRTree_max_components=2, plot_title=""){
  #DDRTree_max_components=2
  
  
  bool_change_cellAttributeFaceWrap=FALSE
  if(!is.na(cell_attribute_facetwrapBy)){
    bool_change_cellAttributeFaceWrap=TRUE
  }
  
  # creation of a cellData for the CellDateSet
  raw_data_transposed=t(raw_data_for_pseudotimeDDRTree[, genes_to_study_v])
  colnames(raw_data_transposed)=raw_data_for_pseudotimeDDRTree$UniqueCellID
  
  # creation of a phenoData for the CellDateSet
  phenoData_data_frame=data.frame(UniqueCellID=raw_data_for_pseudotimeDDRTree[,"UniqueCellID"])
  for(gene in genes_to_study_v){
    phenoData_data_frame[[gene]]=as.vector(raw_data_transposed[gene,])
  }
  for(i in 1:length(cell_attribute_colorBy)){
    if(!cell_attribute_colorBy[i] %in% c("State","Pseudotime")){
      phenoData_data_frame[[cell_attribute_colorBy[i] ]] = raw_data_for_pseudotimeDDRTree[,cell_attribute_colorBy[i]]
    }
  }
  
  row.names(phenoData_data_frame)=colnames(raw_data_transposed)
  phenoData_annotated <- new("AnnotatedDataFrame", data = phenoData_data_frame)
  
  
  #cds=new ("CellDataSet",exprs=raw_data_transposed)
  cds=newCellDataSet(as.matrix(raw_data_transposed),phenoData=phenoData_annotated, expressionFamily=gaussianff())
  
  ggplot_DDRTree_title=paste("DDRTree Pseudotime:",plot_title, "(max CP=",DDRTree_max_components,")")
  
  set.seed(45)
  # -- DDRTree reduce dimension
  cds_reduce_dim_DDRTree=reduceDimension(cds, max_components=DDRTree_max_components, reduction_method="DDRTree",norm_method="none",pseudo_expr=0) #DDRTree
  tree_MST_niaki_DDRTree=orderCells(cds_reduce_dim_DDRTree)
  
  # -- ggplot of the DDRTree
  
  for(cell_attribute_colorBy_value in cell_attribute_colorBy){
    if(bool_change_cellAttributeFaceWrap){
      cell_attribute_facetwrapBy=cell_attribute_colorBy_value
    }
    ggplot_cell_trajectory<-plot_cell_trajectory(tree_MST_niaki_DDRTree,cell_size=1, color_by=cell_attribute_colorBy_value)#+ facet_wrap(~State)#, color_by="SortPheno") #, markers="MKI67")
    ggplot_cell_trajectory<-ggplot_cell_trajectory + ggtitle(ggplot_DDRTree_title)
    ggplot_cell_trajectory<-ggplot_cell_trajectory+ theme(axis.title.x=element_blank(), axis.title.y=element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),line = element_blank())
    
    if(bool_facet_wrap==TRUE && !(identical(cell_attribute_colorBy_value,"Pseudotime")) ){
      ggplot_cell_trajectory_multi<-plot_cell_trajectory(tree_MST_niaki_DDRTree,cell_size=1, show_branch_points = FALSE, color_by=cell_attribute_colorBy_value)+ facet_wrap(cell_attribute_facetwrapBy)#, color_by="SortPheno") #, markers="MKI67")
      ggplot_cell_trajectory_multi<-ggplot_cell_trajectory_multi + ggtitle(ggplot_DDRTree_title)
      ggplot_cell_trajectory_multi<-ggplot_cell_trajectory_multi+ theme(axis.title.x=element_blank(), axis.title.y=element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),line = element_blank())
      print(ggplot_cell_trajectory)
      print(ggplot_cell_trajectory_multi)
    }else{
      print(ggplot_cell_trajectory)
    }
  }
  
  # - recovery of the pseudotime values for each single-cell
  pseudotime_data_Monocle_DDRTree=t(pData(tree_MST_niaki_DDRTree)["Pseudotime"])
  
  # - recovery of the DDRTree state for each single-cell:
  DDRTree_state=t(pData(tree_MST_niaki_DDRTree)["State"])
  
  variables_to_return= list()
  variables_to_return[["DDRTree_pseudotime"]]=pseudotime_data_Monocle_DDRTree
  variables_to_return[["DDRTree_state"]]=DDRTree_state
  variables_to_return[["CellDataSet_DDRTree"]]=tree_MST_niaki_DDRTree
  
  return(variables_to_return)
}


# #####################################################
# #####################################################
# #####################################################
# isOutlier_MAD_function same function as "isOutlier" in scater package 
# #####################################################
# #####################################################
# #####################################################
isOutlier_MAD_function <- function (x_v,  nmads=3, type="both", log=TRUE, title=""){
  vector_v=x_v
  if(log==TRUE){
    vector_v=log10(x_v)
  }
  mad_value = mad(vector_v)
  hist( vector_v, main=title, xlab=ifelse( log == TRUE, paste( "Log", title), title))
  bool_v=c()
  if(type=="both"){
    bool_v = c( (vector_v < (median(vector_v) - nmads*mad_value)) | (vector_v > (median(vector_v) + nmads*mad_value)) )
    abline( v=median(vector_v) - nmads*mad_value, col="blue")
    abline( v=median(vector_v) + nmads*mad_value, col="red")
  } else if(type=="lower"){
    bool_v = c( vector_v < (median(vector_v) - nmads*mad_value) )
    abline( v=median(vector_v) - nmads*mad_value, col="blue")
  } else if(type=="higher"){
    bool_v = c( vector_v > (median(vector_v) + nmads*mad_value) )
    abline( v=median(vector_v) + nmads*mad_value, col="red")
  }
  return(bool_v)
}
