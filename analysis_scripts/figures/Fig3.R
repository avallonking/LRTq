# Fig 3A
# tissue sharing
{
  library(cowplot)
  library(ggplot2)
  library(ggpubr)
  library(ggcorrplot)
  library(reshape2)
  library(Cairo)
  library(dplyr)
  reorder_cormat <- function(cormat) {
    # Use correlation between variables as distance
    na.idx <- apply(cormat, 1, function(x) all(is.na(x)))
    cormat <- cormat[!na.idx, !na.idx]
    cormat[is.na(cormat)] <- 0.20
    
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat <-cormat[hc$order, hc$order]
    return(cormat)
  }
  # Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  plot.sharing.mat <- function(cormat) {
    # cormat <- reorder_cormat(cormat)
    upper_tri = get_lower_tri(cormat)
    melted_cormat <- melt(as.matrix(upper_tri), na.rm = TRUE)
    melted_cormat <- melt(as.matrix(cormat), na.rm = T)
    ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
      geom_tile(color = "white")+ 
      scale_fill_distiller(palette = "YlOrBr", direction = 1) +
      # scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0.20,
      #                      limit = c(0,1), space = "Lab",
      #                      name="Sharing\nfraction") +
      theme_minimal() + # minimal theme
      #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
       #     panel.background = element_blank(), axis.line = element_line(colour = "black"))
      theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                       size = 5, hjust = 1),
            axis.text.y = element_text(vjust = 1, 
                                       size = 5, hjust = 1),
            axis.title.x=element_blank(), axis.title.y=element_blank())+
      coord_fixed()
    return(ggheatmap)
  }
  
  # RV eGenes
  CairoPDF("../materials/tissue.sharing.corr.matrix.pdf")
  corr.matrix <- readRDS("./data/tissue_sharing/tissue_sharing.fraction_shared_egenes.lrt.0.10.rds")
  # plot.sharing.mat(reorder_cormat(corr.matrix))
  plot.sharing.mat(corr.matrix)
  dev.off()
  
  # clustered
  CairoPDF("../materials/tissue.sharing.corr.matrix.clustered.pdf")
  corr.matrix <- readRDS("./data/tissue_sharing/tissue_sharing.fraction_shared_egenes.lrt.0.10.rds")
  plot.sharing.mat(reorder_cormat(corr.matrix))
  # plot.sharing.mat(corr.matrix)
  dev.off()
  
}

# Fig 3B
# tissue sharing - combine RV eGenes and CV eGenes
{
  library(ggpubr)
  library(data.table)
  library(Cairo)
 
  group.limit = 4
  rv = readRDS("./data/tissue_sharing/tissue_sharing.prop_shared_egenes.by.num_tissues.lrt.q.thres.0.05.num.egene.thres.200.rv.egenes.rds")
  cv = readRDS("./data/tissue_sharing/tissue_sharing.prop_shared_egenes.by.num_tissues.lrt.q.thres.0.05.num.egene.thres.200.cv.egenes.rds")
  df3 = rbind(rv[1, ], cv[1, ])
  df3$num.tissues = as.character(df3$num.tissues)
  df3$type = c(rep("RV eGenes", 1), rep("CV eGenes", 1))
  df3 = rbind(df3, c(paste("2 -", group.limit, sep=" "), 
                     sum(rv$prop.shared[2:group.limit]), "RV eGenes"))
  df3 = rbind(df3, c(paste("2 -", group.limit, sep=" "), 
                     sum(cv$prop.shared[2:group.limit]), "CV eGenes"))
  df3 = rbind(df3, c(paste(">", group.limit, sep=" "), 
                     1-sum(rv$prop.shared[1:group.limit]), "RV eGenes"))
  df3 = rbind(df3, c(paste(">", group.limit, sep=" "), 
                     1-sum(cv$prop.shared[1:group.limit]), "CV eGenes"))
  df3$prop.shared = as.numeric(df3$prop.shared)
  CairoPDF("../materials/prop.shared.cv.egenes.rv.egenes.num.tissues.grouped.200.pdf")
  p = ggbarplot(df3, x="num.tissues", y="prop.shared", xlab = "Number of tissues", 
                ylab = "Proportion of shared eGenes", fill = "type", 
                position = position_dodge())
  ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, legend.title = "",
        font.x = 13, font.y = 13, title = "Tissues with more than 200 RV eGenes",
        font.title = 13)
  dev.off()
}
