
library(RColorBrewer)
library(tidyverse)
library(gdata)
library(gplots)

yellowblackpurple <- colorRampPalette(c("yellow","black", "purple"))(n = 100)


### GO heatmap from HVGs ----

# Data import and formatting


GO_HVG_noNA <- read.table("GO_heatmaps/GO_HVG_allZT_summary_smallNames.txt",header=TRUE, sep='\t')

GO_HVG_noNA$numOne <- rowSums(GO_HVG_noNA[2:13]==1)
GO_HVG_noNA_filter <- dplyr::filter(GO_HVG_noNA, numOne<12)

GO_HVG_noNA_filter$pvalueFilter <- rowSums(GO_HVG_noNA_filter[2:13]<=0.005)
GO_HVG_noNA_filter2 <- dplyr::filter(GO_HVG_noNA_filter,pvalueFilter>=1)

GO_HVG_noNA_filter2 <- remove_rownames(GO_HVG_noNA_filter2)
GO_HVG_noNA_filter2 <- column_to_rownames(GO_HVG_noNA_filter2, var="GO")
GO_HVG_noNA_filter2_log <- log10(GO_HVG_noNA_filter2)


expression.pearson <- cor(t(GO_HVG_noNA_filter2_log[,1:12]), method = "pearson")
expression.pearson[is.na(expression.pearson)]<-0
expression.pearson.dist<-as.dist((1-expression.pearson))
expression.pearson.dist[is.na(expression.pearson.dist)]<-0
expression.pearson.dist.clust.complete <- hclust(expression.pearson.dist, method="complete")
expression.pearson.dist.clust.complete.den <- as.dendrogram(expression.pearson.dist.clust.complete)


# Clustering and heatmap plotting

heatmap.2(as.matrix(GO_HVG_noNA_filter2_log[,1:12]), 
          Colv = NA, dendrogram='row', Rowv=expression.pearson.dist.clust.complete.den,
          col= yellowblackpurple,breaks=seq(-3,0.5,length.out=101), trace = 'none' , 
          labRow=rownames(GO_HVG_noNA_filter2_log), labCol=names(GO_HVG_noNA_filter2_log[,1:12]), 
          main="GO, log10(p-value)",
          na.rm=F, na.color="black",margins=c(4,12))



### GO heatmap from LVGs ----

# Data import and formatting


GO_LVG_noNA <- read.table("GO_heatmaps/LVG_allZT_summary_onlySignificants.txt",header=TRUE, sep='\t')

GO_LVG_noNA$numOne <- rowSums(GO_LVG_noNA[2:13]==1)
GO_LVG_noNA_filter <- dplyr::filter(GO_LVG_noNA, numOne<12)

GO_LVG_noNA_filter$pvalueFilter <- rowSums(GO_LVG_noNA_filter[2:13]<=0.005)
GO_LVG_noNA_filter2 <- dplyr::filter(GO_LVG_noNA_filter,pvalueFilter>=1)

GO_LVG_noNA_filter2 <- remove_rownames(GO_LVG_noNA_filter2)
GO_LVG_noNA_filter2 <- column_to_rownames(GO_LVG_noNA_filter2, var="GO")
GO_LVG_noNA_filter2_log <- log10(GO_LVG_noNA_filter2)


expression.pearson <- cor(t(GO_LVG_noNA_filter2_log[,1:12]), method = "pearson")
expression.pearson[is.na(expression.pearson)]<-0
expression.pearson.dist<-as.dist((1-expression.pearson))
expression.pearson.dist[is.na(expression.pearson.dist)]<-0
expression.pearson.dist.clust.complete <- hclust(expression.pearson.dist, method="complete")
expression.pearson.dist.clust.complete.den <- as.dendrogram(expression.pearson.dist.clust.complete)


# Clustering and heatmap plotting

heatmap.2(as.matrix(GO_LVG_noNA_filter2_log[,1:12]), 
          Colv = NA, dendrogram='row', Rowv=expression.pearson.dist.clust.complete.den,
          col= yellowblackpurple,breaks=seq(-3,0.5,length.out=101), trace = 'none' , 
          labRow=rownames(GO_LVG_noNA_filter2_log), labCol=names(GO_LVG_noNA_filter2_log[,1:12]), 
          main="GO, log10(p-value)",
          na.rm=F, na.color="black",margins=c(4,12))




