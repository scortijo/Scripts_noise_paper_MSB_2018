
library(RColorBrewer)
library(tidyverse)
library(gdata)
library(gplots)

purpleblackyellow <- colorRampPalette(c("purple", "black", "yellow"))(n = 100)

ColDay=c("orange", "orange", "orange", "orange", "orange", "red", "black", "black",
         "black", "black", "black", "blue")

### Clustering and heatmap of Corrected CV2 for HVGs ----

# Data import and formatting

AllZT <- read.table("Heatmaps_CV2_expression/BrenneckeMod_allZT_RUV2.txt",
                    header=TRUE, sep='\t')

AllZT_BioVar <- select(AllZT, Gene, contains("BioVar"))

AllZT_BioVar[is.na(AllZT_BioVar)]<-0

HVG <- read.table("Heatmaps_CV2_expression/HVG_allZT_nbZT.txt",header=TRUE,sep='\t')

HVG_allZT_BioVar <- inner_join(AllZT_BioVar,HVG, by=c("Gene"="allZT")) %>%
  select(., Gene, contains("BioVar"))

# Clustering and heatmap plotting

HVG_AllZT_BioVar.pearson <- cor(t(HVG_allZT_BioVar[,c(2:13)]), method = "pearson")
HVG_AllZT_BioVar.pearson[is.na(HVG_AllZT_BioVar.pearson)]<-0
HVG_AllZT_BioVar.pearson.dist<-as.dist((1-HVG_AllZT_BioVar.pearson))
HVG_AllZT_BioVar.pearson.dist[is.na(HVG_AllZT_BioVar.pearson.dist)]<-0
HVG_AllZT_BioVar.pearson.dist.clust.complete <- hclust(HVG_AllZT_BioVar.pearson.dist, 
                                                       method="complete")
HVG_AllZT_BioVar.pearson.dist.clust.complete.den <- as.dendrogram(HVG_AllZT_BioVar.pearson.dist.clust.complete)
HVG_AllZT_BioVar.pearson.dist.clust.complete.col=brewer.pal(12, 'Set3')[cutree(HVG_AllZT_BioVar.pearson.dist.clust.complete, k = 4)]
HVG_allZT_BioVar$cluster_4_pearson <- cutree(HVG_AllZT_BioVar.pearson.dist.clust.complete, 
                                             k = 4)
HVG_allZT_BioVar$cluster_4_pearson_color <- brewer.pal(12, 'Set3')[cutree(HVG_AllZT_BioVar.pearson.dist.clust.complete, k = 4)]


heatmap.2(as.matrix(HVG_allZT_BioVar[,2:13]), Colv=NA, dendrogram='row', 
          Rowv=HVG_AllZT_BioVar.pearson.dist.clust.complete.den, 
          col= purpleblackyellow,breaks=seq(-1.5,4,length.out=101), trace = 'none' ,
          labRow=NA, labCol=names(HVG_allZT_BioVar[,2:13]), 
          main="Highly Variable Genes",ColSideColors=ColDay,
          RowSideColors = HVG_AllZT_BioVar.pearson.dist.clust.complete.col,
          na.rm=F, na.color="black")


# Heatmap for average expression for HVGs, using the hierchical clustering of Corrected CV2

# Data formatting

AllZT_Mean <- select(AllZT, Gene, contains("Mean"))

HVG_AllZT_Mean <- inner_join(AllZT_Mean,HVG, by=c("Gene"="allZT")) %>%
  select(., Gene, contains("Mean"))

HVG_AllZT_Mean_meanNorm <- mutate(HVG_AllZT_Mean,
                                  Mean_ZT2_meanNorm=Mean_ZT2/rowMeans(HVG_AllZT_Mean[,2:13],na.rm=TRUE),
                                  Mean_ZT4_meanNorm=Mean_ZT4/rowMeans(HVG_AllZT_Mean[,2:13],na.rm=TRUE),
                                  Mean_ZT6_meanNorm=Mean_ZT6/rowMeans(HVG_AllZT_Mean[,2:13],na.rm=TRUE),
                                  Mean_ZT8_meanNorm=Mean_ZT8/rowMeans(HVG_AllZT_Mean[,2:13],na.rm=TRUE),
                                  Mean_ZT10_meanNorm=Mean_ZT10/rowMeans(HVG_AllZT_Mean[,2:13],na.rm=TRUE),
                                  Mean_ZT12_meanNorm=Mean_ZT12/rowMeans(HVG_AllZT_Mean[,2:13],na.rm=TRUE),
                                  Mean_ZT14_meanNorm=Mean_ZT14/rowMeans(HVG_AllZT_Mean[,2:13],na.rm=TRUE),
                                  Mean_ZT16_meanNorm=Mean_ZT16/rowMeans(HVG_AllZT_Mean[,2:13],na.rm=TRUE),
                                  Mean_ZT18_meanNorm=Mean_ZT18/rowMeans(HVG_AllZT_Mean[,2:13],na.rm=TRUE),
                                  Mean_ZT20_meanNorm=Mean_ZT20/rowMeans(HVG_AllZT_Mean[,2:13],na.rm=TRUE),
                                  Mean_ZT22_meanNorm=Mean_ZT22/rowMeans(HVG_AllZT_Mean[,2:13],na.rm=TRUE),
                                  Mean_ZT24_meanNorm=Mean_ZT24/rowMeans(HVG_AllZT_Mean[,2:13],na.rm=TRUE))


# Plotting

heatmap.2(as.matrix(HVG_AllZT_Mean_meanNorm[,14:25]), Colv=NA, dendrogram='row', 
          Rowv=HVG_AllZT_BioVar.pearson.dist.clust.complete.den, 
          col= purpleblackyellow,breaks=seq(0.8,1.2,length.out=101), trace = 'none' ,
          labRow=NA, labCol=names(HVG_AllZT_Mean_meanNorm[,2:13]), 
          main="Highly Variable Genes",ColSideColors=ColDay,
          RowSideColors = HVG_AllZT_BioVar.pearson.dist.clust.complete.col,
          na.rm=F, na.color="black")


### Clustering and heatmap of Corrected CV2 for LVGs ----

# Data import and formatting

AllZT <- read.table("Heatmaps_CV2_expression/BrenneckeMod_allZT_RUV2.txt",
                    header=TRUE, sep='\t')

AllZT_BioVar <- select(AllZT, Gene, contains("BioVar"))

AllZT_BioVar[is.na(AllZT_BioVar)]<-0

LVG <- read.table("Heatmaps_CV2_expression/LVG1000_allZT_nbZT.txt",header=TRUE,sep='\t')

LVG_allZT_BioVar <- inner_join(AllZT_BioVar,LVG, by=c("Gene"="allZT")) %>%
  select(., Gene, contains("BioVar"))

# Clustering and heatmap plotting

LVG_AllZT_BioVar.pearson <- cor(t(LVG_allZT_BioVar[,c(2:13)]), method = "pearson")
LVG_AllZT_BioVar.pearson[is.na(LVG_AllZT_BioVar.pearson)]<-0
LVG_AllZT_BioVar.pearson.dist<-as.dist((1-LVG_AllZT_BioVar.pearson))
LVG_AllZT_BioVar.pearson.dist[is.na(LVG_AllZT_BioVar.pearson.dist)]<-0
LVG_AllZT_BioVar.pearson.dist.clust.complete <- hclust(LVG_AllZT_BioVar.pearson.dist, 
                                                       method="complete")
LVG_AllZT_BioVar.pearson.dist.clust.complete.den <- as.dendrogram(LVG_AllZT_BioVar.pearson.dist.clust.complete)
LVG_AllZT_BioVar.pearson.dist.clust.complete.col=brewer.pal(12, 'Set3')[cutree(LVG_AllZT_BioVar.pearson.dist.clust.complete, k = 4)]
LVG_allZT_BioVar$cluster_4_pearson <- cutree(LVG_AllZT_BioVar.pearson.dist.clust.complete, 
                                             k = 4)
LVG_allZT_BioVar$cluster_4_pearson_color <- brewer.pal(12, 'Set3')[cutree(LVG_AllZT_BioVar.pearson.dist.clust.complete, k = 4)]


heatmap.2(as.matrix(LVG_allZT_BioVar[,2:13]), Colv=NA, dendrogram='row', 
          Rowv=LVG_AllZT_BioVar.pearson.dist.clust.complete.den, 
          col= purpleblackyellow,breaks=seq(-3.5,0.1,length.out=101), trace = 'none' ,
          labRow=NA, labCol=names(LVG_allZT_BioVar[,2:13]), 
          main="Lowly Variable Genes",ColSideColors=ColDay,
          RowSideColors = LVG_AllZT_BioVar.pearson.dist.clust.complete.col,
          na.rm=F, na.color="black")


### Clustering and heatmap of Corrected CV2 for une set of random genes ----

# Data import and formatting

AllZT <- read.table("Heatmaps_CV2_expression/BrenneckeMod_allZT_RUV2.txt",
                    header=TRUE, sep='\t')

AllZT_BioVar <- select(AllZT, Gene, contains("BioVar"))

AllZT_BioVar[is.na(AllZT_BioVar)]<-0

random <- read.table("Heatmaps_CV2_expression/Random_HVGnumber_allZT_nbZT.txt",header=TRUE,sep='\t')

random_allZT_BioVar <- inner_join(AllZT_BioVar,random, by=c("Gene"="allZT")) %>%
  select(., Gene, contains("BioVar"))

# Clustering and heatmap plotting

random_AllZT_BioVar.pearson <- cor(t(random_allZT_BioVar[,c(2:13)]), method = "pearson")
random_AllZT_BioVar.pearson[is.na(random_AllZT_BioVar.pearson)]<-0
random_AllZT_BioVar.pearson.dist<-as.dist((1-random_AllZT_BioVar.pearson))
random_AllZT_BioVar.pearson.dist[is.na(random_AllZT_BioVar.pearson.dist)]<-0
random_AllZT_BioVar.pearson.dist.clust.complete <- hclust(random_AllZT_BioVar.pearson.dist, 
                                                       method="complete")
random_AllZT_BioVar.pearson.dist.clust.complete.den <- as.dendrogram(random_AllZT_BioVar.pearson.dist.clust.complete)
random_AllZT_BioVar.pearson.dist.clust.complete.col=brewer.pal(12, 'Set3')[cutree(random_AllZT_BioVar.pearson.dist.clust.complete, k = 4)]
random_allZT_BioVar$cluster_4_pearson <- cutree(random_AllZT_BioVar.pearson.dist.clust.complete, 
                                             k = 4)
random_allZT_BioVar$cluster_4_pearson_color <- brewer.pal(12, 'Set3')[cutree(random_AllZT_BioVar.pearson.dist.clust.complete, k = 4)]


heatmap.2(as.matrix(random_allZT_BioVar[,2:13]), Colv=NA, dendrogram='row', 
          Rowv=random_AllZT_BioVar.pearson.dist.clust.complete.den, 
          col= purpleblackyellow,breaks=seq(-3.5,1.5,length.out=101), trace = 'none' ,
          labRow=NA, labCol=names(random_allZT_BioVar[,2:13]), 
          main="Random Genes",ColSideColors=ColDay,
          RowSideColors = random_AllZT_BioVar.pearson.dist.clust.complete.col,
          na.rm=F, na.color="black")

