
library(tidyverse)
library(gplots)

purpleblackyellow <- colorRampPalette(c("purple","black", "yellow"))(n = 100)
yellowblackpurple <- colorRampPalette(c("yellow","black", "purple"))(n = 100)

# Load expression data ----

tissue_expression <- read.table("Entropy/E-TABM-17_normRatio_noMutant_min7.txt",header=T, sep='\t')

# HVG ----

HVG <- read.table("Entropy/HVG_allZT_nbZT.txt",header=T, sep='\t')

tissue_expression_HVG <- inner_join(tissue_expression, HVG,
                                    by=c("Name"="allZT"))

tissue_expression_HVG_reduced <- select(tissue_expression_HVG, Name,
                                        cauline.leaf.time.21.,
                                        cotyledon.time.7,
                                        flower.time.21.,
                                        hypocotyl.time.7,
                                        inflorescence...shoot.apex.time.14,
                                        internode.time.21.,
                                        leaf...shoot.apex.time.7,
                                        leaf.time.15,
                                        node.time.21.,
                                        petiole.time.17,
                                        root.time.15,
                                        rosette.leaf.time.10,
                                        seed.time.56,
                                        shoot.time.21,
                                        silique.time.56)


names(tissue_expression_HVG_reduced) <- c("Name",
                                          "cauline leaf",
                                          "cotyledon",
                                          "flower",
                                          "hypocotyl",
                                          "inflorescence and shoot apex",
                                          "internode",
                                          "leaf and shoot apex",
                                          "leaf",
                                          "node",
                                          "petiole",
                                          "root",
                                          "rosette leaf",
                                          "seed",
                                          "shoot",
                                          "silique")

expression.pearson <- cor(t(tissue_expression_HVG_reduced[,2:16]), method = "spearman")
expression.pearson[is.na(expression.pearson)]<-0
expression.pearson.dist<-as.dist((1-expression.pearson))
expression.pearson.dist[is.na(expression.pearson.dist)]<-0
expression.pearson.dist.clust.complete <- hclust(expression.pearson.dist, method="complete")
expression.pearson.dist.clust.complete.den <- as.dendrogram(expression.pearson.dist.clust.complete)


heatmap.2(as.matrix(tissue_expression_HVG_reduced[,2:16]), 
          Colv = NA, dendrogram='row', Rowv=expression.pearson.dist.clust.complete.den,
          col= purpleblackyellow,breaks=seq(3,10,length.out=101), trace = 'none' , 
          labRow=NA, labCol=names(tissue_expression_HVG_reduced[,2:16]), main="expression, HVG",
          na.rm=F, na.color="black", cexCol = 2,margins=c(23.5,5))



