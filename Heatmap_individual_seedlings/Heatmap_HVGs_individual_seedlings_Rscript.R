

library(RColorBrewer)
library(tidyverse)
library(gdata)
library(gplots)

purpleblackyellowExtra <- colorRampPalette(c("violet","purple","black", "gold", "yellow"))(n = 100)
purpleblackyellow <- colorRampPalette(c("purple","black", "yellow"))(n = 100)

# Heatmap for individual seedlings on mean normalised expression for highly variable genes -----

# ZT2
Expression_ZT2 <- read.table("Heatmap_individual_seedlings/TPM2_ZT2_genes.txt",
                             header = T, sep='\t')

HVG_ZT2 <- read.table("Heatmap_individual_seedlings/HVG_ZT2.txt", header=T, sep='\t')


HVG_ZT2_expression <- inner_join(HVG_ZT2, Expression_ZT2,
                                 by=c("Gene"="gene"))


HVG_ZT2_expression <- mutate(HVG_ZT2_expression,
                             Seedling1_meanNorm=TPM_WT_ZT2_seedling_1/rowMeans(HVG_ZT2_expression[,9:22],na.rm=TRUE),
                             Seedling3_meanNorm=TPM_WT_ZT2_seedling_3/rowMeans(HVG_ZT2_expression[,9:22],na.rm=TRUE),
                             Seedling4_meanNorm=TPM_WT_ZT2_seedling_4/rowMeans(HVG_ZT2_expression[,9:22],na.rm=TRUE),
                             Seedling5_meanNorm=TPM_WT_ZT2_seedling_5/rowMeans(HVG_ZT2_expression[,9:22],na.rm=TRUE),
                             Seedling6_meanNorm=TPM_WT_ZT2_seedling_6/rowMeans(HVG_ZT2_expression[,9:22],na.rm=TRUE),
                             Seedling7_meanNorm=TPM_WT_ZT2_seedling_7/rowMeans(HVG_ZT2_expression[,9:22],na.rm=TRUE),
                             Seedling8_meanNorm=TPM_WT_ZT2_seedling_8/rowMeans(HVG_ZT2_expression[,9:22],na.rm=TRUE),
                             Seedling9_meanNorm=TPM_WT_ZT2_seedling_9/rowMeans(HVG_ZT2_expression[,9:22],na.rm=TRUE),
                             Seedling10_meanNorm=TPM_WT_ZT2_seedling_10/rowMeans(HVG_ZT2_expression[,9:22],na.rm=TRUE),
                             Seedling11_meanNorm=TPM_WT_ZT2_seedling_11/rowMeans(HVG_ZT2_expression[,9:22],na.rm=TRUE),
                             Seedling12_meanNorm=TPM_WT_ZT2_seedling_12/rowMeans(HVG_ZT2_expression[,9:22],na.rm=TRUE),
                             Seedling14_meanNorm=TPM_WT_ZT2_seedling_14/rowMeans(HVG_ZT2_expression[,9:22],na.rm=TRUE),
                             Seedling15_meanNorm=TPM_WT_ZT2_seedling_15/rowMeans(HVG_ZT2_expression[,9:22],na.rm=TRUE),
                             Seedling16_meanNorm=TPM_WT_ZT2_seedling_16/rowMeans(HVG_ZT2_expression[,9:22],na.rm=TRUE))



heatmap.2(as.matrix(HVG_ZT2_expression[,23:36]), dendrogram='both', 
          Rowv=TRUE, 
          col= purpleblackyellow,
          breaks=seq(0.8,1.2,length.out=101), trace = 'none' ,
          labRow=NA, labCol=names(HVG_ZT2_expression[,9:22]), 
          main="HVG, ZT2",
          na.rm=F, na.color="black", margins=c(6,2))




# ZT4
Expression_ZT4 <- read.table("Heatmap_individual_seedlings/TPM2_ZT4_genes.txt",
                             header = T, sep='\t')

HVG_ZT4 <- read.table("Heatmap_individual_seedlings/HVG_ZT4.txt", header=T, sep='\t')

HVG_ZT4_expression <- inner_join(HVG_ZT4, Expression_ZT4,
                                 by=c("Gene"="gene"))

HVG_ZT4_expression <- mutate(HVG_ZT4_expression,
                             Seedling1_meanNorm=TPM_WT_ZT4_seedling_1/rowMeans(HVG_ZT4_expression[,9:22],na.rm=TRUE),
                             Seedling2_meanNorm=TPM_WT_ZT4_seedling_2/rowMeans(HVG_ZT4_expression[,9:22],na.rm=TRUE),
                             Seedling3_meanNorm=TPM_WT_ZT4_seedling_3/rowMeans(HVG_ZT4_expression[,9:22],na.rm=TRUE),
                             Seedling4_meanNorm=TPM_WT_ZT4_seedling_4/rowMeans(HVG_ZT4_expression[,9:22],na.rm=TRUE),
                             Seedling5_meanNorm=TPM_WT_ZT4_seedling_5/rowMeans(HVG_ZT4_expression[,9:22],na.rm=TRUE),
                             Seedling6_meanNorm=TPM_WT_ZT4_seedling_6/rowMeans(HVG_ZT4_expression[,9:22],na.rm=TRUE),
                             Seedling7_meanNorm=TPM_WT_ZT4_seedling_7/rowMeans(HVG_ZT4_expression[,9:22],na.rm=TRUE),
                             Seedling8_meanNorm=TPM_WT_ZT4_seedling_8/rowMeans(HVG_ZT4_expression[,9:22],na.rm=TRUE),
                             Seedling9_meanNorm=TPM_WT_ZT4_seedling_9/rowMeans(HVG_ZT4_expression[,9:22],na.rm=TRUE),
                             Seedling10_meanNorm=TPM_WT_ZT4_seedling_10/rowMeans(HVG_ZT4_expression[,9:22],na.rm=TRUE),
                             Seedling13_meanNorm=TPM_WT_ZT4_seedling_13/rowMeans(HVG_ZT4_expression[,9:22],na.rm=TRUE),
                             Seedling14_meanNorm=TPM_WT_ZT4_seedling_14/rowMeans(HVG_ZT4_expression[,9:22],na.rm=TRUE),
                             Seedling15_meanNorm=TPM_WT_ZT4_seedling_15/rowMeans(HVG_ZT4_expression[,9:22],na.rm=TRUE),
                             Seedling16_meanNorm=TPM_WT_ZT4_seedling_16/rowMeans(HVG_ZT4_expression[,9:22],na.rm=TRUE))


heatmap.2(as.matrix(HVG_ZT4_expression[,23:36]), dendrogram='both', 
          Rowv=TRUE, 
          col= purpleblackyellow,
          breaks=seq(0.8,1.2,length.out=101), trace = 'none' ,
          labRow=NA, labCol=names(HVG_ZT4_expression[,9:22]), 
          main="HVG, ZT4",
          na.rm=F, na.color="black", margins=c(6,2))




# ZT6
Expression_ZT6 <- read.table("Heatmap_individual_seedlings/TPM2_ZT6_genes.txt",
                             header = T, sep='\t')

HVG_ZT6 <- read.table("Heatmap_individual_seedlings/HVG_ZT6.txt", header=T, sep='\t')

HVG_ZT6_expression <- inner_join(HVG_ZT6, Expression_ZT6,
                                 by=c("Gene"="gene"))

HVG_ZT6_expression <- mutate(HVG_ZT6_expression,
                             Seedling1_meanNorm=TPM_WT_ZT6_seedling_1/rowMeans(HVG_ZT6_expression[,9:22],na.rm=TRUE),
                             Seedling2_meanNorm=TPM_WT_ZT6_seedling_2/rowMeans(HVG_ZT6_expression[,9:22],na.rm=TRUE),
                             Seedling3_meanNorm=TPM_WT_ZT6_seedling_3/rowMeans(HVG_ZT6_expression[,9:22],na.rm=TRUE),
                             Seedling4_meanNorm=TPM_WT_ZT6_seedling_4/rowMeans(HVG_ZT6_expression[,9:22],na.rm=TRUE),
                             Seedling5_meanNorm=TPM_WT_ZT6_seedling_5/rowMeans(HVG_ZT6_expression[,9:22],na.rm=TRUE),
                             Seedling6_meanNorm=TPM_WT_ZT6_seedling_6/rowMeans(HVG_ZT6_expression[,9:22],na.rm=TRUE),
                             Seedling7_meanNorm=TPM_WT_ZT6_seedling_7/rowMeans(HVG_ZT6_expression[,9:22],na.rm=TRUE),
                             Seedling8_meanNorm=TPM_WT_ZT6_seedling_8/rowMeans(HVG_ZT6_expression[,9:22],na.rm=TRUE),
                             Seedling9_meanNorm=TPM_WT_ZT6_seedling_9/rowMeans(HVG_ZT6_expression[,9:22],na.rm=TRUE),
                             Seedling10_meanNorm=TPM_WT_ZT6_seedling_10/rowMeans(HVG_ZT6_expression[,9:22],na.rm=TRUE),
                             Seedling11_meanNorm=TPM_WT_ZT6_seedling_11/rowMeans(HVG_ZT6_expression[,9:22],na.rm=TRUE),
                             Seedling14_meanNorm=TPM_WT_ZT6_seedling_14/rowMeans(HVG_ZT6_expression[,9:22],na.rm=TRUE),
                             Seedling15_meanNorm=TPM_WT_ZT6_seedling_15/rowMeans(HVG_ZT6_expression[,9:22],na.rm=TRUE),
                             Seedling16_meanNorm=TPM_WT_ZT6_seedling_16/rowMeans(HVG_ZT6_expression[,9:22],na.rm=TRUE))



heatmap.2(as.matrix(HVG_ZT6_expression[,23:36]), dendrogram='both', 
          Rowv=TRUE, 
          col= purpleblackyellow,
          breaks=seq(0.8,1.2,length.out=101), trace = 'none' ,
          labRow=NA, labCol=names(HVG_ZT6_expression[,9:22]), 
          main="HVG, ZT6",
          na.rm=F, na.color="black", margins=c(6,2))



# ZT8
Expression_ZT8 <- read.table("Heatmap_individual_seedlings/TPM2_ZT8_genes.txt",
                             header = T, sep='\t')

HVG_ZT8 <- read.table("Heatmap_individual_seedlings/HVG_ZT8.txt", header=T, sep='\t')

HVG_ZT8_expression <- inner_join(HVG_ZT8, Expression_ZT8,
                                 by=c("Gene"="gene"))

HVG_ZT8_expression <- mutate(HVG_ZT8_expression,
                             Seedling1_meanNorm=TPM_WT_ZT8_seedling_1/rowMeans(HVG_ZT8_expression[,9:22],na.rm=TRUE),
                             Seedling2_meanNorm=TPM_WT_ZT8_seedling_2/rowMeans(HVG_ZT8_expression[,9:22],na.rm=TRUE),
                             Seedling3_meanNorm=TPM_WT_ZT8_seedling_3/rowMeans(HVG_ZT8_expression[,9:22],na.rm=TRUE),
                             Seedling4_meanNorm=TPM_WT_ZT8_seedling_4/rowMeans(HVG_ZT8_expression[,9:22],na.rm=TRUE),
                             Seedling5_meanNorm=TPM_WT_ZT8_seedling_5/rowMeans(HVG_ZT8_expression[,9:22],na.rm=TRUE),
                             Seedling6_meanNorm=TPM_WT_ZT8_seedling_6/rowMeans(HVG_ZT8_expression[,9:22],na.rm=TRUE),
                             Seedling7_meanNorm=TPM_WT_ZT8_seedling_7/rowMeans(HVG_ZT8_expression[,9:22],na.rm=TRUE),
                             Seedling9_meanNorm=TPM_WT_ZT8_seedling_9/rowMeans(HVG_ZT8_expression[,9:22],na.rm=TRUE),
                             Seedling10_meanNorm=TPM_WT_ZT8_seedling_10/rowMeans(HVG_ZT8_expression[,9:22],na.rm=TRUE),
                             Seedling11_meanNorm=TPM_WT_ZT8_seedling_11/rowMeans(HVG_ZT8_expression[,9:22],na.rm=TRUE),
                             Seedling12_meanNorm=TPM_WT_ZT8_seedling_12/rowMeans(HVG_ZT8_expression[,9:22],na.rm=TRUE),
                             Seedling13_meanNorm=TPM_WT_ZT8_seedling_13/rowMeans(HVG_ZT8_expression[,9:22],na.rm=TRUE),
                             Seedling14_meanNorm=TPM_WT_ZT8_seedling_14/rowMeans(HVG_ZT8_expression[,9:22],na.rm=TRUE),
                             Seedling16_meanNorm=TPM_WT_ZT8_seedling_16/rowMeans(HVG_ZT8_expression[,9:22],na.rm=TRUE))



heatmap.2(as.matrix(HVG_ZT8_expression[,23:36]), dendrogram='both', 
          Rowv=TRUE, 
          col= purpleblackyellow,
          breaks=seq(0.8,1.2,length.out=101), trace = 'none' ,
          labRow=NA, labCol=names(HVG_ZT8_expression[,9:22]), 
          main="HVG, ZT8",
          na.rm=F, na.color="black", margins=c(6,2))



# ZT10
Expression_ZT10 <- read.table("Heatmap_individual_seedlings/TPM2_ZT10_genes.txt",
                              header = T, sep='\t')

HVG_ZT10 <- read.table("Heatmap_individual_seedlings/HVG_ZT10.txt", header=T, sep='\t')

HVG_ZT10_expression <- inner_join(HVG_ZT10, Expression_ZT10,
                                  by=c("Gene"="gene"))

HVG_ZT10_expression <- mutate(HVG_ZT10_expression,
                              Seedling1_meanNorm=TPM_WT_ZT10_seedling_1/rowMeans(HVG_ZT10_expression[,9:22],na.rm=TRUE),
                              Seedling2_meanNorm=TPM_WT_ZT10_seedling_2/rowMeans(HVG_ZT10_expression[,9:22],na.rm=TRUE),
                              Seedling3_meanNorm=TPM_WT_ZT10_seedling_3/rowMeans(HVG_ZT10_expression[,9:22],na.rm=TRUE),
                              Seedling4_meanNorm=TPM_WT_ZT10_seedling_4/rowMeans(HVG_ZT10_expression[,9:22],na.rm=TRUE),
                              Seedling5_meanNorm=TPM_WT_ZT10_seedling_5/rowMeans(HVG_ZT10_expression[,9:22],na.rm=TRUE),
                              Seedling6_meanNorm=TPM_WT_ZT10_seedling_6/rowMeans(HVG_ZT10_expression[,9:22],na.rm=TRUE),
                              Seedling7_meanNorm=TPM_WT_ZT10_seedling_7/rowMeans(HVG_ZT10_expression[,9:22],na.rm=TRUE),
                              Seedling8_meanNorm=TPM_WT_ZT10_seedling_8/rowMeans(HVG_ZT10_expression[,9:22],na.rm=TRUE),
                              Seedling9_meanNorm=TPM_WT_ZT10_seedling_9/rowMeans(HVG_ZT10_expression[,9:22],na.rm=TRUE),
                              Seedling10_meanNorm=TPM_WT_ZT10_seedling_10/rowMeans(HVG_ZT10_expression[,9:22],na.rm=TRUE),
                              Seedling11_meanNorm=TPM_WT_ZT10_seedling_11/rowMeans(HVG_ZT10_expression[,9:22],na.rm=TRUE),
                              Seedling13_meanNorm=TPM_WT_ZT10_seedling_13/rowMeans(HVG_ZT10_expression[,9:22],na.rm=TRUE),
                              Seedling14_meanNorm=TPM_WT_ZT10_seedling_14/rowMeans(HVG_ZT10_expression[,9:22],na.rm=TRUE),
                              Seedling16_meanNorm=TPM_WT_ZT10_seedling_16/rowMeans(HVG_ZT10_expression[,9:22],na.rm=TRUE))


heatmap.2(as.matrix(HVG_ZT10_expression[,23:36]), dendrogram='both', 
          Rowv=TRUE, 
          col= purpleblackyellow,
          breaks=seq(0.8,1.2,length.out=101), trace = 'none' ,
          labRow=NA, labCol=names(HVG_ZT10_expression[,9:22]), 
          main="HVG, ZT10",
          na.rm=F, na.color="black", margins=c(6,2))



# ZT12
Expression_ZT12 <- read.table("Heatmap_individual_seedlings/TPM2_ZT12_genes.txt",
                              header = T, sep='\t')

HVG_ZT12 <- read.table("Heatmap_individual_seedlings/HVG_ZT12.txt", header=T, sep='\t')

HVG_ZT12_expression <- inner_join(HVG_ZT12, Expression_ZT12,
                                  by=c("Gene"="gene"))

HVG_ZT12_expression <- mutate(HVG_ZT12_expression,
                              Seedling1_meanNorm=TPM_WT_ZT12_seedling_1/rowMeans(HVG_ZT12_expression[,9:22],na.rm=TRUE),
                              Seedling2_meanNorm=TPM_WT_ZT12_seedling_2/rowMeans(HVG_ZT12_expression[,9:22],na.rm=TRUE),
                              Seedling5_meanNorm=TPM_WT_ZT12_seedling_5/rowMeans(HVG_ZT12_expression[,9:22],na.rm=TRUE),
                              Seedling6_meanNorm=TPM_WT_ZT12_seedling_6/rowMeans(HVG_ZT12_expression[,9:22],na.rm=TRUE),
                              Seedling7_meanNorm=TPM_WT_ZT12_seedling_7/rowMeans(HVG_ZT12_expression[,9:22],na.rm=TRUE),
                              Seedling8_meanNorm=TPM_WT_ZT12_seedling_8/rowMeans(HVG_ZT12_expression[,9:22],na.rm=TRUE),
                              Seedling9_meanNorm=TPM_WT_ZT12_seedling_9/rowMeans(HVG_ZT12_expression[,9:22],na.rm=TRUE),
                              Seedling10_meanNorm=TPM_WT_ZT12_seedling_10/rowMeans(HVG_ZT12_expression[,9:22],na.rm=TRUE),
                              Seedling11_meanNorm=TPM_WT_ZT12_seedling_11/rowMeans(HVG_ZT12_expression[,9:22],na.rm=TRUE),
                              Seedling12_meanNorm=TPM_WT_ZT12_seedling_12/rowMeans(HVG_ZT12_expression[,9:22],na.rm=TRUE),
                              Seedling13_meanNorm=TPM_WT_ZT12_seedling_13/rowMeans(HVG_ZT12_expression[,9:22],na.rm=TRUE),
                              Seedling14_meanNorm=TPM_WT_ZT12_seedling_14/rowMeans(HVG_ZT12_expression[,9:22],na.rm=TRUE),
                              Seedling15_meanNorm=TPM_WT_ZT12_seedling_15/rowMeans(HVG_ZT12_expression[,9:22],na.rm=TRUE),
                              Seedling16_meanNorm=TPM_WT_ZT12_seedling_16/rowMeans(HVG_ZT12_expression[,9:22],na.rm=TRUE))



heatmap.2(as.matrix(HVG_ZT12_expression[,23:36]), dendrogram='both', 
          Rowv=TRUE, 
          col= purpleblackyellow,
          breaks=seq(0.8,1.2,length.out=101), trace = 'none' ,
          labRow=NA, labCol=names(HVG_ZT12_expression[,9:22]), 
          main="HVG, ZT12",
          na.rm=F, na.color="black", margins=c(6,2))






# ZT14
Expression_ZT14 <- read.table("Heatmap_individual_seedlings/TPM2_ZT14_genes.txt",
                              header = T, sep='\t')

HVG_ZT14 <- read.table("Heatmap_individual_seedlings/HVG_ZT14.txt", header=T, sep='\t')

HVG_ZT14_expression <- inner_join(HVG_ZT14, Expression_ZT14,
                                  by=c("Gene"="gene"))

HVG_ZT14_expression <- mutate(HVG_ZT14_expression,
                              Seedling1_meanNorm=TPM_WT_ZT14_seedling_1/rowMeans(HVG_ZT14_expression[,9:22],na.rm=TRUE),
                              Seedling2_meanNorm=TPM_WT_ZT14_seedling_2/rowMeans(HVG_ZT14_expression[,9:22],na.rm=TRUE),
                              Seedling4_meanNorm=TPM_WT_ZT14_seedling_4/rowMeans(HVG_ZT14_expression[,9:22],na.rm=TRUE),
                              Seedling5_meanNorm=TPM_WT_ZT14_seedling_5/rowMeans(HVG_ZT14_expression[,9:22],na.rm=TRUE),
                              Seedling6_meanNorm=TPM_WT_ZT14_seedling_6/rowMeans(HVG_ZT14_expression[,9:22],na.rm=TRUE),
                              Seedling7_meanNorm=TPM_WT_ZT14_seedling_7/rowMeans(HVG_ZT14_expression[,9:22],na.rm=TRUE),
                              Seedling8_meanNorm=TPM_WT_ZT14_seedling_8/rowMeans(HVG_ZT14_expression[,9:22],na.rm=TRUE),
                              Seedling9_meanNorm=TPM_WT_ZT14_seedling_9/rowMeans(HVG_ZT14_expression[,9:22],na.rm=TRUE),
                              Seedling11_meanNorm=TPM_WT_ZT14_seedling_11/rowMeans(HVG_ZT14_expression[,9:22],na.rm=TRUE),
                              Seedling12_meanNorm=TPM_WT_ZT14_seedling_12/rowMeans(HVG_ZT14_expression[,9:22],na.rm=TRUE),
                              Seedling13_meanNorm=TPM_WT_ZT14_seedling_13/rowMeans(HVG_ZT14_expression[,9:22],na.rm=TRUE),
                              Seedling14_meanNorm=TPM_WT_ZT14_seedling_14/rowMeans(HVG_ZT14_expression[,9:22],na.rm=TRUE),
                              Seedling15_meanNorm=TPM_WT_ZT14_seedling_15/rowMeans(HVG_ZT14_expression[,9:22],na.rm=TRUE),
                              Seedling16_meanNorm=TPM_WT_ZT14_seedling_16/rowMeans(HVG_ZT14_expression[,9:22],na.rm=TRUE))



heatmap.2(as.matrix(HVG_ZT14_expression[,23:36]), dendrogram='both', 
          Rowv=TRUE, 
          col= purpleblackyellow,
          breaks=seq(0.8,1.2,length.out=101), trace = 'none' ,
          labRow=NA, labCol=names(HVG_ZT14_expression[,9:22]), 
          main="HVG, ZT14",
          na.rm=F, na.color="black", margins=c(6,2))




# ZT16
Expression_ZT16 <- read.table("Heatmap_individual_seedlings/TPM2_ZT16_genes.txt",
                              header = T, sep='\t')

HVG_ZT16 <- read.table("Heatmap_individual_seedlings/HVG_ZT16.txt", header=T, sep='\t')

HVG_ZT16_expression <- inner_join(HVG_ZT16, Expression_ZT16,
                                  by=c("Gene"="gene"))

HVG_ZT16_expression <- mutate(HVG_ZT16_expression,
                              Seedling1_meanNorm=TPM_WT_ZT16_seedling_1/rowMeans(HVG_ZT16_expression[,9:22],na.rm=TRUE),
                              Seedling2_meanNorm=TPM_WT_ZT16_seedling_2/rowMeans(HVG_ZT16_expression[,9:22],na.rm=TRUE),
                              Seedling3_meanNorm=TPM_WT_ZT16_seedling_3/rowMeans(HVG_ZT16_expression[,9:22],na.rm=TRUE),
                              Seedling4_meanNorm=TPM_WT_ZT16_seedling_4/rowMeans(HVG_ZT16_expression[,9:22],na.rm=TRUE),
                              Seedling5_meanNorm=TPM_WT_ZT16_seedling_5/rowMeans(HVG_ZT16_expression[,9:22],na.rm=TRUE),
                              Seedling6_meanNorm=TPM_WT_ZT16_seedling_6/rowMeans(HVG_ZT16_expression[,9:22],na.rm=TRUE),
                              Seedling7_meanNorm=TPM_WT_ZT16_seedling_7/rowMeans(HVG_ZT16_expression[,9:22],na.rm=TRUE),
                              Seedling9_meanNorm=TPM_WT_ZT16_seedling_9/rowMeans(HVG_ZT16_expression[,9:22],na.rm=TRUE),
                              Seedling10_meanNorm=TPM_WT_ZT16_seedling_10/rowMeans(HVG_ZT16_expression[,9:22],na.rm=TRUE),
                              Seedling11_meanNorm=TPM_WT_ZT16_seedling_11/rowMeans(HVG_ZT16_expression[,9:22],na.rm=TRUE),
                              Seedling13_meanNorm=TPM_WT_ZT16_seedling_13/rowMeans(HVG_ZT16_expression[,9:22],na.rm=TRUE),
                              Seedling14_meanNorm=TPM_WT_ZT16_seedling_14/rowMeans(HVG_ZT16_expression[,9:22],na.rm=TRUE),
                              Seedling15_meanNorm=TPM_WT_ZT16_seedling_15/rowMeans(HVG_ZT16_expression[,9:22],na.rm=TRUE),
                              Seedling16_meanNorm=TPM_WT_ZT16_seedling_16/rowMeans(HVG_ZT16_expression[,9:22],na.rm=TRUE))



heatmap.2(as.matrix(HVG_ZT16_expression[,23:36]), dendrogram='both', 
          Rowv=TRUE, 
          col= purpleblackyellow,
          breaks=seq(0.8,1.2,length.out=101), trace = 'none' ,
          labRow=NA, labCol=names(HVG_ZT16_expression[,9:22]), 
          main="HVG, ZT16",
          na.rm=F, na.color="black", margins=c(6,2))




# ZT18
Expression_ZT18 <- read.table("Heatmap_individual_seedlings/TPM2_ZT18_genes.txt",
                              header = T, sep='\t')

HVG_ZT18 <- read.table("Heatmap_individual_seedlings/HVG_ZT18.txt", header=T, sep='\t')

HVG_ZT18_expression <- inner_join(HVG_ZT18, Expression_ZT18,
                                  by=c("Gene"="gene"))

HVG_ZT18_expression <- mutate(HVG_ZT18_expression,
                              Seedling1_meanNorm=TPM_WT_ZT18_seedling_1/rowMeans(HVG_ZT18_expression[,9:22],na.rm=TRUE),
                              Seedling2_meanNorm=TPM_WT_ZT18_seedling_2/rowMeans(HVG_ZT18_expression[,9:22],na.rm=TRUE),
                              Seedling3_meanNorm=TPM_WT_ZT18_seedling_3/rowMeans(HVG_ZT18_expression[,9:22],na.rm=TRUE),
                              Seedling5_meanNorm=TPM_WT_ZT18_seedling_5/rowMeans(HVG_ZT18_expression[,9:22],na.rm=TRUE),
                              Seedling7_meanNorm=TPM_WT_ZT18_seedling_7/rowMeans(HVG_ZT18_expression[,9:22],na.rm=TRUE),
                              Seedling8_meanNorm=TPM_WT_ZT18_seedling_8/rowMeans(HVG_ZT18_expression[,9:22],na.rm=TRUE),
                              Seedling9_meanNorm=TPM_WT_ZT18_seedling_9/rowMeans(HVG_ZT18_expression[,9:22],na.rm=TRUE),
                              Seedling10_meanNorm=TPM_WT_ZT18_seedling_10/rowMeans(HVG_ZT18_expression[,9:22],na.rm=TRUE),
                              Seedling11_meanNorm=TPM_WT_ZT18_seedling_11/rowMeans(HVG_ZT18_expression[,9:22],na.rm=TRUE),
                              Seedling12_meanNorm=TPM_WT_ZT18_seedling_12/rowMeans(HVG_ZT18_expression[,9:22],na.rm=TRUE),
                              Seedling13_meanNorm=TPM_WT_ZT18_seedling_13/rowMeans(HVG_ZT18_expression[,9:22],na.rm=TRUE),
                              Seedling14_meanNorm=TPM_WT_ZT18_seedling_14/rowMeans(HVG_ZT18_expression[,9:22],na.rm=TRUE),
                              Seedling15_meanNorm=TPM_WT_ZT18_seedling_15/rowMeans(HVG_ZT18_expression[,9:22],na.rm=TRUE),
                              Seedling16_meanNorm=TPM_WT_ZT18_seedling_16/rowMeans(HVG_ZT18_expression[,9:22],na.rm=TRUE))



heatmap.2(as.matrix(HVG_ZT18_expression[,23:36]), dendrogram='both', 
          Rowv=TRUE, 
          col= purpleblackyellow,
          breaks=seq(0.8,1.2,length.out=101), trace = 'none' ,
          labRow=NA, labCol=names(HVG_ZT18_expression[,9:22]), 
          main="HVG, ZT18",
          na.rm=F, na.color="black", margins=c(6,2))





# ZT20
Expression_ZT20 <- read.table("Heatmap_individual_seedlings/TPM2_ZT20_genes.txt",
                              header = T, sep='\t')

HVG_ZT20 <- read.table("Heatmap_individual_seedlings/HVG_ZT20.txt", header=T, sep='\t')

HVG_ZT20_expression <- inner_join(HVG_ZT20, Expression_ZT20,
                                  by=c("Gene"="gene"))

HVG_ZT20_expression <- mutate(HVG_ZT20_expression,
                              Seedling1_meanNorm=TPM_WT_ZT20_seedling_1/rowMeans(HVG_ZT20_expression[,9:22],na.rm=TRUE),
                              Seedling3_meanNorm=TPM_WT_ZT20_seedling_3/rowMeans(HVG_ZT20_expression[,9:22],na.rm=TRUE),
                              Seedling4_meanNorm=TPM_WT_ZT20_seedling_4/rowMeans(HVG_ZT20_expression[,9:22],na.rm=TRUE),
                              Seedling5_meanNorm=TPM_WT_ZT20_seedling_5/rowMeans(HVG_ZT20_expression[,9:22],na.rm=TRUE),
                              Seedling7_meanNorm=TPM_WT_ZT20_seedling_7/rowMeans(HVG_ZT20_expression[,9:22],na.rm=TRUE),
                              Seedling8_meanNorm=TPM_WT_ZT20_seedling_8/rowMeans(HVG_ZT20_expression[,9:22],na.rm=TRUE),
                              Seedling9_meanNorm=TPM_WT_ZT20_seedling_9/rowMeans(HVG_ZT20_expression[,9:22],na.rm=TRUE),
                              Seedling10_meanNorm=TPM_WT_ZT20_seedling_10/rowMeans(HVG_ZT20_expression[,9:22],na.rm=TRUE),
                              Seedling11_meanNorm=TPM_WT_ZT20_seedling_11/rowMeans(HVG_ZT20_expression[,9:22],na.rm=TRUE),
                              Seedling12_meanNorm=TPM_WT_ZT20_seedling_12/rowMeans(HVG_ZT20_expression[,9:22],na.rm=TRUE),
                              Seedling13_meanNorm=TPM_WT_ZT20_seedling_13/rowMeans(HVG_ZT20_expression[,9:22],na.rm=TRUE),
                              Seedling14_meanNorm=TPM_WT_ZT20_seedling_14/rowMeans(HVG_ZT20_expression[,9:22],na.rm=TRUE),
                              Seedling15_meanNorm=TPM_WT_ZT20_seedling_15/rowMeans(HVG_ZT20_expression[,9:22],na.rm=TRUE),
                              Seedling16_meanNorm=TPM_WT_ZT20_seedling_16/rowMeans(HVG_ZT20_expression[,9:22],na.rm=TRUE))



heatmap.2(as.matrix(HVG_ZT20_expression[,23:36]), dendrogram='both', 
          Rowv=TRUE, 
          col= purpleblackyellow,
          breaks=seq(0.8,1.2,length.out=101), trace = 'none' ,
          labRow=NA, labCol=names(HVG_ZT20_expression[,9:22]), 
          main="HVG, ZT20",
          na.rm=F, na.color="black", margins=c(6,2))







# ZT22
Expression_ZT22 <- read.table("Heatmap_individual_seedlings/TPM2_ZT22_genes.txt",
                              header = T, sep='\t')

HVG_ZT22 <- read.table("Heatmap_individual_seedlings/HVG_ZT22.txt", header=T, sep='\t')

HVG_ZT22_expression <- inner_join(HVG_ZT22, Expression_ZT22,
                                  by=c("Gene"="gene"))

HVG_ZT22_expression <- mutate(HVG_ZT22_expression,
                              Seedling1_meanNorm=TPM_WT_ZT22_seedling_1/rowMeans(HVG_ZT22_expression[,9:22],na.rm=TRUE),
                              Seedling2_meanNorm=TPM_WT_ZT22_seedling_2/rowMeans(HVG_ZT22_expression[,9:22],na.rm=TRUE),
                              Seedling3_meanNorm=TPM_WT_ZT22_seedling_3/rowMeans(HVG_ZT22_expression[,9:22],na.rm=TRUE),
                              Seedling4_meanNorm=TPM_WT_ZT22_seedling_4/rowMeans(HVG_ZT22_expression[,9:22],na.rm=TRUE),
                              Seedling5_meanNorm=TPM_WT_ZT22_seedling_5/rowMeans(HVG_ZT22_expression[,9:22],na.rm=TRUE),
                              Seedling8_meanNorm=TPM_WT_ZT22_seedling_8/rowMeans(HVG_ZT22_expression[,9:22],na.rm=TRUE),
                              Seedling9_meanNorm=TPM_WT_ZT22_seedling_9/rowMeans(HVG_ZT22_expression[,9:22],na.rm=TRUE),
                              Seedling10_meanNorm=TPM_WT_ZT22_seedling_10/rowMeans(HVG_ZT22_expression[,9:22],na.rm=TRUE),
                              Seedling11_meanNorm=TPM_WT_ZT22_seedling_11/rowMeans(HVG_ZT22_expression[,9:22],na.rm=TRUE),
                              Seedling12_meanNorm=TPM_WT_ZT22_seedling_12/rowMeans(HVG_ZT22_expression[,9:22],na.rm=TRUE),
                              Seedling13_meanNorm=TPM_WT_ZT22_seedling_13/rowMeans(HVG_ZT22_expression[,9:22],na.rm=TRUE),
                              Seedling14_meanNorm=TPM_WT_ZT22_seedling_14/rowMeans(HVG_ZT22_expression[,9:22],na.rm=TRUE),
                              Seedling15_meanNorm=TPM_WT_ZT22_seedling_15/rowMeans(HVG_ZT22_expression[,9:22],na.rm=TRUE),
                              Seedling16_meanNorm=TPM_WT_ZT22_seedling_16/rowMeans(HVG_ZT22_expression[,9:22],na.rm=TRUE))



heatmap.2(as.matrix(HVG_ZT22_expression[,23:36]), dendrogram='both', 
          Rowv=TRUE, 
          col= purpleblackyellow,
          breaks=seq(0.8,1.2,length.out=101), trace = 'none' ,
          labRow=NA, labCol=names(HVG_ZT22_expression[,9:22]), 
          main="HVG, ZT22",
          na.rm=F, na.color="black", margins=c(6,2))







# ZT24
Expression_ZT24 <- read.table("Heatmap_individual_seedlings/TPM2_ZT0_genes.txt",
                              header = T, sep='\t')

HVG_ZT24 <- read.table("Heatmap_individual_seedlings/HVG_ZT0.txt", header=T, sep='\t')

HVG_ZT24_expression <- inner_join(HVG_ZT24, Expression_ZT24,
                                  by=c("Gene"="gene"))

HVG_ZT24_expression <- mutate(HVG_ZT24_expression,
                              Seedling1_meanNorm=TPM_WT_ZT0_seedling_1/rowMeans(HVG_ZT24_expression[,9:22],na.rm=TRUE),
                              Seedling2_meanNorm=TPM_WT_ZT0_seedling_2/rowMeans(HVG_ZT24_expression[,9:22],na.rm=TRUE),
                              Seedling4_meanNorm=TPM_WT_ZT0_seedling_4/rowMeans(HVG_ZT24_expression[,9:22],na.rm=TRUE),
                              Seedling5_meanNorm=TPM_WT_ZT0_seedling_5/rowMeans(HVG_ZT24_expression[,9:22],na.rm=TRUE),
                              Seedling6_meanNorm=TPM_WT_ZT0_seedling_6/rowMeans(HVG_ZT24_expression[,9:22],na.rm=TRUE),
                              Seedling7_meanNorm=TPM_WT_ZT0_seedling_7/rowMeans(HVG_ZT24_expression[,9:22],na.rm=TRUE),
                              Seedling8_meanNorm=TPM_WT_ZT0_seedling_8/rowMeans(HVG_ZT24_expression[,9:22],na.rm=TRUE),
                              Seedling9_meanNorm=TPM_WT_ZT0_seedling_9/rowMeans(HVG_ZT24_expression[,9:22],na.rm=TRUE),
                              Seedling10_meanNorm=TPM_WT_ZT0_seedling_10/rowMeans(HVG_ZT24_expression[,9:22],na.rm=TRUE),
                              Seedling12_meanNorm=TPM_WT_ZT0_seedling_12/rowMeans(HVG_ZT24_expression[,9:22],na.rm=TRUE),
                              Seedling13_meanNorm=TPM_WT_ZT0_seedling_13/rowMeans(HVG_ZT24_expression[,9:22],na.rm=TRUE),
                              Seedling14_meanNorm=TPM_WT_ZT0_seedling_14/rowMeans(HVG_ZT24_expression[,9:22],na.rm=TRUE),
                              Seedling15_meanNorm=TPM_WT_ZT0_seedling_15/rowMeans(HVG_ZT24_expression[,9:22],na.rm=TRUE),
                              Seedling16_meanNorm=TPM_WT_ZT0_seedling_16/rowMeans(HVG_ZT24_expression[,9:22],na.rm=TRUE))



heatmap.2(as.matrix(HVG_ZT24_expression[,23:36]), dendrogram='both', 
          Rowv=TRUE, 
          col= purpleblackyellow,
          breaks=seq(0.8,1.2,length.out=101), trace = 'none' ,
          labRow=NA, labCol=names(HVG_ZT24_expression[,9:22]), 
          main="HVG, ZT24",
          na.rm=F, na.color="black", margins=c(6,2))
























