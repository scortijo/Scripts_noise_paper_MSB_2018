library(tidyverse)
library(M3Drop)
library(dplyr)
library(RUVSeq)
library(statmod)
library(ggpubr)

# Combining all files in raw_data folder, which contains TPM, FPKM and raw reads for each of the 14 seedslings for all 12 time points (168 files in total) ----

filenames <-  list.files(path="Comparison_tophat_hisat2/data", full.names=TRUE)

datalist  <-  lapply(filenames, function(x){read.table(file=x,header=T,sep='\t')})

TPM_FPKM_Raw_14seedling <- Reduce(function(x,y) {merge(x,y)}, datalist)


# Extracting TPM and saving it as TPM_14seedling.txt, which is the basis for all analyses

TPM_ZT0_hisat_14seedling <- select(TPM_FPKM_Raw_14seedling, tracking_id, contains("TPM_"))

write.table(TPM_ZT0_hisat_14seedling, file="Comparison_tophat_hisat2/TPM_ZT0_hisat_14seedling.txt", row.names=FALSE, sep='\t')


# Batch correction and CV2 measurement for hisat2 data ----


### RUV k2

ZT0 <- read.table("Comparison_tophat_hisat2/TPM_ZT0_hisat_14seedling.txt",header=TRUE, sep='\t')

ERCC_ZT0 <- filter(ZT0, grepl("^ERCC" , tracking_id))

data.info_ZT0 <- read.table("Comparison_tophat_hisat2/DataInfo_hisat2_ZT0.txt",header=TRUE, sep='\t')
x_ZT0<-as.factor(data.info_ZT0$ZT)

ZT0$numZero <- rowSums(ZT0[2:15]<0.0001)
ZT0_Max7Zero=filter(ZT0, numZero<5)
ZT0_Max7Zero=ZT0_Max7Zero[,1:15]

ZT0_Max7Zero_ave5 = filter(ZT0_Max7Zero, rowMeans(ZT0_Max7Zero[,2:15])>=5)

Gene_size <- read.table("Extract_HVG_LVG_random/TAIR10_Gene_sizes.txt", header=TRUE, sep='\t')
Gene_size_bis <- Gene_size[,c(1,6)]
Genes_more150bp <- Gene_size_bis[Gene_size_bis$length>150,]

ZT0_Max7Zero_ave5_more150bp <- merge (ZT0_Max7Zero_ave5,Genes_more150bp, 
                                      by.x="tracking_id",by.y="gene_name" )

ZT0_Max7Zero_ave5_more150bp=ZT0_Max7Zero_ave5_more150bp[,1:15]

ZT0_Max7Zero_ave5_more150bp=rbind(ZT0_Max7Zero_ave5_more150bp,ERCC_ZT0)

ZT0_Max7Zero_ave5_more150bp_bis=ZT0_Max7Zero_ave5_more150bp[,-1]
rownames(ZT0_Max7Zero_ave5_more150bp_bis)=ZT0_Max7Zero_ave5_more150bp[,1]

ZT0_Max7Zero_ave5_more150bp_bis <- round(ZT0_Max7Zero_ave5_more150bp_bis)

ZT0_Max7Zero_ave5_more150bp_bis=data.matrix(ZT0_Max7Zero_ave5_more150bp_bis[,1:14])

genes<-rownames(ZT0_Max7Zero_ave5_more150bp_bis)[grep("^AT",rownames(ZT0_Max7Zero_ave5_more150bp_bis))]
ERCC<-rownames(ZT0_Max7Zero_ave5_more150bp_bis)[grep("^ERCC",rownames(ZT0_Max7Zero_ave5_more150bp_bis))]

ZT0_RUV <- newSeqExpressionSet (as.matrix (ZT0_Max7Zero_ave5_more150bp_bis), 
                                phenoData = data.frame (x_ZT0, row.names= colnames(ZT0_Max7Zero_ave5_more150bp_bis)))

RUV2_ZT0 <- RUVg(ZT0_RUV,ERCC,k=2)
RUV2_ZT0_genes <-RUV2_ZT0[genes]

write.table(normCounts(RUV2_ZT0_genes), file="Comparison_tophat_hisat2/RUV2_ZT0_hisat2_genes.txt",sep='\t')


### Data import and formatting for plotting ----

HGV_ZT0_hisat2 <- read.table("Comparison_tophat_hisat2/HVG_BrenneckeMod_ZT0_RUV2_hisat2.txt",
                             header=TRUE, sep='\t')

HGV_ZT0_tophat <- read.table("Comparison_tophat_hisat2/HVG_BrenneckeMod_ZT0_RUV2_tophat.txt",
                             header=TRUE, sep='\t')

HGV_ZT0_hisat2_tophat <- inner_join(HGV_ZT0_hisat2,HGV_ZT0_tophat,
                                    by="Gene")

### Plotting comparison of tophat and hisat2 CV2  ----

ggplot(HGV_ZT0_hisat2_tophat, aes(x=CV2,y=ZT0_CV2)) +
  geom_point() +
  xlab("CV2 hisat2") +
  ylab("CV2 tophat") +
  xlim(0,1.1) +
  ylim(0,1.1) +
  theme(legend.text=element_text(size=20),
        legend.title = element_text(size=20, face="bold"), 
        axis.text=element_text(size=20), 
        axis.title=element_text(size=20),
        plot.title = element_text(size=24,face="bold"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.text.x=element_text(angle=90)) +
  geom_smooth(method='lm', se=FALSE, colour="black") +
  stat_cor(method="spearman", size=5)
