
library(tidyverse)
library(tximport)
library(ensembldb)
library(AnnotationHub)
library(AnnotationDbi)
library(TxDb.Athaliana.BioMart.plantsmart28)
library(M3Drop)
library(RUVSeq)
library(statmod)
library(ggpubr)


### Combining salmon data and summarising transcripts using tximport ----

dir <- "Comparison_tophat_salmon/salmon_raw_data"

samples <- read.table("Comparison_tophat_salmon/data_info.txt", header=T, sep='\t')

files <- file.path(dir,samples$seedling, "quant.sf")
names(files) <- c("seedling1","seedling2","seedling4",
                  "seedling5","seedling6","seedling7",
                  "seedling8","seedling9","seedling10",
                  "seedling12","seedling13","seedling14",
                  "seedling15","seedling16")
all(file.exists(files))

txdb <- TxDb.Athaliana.BioMart.plantsmart28
k <- keys(txdb, keytype = "GENEID")
df <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
tx2gene <- df[, 2:1]

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

TPM_genes_salmon_ZT0_14seedlings <- txi$abundance



### Batch correction and CV2 measurement for salmon data ----

### batch effect correction with RUV (k=2)

ZT0 <- read.table("Comparison_tophat_salmon/TPM_genes_salmon_ZT0_14seedlings_withERCC.txt",
                  header=TRUE, sep='\t')

ERCC_ZT0 <- dplyr::filter(ZT0, grepl("^ERCC" , Name))

data.info_ZT0 <- read.table("Comparison_tophat_salmon/DataInfo_Salmon_ZT0.txt",
                            header=TRUE, sep='\t')
x_ZT0<-as.factor(data.info_ZT0$ZT)

# Selecting genes that are expressed and longer than 150bp
ZT0$numZero <- rowSums(ZT0[2:15]<0.0001)
ZT0_Max7Zero <- dplyr::filter(ZT0, numZero<5)

ZT0_Max7Zero <- ZT0_Max7Zero[,1:15]

ZT0_Max7Zero_ave5 <-  dplyr::filter(ZT0_Max7Zero, rowMeans(ZT0_Max7Zero[,2:15])>=5)


Gene_size <- read.table("Extract_HVG_LVG_random/TAIR10_Gene_sizes.txt", header=TRUE, sep='\t')
Gene_size_bis <- Gene_size[,c(1,6)]
Genes_more150bp <- Gene_size_bis[Gene_size_bis$length>150,]

ZT0_Max7Zero_ave5_more150bp <- merge (ZT0_Max7Zero_ave5,Genes_more150bp, 
                                      by.x="Name",by.y="gene_name" )

ZT0_Max7Zero_ave5_more150bp <- ZT0_Max7Zero_ave5_more150bp[,1:15]

ZT0_Max7Zero_ave5_more150bp <- rbind(ZT0_Max7Zero_ave5_more150bp, ERCC_ZT0)


ZT0_Max7Zero_ave5_more150bp_bis <- ZT0_Max7Zero_ave5_more150bp[,-1]
rownames(ZT0_Max7Zero_ave5_more150bp_bis) <- ZT0_Max7Zero_ave5_more150bp[,1]

ZT0_Max7Zero_ave5_more150bp_bis <- round(ZT0_Max7Zero_ave5_more150bp_bis)

ZT0_Max7Zero_ave5_more150bp_bis <- data.matrix(ZT0_Max7Zero_ave5_more150bp_bis[,1:14])

genes <- rownames(ZT0_Max7Zero_ave5_more150bp_bis)[grep("^AT", rownames(ZT0_Max7Zero_ave5_more150bp_bis))]

ERCC <- rownames(ZT0_Max7Zero_ave5_more150bp_bis)[grep("^ERCC", rownames(ZT0_Max7Zero_ave5_more150bp_bis))]

ZT0_RUV <- newSeqExpressionSet (as.matrix (ZT0_Max7Zero_ave5_more150bp_bis), 
                                phenoData = data.frame (x_ZT0, row.names= colnames(ZT0_Max7Zero_ave5_more150bp_bis)))

RUV2_ZT0 <- RUVg(ZT0_RUV, ERCC, k=2)
RUV2_ZT0_genes <- RUV2_ZT0[genes]

### Data import and formatting for plotting ----

HGV_ZT0_salmon <- read.table("Comparison_tophat_salmon/HVG_BrenneckeMod_ZT0_RUV2_salmon.txt",
                             header=TRUE, sep='\t')

HGV_ZT0_tophat <- read.table("Comparison_tophat_salmon/HVG_BrenneckeMod_ZT0_RUV2.txt",
                             header=TRUE, sep='\t')

HGV_ZT0_salmon_tophat <- inner_join(HGV_ZT0_salmon,HGV_ZT0_tophat,
                                    by="Gene")

### Plotting comparison of tophat and salmon CV2  ----


ggplot(HGV_ZT0_salmon_tophat, aes(x=CV2,y=ZT0_CV2)) +
  geom_point() +
  xlab("CV2 Salmon") +
  ylab("CV2 tophat") +
  xlim(0,1.5) +
  ylim(0,1.5) +
  theme(legend.text=element_text(size=20),
        legend.title = element_text(size=20, face="bold"), 
        axis.text=element_text(size=20), 
        axis.title=element_text(size=20),
        plot.title = element_text(size=24,face="bold"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.text.x=element_text(angle=90)) +
  geom_smooth(method='lm', se=FALSE, colour="black") +
  stat_cor(method="spearman", size=5)

