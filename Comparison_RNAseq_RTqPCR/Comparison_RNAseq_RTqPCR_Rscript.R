
library(tidyverse)
library(ggpubr)

# Plotting for direct comparison of the average of the CV2 of each gene by RNA-seq and RTqPCR (removing time-points in RT-qPCR that are also missing in RNA-seq) ----

BioVar_RNAseq <- read.table("Comparison_RNAseq_RTqPCR/BrenneckeMod_allZT_RUV2.txt",
                            header=TRUE, sep='\t')
BioVar_RNAseq_allZT <- select(BioVar_RNAseq,Gene,  contains("CV2"))

BioVar_RNAseq_allZT <- mutate(BioVar_RNAseq_allZT, 
                              ave_CV2_RNAseq = rowMeans(BioVar_RNAseq_allZT[,2:13], 
                                                        na.rm = TRUE))


CV2_RTqPCR_NA <- read.table("Comparison_RNAseq_RTqPCR/CV2_RTqPCR_10genes_NAnoRNAseq.txt",
                            header=TRUE, sep='\t')

CV2_RTqPCR_NA <- mutate(CV2_RTqPCR_NA, ave_CV2_RT = rowMeans(CV2_RTqPCR_NA[,2:13], 
                                                             na.rm = TRUE))

CV2_RTqPCRna_RNAseq <- inner_join(CV2_RTqPCR_NA, BioVar_RNAseq_allZT,
                                  by="Gene")



ggplot(CV2_RTqPCRna_RNAseq, aes(x=ave_CV2_RNAseq, y=ave_CV2_RT)) +
  geom_point(size=3, aes(colour=Gene)) +
  scale_x_log10(name="log10(CV2) for RNA-seq") +
  scale_y_log10(name="log10(CV2) for RT-qPCR") +
  theme(legend.text=element_text(size=20),
        legend.title = element_text(size=20, face="bold"), 
        axis.text=element_text(size=20), 
        axis.title=element_text(size=20),
        plot.title = element_text(size=24,face="bold"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.text.x=element_text(angle=90)) +
  geom_smooth(method='lm', se=FALSE, colour="black") +
  stat_cor(method="spearman", size=5)


# Data import and formating for plotting of RNA-seq and RTqPCR expression profiles ----

RNAseq <- read.table("Comparison_RNAseq_RTqPCR/TPM2_allSeedlings_allZT_genes.txt",
                     header=TRUE, sep='\t')

RTqPCR <- read.table("Comparison_RNAseq_RTqPCR/Expression_RTqPCR_10genes.txt",
                     header=TRUE, sep='\t')

DataInfo <- read.table("Comparison_RNAseq_RTqPCR/DataInfo_RUV_TPM-RTqPCR_10genes.txt",
                       header=TRUE, sep='\t')

DataInfo$ZT <- as.character(DataInfo$ZT)
DataInfo$ZT <- factor(DataInfo$ZT, levels=unique(DataInfo$ZT))

# Plotting non normalised RNA-seq and RTqPCR expression profiles in separated plots ----


AT4G27410_RNAseq <- filter(RNAseq, gene=="AT4G27410")

gather(AT4G27410_RNAseq, "sample", "Expression", TPM_WT_ZT2_seedling_1:TPM_WT_ZT24_seedling_14, 
       factor_key=TRUE) %>%
  right_join(.,DataInfo, by=c("sample"="Sample"), copy=TRUE) %>%
  ggplot(aes(x=ZT, y=Expression)) + 
  ylab("Expression (RNA-seq)") +
  geom_point(size=4) +
  geom_rect(data=NULL, aes(xmin=6.5, xmax=12.6, ymin=-Inf, ymax=Inf), 
            fill="gray90", alpha=0.1) +
  geom_rect(data=NULL, aes(xmin=0.5, xmax=6.5, ymin=-Inf, ymax=Inf), 
            fill="white", alpha=0.1) +
  geom_point(size=4) +
  ggtitle("AT4G27410") +
  labs(color="") +
  theme(legend.text=element_text(size=20),legend.title = element_text(size=20, face="bold"), 
        axis.text=element_text(size=20), axis.title=element_text(size=20),
        plot.title = element_text(size=24,face="bold"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.text.x=element_text(angle=90))



AT4G27410_RTqPCR <- filter(RTqPCR, Sample=="AT4G27410")

gather(AT4G27410_RTqPCR, "sample", "Expression", RT_WT_ZT2_seedling1:RT_WT_ZT24_seedling16, 
       factor_key=TRUE) %>%
  right_join(.,DataInfo, by=c("sample"="Sample"), copy=TRUE) %>%
  ggplot(aes(x=ZT, y=Expression)) + 
  ylab("Expression (RT-qPCR)") +
  geom_point(size=4, col="grey50") +
  geom_rect(data=NULL, aes(xmin=6.5, xmax=12.6, ymin=-Inf, ymax=Inf), 
            fill="gray90", alpha=0.1) +
  geom_rect(data=NULL, aes(xmin=0.5, xmax=6.5, ymin=-Inf, ymax=Inf), 
            fill="white", alpha=0.1) +
  geom_point(size=4,col="grey50") +
  ggtitle("AT4G27410") +
  labs(color="") +
  theme(legend.text=element_text(size=20),legend.title = element_text(size=20, face="bold"), 
        axis.text=element_text(size=20), axis.title=element_text(size=20),
        plot.title = element_text(size=24,face="bold"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.text.x=element_text(angle=90))




