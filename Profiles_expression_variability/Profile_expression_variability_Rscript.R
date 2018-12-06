
library(tidyverse)
library(tibble)

# Data import and formating for expression data ----

TPM_allZT <- read.table("Profiles_expression_variability/TPM2_allSeedlings_allZT_genes.txt", 
                        header=TRUE, sep='\t')


TPM_Info <-read.table("Profiles_expression_variability/TPM2_allSeedlings_AllZT_DataInfo.txt", 
                      header=TRUE, sep='\t')
TPM_Info$ZT <- as.character(TPM_Info$ZT)
TPM_Info$ZT <- factor(TPM_Info$ZT, levels=unique(TPM_Info$ZT))

# Plotting one gene for expression data ----

AT5G52310  <- TPM_allZT[TPM_allZT$gene=="AT5G52310",]

gather(AT5G52310, "sample", "Expression", c(TPM_WT_ZT2_seedling_1:TPM_WT_ZT24_seedling_9), 
       factor_key=TRUE) %>%
  right_join(.,TPM_Info, by=c("sample"="sample"), copy=TRUE) %>%
  ggplot(aes(x=ZT, y=Expression)) + 
  ylab("TPM") +
  geom_point( size=4) +
  geom_rect(data=NULL, aes(xmin=6.5, xmax=12.6, ymin=-Inf, ymax=Inf), 
            fill="gray90", alpha=0.1) +
  geom_rect(data=NULL, aes(xmin=0.5, xmax=6.5, ymin=-Inf, ymax=Inf), 
            fill="white", alpha=0.1) +
  geom_point( size=4) +
  ggtitle("AT5G52310") +
  labs(color="")+
  theme(legend.text=element_text(size=20),legend.title = element_text(size=20, face="bold"), 
        axis.text=element_text(size=20), axis.title=element_text(size=20),
        plot.title = element_text(size=24,face="bold"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.text.x=element_text(angle=90))

# Data import and formating for variability data ----

BioVar_allZT<- read.table("Profiles_expression_variability/CorCV2_allZT_allGenes.txt",header=T,sep='\t')

DataInfoBioVar<- read.table("Profiles_expression_variability/CorCV2_allZT_allGenes_dataInfo.txt",head=TRUE,sep='\t')
DataInfoBioVar$ZT <- as.character(DataInfoBioVar$ZT)
DataInfoBioVar$ZT <- factor(DataInfoBioVar$ZT, levels=unique(DataInfoBioVar$ZT))

# Plotting one gene for variability data ----

BioVar <- BioVar_allZT[BioVar_allZT$Gene=="AT5G47560",]

gather(BioVar, "sample", "BioVar", BioVar_ZT2:BioVar_ZT24, factor_key=TRUE) %>%
  right_join(.,DataInfoBioVar, by=c("sample"="sample"), copy=TRUE) %>%
  ggplot(aes(x=ZT, y=BioVar, group=Gene)) + 
  geom_line() +
  geom_rect(data=NULL, aes(xmin=6.5, xmax=12.7, ymin=-Inf, ymax=Inf), 
            fill="lightgrey", alpha=0.15) +
  geom_rect(data=NULL, aes(xmin=0, xmax=6.5, ymin=-Inf, ymax=Inf), 
            fill="white", alpha=0.1) +
  geom_line(size=2) +
  ylab("Normalised Variability") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line( size=.3, color="black" ), 
        panel.background = element_rect(fill="white",colour="grey50"), 
        text=element_text(size=18)) +
  ggtitle("AT5G47560")

