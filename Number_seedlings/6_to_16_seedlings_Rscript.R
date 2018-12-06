

library(M3Drop)
library(dplyr)
library(RUVSeq)
library(statmod)
library(tidyverse)


# RUV2 and get HVG function  ----

RUV2_HVG <- function(expr_data, dataInfo, gene_size, spikes=NA, fdr=0.1, minBiolDisp=0.5, fitMeanQuantile=0.8) {
  
  ERCC_TPM = filter(expr_data, grepl("^ERCC",tracking_id))
  expr_data$numZero <- rowSums(expr_data[-1]<0.0001)
  expr_data_Max7Zero=filter(expr_data, numZero<5)
  
  expr_data_Max7Zero=expr_data_Max7Zero[,-ncol(expr_data_Max7Zero)]
  
  expr_data_Max7Zero_ave5 = filter(expr_data_Max7Zero, rowMeans(expr_data_Max7Zero[,-1])>=5)
  
  expr_data_Max7Zero_ave5_more150bp <- merge (expr_data_Max7Zero_ave5,gene_size, 
                                              by.x="tracking_id",by.y="gene_name" )
  
  expr_data_Max7Zero_ave5_more150bp=expr_data_Max7Zero_ave5_more150bp[,-ncol(expr_data_Max7Zero_ave5_more150bp)]
  
  expr_data_Max7Zero_ave5_more150bp=rbind(expr_data_Max7Zero_ave5_more150bp,ERCC_TPM)
  
  expr2_data_Max7Zero_ave5_more150bp=expr_data_Max7Zero_ave5_more150bp[,-1]
  rownames(expr2_data_Max7Zero_ave5_more150bp)=expr_data_Max7Zero_ave5_more150bp[,1]
  
  expr2_data_Max7Zero_ave5_more150bp <- round(expr2_data_Max7Zero_ave5_more150bp)
  
  expr2_data_Max7Zero_ave5_more150bp=data.matrix(expr2_data_Max7Zero_ave5_more150bp[,1:ncol(expr2_data_Max7Zero_ave5_more150bp)])
  
  genes<-rownames(expr2_data_Max7Zero_ave5_more150bp)[grep("^AT",rownames(expr2_data_Max7Zero_ave5_more150bp))]
  
  ERCC<-rownames(expr2_data_Max7Zero_ave5_more150bp)[grep("^ERCC",rownames(expr2_data_Max7Zero_ave5_more150bp))]
  
  
  TPM_RUV <- newSeqExpressionSet (as.matrix (expr2_data_Max7Zero_ave5_more150bp), 
                                  phenoData = data.frame (dataInfo, row.names= colnames(expr2_data_Max7Zero_ave5_more150bp)))
  
  TPM2_RUV <- RUVg(TPM_RUV,ERCC,k=2)
  TPM2_RUV_genes <-TPM2_RUV[genes]
  
  #require(statmod)
  rowVars <- function(x) { unlist(apply(x,1,var, na.rm=TRUE))}
  colGenes <- "black"
  colSp <- "blue"
  fullCountTable <- normCounts(TPM2_RUV_genes);
  if (is.character(spikes)) {
    sp <- rownames(fullCountTable) %in% spikes;
    countsSp <- fullCountTable[sp,];
    countsGenes <- fullCountTable[!sp,];
  } else if (is.numeric(spikes)) {
    countsSp <- fullCountTable[spikes,];
    countsGenes <- fullCountTable[-spikes,];
  } else {
    countsSp <- fullCountTable;
    countsGenes <- fullCountTable;
  }
  meansSp <- rowMeans(countsSp, na.rm=TRUE)
  varsSp <- rowVars(countsSp)
  cv2Sp <- varsSp/meansSp^2
  meansGenes <- rowMeans(countsGenes, na.rm=TRUE)
  varsGenes <- rowVars(countsGenes)
  cv2Genes <- varsGenes/meansGenes^2
  # Fit Model
  minMeanForFit <- unname( quantile( meansSp[ which( cv2Sp > 0.3 ) ], fitMeanQuantile))
  useForFit <- meansSp >= minMeanForFit
  if (sum(useForFit, na.rm=TRUE) < 20) {
    warning("Too few spike-ins exceed minMeanForFit, recomputing using all genes.")
    meansAll <- c(meansGenes, meansSp)
    cv2All <- c(cv2Genes,cv2Sp)
    minMeanForFit <- unname( quantile( meansAll[ which( cv2All > 0.3 ) ], 0.80))
    useForFit <- meansSp >= minMeanForFit
  }
  if (sum(useForFit, na.rm=TRUE) < 30) {warning(paste("Only", sum(useForFit), "spike-ins to be used in fitting, may result in poor fit."))}
  fit <- glmgam.fit( cbind( a0 = 5, a1tilde = 0.4/meansSp[useForFit] ), cv2Sp[useForFit] )
  a0 <- unname( fit$coefficients["a0"] )
  a1 <- unname( fit$coefficients["a1tilde"])
  res <- cv2Genes - (a0 + a1/meansGenes);
  # Test
  psia1theta <- a1
  minBiolDisp <- minBiolDisp^2
  m <- ncol(countsSp);
  cv2th <- a0 + minBiolDisp + a0 * minBiolDisp
  testDenom <- (meansGenes*psia1theta + meansGenes^2*cv2th)/(1+cv2th/m)
  p <- 1-pchisq(varsGenes * (m-1)/testDenom,m-1)
  padj <- p.adjust(p,"BH")
  sig <- padj < fdr
  sig[is.na(sig)] <- FALSE
  TABLE <- data.frame(Gene = names(meansGenes), BioVar=log2(cv2Genes/((a1)/meansGenes+a0)))
  TABLE <- TABLE[order(TABLE[,1]),];
  return(data=TABLE)
}



# 16 seedlings ----
TPM <- read.table("Number_seedlings/TPM_selected_more5Mreads_Max50Zero_ZT6.txt",header=TRUE,sep='\t')
dim(TPM)

data.info_TPM=read.table("Number_seedlings/Data_info_ZT6.txt",header=TRUE,sep='\t',row.names=1)
dim(data.info_TPM)
data_info<-as.factor(data.info_TPM$ZT)

Gene_size <- read.table("Number_seedlings/TAIR10_Gene_sizes.txt", header=TRUE, sep='\t')
Gene_size_bis <- Gene_size[,c(1,6)]
Genes_more150bp <- Gene_size_bis[Gene_size_bis$length>150,]
dim(Genes_more150bp)

Correlations_lessSeedlings <- read.table("Number_seedlings/Correlations_15_14_13_12seedlings.txt",header=TRUE,sep='\t')


BioVar_16seedlings <- RUV2_HVG(TPM,data_info, Genes_more150bp, minBiolDisp = 0.1, fdr=0.01)
names(BioVar_16seedlings) =c("Gene" ,"Allseedlings")

# 15 seedlings ----

BioVar_15seedlings <- BioVar_16seedlings

for (i in 2:17) {
  TPM_15seedlings <- TPM[,-i]
  data.info_TPM_15seedlings <- data.info_TPM[-(i-1),]
  data_info<-as.factor(data.info_TPM_15seedlings$ZT)
  BioVar_15seedlings <- RUV2_HVG(TPM_15seedlings,data_info, Genes_more150bp, minBiolDisp = 0.1, fdr=0.01) %>%
    full_join(BioVar_15seedlings,., by="Gene")
}

names(BioVar_15seedlings) = c("Gene","Allseedlings", paste("BioVar_less_seedlings",1:16, sep=""))

write.table(BioVar_15seedlings, file="Number_seedlings/BioVar_15seedlings.txt", 
            row.names = FALSE, sep='\t')

# loop for correlation
Correlations_lessSeedlings <- read.table("Number_seedlings/correlation_lessSeedlings.txt",header=TRUE,sep='\t')


for (i in 3:ncol(BioVar_15seedlings)) {
  Correlations_lessSeedlings[(i-2),2] <- cor(BioVar_15seedlings[,2], BioVar_15seedlings[,i], 
                                             method="pearson", use="pairwise.complete.obs")
}

names(Correlations_lessSeedlings)=c("position","cor_less_1seedling")

write.table(Correlations_lessSeedlings, file="Number_seedlings/Correlations_15seedlings.txt", 
            row.names = FALSE, sep='\t')

# 14 seedlings ----

BioVar_14seedlings <- BioVar_16seedlings


for (i in 2:17) {
  TPM_15seedlings <- TPM[,-i]
  data.info_TPM_15seedlings <- data.info_TPM[-(i-1),]
  data_info<-as.factor(data.info_TPM_15seedlings$ZT)
  for (j in c(2:16)) {
    TPM_14seedlings <- TPM_15seedlings[,-j]
    data.info_TPM_14seedlings <- data.info_TPM_15seedlings[-(j-1),]
    data_info14<-as.factor(data.info_TPM_14seedlings$ZT)
    BioVar_14seedlings <- RUV2_HVG(TPM_14seedlings,data_info14, Genes_more150bp, minBiolDisp = 0.1, fdr=0.01) %>%
      full_join(BioVar_14seedlings,., by="Gene")
  }
}



write.table(BioVar_14seedlings, file="Number_seedlings/BioVar_14seedlings.txt", 
            row.names = FALSE, sep='\t')


# loop for correlation


for (i in 3:ncol(BioVar_14seedlings)) {
  Correlations_lessSeedlings[(i-2),3] <- cor(BioVar_14seedlings[,2], BioVar_14seedlings[,i], 
                                             method="pearson", use="pairwise.complete.obs")
}

names(Correlations_lessSeedlings)=c("position","cor_less_1seedling","cor_less_2seedling")

write.table(Correlations_lessSeedlings, file="Number_seedlings/Correlations_15_14seedlings.txt", 
            row.names = FALSE, sep='\t')


# 13 seedlings ----

BioVar_13seedlings <- BioVar_16seedlings


for (i in 2:17) {
  TPM_15seedlings <- TPM[,-i]
  data.info_TPM_15seedlings <- data.info_TPM[-(i-1),]
  data_info<-as.factor(data.info_TPM_15seedlings$ZT)
  for (j in c(2:16)) {
    TPM_14seedlings <- TPM_15seedlings[,-j]
    data.info_TPM_14seedlings <- data.info_TPM_15seedlings[-(j-1),]
    data_info14<-as.factor(data.info_TPM_14seedlings$ZT)
    for (k in c(2:15)) {
      TPM_13seedlings <- TPM_14seedlings[,-k]
      data.info_TPM_13seedlings <- data.info_TPM_14seedlings[-(k-1),]
      data_info13<-as.factor(data.info_TPM_13seedlings$ZT)
      BioVar_13seedlings <- RUV2_HVG(TPM_13seedlings,data_info13, Genes_more150bp, minBiolDisp = 0.1, fdr=0.01) %>%
        full_join(BioVar_13seedlings,., by="Gene")
    }
  }
}



write.table(BioVar_13seedlings, file="Number_seedlings/BioVar_13seedlings.txt", 
            row.names = FALSE, sep='\t')


# loop for correlation


for (i in 3:ncol(BioVar_13seedlings)) {
  Correlations_lessSeedlings[(i-2),4] <- cor(BioVar_13seedlings[,2], BioVar_13seedlings[,i], 
                                             method="pearson", use="pairwise.complete.obs")
}

names(Correlations_lessSeedlings)=c("position","cor_less_1seedling","cor_less_2seedling","cor_less_3seedling")

write.table(Correlations_lessSeedlings, file="Number_seedlings/Correlations_15_14_13seedlings.txt", 
            row.names = FALSE, sep='\t')




# 12 seedlings ----


BioVar_12seedlings <- BioVar_16seedlings

TPM_15seedlings <- TPM[,-2]
data.info_TPM_15seedlings <- data.info_TPM[-1,]

for (i in 2:16) {
  TPM_14seedlings <- TPM_15seedlings[,-i]
  data.info_TPM_14seedlings <- data.info_TPM_15seedlings[-(i-1),]
  data_info14<-as.factor(data.info_TPM_14seedlings$ZT)
  for (j in c(2:15)) {
    TPM_13seedlings <- TPM_14seedlings[,-j]
    data.info_TPM_13seedlings <- data.info_TPM_14seedlings[-(j-1),]
    data_info13<-as.factor(data.info_TPM_13seedlings$ZT)
    for (k in c(2:14)) {
      TPM_12seedlings <- TPM_13seedlings[,-k]
      data.info_TPM_12seedlings <- data.info_TPM_13seedlings[-(k-1),]
      data_info12<-as.factor(data.info_TPM_12seedlings$ZT)
      BioVar_12seedlings <- RUV2_HVG(TPM_12seedlings,data_info12, Genes_more150bp, minBiolDisp = 0.1, fdr=0.01) %>%
        full_join(BioVar_12seedlings,., by="Gene")
    }
  }
}



write.table(BioVar_12seedlings, file="Number_seedlings/BioVar_12seedlings.txt", 
            row.names = FALSE, sep='\t')


# loop for correlation


for (i in 3:ncol(BioVar_12seedlings)) {
  Correlations_lessSeedlings[(i-2),5] <- cor(BioVar_12seedlings[,2], BioVar_12seedlings[,i], 
                                             method="pearson", use="pairwise.complete.obs")
}

names(Correlations_lessSeedlings)=c("position","cor_less_1seedling","cor_less_2seedling","cor_less_3seedling","cor_less_4seedling")

write.table(Correlations_lessSeedlings, file="Number_seedlings/Correlations_15_14_13_12seedlings.txt", 
            row.names = FALSE, sep='\t')


# 10 seedlings ----


BioVar_10seedlings <- BioVar_16seedlings

TPM_13seedlings <- TPM[,-c(2:4)]
data.info_TPM_13seedlings <- data.info_TPM[-c(1:3),]

for (i in 2:14) {
  TPM_12seedlings <- TPM_13seedlings[,-i]
  data.info_TPM_12seedlings <- data.info_TPM_13seedlings[-(i-1),]
  data_info12<-as.factor(data.info_TPM_12seedlings$ZT)
  for (j in c(2:13)) {
    TPM_11seedlings <- TPM_12seedlings[,-j]
    data.info_TPM_11seedlings <- data.info_TPM_12seedlings[-(j-1),]
    data_info11<-as.factor(data.info_TPM_11seedlings$ZT)
    for (k in c(2:12)) {
      TPM_10seedlings <- TPM_11seedlings[,-k]
      data.info_TPM_10seedlings <- data.info_TPM_11seedlings[-(k-1),]
      data_info10<-as.factor(data.info_TPM_10seedlings$ZT)
      BioVar_10seedlings <- RUV2_HVG(TPM_10seedlings,data_info10, Genes_more150bp, minBiolDisp = 0.1, fdr=0.01) %>%
        full_join(BioVar_10seedlings,., by="Gene")
    }
  }
}


BioVar_10seedlings[BioVar_10seedlings< (-5)] <- NA

write.table(BioVar_10seedlings, file="Number_seedlings/BioVar_10seedlings.txt", 
            row.names = FALSE, sep='\t')


# loop for correlation


for (i in 3:ncol(BioVar_10seedlings)) {
  Correlations_lessSeedlings[(i-2),6] <- cor(BioVar_10seedlings[,2], BioVar_10seedlings[,i], 
                                             method="pearson", use="pairwise.complete.obs")
}

names(Correlations_lessSeedlings)=c("position","cor_less_1seedling","cor_less_2seedling","cor_less_3seedling","cor_less_4seedling","cor_less_6seedling")

write.table(Correlations_lessSeedlings, file="Number_seedlings/Correlations_15_14_13_12_10seedlings.txt", 
            row.names = FALSE, sep='\t')


# 8 seedlings ----


BioVar_8seedlings <- BioVar_16seedlings

TPM_11seedlings <- TPM[,-c(2:4,8,14)]
data.info_TPM_11seedlings <- data.info_TPM[-c(1:3,7,13),]

for (i in 2:12) {
  TPM_10seedlings <- TPM_11seedlings[,-i]
  data.info_TPM_10seedlings <- data.info_TPM_11seedlings[-(i-1),]
  data_info10<-as.factor(data.info_TPM_10seedlings$ZT)
  for (j in c(2:11)) {
    TPM_9seedlings <- TPM_10seedlings[,-j]
    data.info_TPM_9seedlings <- data.info_TPM_10seedlings[-(j-1),]
    data_info9<-as.factor(data.info_TPM_9seedlings$ZT)
    for (k in c(2:10)) {
      TPM_8seedlings <- TPM_9seedlings[,-k]
      data.info_TPM_8seedlings <- data.info_TPM_9seedlings[-(k-1),]
      data_info8<-as.factor(data.info_TPM_8seedlings$ZT)
      BioVar_8seedlings <- RUV2_HVG(TPM_8seedlings,data_info8, Genes_more150bp, minBiolDisp = 0.1, fdr=0.01) %>%
        full_join(BioVar_8seedlings,., by="Gene")
    }
  }
}


BioVar_8seedlings[BioVar_8seedlings< (-5)] <- NA


write.table(BioVar_8seedlings, file="Number_seedlings/BioVar_8seedlings.txt", 
            row.names = FALSE, sep='\t')


# loop for correlation


for (i in 3:ncol(BioVar_8seedlings)) {
  Correlations_lessSeedlings[(i-2),7] <- cor(BioVar_8seedlings[,2], BioVar_8seedlings[,i], 
                                             method="pearson", use="pairwise.complete.obs")
}

names(Correlations_lessSeedlings)=c("position","cor_less_1seedling","cor_less_2seedling","cor_less_3seedling","cor_less_4seedling","cor_less_6seedling","cor_less_8seedling")

write.table(Correlations_lessSeedlings, file="Number_seedlings/Correlations_15_14_13_12_10_8seedlings.txt", 
            row.names = FALSE, sep='\t')


# 6 seedlings ----


BioVar_6seedlings <- BioVar_16seedlings

TPM_10seedlings <- TPM[,-c(2:4,8,14,16)]
data.info_TPM_10seedlings <- data.info_TPM[-c(1:3,7,13,15),]

for (i in 2:11) {
  TPM_9seedlings <- TPM_10seedlings[,-i]
  data.info_TPM_9seedlings <- data.info_TPM_10seedlings[-(i-1),]
  data_info9<-as.factor(data.info_TPM_9seedlings$ZT)
  for (j in c(2:10)) {
    TPM_8seedlings <- TPM_9seedlings[,-j]
    data.info_TPM_8seedlings <- data.info_TPM_9seedlings[-(j-1),]
    data_info8<-as.factor(data.info_TPM_8seedlings$ZT)
    for (j in c(2:9)) {
      TPM_7seedlings <- TPM_8seedlings[,-j]
      data.info_TPM_7seedlings <- data.info_TPM_8seedlings[-(j-1),]
      data_info7<-as.factor(data.info_TPM_7seedlings$ZT)
      for (k in c(2:8)) {
        TPM_6seedlings <- TPM_7seedlings[,-k]
        data.info_TPM_6seedlings <- data.info_TPM_7seedlings[-(k-1),]
        data_info6<-as.factor(data.info_TPM_6seedlings$ZT)
        BioVar_6seedlings <- RUV2_HVG(TPM_6seedlings,data_info6, Genes_more150bp, minBiolDisp = 0.1, fdr=0.01) %>%
          full_join(BioVar_6seedlings,., by="Gene")
      }
    }
  }
}

BioVar_6seedlings[BioVar_6seedlings< (-5)] <- NA


write.table(BioVar_6seedlings, file="Number_seedlings/BioVar_6seedlings.txt", 
            row.names = FALSE, sep='\t')


# loop for correlation


for (i in 3:ncol(BioVar_6seedlings)) {
  Correlations_lessSeedlings[(i-2),8] <- cor(BioVar_6seedlings[,2], BioVar_6seedlings[,i], 
                                             method="pearson", use="pairwise.complete.obs")
}

names(Correlations_lessSeedlings)=c("position","cor_less_1seedling",
                                    "cor_less_2seedling","cor_less_3seedling",
                                    "cor_less_4seedling","cor_less_6seedling",
                                    "cor_less_8seedling","cor_less_10seedling")

write.table(Correlations_lessSeedlings, file="Number_seedlings/Correlations_15_14_13_12_10_8_6seedlings.txt", 
            row.names = FALSE, sep='\t')


# boxplot  ----

Correlations_lessSeedlings <- read.table("Number_seedlings/Correlations_15_14_13_12_10_8_6seedlings.txt",
                                         header=TRUE,sep='\t')

pdf("correlation_6to15_seedlings.pdf", 
    width= 10, height= 7)
boxplot(Correlations_lessSeedlings$cor_less_10seedling,
        Correlations_lessSeedlings$cor_less_8seedling,
        Correlations_lessSeedlings$cor_less_6seedling,
        Correlations_lessSeedlings$cor_less_4seedling,
        Correlations_lessSeedlings$cor_less_3seedling,
        Correlations_lessSeedlings$cor_less_2seedling,
        Correlations_lessSeedlings$cor_less_1seedling,names=c("6 seedlings", "8 seedlings", "10 seedlings",
                                                              "12 seedlings", "13 seedlings",
                                                              "14 seedlings", "15 seedlings"),
        ylab="Pearson correlation",
        cex.lab=1.5, cex.axis=1.3)
dev.off()

gather(Correlations_lessSeedlings, number_seedlings, correlation,
       cor_less_10seedling:cor_less_1seedling, factor_key = TRUE) %>%
  ggplot(.,aes(number_seedlings, correlation)) +
  geom_boxplot() +
  scale_x_discrete(labels=c("6 seedlings","8 seedlings", "10 seedlings",
                            "12 seedlings", "13 seedlings",
                            "14 seedlings", "15 seedlings")) +
  xlab("") +
  ylab("Pearson correlation") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line( size=.1, color="black" ),
        panel.background = element_rect(fill="white",colour="grey50"), 
        text=element_text(size=22),
        axis.text.x = element_text(angle=45, hjust=1)) 

