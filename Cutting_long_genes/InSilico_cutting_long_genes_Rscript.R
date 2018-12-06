
library(tidyverse)
library(M3Drop)
library(RUVSeq)
library(statmod)

### Combine all TPM tables in one file ----

filenames=list.files(path="Cutting_long_genes/data", full.names=TRUE)
datalist = lapply(filenames, function(x){read.table(file=x,header=T,sep='\t')})
TPM_FPKM_Raw_ZT0_cutGene=Reduce(function(x,y) {merge(x,y)}, datalist)

TPM_ZT0_cutGene <- select (TPM_FPKM_Raw_ZT0_cutGene, tracking_id, contains("TPM_"))

write.table(TPM_ZT0_cutGene, file="Cutting_long_genes/TPM_ZT0_cutGene.txt", row.names=FALSE, sep='\t')

### function to identify HVGs and LVGs ----

BrenneckeGetVariableGenesMod <- function(expr_mat, spikes=NA, suppress.plot=FALSE, fdr=0.1, minBiolDisp=0.5, fitMeanQuantile=0.8) {
  #require(statmod)
  rowVars <- function(x) { unlist(apply(x,1,var, na.rm=TRUE))}
  colGenes <- "black"
  colSp <- "blue"
  fullCountTable <- expr_mat;
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
  if (!suppress.plot) {
    plot( meansGenes,cv2Genes, xaxt="n", yaxt="n", log="xy",
          xlab = "average normalized read count",
          ylab = "squared coefficient of variation (CV^2)", col="white")
    axis( 1, 10^(-2:5), c( "0.01", "0.1", "1", "10", "100", "1000",
                           expression(10^4), expression(10^5) ) )
    axis( 2, 10^(-2:3), c( "0.01", "0.1", "1", "10", "100","1000" ), las=2 )
    abline( h=10^(-2:1), v=10^(-1:5), col="#D0D0D0", lwd=2 )
    # Plot the genes, use a different color if they are highly variable
    points( meansGenes, cv2Genes, pch=19, cex=.6,
            col = ifelse( padj < .1, "blue", colGenes ) )
    # Plot/highlight the spike-ins if they are different from the genes
    if (length(meansSp) < length(meansGenes)) {
      points(meansSp, cv2Sp, pch=20, cex=.5, col=colSp)
    }
    # Add the technical noise fit
    xg <- 10^seq( -2, 6, length.out=1000 )
    lines( xg, (a1)/xg + a0, col="red", lwd=3 )
    # Add a curve showing the expectation for the chosen biological CV^2 thershold
    lines( xg, psia1theta/xg + a0 + minBiolDisp, lty="dashed", col="purple", lwd=3)
  }
  TABLE <- data.frame(Gene = names(meansGenes), effect.size=res, p.value = p, q.value= padj, Mean=meansGenes, CV2=cv2Genes , 
                      TechNoise=(a1)/meansGenes+a0, BioVar=cv2Genes-((a1)/meansGenes+a0), Stable=cv2Genes-((a1/10.5)/meansGenes+(a0/2)))
  TABLE <- TABLE[order(-TABLE[,2]),];
  return(list(data=TABLE,A0=a0,A1=a1))
}

### Identify HVG and LVGs ----

### RUV k2

TPM_ZT0 <- read.table("Cutting_long_genes/TPM_ZT0_cutGene.txt",header=T, sep='\t')

ERCC_ZT0 <- filter(TPM_ZT0, grepl("^ERCC" , tracking_id))

Cut_genes <- filter(TPM_ZT0, grepl("_" , tracking_id))

data.info_TPM_ZT0 <- read.table("Cutting_long_genes/DataInfo_ZT0.txt",header=T,sep='\t')
x_ZT0<-as.factor(data.info_TPM_ZT0$ZT)

TPM_ZT0$numZero <- rowSums(TPM_ZT0[2:15]<0.0001)
TPM_ZT0_Max7Zero=filter(TPM_ZT0, numZero<5)

TPM_ZT0_Max7Zero=TPM_ZT0_Max7Zero[,1:15]

TPM_ZT0_Max7Zero_ave5 = filter(TPM_ZT0_Max7Zero, rowMeans(TPM_ZT0_Max7Zero[,2:15])>=5)

Gene_size <- read.table("Extract_HVG_LVG_random/TAIR10_Gene_sizes.txt", header=TRUE, sep='\t')
Gene_size_bis <- Gene_size[,c(1,6)]
Genes_more150bp <- Gene_size_bis[Gene_size_bis$length>150,]

TPM_ZT0_Max7Zero_ave5_more150bp <- merge (TPM_ZT0_Max7Zero_ave5,Genes_more150bp, 
                                          by.x="tracking_id",by.y="gene_name" )

TPM_ZT0_Max7Zero_ave5_more150bp=TPM_ZT0_Max7Zero_ave5_more150bp[,1:15]

TPM_ZT0_Max7Zero_ave5_more150bp=rbind(TPM_ZT0_Max7Zero_ave5_more150bp,ERCC_ZT0)
TPM_ZT0_Max7Zero_ave5_more150bp=rbind(TPM_ZT0_Max7Zero_ave5_more150bp,Cut_genes)

TPM2_ZT0_Max7Zero_ave5_more150bp=TPM_ZT0_Max7Zero_ave5_more150bp[,-1]
rownames(TPM2_ZT0_Max7Zero_ave5_more150bp)=TPM_ZT0_Max7Zero_ave5_more150bp[,1]

TPM2_ZT0_Max7Zero_ave5_more150bp <- round(TPM2_ZT0_Max7Zero_ave5_more150bp)

TPM2_ZT0_Max7Zero_ave5_more150bp=data.matrix(TPM2_ZT0_Max7Zero_ave5_more150bp[,1:14])

genes<-rownames(TPM2_ZT0_Max7Zero_ave5_more150bp)[grep("^AT",rownames(TPM2_ZT0_Max7Zero_ave5_more150bp))]
ERCC<-rownames(TPM2_ZT0_Max7Zero_ave5_more150bp)[grep("^ERCC",rownames(TPM2_ZT0_Max7Zero_ave5_more150bp))]

TPM_ZT0_RUV <- newSeqExpressionSet (as.matrix (TPM2_ZT0_Max7Zero_ave5_more150bp), 
                                    phenoData = data.frame (x_ZT0, row.names= colnames(TPM2_ZT0_Max7Zero_ave5_more150bp)))

TPM2_ZT0 <- RUVg(TPM_ZT0_RUV,ERCC,k=2)
TPM2_ZT0_genes <-TPM2_ZT0[genes]

write.table(normCounts(TPM2_ZT0_genes), file="Cutting_long_genes/TPM2_ZT0_cutGenes.txt",sep='\t')

## Identify stable genes with Brennecke method

HVG_ZT0_genes <- BrenneckeGetVariableGenesMod(normCounts(TPM2_ZT0_genes), minBiolDisp = 0.1, fdr=0.01)

ZT0_a0 <- HVG_ZT0_genes$A0 
ZT0_a1 <- HVG_ZT0_genes$A1 

write.table(HVG_ZT0_genes$data,file="Cutting_long_genes/HVG_BrenneckeMod_ZT0_RUV2_cut_gene.txt",sep='\t')

data_ZT0<- HVG_ZT0_genes$data

ZT0_highNoise <- filter(data_ZT0, q.value<0.01)

write.table(ZT0_highNoise, file="Cutting_long_genes/Variable_BrenneckeMod_ZT0_RUV2.txt",sep='\t')

ZT0_stable <- filter(data_ZT0, Stable<0)

write.table(ZT0_stable, file="Cutting_long_genes/Stable_BrenneckeMod_ZT0_RUV2.txt",sep='\t')

### Plot of Biovar for full genes and fragments ----

Fragments <- read.table("Cutting_long_genes/BioVar_gene_fragments.txt",header=T, sep='\t')
Fragments$Gene <- factor(Fragments$Gene,
                         levels=c("AT2G41480","AT5G47500","AT3G14210",
                                  "AT4G01070","AT5G46690","AT1G25260",
                                  "AT5G64380","AT3G07310","AT4G11150"))
Fragments$Stable_variable <- factor(Fragments$Stable_variable,
                                    levels=c("variable","trend","stable"))

ggplot(Fragments, aes(x=Gene, y=BioVar, colour=type, shape=Stable_variable)) +
  geom_point(size=3)+
  theme(axis.text.x = element_text(angle=90))





