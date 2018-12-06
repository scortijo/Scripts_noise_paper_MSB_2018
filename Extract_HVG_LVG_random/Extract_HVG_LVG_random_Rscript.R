
library(M3Drop)
library(tidyverse)
library(RUVSeq)
library(statmod)


# Get HVG function  ----

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
                      TechNoise=(a1)/meansGenes+a0, BioVar=log2(cv2Genes/((a1)/meansGenes+a0)))
  TABLE <- TABLE[order(-TABLE[,2]),];
  return(list(data=TABLE,A0=a0,A1=a1))
}


# all ZT data preparation   ----

TPM <- read.table("2.combining_data/TPM_14seedling.txt", header=TRUE, sep='\t')

data.info_TPM <- read.table("Extract_HVG_LVG_random/DataInfo_TPM_14seedling.txt",
                         header=TRUE,sep='\t',row.names=1)

# Keeping Spike-ins on a side
ERCC_TPM = filter(TPM, grepl("^ERCC",tracking_id))

# Selecting genes that are expressed
TPM_ave5 <-  filter(TPM, rowMeans(TPM[,2:169])>=5)

# Selecting genes that are longer than 150bp
Gene_size <- read.table("Extract_HVG_LVG_random/TAIR10_Gene_sizes.txt", header=TRUE, sep='\t')
Gene_size_bis <- Gene_size[,c(1,6)]
Genes_more150bp <- Gene_size_bis[Gene_size_bis$length>150,]




# ZT2   ----

### batch effect correction with RUV (k=2)

ERCC_ZT2 <- select (ERCC_TPM,tracking_id,contains("ZT2_"))

TPM_ZT2 <- select (TPM_ave5,tracking_id,contains("ZT2_"))

data.info_TPM_ZT2 <- data.info_TPM[data.info_TPM$ZT=="ZT2",]
x_ZT2<-as.factor(data.info_TPM_ZT2$ZT)

# Selecting genes that are expressed and longer than 150bp
TPM_ZT2$numZero <- rowSums(TPM_ZT2[2:15]<0.0001)
TPM_ZT2_Max7Zero <- filter(TPM_ZT2, numZero<5)

TPM_ZT2_Max7Zero <- TPM_ZT2_Max7Zero[,1:15]

TPM_ZT2_Max7Zero_ave5 <- filter(TPM_ZT2_Max7Zero, rowMeans(TPM_ZT2_Max7Zero[,2:15])>=5)

TPM_ZT2_Max7Zero_ave5_more150bp <- merge (TPM_ZT2_Max7Zero_ave5, Genes_more150bp, 
                                          by.x="tracking_id", by.y="gene_name" )

TPM_ZT2_Max7Zero_ave5_more150bp <- TPM_ZT2_Max7Zero_ave5_more150bp[,1:15]

TPM_ZT2_Max7Zero_ave5_more150bp <- rbind(TPM_ZT2_Max7Zero_ave5_more150bp, ERCC_ZT2)

TPM2_ZT2_Max7Zero_ave5_more150bp <- TPM_ZT2_Max7Zero_ave5_more150bp[,-1]

rownames(TPM2_ZT2_Max7Zero_ave5_more150bp) <- TPM_ZT2_Max7Zero_ave5_more150bp[,1]

TPM2_ZT2_Max7Zero_ave5_more150bp <- round(TPM2_ZT2_Max7Zero_ave5_more150bp)

TPM2_ZT2_Max7Zero_ave5_more150bp <- data.matrix(TPM2_ZT2_Max7Zero_ave5_more150bp[,1:14])

# Organising data for RUV
genes <- rownames(TPM2_ZT2_Max7Zero_ave5_more150bp)[grep("^AT",rownames(TPM2_ZT2_Max7Zero_ave5_more150bp))]

ERCC <- rownames(TPM2_ZT2_Max7Zero_ave5_more150bp)[grep("^ERCC",rownames(TPM2_ZT2_Max7Zero_ave5_more150bp))]


TPM_ZT2_RUV <- newSeqExpressionSet (as.matrix (TPM2_ZT2_Max7Zero_ave5_more150bp), 
                                    phenoData = data.frame (x_ZT2, row.names= colnames(TPM2_ZT2_Max7Zero_ave5_more150bp)))

# Running RUV (K=2) and extracting genes
TPM2_ZT2 <- RUVg(TPM_ZT2_RUV, ERCC, k=2)
TPM2_ZT2_genes <-TPM2_ZT2[genes]

## Identifying highly variable genes

HVG_ZT2_genes <- BrenneckeGetVariableGenesMod(normCounts(TPM2_ZT2_genes), minBiolDisp = 0.1, fdr=0.01)

ZT2_a0 <- HVG_ZT2_genes$A0 
ZT2_a1 <- HVG_ZT2_genes$A1 

data_ZT2 <- HVG_ZT2_genes$data

ZT2_highNoise <- filter(data_ZT2, q.value<0.01)

## Extracting random genes

ZT2_random <- data_ZT2[sample(nrow(data_ZT2),nrow(ZT2_highNoise)),]

## Identifying lowly variable genes: 1000 least variable genes

data_ZT2_orderedBioVar <- arrange(data_ZT2, BioVar)

ZT2_1000LVG <- data_ZT2_orderedBioVar[1:1000,]

## plot

plot( data_ZT2$Mean,data_ZT2$CV2, xaxt="n", yaxt="n", log="xy",
      xlab = "average normalized read count",
      ylab = "", pch=19, cex=.6, cex.lab=1.5)
axis( 1, 10^(-2:5), c( "0.01", "0.1", "1", "10", "100", "1000",
                       expression(10^4), expression(10^5) ) , cex.axis=1.1)
axis( 2, 10^(-2:3), c( "0.01", "0.1", "1", "10", "100","1000" ), las=2 , cex.axis=1.1)
points(ZT2_highNoise$Mean,ZT2_highNoise$CV2, col="blue", pch=19, cex=.6)
points(ZT2_1000LVG$Mean,ZT2_1000LVG$CV2, 
       col="green3", pch=19, cex=.6)
points(ZT2_random$Mean,ZT2_random$CV2, 
       col="grey", pch=19, cex=.6)
xg <- 10^seq( -2, 6, length.out=1000 )
lines( xg, (ZT2_a1)/xg + ZT2_a0, col="red", lwd=2 , lty=1)
mtext(expression(paste("squared coefficient of variation (CV^2)")), side=2, 
      cex=1.5, line=2.5)

# ZT4   ----

### batch effect correction with RUV (k=2)

ERCC_ZT4 <- select (ERCC_TPM,tracking_id,contains("ZT4_"))

TPM_ZT4 <- select (TPM_ave5,tracking_id,contains("ZT4_"))

data.info_TPM_ZT4 <- data.info_TPM[data.info_TPM$ZT=="ZT4",]
x_ZT4<-as.factor(data.info_TPM_ZT4$ZT)

# Selecting genes that are expressed and longer than 150bp
TPM_ZT4$numZero <- rowSums(TPM_ZT4[2:15]<0.0001)
TPM_ZT4_Max7Zero <- filter(TPM_ZT4, numZero<5)

TPM_ZT4_Max7Zero <- TPM_ZT4_Max7Zero[,1:15]

TPM_ZT4_Max7Zero_ave5 <- filter(TPM_ZT4_Max7Zero, rowMeans(TPM_ZT4_Max7Zero[,2:15])>=5)

TPM_ZT4_Max7Zero_ave5_more150bp <- merge (TPM_ZT4_Max7Zero_ave5, Genes_more150bp, 
                                          by.x="tracking_id", by.y="gene_name" )

TPM_ZT4_Max7Zero_ave5_more150bp <- TPM_ZT4_Max7Zero_ave5_more150bp[,1:15]

TPM_ZT4_Max7Zero_ave5_more150bp <- rbind(TPM_ZT4_Max7Zero_ave5_more150bp, ERCC_ZT4)

TPM2_ZT4_Max7Zero_ave5_more150bp <- TPM_ZT4_Max7Zero_ave5_more150bp[,-1]

rownames(TPM2_ZT4_Max7Zero_ave5_more150bp) <- TPM_ZT4_Max7Zero_ave5_more150bp[,1]

TPM2_ZT4_Max7Zero_ave5_more150bp <- round(TPM2_ZT4_Max7Zero_ave5_more150bp)

TPM2_ZT4_Max7Zero_ave5_more150bp <- data.matrix(TPM2_ZT4_Max7Zero_ave5_more150bp[,1:14])

# Organising data for RUV
genes <- rownames(TPM2_ZT4_Max7Zero_ave5_more150bp)[grep("^AT",rownames(TPM2_ZT4_Max7Zero_ave5_more150bp))]

ERCC <- rownames(TPM2_ZT4_Max7Zero_ave5_more150bp)[grep("^ERCC",rownames(TPM2_ZT4_Max7Zero_ave5_more150bp))]


TPM_ZT4_RUV <- newSeqExpressionSet (as.matrix (TPM2_ZT4_Max7Zero_ave5_more150bp), 
                                    phenoData = data.frame (x_ZT4, row.names= colnames(TPM2_ZT4_Max7Zero_ave5_more150bp)))

# Running RUV (K=2) and extracting genes
TPM2_ZT4 <- RUVg(TPM_ZT4_RUV, ERCC, k=2)
TPM2_ZT4_genes <-TPM2_ZT4[genes]

## Identifying highly variable genes

HVG_ZT4_genes <- BrenneckeGetVariableGenesMod(normCounts(TPM2_ZT4_genes), minBiolDisp = 0.1, fdr=0.01)

ZT4_a0 <- HVG_ZT4_genes$A0 
ZT4_a1 <- HVG_ZT4_genes$A1 

data_ZT4 <- HVG_ZT4_genes$data

ZT4_highNoise <- filter(data_ZT4, q.value<0.01)

## Extracting random genes

ZT4_random <- data_ZT4[sample(nrow(data_ZT4),nrow(ZT4_highNoise)),]

## Identifying lowly variable genes: 1000 least variable genes

data_ZT4_orderedBioVar <- arrange(data_ZT4, BioVar)

ZT4_1000LVG <- data_ZT4_orderedBioVar[1:1000,]

## plot

plot( data_ZT4$Mean,data_ZT4$CV2, xaxt="n", yaxt="n", log="xy",
      xlab = "average normalized read count",
      ylab = "", pch=19, cex=.6, cex.lab=1.5)
axis( 1, 10^(-2:5), c( "0.01", "0.1", "1", "10", "100", "1000",
                       expression(10^4), expression(10^5) ) , cex.axis=1.1)
axis( 2, 10^(-2:3), c( "0.01", "0.1", "1", "10", "100","1000" ), las=2 , cex.axis=1.1)
points(ZT4_highNoise$Mean,ZT4_highNoise$CV2, col="blue", pch=19, cex=.6)
points(ZT4_1000LVG$Mean,ZT4_1000LVG$CV2, 
       col="green3", pch=19, cex=.6)
points(ZT4_random$Mean,ZT4_random$CV2, 
       col="grey", pch=19, cex=.6)
xg <- 10^seq( -2, 6, length.out=1000 )
lines( xg, (ZT4_a1)/xg + ZT4_a0, col="red", lwd=2 , lty=1)
mtext(expression(paste("squared coefficient of variation (CV^2)")), side=2, 
      cex=1.5, line=2.5)


# ZT6   ----

### batch effect correction with RUV (k=2)

ERCC_ZT6 <- select (ERCC_TPM,tracking_id,contains("ZT6_"))

TPM_ZT6 <- select (TPM_ave5,tracking_id,contains("ZT6_"))

data.info_TPM_ZT6 <- data.info_TPM[data.info_TPM$ZT=="ZT6",]
x_ZT6<-as.factor(data.info_TPM_ZT6$ZT)

# Selecting genes that are expressed and longer than 150bp
TPM_ZT6$numZero <- rowSums(TPM_ZT6[2:15]<0.0001)
TPM_ZT6_Max7Zero <- filter(TPM_ZT6, numZero<5)

TPM_ZT6_Max7Zero <- TPM_ZT6_Max7Zero[,1:15]

TPM_ZT6_Max7Zero_ave5 <- filter(TPM_ZT6_Max7Zero, rowMeans(TPM_ZT6_Max7Zero[,2:15])>=5)

TPM_ZT6_Max7Zero_ave5_more150bp <- merge (TPM_ZT6_Max7Zero_ave5, Genes_more150bp, 
                                          by.x="tracking_id", by.y="gene_name" )

TPM_ZT6_Max7Zero_ave5_more150bp <- TPM_ZT6_Max7Zero_ave5_more150bp[,1:15]

TPM_ZT6_Max7Zero_ave5_more150bp <- rbind(TPM_ZT6_Max7Zero_ave5_more150bp, ERCC_ZT6)

TPM2_ZT6_Max7Zero_ave5_more150bp <- TPM_ZT6_Max7Zero_ave5_more150bp[,-1]

rownames(TPM2_ZT6_Max7Zero_ave5_more150bp) <- TPM_ZT6_Max7Zero_ave5_more150bp[,1]

TPM2_ZT6_Max7Zero_ave5_more150bp <- round(TPM2_ZT6_Max7Zero_ave5_more150bp)

TPM2_ZT6_Max7Zero_ave5_more150bp <- data.matrix(TPM2_ZT6_Max7Zero_ave5_more150bp[,1:14])

# Organising data for RUV
genes <- rownames(TPM2_ZT6_Max7Zero_ave5_more150bp)[grep("^AT",rownames(TPM2_ZT6_Max7Zero_ave5_more150bp))]

ERCC <- rownames(TPM2_ZT6_Max7Zero_ave5_more150bp)[grep("^ERCC",rownames(TPM2_ZT6_Max7Zero_ave5_more150bp))]


TPM_ZT6_RUV <- newSeqExpressionSet (as.matrix (TPM2_ZT6_Max7Zero_ave5_more150bp), 
                                    phenoData = data.frame (x_ZT6, row.names= colnames(TPM2_ZT6_Max7Zero_ave5_more150bp)))

# Running RUV (K=2) and extracting genes
TPM2_ZT6 <- RUVg(TPM_ZT6_RUV, ERCC, k=2)
TPM2_ZT6_genes <-TPM2_ZT6[genes]

## Identifying highly variable genes

HVG_ZT6_genes <- BrenneckeGetVariableGenesMod(normCounts(TPM2_ZT6_genes), minBiolDisp = 0.1, fdr=0.01)

ZT6_a0 <- HVG_ZT6_genes$A0 
ZT6_a1 <- HVG_ZT6_genes$A1 

data_ZT6 <- HVG_ZT6_genes$data

ZT6_highNoise <- filter(data_ZT6, q.value<0.01)

## Extracting random genes

ZT6_random <- data_ZT6[sample(nrow(data_ZT6),nrow(ZT6_highNoise)),]

## Identifying lowly variable genes: 1000 least variable genes

data_ZT6_orderedBioVar <- arrange(data_ZT6, BioVar)

ZT6_1000LVG <- data_ZT6_orderedBioVar[1:1000,]

## plot

plot( data_ZT6$Mean,data_ZT6$CV2, xaxt="n", yaxt="n", log="xy",
      xlab = "average normalized read count",
      ylab = "", pch=19, cex=.6, cex.lab=1.5)
axis( 1, 10^(-2:5), c( "0.01", "0.1", "1", "10", "100", "1000",
                       expression(10^4), expression(10^5) ) , cex.axis=1.1)
axis( 2, 10^(-2:3), c( "0.01", "0.1", "1", "10", "100","1000" ), las=2 , cex.axis=1.1)
points(ZT6_highNoise$Mean,ZT6_highNoise$CV2, col="blue", pch=19, cex=.6)
points(ZT6_1000LVG$Mean,ZT6_1000LVG$CV2, 
       col="green3", pch=19, cex=.6)
points(ZT6_random$Mean,ZT6_random$CV2, 
       col="grey", pch=19, cex=.6)
xg <- 10^seq( -2, 6, length.out=1000 )
lines( xg, (ZT6_a1)/xg + ZT6_a0, col="red", lwd=2 , lty=1)
mtext(expression(paste("squared coefficient of variation (CV^2)")), side=2, 
      cex=1.5, line=2.5)



# ZT8   ----

### batch effect correction with RUV (k=2)

ERCC_ZT8 <- select (ERCC_TPM,tracking_id,contains("ZT8_"))

TPM_ZT8 <- select (TPM_ave5,tracking_id,contains("ZT8_"))

data.info_TPM_ZT8 <- data.info_TPM[data.info_TPM$ZT=="ZT8",]
x_ZT8<-as.factor(data.info_TPM_ZT8$ZT)

# Selecting genes that are expressed and longer than 150bp
TPM_ZT8$numZero <- rowSums(TPM_ZT8[2:15]<0.0001)
TPM_ZT8_Max7Zero <- filter(TPM_ZT8, numZero<5)

TPM_ZT8_Max7Zero <- TPM_ZT8_Max7Zero[,1:15]

TPM_ZT8_Max7Zero_ave5 <- filter(TPM_ZT8_Max7Zero, rowMeans(TPM_ZT8_Max7Zero[,2:15])>=5)

TPM_ZT8_Max7Zero_ave5_more150bp <- merge (TPM_ZT8_Max7Zero_ave5, Genes_more150bp, 
                                          by.x="tracking_id", by.y="gene_name" )

TPM_ZT8_Max7Zero_ave5_more150bp <- TPM_ZT8_Max7Zero_ave5_more150bp[,1:15]

TPM_ZT8_Max7Zero_ave5_more150bp <- rbind(TPM_ZT8_Max7Zero_ave5_more150bp, ERCC_ZT8)

TPM2_ZT8_Max7Zero_ave5_more150bp <- TPM_ZT8_Max7Zero_ave5_more150bp[,-1]

rownames(TPM2_ZT8_Max7Zero_ave5_more150bp) <- TPM_ZT8_Max7Zero_ave5_more150bp[,1]

TPM2_ZT8_Max7Zero_ave5_more150bp <- round(TPM2_ZT8_Max7Zero_ave5_more150bp)

TPM2_ZT8_Max7Zero_ave5_more150bp <- data.matrix(TPM2_ZT8_Max7Zero_ave5_more150bp[,1:14])

# Organising data for RUV
genes <- rownames(TPM2_ZT8_Max7Zero_ave5_more150bp)[grep("^AT",rownames(TPM2_ZT8_Max7Zero_ave5_more150bp))]

ERCC <- rownames(TPM2_ZT8_Max7Zero_ave5_more150bp)[grep("^ERCC",rownames(TPM2_ZT8_Max7Zero_ave5_more150bp))]


TPM_ZT8_RUV <- newSeqExpressionSet (as.matrix (TPM2_ZT8_Max7Zero_ave5_more150bp), 
                                    phenoData = data.frame (x_ZT8, row.names= colnames(TPM2_ZT8_Max7Zero_ave5_more150bp)))

# Running RUV (K=2) and extracting genes
TPM2_ZT8 <- RUVg(TPM_ZT8_RUV, ERCC, k=2)
TPM2_ZT8_genes <-TPM2_ZT8[genes]

## Identifying highly variable genes

HVG_ZT8_genes <- BrenneckeGetVariableGenesMod(normCounts(TPM2_ZT8_genes), minBiolDisp = 0.1, fdr=0.01)

ZT8_a0 <- HVG_ZT8_genes$A0 
ZT8_a1 <- HVG_ZT8_genes$A1 

data_ZT8 <- HVG_ZT8_genes$data

ZT8_highNoise <- filter(data_ZT8, q.value<0.01)

## Extracting random genes

ZT8_random <- data_ZT8[sample(nrow(data_ZT8),nrow(ZT8_highNoise)),]

## Identifying lowly variable genes: 1000 least variable genes

data_ZT8_orderedBioVar <- arrange(data_ZT8, BioVar)

ZT8_1000LVG <- data_ZT8_orderedBioVar[1:1000,]

## plot

plot( data_ZT8$Mean,data_ZT8$CV2, xaxt="n", yaxt="n", log="xy",
      xlab = "average normalized read count",
      ylab = "", pch=19, cex=.6, cex.lab=1.5)
axis( 1, 10^(-2:5), c( "0.01", "0.1", "1", "10", "100", "1000",
                       expression(10^4), expression(10^5) ) , cex.axis=1.1)
axis( 2, 10^(-2:3), c( "0.01", "0.1", "1", "10", "100","1000" ), las=2 , cex.axis=1.1)
points(ZT8_highNoise$Mean,ZT8_highNoise$CV2, col="blue", pch=19, cex=.6)
points(ZT8_1000LVG$Mean,ZT8_1000LVG$CV2, 
       col="green3", pch=19, cex=.6)
points(ZT8_random$Mean,ZT8_random$CV2, 
       col="grey", pch=19, cex=.6)
xg <- 10^seq( -2, 6, length.out=1000 )
lines( xg, (ZT8_a1)/xg + ZT8_a0, col="red", lwd=2 , lty=1)
mtext(expression(paste("squared coefficient of variation (CV^2)")), side=2, 
      cex=1.5, line=2.5)

# ZT10   ----

### batch effect correction with RUV (k=2)

ERCC_ZT10 <- select (ERCC_TPM,tracking_id,contains("ZT10_"))

TPM_ZT10 <- select (TPM_ave5,tracking_id,contains("ZT10_"))

data.info_TPM_ZT10 <- data.info_TPM[data.info_TPM$ZT=="ZT10",]
x_ZT10<-as.factor(data.info_TPM_ZT10$ZT)

# Selecting genes that are expressed and longer than 150bp
TPM_ZT10$numZero <- rowSums(TPM_ZT10[2:15]<0.0001)
TPM_ZT10_Max7Zero <- filter(TPM_ZT10, numZero<5)

TPM_ZT10_Max7Zero <- TPM_ZT10_Max7Zero[,1:15]

TPM_ZT10_Max7Zero_ave5 <- filter(TPM_ZT10_Max7Zero, rowMeans(TPM_ZT10_Max7Zero[,2:15])>=5)

TPM_ZT10_Max7Zero_ave5_more150bp <- merge (TPM_ZT10_Max7Zero_ave5, Genes_more150bp, 
                                          by.x="tracking_id", by.y="gene_name" )

TPM_ZT10_Max7Zero_ave5_more150bp <- TPM_ZT10_Max7Zero_ave5_more150bp[,1:15]

TPM_ZT10_Max7Zero_ave5_more150bp <- rbind(TPM_ZT10_Max7Zero_ave5_more150bp, ERCC_ZT10)

TPM2_ZT10_Max7Zero_ave5_more150bp <- TPM_ZT10_Max7Zero_ave5_more150bp[,-1]

rownames(TPM2_ZT10_Max7Zero_ave5_more150bp) <- TPM_ZT10_Max7Zero_ave5_more150bp[,1]

TPM2_ZT10_Max7Zero_ave5_more150bp <- round(TPM2_ZT10_Max7Zero_ave5_more150bp)

TPM2_ZT10_Max7Zero_ave5_more150bp <- data.matrix(TPM2_ZT10_Max7Zero_ave5_more150bp[,1:14])

# Organising data for RUV
genes <- rownames(TPM2_ZT10_Max7Zero_ave5_more150bp)[grep("^AT",rownames(TPM2_ZT10_Max7Zero_ave5_more150bp))]

ERCC <- rownames(TPM2_ZT10_Max7Zero_ave5_more150bp)[grep("^ERCC",rownames(TPM2_ZT10_Max7Zero_ave5_more150bp))]


TPM_ZT10_RUV <- newSeqExpressionSet (as.matrix (TPM2_ZT10_Max7Zero_ave5_more150bp), 
                                    phenoData = data.frame (x_ZT10, row.names= colnames(TPM2_ZT10_Max7Zero_ave5_more150bp)))

# Running RUV (K=2) and extracting genes
TPM2_ZT10 <- RUVg(TPM_ZT10_RUV, ERCC, k=2)
TPM2_ZT10_genes <-TPM2_ZT10[genes]

## Identifying highly variable genes

HVG_ZT10_genes <- BrenneckeGetVariableGenesMod(normCounts(TPM2_ZT10_genes), minBiolDisp = 0.1, fdr=0.01)

ZT10_a0 <- HVG_ZT10_genes$A0 
ZT10_a1 <- HVG_ZT10_genes$A1 

data_ZT10 <- HVG_ZT10_genes$data

ZT10_highNoise <- filter(data_ZT10, q.value<0.01)

## Extracting random genes

ZT10_random <- data_ZT10[sample(nrow(data_ZT10),nrow(ZT10_highNoise)),]

## Identifying lowly variable genes: 1000 least variable genes

data_ZT10_orderedBioVar <- arrange(data_ZT10, BioVar)

ZT10_1000LVG <- data_ZT10_orderedBioVar[1:1000,]

## plot

plot( data_ZT10$Mean,data_ZT10$CV2, xaxt="n", yaxt="n", log="xy",
      xlab = "average normalized read count",
      ylab = "", pch=19, cex=.6, cex.lab=1.5)
axis( 1, 10^(-2:5), c( "0.01", "0.1", "1", "10", "100", "1000",
                       expression(10^4), expression(10^5) ) , cex.axis=1.1)
axis( 2, 10^(-2:3), c( "0.01", "0.1", "1", "10", "100","1000" ), las=2 , cex.axis=1.1)
points(ZT10_highNoise$Mean,ZT10_highNoise$CV2, col="blue", pch=19, cex=.6)
points(ZT10_1000LVG$Mean,ZT10_1000LVG$CV2, 
       col="green3", pch=19, cex=.6)
points(ZT10_random$Mean,ZT10_random$CV2, 
       col="grey", pch=19, cex=.6)
xg <- 10^seq( -2, 6, length.out=1000 )
lines( xg, (ZT10_a1)/xg + ZT10_a0, col="red", lwd=2 , lty=1)
mtext(expression(paste("squared coefficient of variation (CV^2)")), side=2, 
      cex=1.5, line=2.5)


# ZT12   ----

### batch effect correction with RUV (k=2)

ERCC_ZT12 <- select (ERCC_TPM,tracking_id,contains("ZT12_"))

TPM_ZT12 <- select (TPM_ave5,tracking_id,contains("ZT12_"))

data.info_TPM_ZT12 <- data.info_TPM[data.info_TPM$ZT=="ZT12",]
x_ZT12<-as.factor(data.info_TPM_ZT12$ZT)

# Selecting genes that are expressed and longer than 150bp
TPM_ZT12$numZero <- rowSums(TPM_ZT12[2:15]<0.0001)
TPM_ZT12_Max7Zero <- filter(TPM_ZT12, numZero<5)

TPM_ZT12_Max7Zero <- TPM_ZT12_Max7Zero[,1:15]

TPM_ZT12_Max7Zero_ave5 <- filter(TPM_ZT12_Max7Zero, rowMeans(TPM_ZT12_Max7Zero[,2:15])>=5)

TPM_ZT12_Max7Zero_ave5_more150bp <- merge (TPM_ZT12_Max7Zero_ave5, Genes_more150bp, 
                                           by.x="tracking_id", by.y="gene_name" )

TPM_ZT12_Max7Zero_ave5_more150bp <- TPM_ZT12_Max7Zero_ave5_more150bp[,1:15]

TPM_ZT12_Max7Zero_ave5_more150bp <- rbind(TPM_ZT12_Max7Zero_ave5_more150bp, ERCC_ZT12)

TPM2_ZT12_Max7Zero_ave5_more150bp <- TPM_ZT12_Max7Zero_ave5_more150bp[,-1]

rownames(TPM2_ZT12_Max7Zero_ave5_more150bp) <- TPM_ZT12_Max7Zero_ave5_more150bp[,1]

TPM2_ZT12_Max7Zero_ave5_more150bp <- round(TPM2_ZT12_Max7Zero_ave5_more150bp)

TPM2_ZT12_Max7Zero_ave5_more150bp <- data.matrix(TPM2_ZT12_Max7Zero_ave5_more150bp[,1:14])

# Organising data for RUV
genes <- rownames(TPM2_ZT12_Max7Zero_ave5_more150bp)[grep("^AT",rownames(TPM2_ZT12_Max7Zero_ave5_more150bp))]

ERCC <- rownames(TPM2_ZT12_Max7Zero_ave5_more150bp)[grep("^ERCC",rownames(TPM2_ZT12_Max7Zero_ave5_more150bp))]


TPM_ZT12_RUV <- newSeqExpressionSet (as.matrix (TPM2_ZT12_Max7Zero_ave5_more150bp), 
                                     phenoData = data.frame (x_ZT12, row.names= colnames(TPM2_ZT12_Max7Zero_ave5_more150bp)))

# Running RUV (K=2) and extracting genes
TPM2_ZT12 <- RUVg(TPM_ZT12_RUV, ERCC, k=2)
TPM2_ZT12_genes <-TPM2_ZT12[genes]

## Identifying highly variable genes

HVG_ZT12_genes <- BrenneckeGetVariableGenesMod(normCounts(TPM2_ZT12_genes), minBiolDisp = 0.1, fdr=0.01)

ZT12_a0 <- HVG_ZT12_genes$A0 
ZT12_a1 <- HVG_ZT12_genes$A1 

data_ZT12 <- HVG_ZT12_genes$data

ZT12_highNoise <- filter(data_ZT12, q.value<0.01)

## Extracting random genes

ZT12_random <- data_ZT12[sample(nrow(data_ZT12),nrow(ZT12_highNoise)),]

## Identifying lowly variable genes: 1000 least variable genes

data_ZT12_orderedBioVar <- arrange(data_ZT12, BioVar)

ZT12_1000LVG <- data_ZT12_orderedBioVar[1:1000,]

## plot

plot( data_ZT12$Mean,data_ZT12$CV2, xaxt="n", yaxt="n", log="xy",
      xlab = "average normalized read count",
      ylab = "", pch=19, cex=.6, cex.lab=1.5)
axis( 1, 10^(-2:5), c( "0.01", "0.1", "1", "10", "100", "1000",
                       expression(10^4), expression(10^5) ) , cex.axis=1.1)
axis( 2, 10^(-2:3), c( "0.01", "0.1", "1", "10", "100","1000" ), las=2 , cex.axis=1.1)
points(ZT12_highNoise$Mean,ZT12_highNoise$CV2, col="blue", pch=19, cex=.6)
points(ZT12_1000LVG$Mean,ZT12_1000LVG$CV2, 
       col="green3", pch=19, cex=.6)
points(ZT12_random$Mean,ZT12_random$CV2, 
       col="grey", pch=19, cex=.6)
xg <- 10^seq( -2, 6, length.out=1000 )
lines( xg, (ZT12_a1)/xg + ZT12_a0, col="red", lwd=2 , lty=1)
mtext(expression(paste("squared coefficient of variation (CV^2)")), side=2, 
      cex=1.5, line=2.5)

# ZT14   ----

### batch effect correction with RUV (k=2)

ERCC_ZT14 <- select (ERCC_TPM,tracking_id,contains("ZT14_"))

TPM_ZT14 <- select (TPM_ave5,tracking_id,contains("ZT14_"))

data.info_TPM_ZT14 <- data.info_TPM[data.info_TPM$ZT=="ZT14",]
x_ZT14<-as.factor(data.info_TPM_ZT14$ZT)

# Selecting genes that are expressed and longer than 150bp
TPM_ZT14$numZero <- rowSums(TPM_ZT14[2:15]<0.0001)
TPM_ZT14_Max7Zero <- filter(TPM_ZT14, numZero<5)

TPM_ZT14_Max7Zero <- TPM_ZT14_Max7Zero[,1:15]

TPM_ZT14_Max7Zero_ave5 <- filter(TPM_ZT14_Max7Zero, rowMeans(TPM_ZT14_Max7Zero[,2:15])>=5)

TPM_ZT14_Max7Zero_ave5_more150bp <- merge (TPM_ZT14_Max7Zero_ave5, Genes_more150bp, 
                                           by.x="tracking_id", by.y="gene_name" )

TPM_ZT14_Max7Zero_ave5_more150bp <- TPM_ZT14_Max7Zero_ave5_more150bp[,1:15]

TPM_ZT14_Max7Zero_ave5_more150bp <- rbind(TPM_ZT14_Max7Zero_ave5_more150bp, ERCC_ZT14)

TPM2_ZT14_Max7Zero_ave5_more150bp <- TPM_ZT14_Max7Zero_ave5_more150bp[,-1]

rownames(TPM2_ZT14_Max7Zero_ave5_more150bp) <- TPM_ZT14_Max7Zero_ave5_more150bp[,1]

TPM2_ZT14_Max7Zero_ave5_more150bp <- round(TPM2_ZT14_Max7Zero_ave5_more150bp)

TPM2_ZT14_Max7Zero_ave5_more150bp <- data.matrix(TPM2_ZT14_Max7Zero_ave5_more150bp[,1:14])

# Organising data for RUV
genes <- rownames(TPM2_ZT14_Max7Zero_ave5_more150bp)[grep("^AT",rownames(TPM2_ZT14_Max7Zero_ave5_more150bp))]

ERCC <- rownames(TPM2_ZT14_Max7Zero_ave5_more150bp)[grep("^ERCC",rownames(TPM2_ZT14_Max7Zero_ave5_more150bp))]


TPM_ZT14_RUV <- newSeqExpressionSet (as.matrix (TPM2_ZT14_Max7Zero_ave5_more150bp), 
                                     phenoData = data.frame (x_ZT14, row.names= colnames(TPM2_ZT14_Max7Zero_ave5_more150bp)))

# Running RUV (K=2) and extracting genes
TPM2_ZT14 <- RUVg(TPM_ZT14_RUV, ERCC, k=2)
TPM2_ZT14_genes <-TPM2_ZT14[genes]

## Identifying highly variable genes

HVG_ZT14_genes <- BrenneckeGetVariableGenesMod(normCounts(TPM2_ZT14_genes), minBiolDisp = 0.1, fdr=0.01)

ZT14_a0 <- HVG_ZT14_genes$A0 
ZT14_a1 <- HVG_ZT14_genes$A1 

data_ZT14 <- HVG_ZT14_genes$data

ZT14_highNoise <- filter(data_ZT14, q.value<0.01)

## Extracting random genes

ZT14_random <- data_ZT14[sample(nrow(data_ZT14),nrow(ZT14_highNoise)),]

## Identifying lowly variable genes: 1000 least variable genes

data_ZT14_orderedBioVar <- arrange(data_ZT14, BioVar)

ZT14_1000LVG <- data_ZT14_orderedBioVar[1:1000,]

## plot

plot( data_ZT14$Mean,data_ZT14$CV2, xaxt="n", yaxt="n", log="xy",
      xlab = "average normalized read count",
      ylab = "", pch=19, cex=.6, cex.lab=1.5)
axis( 1, 10^(-2:5), c( "0.01", "0.1", "1", "10", "100", "1000",
                       expression(10^4), expression(10^5) ) , cex.axis=1.1)
axis( 2, 10^(-2:3), c( "0.01", "0.1", "1", "10", "100","1000" ), las=2 , cex.axis=1.1)
points(ZT14_highNoise$Mean,ZT14_highNoise$CV2, col="blue", pch=19, cex=.6)
points(ZT14_1000LVG$Mean,ZT14_1000LVG$CV2, 
       col="green3", pch=19, cex=.6)
points(ZT14_random$Mean,ZT14_random$CV2, 
       col="grey", pch=19, cex=.6)
xg <- 10^seq( -2, 6, length.out=1000 )
lines( xg, (ZT14_a1)/xg + ZT14_a0, col="red", lwd=2 , lty=1)
mtext(expression(paste("squared coefficient of variation (CV^2)")), side=2, 
      cex=1.5, line=2.5)



# ZT16   ----

### batch effect correction with RUV (k=2)

ERCC_ZT16 <- select (ERCC_TPM,tracking_id,contains("ZT16_"))

TPM_ZT16 <- select (TPM_ave5,tracking_id,contains("ZT16_"))

data.info_TPM_ZT16 <- data.info_TPM[data.info_TPM$ZT=="ZT16",]
x_ZT16<-as.factor(data.info_TPM_ZT16$ZT)

# Selecting genes that are expressed and longer than 150bp
TPM_ZT16$numZero <- rowSums(TPM_ZT16[2:15]<0.0001)
TPM_ZT16_Max7Zero <- filter(TPM_ZT16, numZero<5)

TPM_ZT16_Max7Zero <- TPM_ZT16_Max7Zero[,1:15]

TPM_ZT16_Max7Zero_ave5 <- filter(TPM_ZT16_Max7Zero, rowMeans(TPM_ZT16_Max7Zero[,2:15])>=5)

TPM_ZT16_Max7Zero_ave5_more150bp <- merge (TPM_ZT16_Max7Zero_ave5, Genes_more150bp, 
                                           by.x="tracking_id", by.y="gene_name" )

TPM_ZT16_Max7Zero_ave5_more150bp <- TPM_ZT16_Max7Zero_ave5_more150bp[,1:15]

TPM_ZT16_Max7Zero_ave5_more150bp <- rbind(TPM_ZT16_Max7Zero_ave5_more150bp, ERCC_ZT16)

TPM2_ZT16_Max7Zero_ave5_more150bp <- TPM_ZT16_Max7Zero_ave5_more150bp[,-1]

rownames(TPM2_ZT16_Max7Zero_ave5_more150bp) <- TPM_ZT16_Max7Zero_ave5_more150bp[,1]

TPM2_ZT16_Max7Zero_ave5_more150bp <- round(TPM2_ZT16_Max7Zero_ave5_more150bp)

TPM2_ZT16_Max7Zero_ave5_more150bp <- data.matrix(TPM2_ZT16_Max7Zero_ave5_more150bp[,1:14])

# Organising data for RUV
genes <- rownames(TPM2_ZT16_Max7Zero_ave5_more150bp)[grep("^AT",rownames(TPM2_ZT16_Max7Zero_ave5_more150bp))]

ERCC <- rownames(TPM2_ZT16_Max7Zero_ave5_more150bp)[grep("^ERCC",rownames(TPM2_ZT16_Max7Zero_ave5_more150bp))]


TPM_ZT16_RUV <- newSeqExpressionSet (as.matrix (TPM2_ZT16_Max7Zero_ave5_more150bp), 
                                     phenoData = data.frame (x_ZT16, row.names= colnames(TPM2_ZT16_Max7Zero_ave5_more150bp)))

# Running RUV (K=2) and extracting genes
TPM2_ZT16 <- RUVg(TPM_ZT16_RUV, ERCC, k=2)
TPM2_ZT16_genes <-TPM2_ZT16[genes]

## Identifying highly variable genes

HVG_ZT16_genes <- BrenneckeGetVariableGenesMod(normCounts(TPM2_ZT16_genes), minBiolDisp = 0.1, fdr=0.01)

ZT16_a0 <- HVG_ZT16_genes$A0 
ZT16_a1 <- HVG_ZT16_genes$A1 

data_ZT16 <- HVG_ZT16_genes$data

ZT16_highNoise <- filter(data_ZT16, q.value<0.01)

## Extracting random genes

ZT16_random <- data_ZT16[sample(nrow(data_ZT16),nrow(ZT16_highNoise)),]

## Identifying lowly variable genes: 1000 least variable genes

data_ZT16_orderedBioVar <- arrange(data_ZT16, BioVar)

ZT16_1000LVG <- data_ZT16_orderedBioVar[1:1000,]

## plot

plot( data_ZT16$Mean,data_ZT16$CV2, xaxt="n", yaxt="n", log="xy",
      xlab = "average normalized read count",
      ylab = "", pch=19, cex=.6, cex.lab=1.5)
axis( 1, 10^(-2:5), c( "0.01", "0.1", "1", "10", "100", "1000",
                       expression(10^4), expression(10^5) ) , cex.axis=1.1)
axis( 2, 10^(-2:3), c( "0.01", "0.1", "1", "10", "100","1000" ), las=2 , cex.axis=1.1)
points(ZT16_highNoise$Mean,ZT16_highNoise$CV2, col="blue", pch=19, cex=.6)
points(ZT16_1000LVG$Mean,ZT16_1000LVG$CV2, 
       col="green3", pch=19, cex=.6)
points(ZT16_random$Mean,ZT16_random$CV2, 
       col="grey", pch=19, cex=.6)
xg <- 10^seq( -2, 6, length.out=1000 )
lines( xg, (ZT16_a1)/xg + ZT16_a0, col="red", lwd=2 , lty=1)
mtext(expression(paste("squared coefficient of variation (CV^2)")), side=2, 
      cex=1.5, line=2.5)


# ZT18   ----

### batch effect correction with RUV (k=2)

ERCC_ZT18 <- select (ERCC_TPM,tracking_id,contains("ZT18_"))

TPM_ZT18 <- select (TPM_ave5,tracking_id,contains("ZT18_"))

data.info_TPM_ZT18 <- data.info_TPM[data.info_TPM$ZT=="ZT18",]
x_ZT18 <- as.factor(data.info_TPM_ZT18$ZT)

# Selecting genes that are expressed and longer than 150bp
TPM_ZT18$numZero <- rowSums(TPM_ZT18[2:15]<0.0001)
TPM_ZT18_Max7Zero <- filter(TPM_ZT18, numZero<5)

TPM_ZT18_Max7Zero <- TPM_ZT18_Max7Zero[,1:15]

TPM_ZT18_Max7Zero_ave5 <- filter(TPM_ZT18_Max7Zero, rowMeans(TPM_ZT18_Max7Zero[,2:15])>=5)

TPM_ZT18_Max7Zero_ave5_more150bp <- merge (TPM_ZT18_Max7Zero_ave5, Genes_more150bp, 
                                           by.x="tracking_id", by.y="gene_name" )

TPM_ZT18_Max7Zero_ave5_more150bp <- TPM_ZT18_Max7Zero_ave5_more150bp[,1:15]

TPM_ZT18_Max7Zero_ave5_more150bp <- rbind(TPM_ZT18_Max7Zero_ave5_more150bp, ERCC_ZT18)

TPM2_ZT18_Max7Zero_ave5_more150bp <- TPM_ZT18_Max7Zero_ave5_more150bp[,-1]

rownames(TPM2_ZT18_Max7Zero_ave5_more150bp) <- TPM_ZT18_Max7Zero_ave5_more150bp[,1]

TPM2_ZT18_Max7Zero_ave5_more150bp <- round(TPM2_ZT18_Max7Zero_ave5_more150bp)

TPM2_ZT18_Max7Zero_ave5_more150bp <- data.matrix(TPM2_ZT18_Max7Zero_ave5_more150bp[,1:14])

# Organising data for RUV
genes <- rownames(TPM2_ZT18_Max7Zero_ave5_more150bp)[grep("^AT",rownames(TPM2_ZT18_Max7Zero_ave5_more150bp))]

ERCC <- rownames(TPM2_ZT18_Max7Zero_ave5_more150bp)[grep("^ERCC",rownames(TPM2_ZT18_Max7Zero_ave5_more150bp))]


TPM_ZT18_RUV <- newSeqExpressionSet (as.matrix (TPM2_ZT18_Max7Zero_ave5_more150bp), 
                                     phenoData = data.frame (x_ZT18, row.names= colnames(TPM2_ZT18_Max7Zero_ave5_more150bp)))

# Running RUV (K=2) and extracting genes
TPM2_ZT18 <- RUVg(TPM_ZT18_RUV, ERCC, k=2)
TPM2_ZT18_genes <-TPM2_ZT18[genes]

## Identifying highly variable genes

HVG_ZT18_genes <- BrenneckeGetVariableGenesMod(normCounts(TPM2_ZT18_genes), minBiolDisp = 0.1, fdr=0.01)

ZT18_a0 <- HVG_ZT18_genes$A0 
ZT18_a1 <- HVG_ZT18_genes$A1 

data_ZT18 <- HVG_ZT18_genes$data

ZT18_highNoise <- filter(data_ZT18, q.value<0.01)

## Extracting random genes

ZT18_random <- data_ZT18[sample(nrow(data_ZT18),nrow(ZT18_highNoise)),]

## Identifying lowly variable genes: 1000 least variable genes

data_ZT18_orderedBioVar <- arrange(data_ZT18, BioVar)

ZT18_1000LVG <- data_ZT18_orderedBioVar[1:1000,]

## plot

plot( data_ZT18$Mean,data_ZT18$CV2, xaxt="n", yaxt="n", log="xy",
      xlab = "average normalized read count",
      ylab = "", pch=19, cex=.6, cex.lab=1.5)
axis( 1, 10^(-2:5), c( "0.01", "0.1", "1", "10", "100", "1000",
                       expression(10^4), expression(10^5) ) , cex.axis=1.1)
axis( 2, 10^(-2:3), c( "0.01", "0.1", "1", "10", "100","1000" ), las=2 , cex.axis=1.1)
points(ZT18_highNoise$Mean,ZT18_highNoise$CV2, col="blue", pch=19, cex=.6)
points(ZT18_1000LVG$Mean,ZT18_1000LVG$CV2, 
       col="green3", pch=19, cex=.6)
points(ZT18_random$Mean,ZT18_random$CV2, 
       col="grey", pch=19, cex=.6)
xg <- 10^seq( -2, 6, length.out=1000 )
lines( xg, (ZT18_a1)/xg + ZT18_a0, col="red", lwd=2 , lty=1)
mtext(expression(paste("squared coefficient of variation (CV^2)")), side=2, 
      cex=1.5, line=2.5)


# ZT20   ----

### batch effect correction with RUV (k=2)

ERCC_ZT20 <- select (ERCC_TPM,tracking_id,contains("ZT20_"))

TPM_ZT20 <- select (TPM_ave5,tracking_id,contains("ZT20_"))

data.info_TPM_ZT20 <- data.info_TPM[data.info_TPM$ZT=="ZT20",]
x_ZT20 <- as.factor(data.info_TPM_ZT20$ZT)

# Selecting genes that are expressed and longer than 150bp
TPM_ZT20$numZero <- rowSums(TPM_ZT20[2:15]<0.0001)
TPM_ZT20_Max7Zero <- filter(TPM_ZT20, numZero<5)

TPM_ZT20_Max7Zero <- TPM_ZT20_Max7Zero[,1:15]

TPM_ZT20_Max7Zero_ave5 <- filter(TPM_ZT20_Max7Zero, rowMeans(TPM_ZT20_Max7Zero[,2:15])>=5)

TPM_ZT20_Max7Zero_ave5_more150bp <- merge (TPM_ZT20_Max7Zero_ave5, Genes_more150bp, 
                                           by.x="tracking_id", by.y="gene_name" )

TPM_ZT20_Max7Zero_ave5_more150bp <- TPM_ZT20_Max7Zero_ave5_more150bp[,1:15]

TPM_ZT20_Max7Zero_ave5_more150bp <- rbind(TPM_ZT20_Max7Zero_ave5_more150bp, ERCC_ZT20)

TPM2_ZT20_Max7Zero_ave5_more150bp <- TPM_ZT20_Max7Zero_ave5_more150bp[,-1]

rownames(TPM2_ZT20_Max7Zero_ave5_more150bp) <- TPM_ZT20_Max7Zero_ave5_more150bp[,1]

TPM2_ZT20_Max7Zero_ave5_more150bp <- round(TPM2_ZT20_Max7Zero_ave5_more150bp)

TPM2_ZT20_Max7Zero_ave5_more150bp <- data.matrix(TPM2_ZT20_Max7Zero_ave5_more150bp[,1:14])

# Organising data for RUV
genes <- rownames(TPM2_ZT20_Max7Zero_ave5_more150bp)[grep("^AT",rownames(TPM2_ZT20_Max7Zero_ave5_more150bp))]

ERCC <- rownames(TPM2_ZT20_Max7Zero_ave5_more150bp)[grep("^ERCC",rownames(TPM2_ZT20_Max7Zero_ave5_more150bp))]


TPM_ZT20_RUV <- newSeqExpressionSet (as.matrix (TPM2_ZT20_Max7Zero_ave5_more150bp), 
                                     phenoData = data.frame (x_ZT20, row.names= colnames(TPM2_ZT20_Max7Zero_ave5_more150bp)))

# Running RUV (K=2) and extracting genes
TPM2_ZT20 <- RUVg(TPM_ZT20_RUV, ERCC, k=2)
TPM2_ZT20_genes <-TPM2_ZT20[genes]

## Identifying highly variable genes

HVG_ZT20_genes <- BrenneckeGetVariableGenesMod(normCounts(TPM2_ZT20_genes), minBiolDisp = 0.1, fdr=0.01)

ZT20_a0 <- HVG_ZT20_genes$A0 
ZT20_a1 <- HVG_ZT20_genes$A1 

data_ZT20 <- HVG_ZT20_genes$data

ZT20_highNoise <- filter(data_ZT20, q.value<0.01)

## Extracting random genes

ZT20_random <- data_ZT20[sample(nrow(data_ZT20),nrow(ZT20_highNoise)),]

## Identifying lowly variable genes: 1000 least variable genes

data_ZT20_orderedBioVar <- arrange(data_ZT20, BioVar)

ZT20_1000LVG <- data_ZT20_orderedBioVar[1:1000,]

## plot

plot( data_ZT20$Mean,data_ZT20$CV2, xaxt="n", yaxt="n", log="xy",
      xlab = "average normalized read count",
      ylab = "", pch=19, cex=.6, cex.lab=1.5)
axis( 1, 10^(-2:5), c( "0.01", "0.1", "1", "10", "100", "1000",
                       expression(10^4), expression(10^5) ) , cex.axis=1.1)
axis( 2, 10^(-2:3), c( "0.01", "0.1", "1", "10", "100","1000" ), las=2 , cex.axis=1.1)
points(ZT20_highNoise$Mean,ZT20_highNoise$CV2, col="blue", pch=19, cex=.6)
points(ZT20_1000LVG$Mean,ZT20_1000LVG$CV2, 
       col="green3", pch=19, cex=.6)
points(ZT20_random$Mean,ZT20_random$CV2, 
       col="grey", pch=19, cex=.6)
xg <- 10^seq( -2, 6, length.out=1000 )
lines( xg, (ZT20_a1)/xg + ZT20_a0, col="red", lwd=2 , lty=1)
mtext(expression(paste("squared coefficient of variation (CV^2)")), side=2, 
      cex=1.5, line=2.5)


# ZT22   ----

### batch effect correction with RUV (k=2)

ERCC_ZT22 <- select (ERCC_TPM,tracking_id,contains("ZT22_"))

TPM_ZT22 <- select (TPM_ave5,tracking_id,contains("ZT22_"))

data.info_TPM_ZT22 <- data.info_TPM[data.info_TPM$ZT=="ZT22",]
x_ZT22 <- as.factor(data.info_TPM_ZT22$ZT)

# Selecting genes that are expressed and longer than 150bp
TPM_ZT22$numZero <- rowSums(TPM_ZT22[2:15]<0.0001)
TPM_ZT22_Max7Zero <- filter(TPM_ZT22, numZero<5)

TPM_ZT22_Max7Zero <- TPM_ZT22_Max7Zero[,1:15]

TPM_ZT22_Max7Zero_ave5 <- filter(TPM_ZT22_Max7Zero, rowMeans(TPM_ZT22_Max7Zero[,2:15])>=5)

TPM_ZT22_Max7Zero_ave5_more150bp <- merge (TPM_ZT22_Max7Zero_ave5, Genes_more150bp, 
                                           by.x="tracking_id", by.y="gene_name" )

TPM_ZT22_Max7Zero_ave5_more150bp <- TPM_ZT22_Max7Zero_ave5_more150bp[,1:15]

TPM_ZT22_Max7Zero_ave5_more150bp <- rbind(TPM_ZT22_Max7Zero_ave5_more150bp, ERCC_ZT22)

TPM2_ZT22_Max7Zero_ave5_more150bp <- TPM_ZT22_Max7Zero_ave5_more150bp[,-1]

rownames(TPM2_ZT22_Max7Zero_ave5_more150bp) <- TPM_ZT22_Max7Zero_ave5_more150bp[,1]

TPM2_ZT22_Max7Zero_ave5_more150bp <- round(TPM2_ZT22_Max7Zero_ave5_more150bp)

TPM2_ZT22_Max7Zero_ave5_more150bp <- data.matrix(TPM2_ZT22_Max7Zero_ave5_more150bp[,1:14])

# Organising data for RUV
genes <- rownames(TPM2_ZT22_Max7Zero_ave5_more150bp)[grep("^AT",rownames(TPM2_ZT22_Max7Zero_ave5_more150bp))]

ERCC <- rownames(TPM2_ZT22_Max7Zero_ave5_more150bp)[grep("^ERCC",rownames(TPM2_ZT22_Max7Zero_ave5_more150bp))]


TPM_ZT22_RUV <- newSeqExpressionSet (as.matrix (TPM2_ZT22_Max7Zero_ave5_more150bp), 
                                     phenoData = data.frame (x_ZT22, row.names= colnames(TPM2_ZT22_Max7Zero_ave5_more150bp)))

# Running RUV (K=2) and extracting genes
TPM2_ZT22 <- RUVg(TPM_ZT22_RUV, ERCC, k=2)
TPM2_ZT22_genes <-TPM2_ZT22[genes]

## Identifying highly variable genes

HVG_ZT22_genes <- BrenneckeGetVariableGenesMod(normCounts(TPM2_ZT22_genes), minBiolDisp = 0.1, fdr=0.01)

ZT22_a0 <- HVG_ZT22_genes$A0 
ZT22_a1 <- HVG_ZT22_genes$A1 

data_ZT22 <- HVG_ZT22_genes$data

ZT22_highNoise <- filter(data_ZT22, q.value<0.01)

## Extracting random genes

ZT22_random <- data_ZT22[sample(nrow(data_ZT22),nrow(ZT22_highNoise)),]

## Identifying lowly variable genes: 1000 least variable genes

data_ZT22_orderedBioVar <- arrange(data_ZT22, BioVar)

ZT22_1000LVG <- data_ZT22_orderedBioVar[1:1000,]

## plot

plot( data_ZT22$Mean,data_ZT22$CV2, xaxt="n", yaxt="n", log="xy",
      xlab = "average normalized read count",
      ylab = "", pch=19, cex=.6, cex.lab=1.5)
axis( 1, 10^(-2:5), c( "0.01", "0.1", "1", "10", "100", "1000",
                       expression(10^4), expression(10^5) ) , cex.axis=1.1)
axis( 2, 10^(-2:3), c( "0.01", "0.1", "1", "10", "100","1000" ), las=2 , cex.axis=1.1)
points(ZT22_highNoise$Mean,ZT22_highNoise$CV2, col="blue", pch=19, cex=.6)
points(ZT22_1000LVG$Mean,ZT22_1000LVG$CV2, 
       col="green3", pch=19, cex=.6)
points(ZT22_random$Mean,ZT22_random$CV2, 
       col="grey", pch=19, cex=.6)
xg <- 10^seq( -2, 6, length.out=1000 )
lines( xg, (ZT22_a1)/xg + ZT22_a0, col="red", lwd=2 , lty=1)
mtext(expression(paste("squared coefficient of variation (CV^2)")), side=2, 
      cex=1.5, line=2.5)


# ZT24   ----

### batch effect correction with RUV (k=2)

ERCC_ZT24 <- select (ERCC_TPM,tracking_id,contains("ZT24_"))

TPM_ZT24 <- select (TPM_ave5,tracking_id,contains("ZT24_"))

data.info_TPM_ZT24 <- data.info_TPM[data.info_TPM$ZT=="ZT24",]
x_ZT24<-as.factor(data.info_TPM_ZT24$ZT)

# Selecting genes that are expressed and longer than 150bp
TPM_ZT24$numZero <- rowSums(TPM_ZT24[2:15]<0.0001)
TPM_ZT24_Max7Zero <- filter(TPM_ZT24, numZero<5)

TPM_ZT24_Max7Zero <- TPM_ZT24_Max7Zero[,1:15]

TPM_ZT24_Max7Zero_ave5 <- filter(TPM_ZT24_Max7Zero, rowMeans(TPM_ZT24_Max7Zero[,2:15])>=5)

TPM_ZT24_Max7Zero_ave5_more150bp <- merge (TPM_ZT24_Max7Zero_ave5, Genes_more150bp, 
                                          by.x="tracking_id", by.y="gene_name" )

TPM_ZT24_Max7Zero_ave5_more150bp <- TPM_ZT24_Max7Zero_ave5_more150bp[,1:15]

TPM_ZT24_Max7Zero_ave5_more150bp <- rbind(TPM_ZT24_Max7Zero_ave5_more150bp, ERCC_ZT24)

TPM2_ZT24_Max7Zero_ave5_more150bp <- TPM_ZT24_Max7Zero_ave5_more150bp[,-1]

rownames(TPM2_ZT24_Max7Zero_ave5_more150bp) <- TPM_ZT24_Max7Zero_ave5_more150bp[,1]

TPM2_ZT24_Max7Zero_ave5_more150bp <- round(TPM2_ZT24_Max7Zero_ave5_more150bp)

TPM2_ZT24_Max7Zero_ave5_more150bp <- data.matrix(TPM2_ZT24_Max7Zero_ave5_more150bp[,1:14])

# Organising data for RUV
genes <- rownames(TPM2_ZT24_Max7Zero_ave5_more150bp)[grep("^AT",rownames(TPM2_ZT24_Max7Zero_ave5_more150bp))]

ERCC <- rownames(TPM2_ZT24_Max7Zero_ave5_more150bp)[grep("^ERCC",rownames(TPM2_ZT24_Max7Zero_ave5_more150bp))]


TPM_ZT24_RUV <- newSeqExpressionSet (as.matrix (TPM2_ZT24_Max7Zero_ave5_more150bp), 
                                    phenoData = data.frame (x_ZT24, row.names= colnames(TPM2_ZT24_Max7Zero_ave5_more150bp)))

# Running RUV (K=2) and extracting genes
TPM2_ZT24 <- RUVg(TPM_ZT24_RUV, ERCC, k=2)
TPM2_ZT24_genes <-TPM2_ZT24[genes]

## Identifying highly variable genes

HVG_ZT24_genes <- BrenneckeGetVariableGenesMod(normCounts(TPM2_ZT24_genes), minBiolDisp = 0.1, fdr=0.01)

ZT24_a0 <- HVG_ZT24_genes$A0 
ZT24_a1 <- HVG_ZT24_genes$A1 

data_ZT24 <- HVG_ZT24_genes$data

ZT24_highNoise <- filter(data_ZT24, q.value<0.01)

## Extracting random genes

ZT24_random <- data_ZT24[sample(nrow(data_ZT24),nrow(ZT24_highNoise)),]

## Identifying lowly variable genes: 1000 least variable genes

data_ZT24_orderedBioVar <- arrange(data_ZT24, BioVar)

ZT24_1000LVG <- data_ZT24_orderedBioVar[1:1000,]

## plot

plot( data_ZT24$Mean,data_ZT24$CV2, xaxt="n", yaxt="n", log="xy",
      xlab = "average normalized read count",
      ylab = "", pch=19, cex=.6, cex.lab=1.5)
axis( 1, 10^(-2:5), c( "0.01", "0.1", "1", "10", "100", "1000",
                       expression(10^4), expression(10^5) ) , cex.axis=1.1)
axis( 2, 10^(-2:3), c( "0.01", "0.1", "1", "10", "100","1000" ), las=2 , cex.axis=1.1)
points(ZT24_highNoise$Mean,ZT24_highNoise$CV2, col="blue", pch=19, cex=.6)
points(ZT24_1000LVG$Mean,ZT24_1000LVG$CV2, 
       col="green3", pch=19, cex=.6)
points(ZT24_random$Mean,ZT24_random$CV2, 
       col="grey", pch=19, cex=.6)
xg <- 10^seq( -2, 6, length.out=1000 )
lines( xg, (ZT24_a1)/xg + ZT24_a0, col="red", lwd=2 , lty=1)
mtext(expression(paste("squared coefficient of variation (CV^2)")), side=2, 
      cex=1.5, line=2.5)


# Global trends ----

par(mar=c(5,6,2,2))
xg <- 10^seq( -2, 6, length.out=1000 )
plot( xg, (ZT24_a1)/xg + ZT24_a0, col="#a6cee3", lwd=2 , lty=1,type = "l",log="xy",xlim=c(5,10000),ylim=c(0.001,2),xaxt="n", yaxt="n", xlab = "Average normalized read count",
      ylab = "",cex.lab=1.5)
mtext(expression(CV^2),2,cex=1.9, line=4)
axis( 1, 10^(0:5), c( "1", "10", "100", "1000",
                      expression(10^4), expression(10^5) ) , cex.axis=1.8)
axis( 2, 10^(-3:1), c("0.001", "0.01", "0.1", "1", "10" ), las=2 , cex.axis=1.8)
lines( xg, (ZT2_a1)/xg + ZT2_a0, col="#1f78b4", lwd=2 , lty=1)
lines( xg, (ZT4_a1)/xg + ZT4_a0, col="#b2df8a", lwd=2 , lty=1)
lines( xg, (ZT6_a1)/xg + ZT6_a0, col="#33a02c", lwd=2 , lty=1)
lines( xg, (ZT8_a1)/xg + ZT8_a0, col="#fb9a99", lwd=2 , lty=1)
lines( xg, (ZT10_a1)/xg + ZT10_a0, col="#e31a1c", lwd=2 , lty=1)
lines( xg, (ZT12_a1)/xg + ZT12_a0, col="#fdbf6f", lwd=2 , lty=1)
lines( xg, (ZT14_a1)/xg + ZT14_a0, col="#ff7f00", lwd=2 , lty=1)
lines( xg, (ZT16_a1)/xg + ZT16_a0, col="#cab2d6", lwd=2, lty=1 )
lines( xg, (ZT18_a1)/xg + ZT18_a0, col="#6a3d9a", lwd=2 , lty=1)
lines( xg, (ZT20_a1)/xg + ZT20_a0, col="#ffff99", lwd=2 , lty=1)
lines( xg, (ZT22_a1)/xg + ZT22_a0, col="#b15928", lwd=2 , lty=1)
legend("topright",legend=c("ZT2","ZT4","ZT6","ZT8","ZT10","ZT12","ZT14","ZT16","ZT18","ZT20","ZT22", "ZT24"), 
       col=c("#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f",
             "#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928", "#a6cee3"), lty=1, cex=1.2,lwd=3)



