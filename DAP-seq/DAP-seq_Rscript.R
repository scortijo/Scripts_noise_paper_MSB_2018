
library(tidyverse)

#### Plotting of the number of TFs with enriched targets in each family ----

DAPseq_counts <- read.table("DAP-seq/DAP-seq_family_nbTF_HVG_LVG_random_background.txt",
                            header=T, sep='\t')

ggplot(DAPseq_counts, aes(x=total.number, y=HVG, colour=status)) +
  geom_point(size=3) +
  geom_text(aes(label=family_short),hjust=1, vjust=-0.5, size=6) +
  #geom_smooth(method='lm',formula=y~x)  +
  xlab("DAP-seq dataset") +
  ylim(0,15) +
  xlim(-2,53) +
  scale_color_manual(values=c("red","green3","black")) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line( size=.1, color="black" ),
        panel.background = element_rect(fill="white",colour="grey50"), 
        text=element_text(size=22)) 


#### Plotting of the number of TFs targetting each gene for HVGs, LVGs and a thousand set of random genes ----

#data preparation


DAPseq <- read.table("DAP-seq/DAPseq_summary.txt",header=T, sep='\t')

DAPseq_nbTF <- dplyr::mutate(DAPseq, nbTF=rowSums(DAPseq[,-1])) %>%
  select(., gene_name, nbTF)

Background <- read.table("DAP-seq/All_genes_background.txt", header=TRUE, sep='\t')

Background_gene_nbTF <- inner_join(Background, DAPseq_nbTF, by=c("Gene"="gene_name")) %>%
  select(Gene, nbTF)


# extract 1000 sets of random genes ----

random_distribution_nbTF_1000=matrix(nrow=512,ncol=1000)
random_distribution_nbTF_1000=data.frame(random_distribution_nbTF_1000)

random_nbTF_1000=matrix(nrow=1358,ncol=1000)
random_nbTF_1000=data.frame(random_nbTF_1000)


for (i in 1:1000) {
  random_set <- Background_gene_nbTF[sample(nrow(Background_gene_nbTF), 1358),]
  random_nbTF_1000[,i] <- random_set[,2]
  d <- density(random_set[,2], from=0, to=253)
  random_distribution_nbTF_1000[,i] <- d$y
}

colnames(random_nbTF_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(Background_gene_nbTF$nbTF, from=0, to=253)
random_distribution_nbTF_1000 <- cbind(d_background$x, random_distribution_nbTF_1000)
colnames(random_distribution_nbTF_1000) <- c("nbTF",paste("random_set",1:100, sep="_"))


HVG <- read.table("DAP-seq/HVG_allZT_background.txt", header=TRUE, sep='\t')

Variable_nbTF <- inner_join(HVG, DAPseq_nbTF, by=c("allZT"="gene_name")) %>%
  select(allZT, nbTF) 

LVG <- read.table("DAP-seq/LVG1000_allZT_background.txt", header=TRUE, sep='\t')

Stable_nbTF <- inner_join(LVG, DAPseq_nbTF, by=c("allZT"="gene_name")) %>%
  select(allZT, nbTF) 


HVG_LVG_distribution_nbTF=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_nbTF=data.frame(HVG_LVG_distribution_nbTF)

d_background <-  density(Background_gene_nbTF$nbTF, from=0, to=253)
HVG_LVG_distribution_nbTF[,1] <- d_background$x

d_HVG <- density(Variable_nbTF$nbTF, from=0, to=253)
d_LVG <- density(Stable_nbTF$nbTF, from=0, to=253)

HVG_LVG_distribution_nbTF[,2] <- d_HVG$y
HVG_LVG_distribution_nbTF[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_nbTF) <- c("nbTF", "HVG", "LVG")
HVG_LVG_distribution_nbTF_long <- gather(HVG_LVG_distribution_nbTF, type, distr, HVG:LVG)


# graph

HVG_LVG_distribution_nbTF <- mutate(HVG_LVG_distribution_nbTF, mean_bootstrat=rowMeans(random_distribution_nbTF_1000[,2:1001]),
                                    bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_nbTF_1000[,2:1001])),
                                    top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_nbTF_1000[,2:1001])))



gather(HVG_LVG_distribution_nbTF, type, distr, HVG:top_CI, factor_key=TRUE) %>%
  ggplot(aes(x=nbTF, y=distr, color=type, linetype=type)) +
  geom_line(size=1) +
  scale_x_continuous(limits = c(0, 250)) +
  xlab("Number of TFs") +
  ylab("Density") +
  scale_color_manual(values=c("blue","green4","black","grey30","grey30")) +
  scale_linetype_manual(values=c("solid","solid","solid", "dashed", "dashed")) +
  theme(panel.grid.major.x = element_blank(),panel.grid.major.y = element_line( size=.1, color="black" ) 
        , panel.background = element_rect(fill="white",colour="grey50"), text=element_text(size=22))


# wilcoxon test

Wilcox_stat=matrix(nrow=1000,ncol=1)
Wilcox_stat=data.frame(Wilcox_stat)


for (i in 1:1000) {
  Wilcox_stat[i,] <- wilcox.test(Variable_nbTF$nbTF,random_nbTF_1000[,i])$statistic
}

names(Wilcox_stat) <- "result"


Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(Variable_nbTF$nbTF,random_nbTF_1000[,i])$p.value
}

names(Wilcox_pvalues) <- "result"



