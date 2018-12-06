
library(tidyverse)

### Plotting of the gene entropy (tissue specificity in gene expression) for HVGs, LVGs and random genes ----


#data preparation


Entropy <- read.table("Entropy/name_entropylog2_FR.txt", header=TRUE, sep='\t')

Background <- read.table("Entropy/All_genes_background.txt", header=TRUE, sep='\t')


Background_Entropy <- inner_join(Background, Entropy, by=c("Gene"="Name")) 


# extract 1000 sets of random genes ----


random_distribution_Entropy_1000=matrix(nrow=512,ncol=1000)
random_distribution_Entropy_1000=data.frame(random_distribution_Entropy_1000)

random_Entropy_1000=matrix(nrow=1358,ncol=1000)
random_Entropy_1000=data.frame(random_Entropy_1000)


for (i in 1:1000) {
  random_set <- Background_Entropy[sample(nrow(Background_Entropy), 1358),]
  random_Entropy_1000[,i] <- random_set[,2]
  d <- density(random_set[,2], from=2.6, to=6)
  random_distribution_Entropy_1000[,i] <- d$y
}

colnames(random_Entropy_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(Background_Entropy$entropyloG2, from=2.6, to=6)
random_distribution_Entropy_1000 <- cbind(d_background$x, random_distribution_Entropy_1000)
colnames(random_distribution_Entropy_1000) <- c("entropy",paste("random_set",1:100, sep="_"))


HVG <- read.table("Entropy/HVG_allZT_background.txt", header=TRUE, sep='\t')

Variable_Entropy <- inner_join(HVG, Entropy, by=c("allZT"="Name")) 


LVG <- read.table("Entropy/LVG1000_allZT_background.txt", header=TRUE, sep='\t')

Stable_Entropy <- inner_join(LVG, Entropy, by=c("allZT"="Name")) 


HVG_LVG_distribution_Entropy=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_Entropy=data.frame(HVG_LVG_distribution_Entropy)

d_background <-  density(Background_Entropy$entropyloG2, from=2.6, to=6)
HVG_LVG_distribution_Entropy[,1] <- d_background$x

d_HVG <- density(Variable_Entropy$entropyloG2, from=2.6, to=6)
d_LVG <- density(Stable_Entropy$entropyloG2, from=2.6, to=6)

HVG_LVG_distribution_Entropy[,2] <- d_HVG$y
HVG_LVG_distribution_Entropy[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_Entropy) <- c("nbTF", "HVG", "LVG")
HVG_LVG_distribution_Entropy_long <- gather(HVG_LVG_distribution_Entropy, type, distr, HVG:LVG)


# graph

HVG_LVG_distribution_Entropy <- mutate(HVG_LVG_distribution_Entropy, mean_bootstrat=rowMeans(random_distribution_Entropy_1000[,2:1001]),
                                       bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_Entropy_1000[,2:1001])),
                                       top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_Entropy_1000[,2:1001])))



gather(HVG_LVG_distribution_Entropy, type, distr, HVG:top_CI, factor_key=TRUE) %>%
  ggplot(aes(x=nbTF, y=distr, color=type, linetype=type)) +
  geom_line(size=1) +
  scale_x_continuous(limits = c(2.6, 6)) +
  xlab("Entropy") +
  ylab("Density") +
  scale_color_manual(values=c("blue","green4","black","grey30","grey30")) +
  scale_linetype_manual(values=c("solid","solid","solid", "dashed", "dashed")) +
  theme(panel.grid.major.x = element_blank(),panel.grid.major.y = element_line( size=.1, color="black" ) 
        , panel.background = element_rect(fill="white",colour="grey50"), text=element_text(size=22))


# wilcoxon test

Wilcox_stat=matrix(nrow=1000,ncol=1)
Wilcox_stat=data.frame(Wilcox_stat)


for (i in 1:1000) {
  Wilcox_stat[i,] <- wilcox.test(Variable_Entropy$entropyloG2,random_Entropy_1000[,i])$statistic
}

names(Wilcox_stat) <- "result"


Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(Variable_Entropy$entropyloG2,random_Entropy_1000[,i])$p.value
}

names(Wilcox_pvalues) <- "result"

