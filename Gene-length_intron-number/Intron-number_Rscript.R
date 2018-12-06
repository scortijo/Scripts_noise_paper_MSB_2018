

### Plotting of the number of introns for HVGs, LVGs and a 1000 sets of random genes ----

#data preparation

Gene_nbIntrons <- read.table("Gene-length_intron-number/TAIR10_Gene_nbIntrons.txt", header=TRUE, sep='\t')
dim(Gene_nbIntrons)

Background <- read.table("Gene-length_intron-number/All_genes_background.txt", header=TRUE, sep='\t')
dim(Background)

Background_gene_nbIntrons <- inner_join(Background, Gene_nbIntrons, by=c("Gene"="gene_name")) %>%
  select(Gene, Number_of_introns)


# extract 1000 sets of random genes ----


random_distribution_nbIntrons_1000=matrix(nrow=512,ncol=1000)
random_distribution_nbIntrons_1000=data.frame(random_distribution_nbIntrons_1000)

random_nbIntrons_1000=matrix(nrow=1358,ncol=1000)
random_nbIntrons_1000=data.frame(random_nbIntrons_1000)


for (i in 1:1000) {
  random_set <- Background_gene_nbIntrons[sample(nrow(Background_gene_nbIntrons), 1358),]
  random_nbIntrons_1000[,i] <- random_set[,2]
  d <- density(random_set[,2], from=0, to=196)
  random_distribution_nbIntrons_1000[,i] <- d$y
}

colnames(random_nbIntrons_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(Background_gene_nbIntrons$Number_of_introns, from=0, to=196)
random_distribution_nbIntrons_1000 <- cbind(d_background$x, random_distribution_nbIntrons_1000)
colnames(random_distribution_nbIntrons_1000) <- c("Number_of_introns",paste("random_set",1:100, sep="_"))


HVG <- read.table("Gene-length_intron-number/HVG_allZT_background.txt", header=TRUE, sep='\t')
dim(HVG)
[1] 1358    2

Variable_nbIntrons <- inner_join(HVG, Gene_nbIntrons, by=c("allZT"="gene_name")) %>%
  select(allZT, Number_of_introns) 


LVG <- read.table("Gene-length_intron-number/LVG1000_allZT_background.txt", header=TRUE, sep='\t')
dim(LVG)
[1] 5727    2

Stable_nbIntrons <- inner_join(LVG, Gene_nbIntrons, by=c("allZT"="gene_name")) %>%
  select(allZT, Number_of_introns) 



HVG_LVG_distribution_nbIntrons=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_nbIntrons=data.frame(HVG_LVG_distribution_nbIntrons)

d_background <-  density(Background_gene_nbIntrons$Number_of_introns, from=0, to=196)
HVG_LVG_distribution_nbIntrons[,1] <- d_background$x

d_HVG <- density(Variable_nbIntrons$Number_of_introns, from=0, to=196)
d_LVG <- density(Stable_nbIntrons$Number_of_introns, from=0, to=196)

HVG_LVG_distribution_nbIntrons[,2] <- d_HVG$y
HVG_LVG_distribution_nbIntrons[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_nbIntrons) <- c("Number_of_introns", "HVG", "LVG")
HVG_LVG_distribution_nbIntrons_long <- gather(HVG_LVG_distribution_nbIntrons, type, distr, HVG:LVG)

# graph

HVG_LVG_distribution_nbIntrons <- mutate(HVG_LVG_distribution_nbIntrons, mean_bootstrat=rowMeans(random_distribution_nbIntrons_1000[,2:1001]),
                                         bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_nbIntrons_1000[,2:1001])),
                                         top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_nbIntrons_1000[,2:1001])))



gather(HVG_LVG_distribution_nbIntrons, type, distr, HVG:top_CI, factor_key=TRUE) %>%
  ggplot(aes(x=Number_of_introns, y=distr, color=type, linetype=type)) +
  geom_line(size=1) +
  scale_x_continuous(limits = c(0, 50)) +
  xlab("Number of introns") +
  ylab("Density") +
  scale_color_manual(values=c("blue","green4","black","grey30","grey30")) +
  scale_linetype_manual(values=c("solid","solid","solid", "dashed", "dashed")) +
  theme(panel.grid.major.x = element_blank(),panel.grid.major.y = element_line( size=.1, color="black" ) 
        , panel.background = element_rect(fill="white",colour="grey50"), text=element_text(size=22))

# wilcoxon test

Wilcox_stat=matrix(nrow=1000,ncol=1)
Wilcox_stat=data.frame(Wilcox_stat)


for (i in 1:1000) {
  Wilcox_stat[i,] <- wilcox.test(Variable_nbIntrons$Number_of_introns,random_nbIntrons_1000[,i])$statistic
}

names(Wilcox_stat) <- "result"


Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(Variable_nbIntrons$Number_of_introns,random_nbIntrons_1000[,i])$p.value
}

names(Wilcox_pvalues) <- "result"


#All values between 3.852037e-38 and 1.002191e-64 #




### Plotting of the number of introns of HVGs, LVGs and a 1000 sets of random genes for each ZT ----


# data preparation

Gene_nbIntrons <- read.table("Gene-length_intron-number/TAIR10_Gene_nbIntrons.txt", header=TRUE, sep='\t')
dim(Gene_nbIntrons)

Background <- read.table("Gene-length_intron-number/All_genes_background.txt", header=TRUE, sep='\t')
dim(Background)

Gene_nbIntrons <- inner_join(Gene_nbIntrons, Background, 
                             by=c("gene_name"= "Gene"))


# bootstrap for each ZT


# ZT2 ----
random_distribution_nbIntron_ZT2_1000=matrix(nrow=512,ncol=1000)
random_distribution_nbIntron_ZT2_1000=data.frame(random_distribution_nbIntron_ZT2_1000)

random_nbIntron_ZT2_1000=matrix(nrow=313,ncol=1000)
random_nbIntron_ZT2_1000=data.frame(random_nbIntron_ZT2_1000)


for (i in 1:1000) {
  random_set <- Gene_nbIntrons[sample(nrow(Gene_nbIntrons), 313),]
  random_nbIntron_ZT2_1000[,i] <- random_set[,2]
  d <- density(random_set[,2], from=0, to=196, na.rm=TRUE)
  random_distribution_nbIntron_ZT2_1000[,i] <- d$y
}

colnames(random_nbIntron_ZT2_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(Gene_nbIntrons$Number_of_introns, from=0, to=196, na.rm=TRUE)
random_distribution_nbIntron_ZT2_1000 <- cbind(d_background$x, random_distribution_nbIntron_ZT2_1000)
colnames(random_distribution_nbIntron_ZT2_1000) <- c("nbIntron",paste("random_set",1:100, sep="_"))


ZT2 <- read.table("Gene-length_intron-number/ZT2_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT2)

ZT2_nbIntron_HVG <- merge(ZT2,Gene_nbIntrons,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,Number_of_introns) %>%
  filter(TypeV=="Variable")

ZT2_nbIntron_LVG <- merge(ZT2,Gene_nbIntrons,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,Number_of_introns) %>%
  filter(TypeV=="Stable")



HVG_LVG_distribution_nbIntron_ZT2=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_nbIntron_ZT2=data.frame(HVG_LVG_distribution_nbIntron_ZT2)

d_background <-  density(Gene_nbIntrons$Number_of_introns, from=0, to=196, na.rm=TRUE)
HVG_LVG_distribution_nbIntron_ZT2[,1] <- d_background$x

d_HVG <- density(ZT2_nbIntron_HVG$Number_of_introns, from=0, to=196, na.rm=TRUE)
d_LVG <- density(ZT2_nbIntron_LVG$Number_of_introns, from=0, to=196, na.rm=TRUE)

HVG_LVG_distribution_nbIntron_ZT2[,2] <- d_HVG$y
HVG_LVG_distribution_nbIntron_ZT2[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_nbIntron_ZT2) <- c("nbIntron", "HVG", "LVG")
HVG_LVG_distribution_nbIntron_ZT2_long <- gather(HVG_LVG_distribution_nbIntron_ZT2, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_nbIntron_ZT2 <- mutate(HVG_LVG_distribution_nbIntron_ZT2, mean_bootstrat=rowMeans(random_distribution_nbIntron_ZT2_1000[,2:1001]),
                                            bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_nbIntron_ZT2_1000[,2:1001])),
                                            top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_nbIntron_ZT2_1000[,2:1001])),
                                            ZT= "ZT2")


# wilcoxon test

Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(ZT2_nbIntron_HVG$Number_of_introns,random_nbIntron_ZT2_1000[,i])$p.value
}

# All values between 2.281532e-14 and 4.119118e-34 #



# ZT4 ----
random_distribution_nbIntron_ZT4_1000=matrix(nrow=512,ncol=1000)
random_distribution_nbIntron_ZT4_1000=data.frame(random_distribution_nbIntron_ZT4_1000)

random_nbIntron_ZT4_1000=matrix(nrow=257,ncol=1000)
random_nbIntron_ZT4_1000=data.frame(random_nbIntron_ZT4_1000)


for (i in 1:1000) {
  random_set <- Gene_nbIntrons[sample(nrow(Gene_nbIntrons), 257),]
  random_nbIntron_ZT4_1000[,i] <- random_set[,2]
  d <- density(random_set[,2], from=0, to=196, na.rm=TRUE)
  random_distribution_nbIntron_ZT4_1000[,i] <- d$y
}

colnames(random_nbIntron_ZT4_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(Gene_nbIntrons$Number_of_introns, from=0, to=196, na.rm=TRUE)
random_distribution_nbIntron_ZT4_1000 <- cbind(d_background$x, random_distribution_nbIntron_ZT4_1000)
colnames(random_distribution_nbIntron_ZT4_1000) <- c("nbIntron",paste("random_set",1:100, sep="_"))


ZT4 <- read.table("Gene-length_intron-number/ZT4_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT4)

ZT4_nbIntron_HVG <- merge(ZT4,Gene_nbIntrons,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,Number_of_introns) %>%
  filter(TypeV=="Variable")

ZT4_nbIntron_LVG <- merge(ZT4,Gene_nbIntrons,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,Number_of_introns) %>%
  filter(TypeV=="Stable")



HVG_LVG_distribution_nbIntron_ZT4=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_nbIntron_ZT4=data.frame(HVG_LVG_distribution_nbIntron_ZT4)

d_background <-  density(Gene_nbIntrons$Number_of_introns, from=0, to=196, na.rm=TRUE)
HVG_LVG_distribution_nbIntron_ZT4[,1] <- d_background$x

d_HVG <- density(ZT4_nbIntron_HVG$Number_of_introns, from=0, to=196, na.rm=TRUE)
d_LVG <- density(ZT4_nbIntron_LVG$Number_of_introns, from=0, to=196, na.rm=TRUE)

HVG_LVG_distribution_nbIntron_ZT4[,2] <- d_HVG$y
HVG_LVG_distribution_nbIntron_ZT4[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_nbIntron_ZT4) <- c("nbIntron", "HVG", "LVG")
HVG_LVG_distribution_nbIntron_ZT4_long <- gather(HVG_LVG_distribution_nbIntron_ZT4, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_nbIntron_ZT4 <- mutate(HVG_LVG_distribution_nbIntron_ZT4, mean_bootstrat=rowMeans(random_distribution_nbIntron_ZT4_1000[,2:1001]),
                                            bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_nbIntron_ZT4_1000[,2:1001])),
                                            top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_nbIntron_ZT4_1000[,2:1001])),
                                            ZT= "ZT4")


# wilcoxon test

Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(ZT4_nbIntron_HVG$Number_of_introns,random_nbIntron_ZT4_1000[,i])$p.value
}

# All values between 9.453940e-11 and 6.608813e-26 #



# ZT6 ----
random_distribution_nbIntron_ZT6_1000=matrix(nrow=512,ncol=1000)
random_distribution_nbIntron_ZT6_1000=data.frame(random_distribution_nbIntron_ZT6_1000)

random_nbIntron_ZT6_1000=matrix(nrow=374,ncol=1000)
random_nbIntron_ZT6_1000=data.frame(random_nbIntron_ZT6_1000)


for (i in 1:1000) {
  random_set <- Gene_nbIntrons[sample(nrow(Gene_nbIntrons), 374),]
  random_nbIntron_ZT6_1000[,i] <- random_set[,2]
  d <- density(random_set[,2], from=0, to=196, na.rm=TRUE)
  random_distribution_nbIntron_ZT6_1000[,i] <- d$y
}

colnames(random_nbIntron_ZT6_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(Gene_nbIntrons$Number_of_introns, from=0, to=196, na.rm=TRUE)
random_distribution_nbIntron_ZT6_1000 <- cbind(d_background$x, random_distribution_nbIntron_ZT6_1000)
colnames(random_distribution_nbIntron_ZT6_1000) <- c("nbIntron",paste("random_set",1:100, sep="_"))


ZT6 <- read.table("Gene-length_intron-number/ZT6_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT6)

ZT6_nbIntron_HVG <- merge(ZT6,Gene_nbIntrons,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,Number_of_introns) %>%
  filter(TypeV=="Variable")

ZT6_nbIntron_LVG <- merge(ZT6,Gene_nbIntrons,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,Number_of_introns) %>%
  filter(TypeV=="Stable")



HVG_LVG_distribution_nbIntron_ZT6=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_nbIntron_ZT6=data.frame(HVG_LVG_distribution_nbIntron_ZT6)

d_background <-  density(Gene_nbIntrons$Number_of_introns, from=0, to=196, na.rm=TRUE)
HVG_LVG_distribution_nbIntron_ZT6[,1] <- d_background$x

d_HVG <- density(ZT6_nbIntron_HVG$Number_of_introns, from=0, to=196, na.rm=TRUE)
d_LVG <- density(ZT6_nbIntron_LVG$Number_of_introns, from=0, to=196, na.rm=TRUE)

HVG_LVG_distribution_nbIntron_ZT6[,2] <- d_HVG$y
HVG_LVG_distribution_nbIntron_ZT6[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_nbIntron_ZT6) <- c("nbIntron", "HVG", "LVG")
HVG_LVG_distribution_nbIntron_ZT6_long <- gather(HVG_LVG_distribution_nbIntron_ZT6, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_nbIntron_ZT6 <- mutate(HVG_LVG_distribution_nbIntron_ZT6, mean_bootstrat=rowMeans(random_distribution_nbIntron_ZT6_1000[,2:1001]),
                                            bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_nbIntron_ZT6_1000[,2:1001])),
                                            top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_nbIntron_ZT6_1000[,2:1001])),
                                            ZT= "ZT6")


# wilcoxon test

Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(ZT6_nbIntron_HVG$Number_of_introns,random_nbIntron_ZT6_1000[,i])$p.value
}

# All values between 1.764752e-13 and 2.000417e-33 #


# ZT8 ----
random_distribution_nbIntron_ZT8_1000=matrix(nrow=512,ncol=1000)
random_distribution_nbIntron_ZT8_1000=data.frame(random_distribution_nbIntron_ZT8_1000)

random_nbIntron_ZT8_1000=matrix(nrow=327,ncol=1000)
random_nbIntron_ZT8_1000=data.frame(random_nbIntron_ZT8_1000)


for (i in 1:1000) {
  random_set <- Gene_nbIntrons[sample(nrow(Gene_nbIntrons), 327),]
  random_nbIntron_ZT8_1000[,i] <- random_set[,2]
  d <- density(random_set[,2], from=0, to=196, na.rm=TRUE)
  random_distribution_nbIntron_ZT8_1000[,i] <- d$y
}

colnames(random_nbIntron_ZT8_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(Gene_nbIntrons$Number_of_introns, from=0, to=196, na.rm=TRUE)
random_distribution_nbIntron_ZT8_1000 <- cbind(d_background$x, random_distribution_nbIntron_ZT8_1000)
colnames(random_distribution_nbIntron_ZT8_1000) <- c("nbIntron",paste("random_set",1:100, sep="_"))


ZT8 <- read.table("Gene-length_intron-number/ZT8_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT8)

ZT8_nbIntron_HVG <- merge(ZT8,Gene_nbIntrons,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,Number_of_introns) %>%
  filter(TypeV=="Variable")

ZT8_nbIntron_LVG <- merge(ZT8,Gene_nbIntrons,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,Number_of_introns) %>%
  filter(TypeV=="Stable")



HVG_LVG_distribution_nbIntron_ZT8=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_nbIntron_ZT8=data.frame(HVG_LVG_distribution_nbIntron_ZT8)

d_background <-  density(Gene_nbIntrons$Number_of_introns, from=0, to=196, na.rm=TRUE)
HVG_LVG_distribution_nbIntron_ZT8[,1] <- d_background$x

d_HVG <- density(ZT8_nbIntron_HVG$Number_of_introns, from=0, to=196, na.rm=TRUE)
d_LVG <- density(ZT8_nbIntron_LVG$Number_of_introns, from=0, to=196, na.rm=TRUE)

HVG_LVG_distribution_nbIntron_ZT8[,2] <- d_HVG$y
HVG_LVG_distribution_nbIntron_ZT8[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_nbIntron_ZT8) <- c("nbIntron", "HVG", "LVG")
HVG_LVG_distribution_nbIntron_ZT8_long <- gather(HVG_LVG_distribution_nbIntron_ZT8, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_nbIntron_ZT8 <- mutate(HVG_LVG_distribution_nbIntron_ZT8, mean_bootstrat=rowMeans(random_distribution_nbIntron_ZT8_1000[,2:1001]),
                                            bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_nbIntron_ZT8_1000[,2:1001])),
                                            top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_nbIntron_ZT8_1000[,2:1001])),
                                            ZT= "ZT8")


# wilcoxon test

Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(ZT8_nbIntron_HVG$Number_of_introns,random_nbIntron_ZT8_1000[,i])$p.value
}

# All values between 2.848094e-11 and 5.415106e-28 #



# ZT10 ----
random_distribution_nbIntron_ZT10_1000=matrix(nrow=512,ncol=1000)
random_distribution_nbIntron_ZT10_1000=data.frame(random_distribution_nbIntron_ZT10_1000)

random_nbIntron_ZT10_1000=matrix(nrow=310,ncol=1000)
random_nbIntron_ZT10_1000=data.frame(random_nbIntron_ZT10_1000)


for (i in 1:1000) {
  random_set <- Gene_nbIntrons[sample(nrow(Gene_nbIntrons), 310),]
  random_nbIntron_ZT10_1000[,i] <- random_set[,2]
  d <- density(random_set[,2], from=0, to=196, na.rm=TRUE)
  random_distribution_nbIntron_ZT10_1000[,i] <- d$y
}

colnames(random_nbIntron_ZT10_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(Gene_nbIntrons$Number_of_introns, from=0, to=196, na.rm=TRUE)
random_distribution_nbIntron_ZT10_1000 <- cbind(d_background$x, random_distribution_nbIntron_ZT10_1000)
colnames(random_distribution_nbIntron_ZT10_1000) <- c("nbIntron",paste("random_set",1:100, sep="_"))


ZT10 <- read.table("Gene-length_intron-number/ZT10_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT10)

ZT10_nbIntron_HVG <- merge(ZT10,Gene_nbIntrons,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,Number_of_introns) %>%
  filter(TypeV=="Variable")

ZT10_nbIntron_LVG <- merge(ZT10,Gene_nbIntrons,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,Number_of_introns) %>%
  filter(TypeV=="Stable")



HVG_LVG_distribution_nbIntron_ZT10=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_nbIntron_ZT10=data.frame(HVG_LVG_distribution_nbIntron_ZT10)

d_background <-  density(Gene_nbIntrons$Number_of_introns, from=0, to=196, na.rm=TRUE)
HVG_LVG_distribution_nbIntron_ZT10[,1] <- d_background$x

d_HVG <- density(ZT10_nbIntron_HVG$Number_of_introns, from=0, to=196, na.rm=TRUE)
d_LVG <- density(ZT10_nbIntron_LVG$Number_of_introns, from=0, to=196, na.rm=TRUE)

HVG_LVG_distribution_nbIntron_ZT10[,2] <- d_HVG$y
HVG_LVG_distribution_nbIntron_ZT10[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_nbIntron_ZT10) <- c("nbIntron", "HVG", "LVG")
HVG_LVG_distribution_nbIntron_ZT10_long <- gather(HVG_LVG_distribution_nbIntron_ZT10, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_nbIntron_ZT10 <- mutate(HVG_LVG_distribution_nbIntron_ZT10, mean_bootstrat=rowMeans(random_distribution_nbIntron_ZT10_1000[,2:1001]),
                                             bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_nbIntron_ZT10_1000[,2:1001])),
                                             top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_nbIntron_ZT10_1000[,2:1001])),
                                             ZT= "ZT10")


# wilcoxon test

Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(ZT10_nbIntron_HVG$Number_of_introns,random_nbIntron_ZT10_1000[,i])$p.value
}

# All values between 2.848094e-11 and 5.415106e-28 #



# ZT12 ----
random_distribution_nbIntron_ZT12_1000=matrix(nrow=512,ncol=1000)
random_distribution_nbIntron_ZT12_1000=data.frame(random_distribution_nbIntron_ZT12_1000)

random_nbIntron_ZT12_1000=matrix(nrow=636,ncol=1000)
random_nbIntron_ZT12_1000=data.frame(random_nbIntron_ZT12_1000)


for (i in 1:1000) {
  random_set <- Gene_nbIntrons[sample(nrow(Gene_nbIntrons), 636),]
  random_nbIntron_ZT12_1000[,i] <- random_set[,2]
  d <- density(random_set[,2], from=0, to=196, na.rm=TRUE)
  random_distribution_nbIntron_ZT12_1000[,i] <- d$y
}

colnames(random_nbIntron_ZT12_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(Gene_nbIntrons$Number_of_introns, from=0, to=196, na.rm=TRUE)
random_distribution_nbIntron_ZT12_1000 <- cbind(d_background$x, random_distribution_nbIntron_ZT12_1000)
colnames(random_distribution_nbIntron_ZT12_1000) <- c("nbIntron",paste("random_set",1:100, sep="_"))


ZT12 <- read.table("Gene-length_intron-number/ZT12_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT12)

ZT12_nbIntron_HVG <- merge(ZT12,Gene_nbIntrons,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,Number_of_introns) %>%
  filter(TypeV=="Variable")

ZT12_nbIntron_LVG <- merge(ZT12,Gene_nbIntrons,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,Number_of_introns) %>%
  filter(TypeV=="Stable")



HVG_LVG_distribution_nbIntron_ZT12=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_nbIntron_ZT12=data.frame(HVG_LVG_distribution_nbIntron_ZT12)

d_background <-  density(Gene_nbIntrons$Number_of_introns, from=0, to=196, na.rm=TRUE)
HVG_LVG_distribution_nbIntron_ZT12[,1] <- d_background$x

d_HVG <- density(ZT12_nbIntron_HVG$Number_of_introns, from=0, to=196, na.rm=TRUE)
d_LVG <- density(ZT12_nbIntron_LVG$Number_of_introns, from=0, to=196, na.rm=TRUE)

HVG_LVG_distribution_nbIntron_ZT12[,2] <- d_HVG$y
HVG_LVG_distribution_nbIntron_ZT12[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_nbIntron_ZT12) <- c("nbIntron", "HVG", "LVG")
HVG_LVG_distribution_nbIntron_ZT12_long <- gather(HVG_LVG_distribution_nbIntron_ZT12, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_nbIntron_ZT12 <- mutate(HVG_LVG_distribution_nbIntron_ZT12, mean_bootstrat=rowMeans(random_distribution_nbIntron_ZT12_1000[,2:1001]),
                                             bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_nbIntron_ZT12_1000[,2:1001])),
                                             top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_nbIntron_ZT12_1000[,2:1001])),
                                             ZT= "ZT12")


# wilcoxon test

Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(ZT12_nbIntron_HVG$Number_of_introns,random_nbIntron_ZT12_1000[,i])$p.value
}

# All values between 9.616718e-20 and 2.302577e-40 #


# ZT14 ----
random_distribution_nbIntron_ZT14_1000=matrix(nrow=512,ncol=1000)
random_distribution_nbIntron_ZT14_1000=data.frame(random_distribution_nbIntron_ZT14_1000)

random_nbIntron_ZT14_1000=matrix(nrow=418,ncol=1000)
random_nbIntron_ZT14_1000=data.frame(random_nbIntron_ZT14_1000)


for (i in 1:1000) {
  random_set <- Gene_nbIntrons[sample(nrow(Gene_nbIntrons), 418),]
  random_nbIntron_ZT14_1000[,i] <- random_set[,2]
  d <- density(random_set[,2], from=0, to=196, na.rm=TRUE)
  random_distribution_nbIntron_ZT14_1000[,i] <- d$y
}

colnames(random_nbIntron_ZT14_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(Gene_nbIntrons$Number_of_introns, from=0, to=196, na.rm=TRUE)
random_distribution_nbIntron_ZT14_1000 <- cbind(d_background$x, random_distribution_nbIntron_ZT14_1000)
colnames(random_distribution_nbIntron_ZT14_1000) <- c("nbIntron",paste("random_set",1:100, sep="_"))


ZT14 <- read.table("Gene-length_intron-number/ZT14_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT14)

ZT14_nbIntron_HVG <- merge(ZT14,Gene_nbIntrons,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,Number_of_introns) %>%
  filter(TypeV=="Variable")

ZT14_nbIntron_LVG <- merge(ZT14,Gene_nbIntrons,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,Number_of_introns) %>%
  filter(TypeV=="Stable")



HVG_LVG_distribution_nbIntron_ZT14=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_nbIntron_ZT14=data.frame(HVG_LVG_distribution_nbIntron_ZT14)

d_background <-  density(Gene_nbIntrons$Number_of_introns, from=0, to=196, na.rm=TRUE)
HVG_LVG_distribution_nbIntron_ZT14[,1] <- d_background$x

d_HVG <- density(ZT14_nbIntron_HVG$Number_of_introns, from=0, to=196, na.rm=TRUE)
d_LVG <- density(ZT14_nbIntron_LVG$Number_of_introns, from=0, to=196, na.rm=TRUE)

HVG_LVG_distribution_nbIntron_ZT14[,2] <- d_HVG$y
HVG_LVG_distribution_nbIntron_ZT14[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_nbIntron_ZT14) <- c("nbIntron", "HVG", "LVG")
HVG_LVG_distribution_nbIntron_ZT14_long <- gather(HVG_LVG_distribution_nbIntron_ZT14, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_nbIntron_ZT14 <- mutate(HVG_LVG_distribution_nbIntron_ZT14, mean_bootstrat=rowMeans(random_distribution_nbIntron_ZT14_1000[,2:1001]),
                                             bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_nbIntron_ZT14_1000[,2:1001])),
                                             top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_nbIntron_ZT14_1000[,2:1001])),
                                             ZT= "ZT14")


# wilcoxon test

Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(ZT14_nbIntron_HVG$Number_of_introns,random_nbIntron_ZT14_1000[,i])$p.value
}

# All values between 5.637005e-11 and 3.400728e-28 #



# ZT16 ----
random_distribution_nbIntron_ZT16_1000=matrix(nrow=512,ncol=1000)
random_distribution_nbIntron_ZT16_1000=data.frame(random_distribution_nbIntron_ZT16_1000)

random_nbIntron_ZT16_1000=matrix(nrow=630,ncol=1000)
random_nbIntron_ZT16_1000=data.frame(random_nbIntron_ZT16_1000)


for (i in 1:1000) {
  random_set <- Gene_nbIntrons[sample(nrow(Gene_nbIntrons), 630),]
  random_nbIntron_ZT16_1000[,i] <- random_set[,2]
  d <- density(random_set[,2], from=0, to=196, na.rm=TRUE)
  random_distribution_nbIntron_ZT16_1000[,i] <- d$y
}

colnames(random_nbIntron_ZT16_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(Gene_nbIntrons$Number_of_introns, from=0, to=196, na.rm=TRUE)
random_distribution_nbIntron_ZT16_1000 <- cbind(d_background$x, random_distribution_nbIntron_ZT16_1000)
colnames(random_distribution_nbIntron_ZT16_1000) <- c("nbIntron",paste("random_set",1:100, sep="_"))


ZT16 <- read.table("Gene-length_intron-number/ZT16_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT16)

ZT16_nbIntron_HVG <- merge(ZT16,Gene_nbIntrons,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,Number_of_introns) %>%
  filter(TypeV=="Variable")

ZT16_nbIntron_LVG <- merge(ZT16,Gene_nbIntrons,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,Number_of_introns) %>%
  filter(TypeV=="Stable")



HVG_LVG_distribution_nbIntron_ZT16=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_nbIntron_ZT16=data.frame(HVG_LVG_distribution_nbIntron_ZT16)

d_background <-  density(Gene_nbIntrons$Number_of_introns, from=0, to=196, na.rm=TRUE)
HVG_LVG_distribution_nbIntron_ZT16[,1] <- d_background$x

d_HVG <- density(ZT16_nbIntron_HVG$Number_of_introns, from=0, to=196, na.rm=TRUE)
d_LVG <- density(ZT16_nbIntron_LVG$Number_of_introns, from=0, to=196, na.rm=TRUE)

HVG_LVG_distribution_nbIntron_ZT16[,2] <- d_HVG$y
HVG_LVG_distribution_nbIntron_ZT16[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_nbIntron_ZT16) <- c("nbIntron", "HVG", "LVG")
HVG_LVG_distribution_nbIntron_ZT16_long <- gather(HVG_LVG_distribution_nbIntron_ZT16, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_nbIntron_ZT16 <- mutate(HVG_LVG_distribution_nbIntron_ZT16, mean_bootstrat=rowMeans(random_distribution_nbIntron_ZT16_1000[,2:1001]),
                                             bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_nbIntron_ZT16_1000[,2:1001])),
                                             top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_nbIntron_ZT16_1000[,2:1001])),
                                             ZT= "ZT16")


# wilcoxon test

Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(ZT16_nbIntron_HVG$Number_of_introns,random_nbIntron_ZT16_1000[,i])$p.value
}

# All values between 4.302776e-14 and 7.273272e-31 #




# ZT18 ----
random_distribution_nbIntron_ZT18_1000=matrix(nrow=512,ncol=1000)
random_distribution_nbIntron_ZT18_1000=data.frame(random_distribution_nbIntron_ZT18_1000)

random_nbIntron_ZT18_1000=matrix(nrow=408,ncol=1000)
random_nbIntron_ZT18_1000=data.frame(random_nbIntron_ZT18_1000)


for (i in 1:1000) {
  random_set <- Gene_nbIntrons[sample(nrow(Gene_nbIntrons), 408),]
  random_nbIntron_ZT18_1000[,i] <- random_set[,2]
  d <- density(random_set[,2], from=0, to=196, na.rm=TRUE)
  random_distribution_nbIntron_ZT18_1000[,i] <- d$y
}

colnames(random_nbIntron_ZT18_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(Gene_nbIntrons$Number_of_introns, from=0, to=196, na.rm=TRUE)
random_distribution_nbIntron_ZT18_1000 <- cbind(d_background$x, random_distribution_nbIntron_ZT18_1000)
colnames(random_distribution_nbIntron_ZT18_1000) <- c("nbIntron",paste("random_set",1:100, sep="_"))


ZT18 <- read.table("Gene-length_intron-number/ZT18_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT18)

ZT18_nbIntron_HVG <- merge(ZT18,Gene_nbIntrons,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,Number_of_introns) %>%
  filter(TypeV=="Variable")

ZT18_nbIntron_LVG <- merge(ZT18,Gene_nbIntrons,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,Number_of_introns) %>%
  filter(TypeV=="Stable")



HVG_LVG_distribution_nbIntron_ZT18=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_nbIntron_ZT18=data.frame(HVG_LVG_distribution_nbIntron_ZT18)

d_background <-  density(Gene_nbIntrons$Number_of_introns, from=0, to=196, na.rm=TRUE)
HVG_LVG_distribution_nbIntron_ZT18[,1] <- d_background$x

d_HVG <- density(ZT18_nbIntron_HVG$Number_of_introns, from=0, to=196, na.rm=TRUE)
d_LVG <- density(ZT18_nbIntron_LVG$Number_of_introns, from=0, to=196, na.rm=TRUE)

HVG_LVG_distribution_nbIntron_ZT18[,2] <- d_HVG$y
HVG_LVG_distribution_nbIntron_ZT18[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_nbIntron_ZT18) <- c("nbIntron", "HVG", "LVG")
HVG_LVG_distribution_nbIntron_ZT18_long <- gather(HVG_LVG_distribution_nbIntron_ZT18, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_nbIntron_ZT18 <- mutate(HVG_LVG_distribution_nbIntron_ZT18, mean_bootstrat=rowMeans(random_distribution_nbIntron_ZT18_1000[,2:1001]),
                                             bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_nbIntron_ZT18_1000[,2:1001])),
                                             top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_nbIntron_ZT18_1000[,2:1001])),
                                             ZT= "ZT18")


# wilcoxon test

Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(ZT18_nbIntron_HVG$Number_of_introns,random_nbIntron_ZT18_1000[,i])$p.value
}

# All values between 5.839084e-12 and 4.066486e-28 #


# ZT20 ----
random_distribution_nbIntron_ZT20_1000=matrix(nrow=512,ncol=1000)
random_distribution_nbIntron_ZT20_1000=data.frame(random_distribution_nbIntron_ZT20_1000)

random_nbIntron_ZT20_1000=matrix(nrow=716,ncol=1000)
random_nbIntron_ZT20_1000=data.frame(random_nbIntron_ZT20_1000)


for (i in 1:1000) {
  random_set <- Gene_nbIntrons[sample(nrow(Gene_nbIntrons), 716),]
  random_nbIntron_ZT20_1000[,i] <- random_set[,2]
  d <- density(random_set[,2], from=0, to=196, na.rm=TRUE)
  random_distribution_nbIntron_ZT20_1000[,i] <- d$y
}

colnames(random_nbIntron_ZT20_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(Gene_nbIntrons$Number_of_introns, from=0, to=196, na.rm=TRUE)
random_distribution_nbIntron_ZT20_1000 <- cbind(d_background$x, random_distribution_nbIntron_ZT20_1000)
colnames(random_distribution_nbIntron_ZT20_1000) <- c("nbIntron",paste("random_set",1:100, sep="_"))


ZT20 <- read.table("Gene-length_intron-number/ZT20_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT20)

ZT20_nbIntron_HVG <- merge(ZT20,Gene_nbIntrons,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,Number_of_introns) %>%
  filter(TypeV=="Variable")

ZT20_nbIntron_LVG <- merge(ZT20,Gene_nbIntrons,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,Number_of_introns) %>%
  filter(TypeV=="Stable")



HVG_LVG_distribution_nbIntron_ZT20=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_nbIntron_ZT20=data.frame(HVG_LVG_distribution_nbIntron_ZT20)

d_background <-  density(Gene_nbIntrons$Number_of_introns, from=0, to=196, na.rm=TRUE)
HVG_LVG_distribution_nbIntron_ZT20[,1] <- d_background$x

d_HVG <- density(ZT20_nbIntron_HVG$Number_of_introns, from=0, to=196, na.rm=TRUE)
d_LVG <- density(ZT20_nbIntron_LVG$Number_of_introns, from=0, to=196, na.rm=TRUE)

HVG_LVG_distribution_nbIntron_ZT20[,2] <- d_HVG$y
HVG_LVG_distribution_nbIntron_ZT20[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_nbIntron_ZT20) <- c("nbIntron", "HVG", "LVG")
HVG_LVG_distribution_nbIntron_ZT20_long <- gather(HVG_LVG_distribution_nbIntron_ZT20, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_nbIntron_ZT20 <- mutate(HVG_LVG_distribution_nbIntron_ZT20, mean_bootstrat=rowMeans(random_distribution_nbIntron_ZT20_1000[,2:1001]),
                                             bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_nbIntron_ZT20_1000[,2:1001])),
                                             top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_nbIntron_ZT20_1000[,2:1001])),
                                             ZT= "ZT20")


# wilcoxon test

Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(ZT20_nbIntron_HVG$Number_of_introns,random_nbIntron_ZT20_1000[,i])$p.value
}

# All values between 1.765287e-19 and 1.289474e-42 #



# ZT22 ----
random_distribution_nbIntron_ZT22_1000=matrix(nrow=512,ncol=1000)
random_distribution_nbIntron_ZT22_1000=data.frame(random_distribution_nbIntron_ZT22_1000)

random_nbIntron_ZT22_1000=matrix(nrow=566,ncol=1000)
random_nbIntron_ZT22_1000=data.frame(random_nbIntron_ZT22_1000)


for (i in 1:1000) {
  random_set <- Gene_nbIntrons[sample(nrow(Gene_nbIntrons), 566),]
  random_nbIntron_ZT22_1000[,i] <- random_set[,2]
  d <- density(random_set[,2], from=0, to=196, na.rm=TRUE)
  random_distribution_nbIntron_ZT22_1000[,i] <- d$y
}

colnames(random_nbIntron_ZT22_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(Gene_nbIntrons$Number_of_introns, from=0, to=196, na.rm=TRUE)
random_distribution_nbIntron_ZT22_1000 <- cbind(d_background$x, random_distribution_nbIntron_ZT22_1000)
colnames(random_distribution_nbIntron_ZT22_1000) <- c("nbIntron",paste("random_set",1:100, sep="_"))


ZT22 <- read.table("Gene-length_intron-number/ZT22_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT22)

ZT22_nbIntron_HVG <- merge(ZT22,Gene_nbIntrons,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,Number_of_introns) %>%
  filter(TypeV=="Variable")

ZT22_nbIntron_LVG <- merge(ZT22,Gene_nbIntrons,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,Number_of_introns) %>%
  filter(TypeV=="Stable")



HVG_LVG_distribution_nbIntron_ZT22=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_nbIntron_ZT22=data.frame(HVG_LVG_distribution_nbIntron_ZT22)

d_background <-  density(Gene_nbIntrons$Number_of_introns, from=0, to=196, na.rm=TRUE)
HVG_LVG_distribution_nbIntron_ZT22[,1] <- d_background$x

d_HVG <- density(ZT22_nbIntron_HVG$Number_of_introns, from=0, to=196, na.rm=TRUE)
d_LVG <- density(ZT22_nbIntron_LVG$Number_of_introns, from=0, to=196, na.rm=TRUE)

HVG_LVG_distribution_nbIntron_ZT22[,2] <- d_HVG$y
HVG_LVG_distribution_nbIntron_ZT22[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_nbIntron_ZT22) <- c("nbIntron", "HVG", "LVG")
HVG_LVG_distribution_nbIntron_ZT22_long <- gather(HVG_LVG_distribution_nbIntron_ZT22, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_nbIntron_ZT22 <- mutate(HVG_LVG_distribution_nbIntron_ZT22, mean_bootstrat=rowMeans(random_distribution_nbIntron_ZT22_1000[,2:1001]),
                                             bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_nbIntron_ZT22_1000[,2:1001])),
                                             top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_nbIntron_ZT22_1000[,2:1001])),
                                             ZT= "ZT22")


# wilcoxon test

Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(ZT22_nbIntron_HVG$Number_of_introns,random_nbIntron_ZT22_1000[,i])$p.value
}

# All values between 2.302640e-21 and 2.301209e-41 #



# ZT24 ----
random_distribution_nbIntron_ZT24_1000=matrix(nrow=512,ncol=1000)
random_distribution_nbIntron_ZT24_1000=data.frame(random_distribution_nbIntron_ZT24_1000)

random_nbIntron_ZT24_1000=matrix(nrow=377,ncol=1000)
random_nbIntron_ZT24_1000=data.frame(random_nbIntron_ZT24_1000)


for (i in 1:1000) {
  random_set <- Gene_nbIntrons[sample(nrow(Gene_nbIntrons), 377),]
  random_nbIntron_ZT24_1000[,i] <- random_set[,2]
  d <- density(random_set[,2], from=0, to=196, na.rm=TRUE)
  random_distribution_nbIntron_ZT24_1000[,i] <- d$y
}

colnames(random_nbIntron_ZT24_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(Gene_nbIntrons$Number_of_introns, from=0, to=196, na.rm=TRUE)
random_distribution_nbIntron_ZT24_1000 <- cbind(d_background$x, random_distribution_nbIntron_ZT24_1000)
colnames(random_distribution_nbIntron_ZT24_1000) <- c("nbIntron",paste("random_set",1:100, sep="_"))


ZT24 <- read.table("Gene-length_intron-number/ZT0_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT24)

ZT24_nbIntron_HVG <- merge(ZT24,Gene_nbIntrons,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,Number_of_introns) %>%
  filter(TypeV=="Variable")

ZT24_nbIntron_LVG <- merge(ZT24,Gene_nbIntrons,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,Number_of_introns) %>%
  filter(TypeV=="Stable")



HVG_LVG_distribution_nbIntron_ZT24=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_nbIntron_ZT24=data.frame(HVG_LVG_distribution_nbIntron_ZT24)

d_background <-  density(Gene_nbIntrons$Number_of_introns, from=0, to=196, na.rm=TRUE)
HVG_LVG_distribution_nbIntron_ZT24[,1] <- d_background$x

d_HVG <- density(ZT24_nbIntron_HVG$Number_of_introns, from=0, to=196, na.rm=TRUE)
d_LVG <- density(ZT24_nbIntron_LVG$Number_of_introns, from=0, to=196, na.rm=TRUE)

HVG_LVG_distribution_nbIntron_ZT24[,2] <- d_HVG$y
HVG_LVG_distribution_nbIntron_ZT24[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_nbIntron_ZT24) <- c("nbIntron", "HVG", "LVG")
HVG_LVG_distribution_nbIntron_ZT24_long <- gather(HVG_LVG_distribution_nbIntron_ZT24, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_nbIntron_ZT24 <- mutate(HVG_LVG_distribution_nbIntron_ZT24, mean_bootstrat=rowMeans(random_distribution_nbIntron_ZT24_1000[,2:1001]),
                                             bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_nbIntron_ZT24_1000[,2:1001])),
                                             top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_nbIntron_ZT24_1000[,2:1001])),
                                             ZT= "ZT24")


# wilcoxon test

Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(ZT24_nbIntron_HVG$Number_of_introns,random_nbIntron_ZT24_1000[,i])$p.value
}

# All values between 1.811546e-18 and 2.181257e-40 #




# all ZT combined and graph


allZT_boostrap_nbIntron <- bind_rows(HVG_LVG_distribution_nbIntron_ZT2,
                                     HVG_LVG_distribution_nbIntron_ZT4,
                                     HVG_LVG_distribution_nbIntron_ZT6,
                                     HVG_LVG_distribution_nbIntron_ZT8,
                                     HVG_LVG_distribution_nbIntron_ZT10,
                                     HVG_LVG_distribution_nbIntron_ZT12,
                                     HVG_LVG_distribution_nbIntron_ZT14,
                                     HVG_LVG_distribution_nbIntron_ZT16,
                                     HVG_LVG_distribution_nbIntron_ZT18,
                                     HVG_LVG_distribution_nbIntron_ZT20,
                                     HVG_LVG_distribution_nbIntron_ZT22,
                                     HVG_LVG_distribution_nbIntron_ZT24)



gather(allZT_boostrap_nbIntron, type, distr, HVG:top_CI, factor_key=TRUE) %>%
  ggplot(aes(x=nbIntron, y=distr, color=type, linetype=type)) +
  geom_line(size=1) +
  scale_x_continuous(limits = c(0, 50)) +
  facet_wrap(~ZT) +
  xlab("Intron number") +
  ylab("Density") +
  scale_color_manual(values=c("blue","green4","black","grey30","grey30")) +
  scale_linetype_manual(values=c("solid","solid","solid", "dashed", "dashed")) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line( size=.1, color="black" ),
        panel.background = element_rect(fill="white",colour="grey50"), 
        text=element_text(size=22),
        axis.text.x = element_text(angle=90))



### Quantitative comparison of gene length and Corrrected CV2 ----


Gene_nbIntrons <- read.table("Gene-length_intron-number/TAIR10_Gene_nbIntrons.txt", 
                             header=TRUE, sep='\t')

CorCV <- read.table("Gene-length_intron-number/BrenneckeMod_allZT_RUV2.txt",
                    header=T,sep='\t') %>%
  select(Gene, contains("Mean"), contains("BioVar"))


CorCV_gene_nbIntrons <- merge( CorCV, Gene_nbIntrons,
                               by.x="Gene",by.y="gene_name")


CorCV_gene_nbIntrons %>%
  gather("sample", "BioVar", c(BioVar_ZT2:BioVar_ZT0),
         factor_key=TRUE) %>%
  ggplot(aes(x=BioVar, y=Number_of_introns)) + 
  geom_point() + 
  facet_wrap(~ sample) +
  stat_bin2d(bins=40) +
  scale_fill_gradient(limits=c(0,100)) +
  xlim(-5,6) +
  ylim(0,90)

