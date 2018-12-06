
library(tidyverse)

### Plotting of the gene length for HVGs, LVGs and a 1000 sets of random genes ----

#data preparation

Gene_size <- read.table("Gene-length_intron-number/TAIR10_Gene_sizes.txt", header=TRUE, sep='\t')
dim(Gene_size)

Background <- read.table("Gene-length_intron-number/All_genes_background.txt", header=TRUE, sep='\t')
dim(Background)

Background_gene_size <- inner_join(Background, Gene_size, by=c("Gene"="gene_name")) %>%
  select(Gene, length)


# extract 1000 sets of random genes ----


random_distribution_length_1000=matrix(nrow=512,ncol=1000)
random_distribution_length_1000=data.frame(random_distribution_length_1000)

random_length_1000=matrix(nrow=1358,ncol=1000)
random_length_1000=data.frame(random_length_1000)


for (i in 1:1000) {
  random_set <- Background_gene_size[sample(nrow(Background_gene_size), 1358),]
  random_length_1000[,i] <- random_set[,2]
  d <- density(random_set[,2], from=153, to=26435)
  random_distribution_length_1000[,i] <- d$y
}

colnames(random_length_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(Background_gene_size$length, from=153, to=26435)
random_distribution_length_1000 <- cbind(d_background$x, random_distribution_length_1000)
colnames(random_distribution_length_1000) <- c("length",paste("random_set",1:100, sep="_"))


HVG <- read.table("Gene-length_intron-number/HVG_allZT_background.txt", header=TRUE, sep='\t')

Variable_length <- inner_join(HVG, Gene_size, by=c("allZT"="gene_name")) %>%
  select(allZT, length) 


LVG <- read.table("Gene-length_intron-number/LVG1000_allZT_background.txt", header=TRUE, sep='\t')

Stable_length <- inner_join(LVG, Gene_size, by=c("allZT"="gene_name")) %>%
  select(allZT, length) 



HVG_LVG_distribution_length=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_length=data.frame(HVG_LVG_distribution_length)

d_background <-  density(Background_gene_size$length, from=153, to=26435)
HVG_LVG_distribution_length[,1] <- d_background$x

d_HVG <- density(Variable_length$length, from=153, to=26435)
d_LVG <- density(Stable_length$length, from=153, to=26435)

HVG_LVG_distribution_length[,2] <- d_HVG$y
HVG_LVG_distribution_length[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_length) <- c("length", "HVG", "LVG")
HVG_LVG_distribution_length_long <- gather(HVG_LVG_distribution_length, type, distr, HVG:LVG)

# graph

HVG_LVG_distribution_length <- mutate(HVG_LVG_distribution_length, mean_bootstrat=rowMeans(random_distribution_length_1000[,2:1001]),
                                      bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_length_1000[,2:1001])),
                                      top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_length_1000[,2:1001])))



gather(HVG_LVG_distribution_length, type, distr, HVG:top_CI, factor_key=TRUE) %>%
  ggplot(aes(x=length, y=distr, color=type, linetype=type)) +
  geom_line(size=1) +
  scale_x_continuous(limits = c(0, 8000)) +
  xlab("Gene length (bp)") +
  ylab("Density") +
  scale_color_manual(values=c("blue","green4","black","grey30","grey30")) +
  scale_linetype_manual(values=c("solid","solid","solid", "dashed", "dashed")) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line( size=.1, color="black" ), 
        panel.background = element_rect(fill="white",colour="grey50"), 
        text=element_text(size=22))


# wilcoxon test

Wilcox_stat=matrix(nrow=1000,ncol=1)
Wilcox_stat=data.frame(Wilcox_stat)


for (i in 1:1000) {
  Wilcox_stat[i,] <- wilcox.test(Variable_length$length,random_length_1000[,i])$statistic
}

names(Wilcox_stat) <- "result"

Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(Variable_length$length,random_length_1000[,i])$p.value
}

names(Wilcox_pvalues) <- "result"

# All values between 6.344287e-59 and 9.368690e-92 #



### Plotting of the gene length of HVGs, LVGs and a 1000 sets of random genes for each ZT ----


# data preparation

Gene_size <- read.table("Gene-length_intron-number/TAIR10_Gene_sizes.txt", header=TRUE, sep='\t')
dim(Gene_size)


Background <- read.table("Gene-length_intron-number/All_genes_background.txt", header=TRUE, sep='\t')
dim(Background)

Gene_size <- inner_join(Gene_size, Background, by=c("gene_name"="Gene"))



# bootstrap for each ZT

# ZT2 ----
random_distribution_length_ZT2_1000=matrix(nrow=512,ncol=1000)
random_distribution_length_ZT2_1000=data.frame(random_distribution_length_ZT2_1000)

random_length_ZT2_1000=matrix(nrow=313,ncol=1000)
random_length_ZT2_1000=data.frame(random_length_ZT2_1000)


for (i in 1:1000) {
  random_set <- Gene_size[sample(nrow(Gene_size), 313),]
  random_length_ZT2_1000[,i] <- random_set[,6]
  d <- density(random_set[,6], from=153, to=26435, na.rm=TRUE)
  random_distribution_length_ZT2_1000[,i] <- d$y
}

colnames(random_length_ZT2_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(Gene_size$length, from=153, to=26435, na.rm=TRUE)
random_distribution_length_ZT2_1000 <- cbind(d_background$x, random_distribution_length_ZT2_1000)
colnames(random_distribution_length_ZT2_1000) <- c("length",paste("random_set",1:100, sep="_"))


ZT2 <- read.table("Gene-length_intron-number/ZT2_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT2)

ZT2_length_HVG <- merge(ZT2,Gene_size,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,length) %>%
  filter(TypeV=="Variable")

ZT2_length_LVG <- merge(ZT2,Gene_size,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,length) %>%
  filter(TypeV=="Stable")



HVG_LVG_distribution_length_ZT2=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_length_ZT2=data.frame(HVG_LVG_distribution_length_ZT2)

d_background <-  density(Gene_size$length, from=153, to=26435, na.rm=TRUE)
HVG_LVG_distribution_length_ZT2[,1] <- d_background$x

d_HVG <- density(ZT2_length_HVG$length, from=153, to=26435, na.rm=TRUE)
d_LVG <- density(ZT2_length_LVG$length, from=153, to=26435, na.rm=TRUE)

HVG_LVG_distribution_length_ZT2[,2] <- d_HVG$y
HVG_LVG_distribution_length_ZT2[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_length_ZT2) <- c("length", "HVG", "LVG")
HVG_LVG_distribution_length_ZT2_long <- gather(HVG_LVG_distribution_length_ZT2, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_length_ZT2 <- mutate(HVG_LVG_distribution_length_ZT2, mean_bootstrat=rowMeans(random_distribution_length_ZT2_1000[,2:1001]),
                                          bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_length_ZT2_1000[,2:1001])),
                                          top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_length_ZT2_1000[,2:1001])),
                                          ZT= "ZT2")


# wilcoxon test

Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(ZT2_length_HVG$length,random_length_ZT2_1000[,i])$p.value
}

# All values between 1.701999e-21 and 1.157552e-37 #




# ZT4 ----
random_distribution_length_ZT4_1000=matrix(nrow=512,ncol=1000)
random_distribution_length_ZT4_1000=data.frame(random_distribution_length_ZT4_1000)

random_length_ZT4_1000=matrix(nrow=257,ncol=1000)
random_length_ZT4_1000=data.frame(random_length_ZT4_1000)


for (i in 1:1000) {
  random_set <- Gene_size[sample(nrow(Gene_size), 257),]
  random_length_ZT4_1000[,i] <- random_set[,6]
  d <- density(random_set[,6], from=153, to=26435, na.rm=TRUE)
  random_distribution_length_ZT4_1000[,i] <- d$y
}

colnames(random_length_ZT4_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(Gene_size$length, from=153, to=26435, na.rm=TRUE)
random_distribution_length_ZT4_1000 <- cbind(d_background$x, random_distribution_length_ZT4_1000)
colnames(random_distribution_length_ZT4_1000) <- c("length",paste("random_set",1:100, sep="_"))


ZT4 <- read.table("Gene-length_intron-number/ZT4_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT4)

ZT4_length_HVG <- merge(ZT4,Gene_size,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,length) %>%
  filter(TypeV=="Variable")

ZT4_length_LVG <- merge(ZT4,Gene_size,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,length) %>%
  filter(TypeV=="Stable")



HVG_LVG_distribution_length_ZT4=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_length_ZT4=data.frame(HVG_LVG_distribution_length_ZT4)

d_background <-  density(Gene_size$length, from=153, to=26435, na.rm=TRUE)
HVG_LVG_distribution_length_ZT4[,1] <- d_background$x

d_HVG <- density(ZT4_length_HVG$length, from=153, to=26435, na.rm=TRUE)
d_LVG <- density(ZT4_length_LVG$length, from=153, to=26435, na.rm=TRUE)

HVG_LVG_distribution_length_ZT4[,2] <- d_HVG$y
HVG_LVG_distribution_length_ZT4[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_length_ZT4) <- c("length", "HVG", "LVG")
HVG_LVG_distribution_length_ZT4_long <- gather(HVG_LVG_distribution_length_ZT4, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_length_ZT4 <- mutate(HVG_LVG_distribution_length_ZT4, mean_bootstrat=rowMeans(random_distribution_length_ZT4_1000[,2:1001]),
                                          bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_length_ZT4_1000[,2:1001])),
                                          top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_length_ZT4_1000[,2:1001])),
                                          ZT= "ZT4")


# wilcoxon test

Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(ZT4_length_HVG$length,random_length_ZT4_1000[,i])$p.value
}

# All values between 1.527582e-14 and 2.505410e-32 #



# ZT6 ----
random_distribution_length_ZT6_1000=matrix(nrow=512,ncol=1000)
random_distribution_length_ZT6_1000=data.frame(random_distribution_length_ZT6_1000)

random_length_ZT6_1000=matrix(nrow=374,ncol=1000)
random_length_ZT6_1000=data.frame(random_length_ZT6_1000)


for (i in 1:1000) {
  random_set <- Gene_size[sample(nrow(Gene_size), 374),]
  random_length_ZT6_1000[,i] <- random_set[,6]
  d <- density(random_set[,6], from=153, to=26435, na.rm=TRUE)
  random_distribution_length_ZT6_1000[,i] <- d$y
}

colnames(random_length_ZT6_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(Gene_size$length, from=153, to=26435, na.rm=TRUE)
random_distribution_length_ZT6_1000 <- cbind(d_background$x, random_distribution_length_ZT6_1000)
colnames(random_distribution_length_ZT6_1000) <- c("length",paste("random_set",1:100, sep="_"))


ZT6 <- read.table("Gene-length_intron-number/ZT6_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT6)

ZT6_length_HVG <- merge(ZT6,Gene_size,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,length) %>%
  filter(TypeV=="Variable")

ZT6_length_LVG <- merge(ZT6,Gene_size,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,length) %>%
  filter(TypeV=="Stable")



HVG_LVG_distribution_length_ZT6=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_length_ZT6=data.frame(HVG_LVG_distribution_length_ZT6)

d_background <-  density(Gene_size$length, from=153, to=26435, na.rm=TRUE)
HVG_LVG_distribution_length_ZT6[,1] <- d_background$x

d_HVG <- density(ZT6_length_HVG$length, from=153, to=26435, na.rm=TRUE)
d_LVG <- density(ZT6_length_LVG$length, from=153, to=26435, na.rm=TRUE)

HVG_LVG_distribution_length_ZT6[,2] <- d_HVG$y
HVG_LVG_distribution_length_ZT6[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_length_ZT6) <- c("length", "HVG", "LVG")
HVG_LVG_distribution_length_ZT6_long <- gather(HVG_LVG_distribution_length_ZT6, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_length_ZT6 <- mutate(HVG_LVG_distribution_length_ZT6, mean_bootstrat=rowMeans(random_distribution_length_ZT6_1000[,2:1001]),
                                          bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_length_ZT6_1000[,2:1001])),
                                          top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_length_ZT6_1000[,2:1001])),
                                          ZT= "ZT6")


# wilcoxon test

Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(ZT6_length_HVG$length,random_length_ZT6_1000[,i])$p.value
}

# All values between 4.724373e-20 and 8.814677e-38 #



# ZT8 ----
random_distribution_length_ZT8_1000=matrix(nrow=512,ncol=1000)
random_distribution_length_ZT8_1000=data.frame(random_distribution_length_ZT8_1000)

random_length_ZT8_1000=matrix(nrow=327,ncol=1000)
random_length_ZT8_1000=data.frame(random_length_ZT8_1000)


for (i in 1:1000) {
  random_set <- Gene_size[sample(nrow(Gene_size), 327),]
  random_length_ZT8_1000[,i] <- random_set[,6]
  d <- density(random_set[,6], from=153, to=26435, na.rm=TRUE)
  random_distribution_length_ZT8_1000[,i] <- d$y
}

colnames(random_length_ZT8_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(Gene_size$length, from=153, to=26435, na.rm=TRUE)
random_distribution_length_ZT8_1000 <- cbind(d_background$x, random_distribution_length_ZT8_1000)
colnames(random_distribution_length_ZT8_1000) <- c("length",paste("random_set",1:100, sep="_"))


ZT8 <- read.table("Gene-length_intron-number/ZT8_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT8)

ZT8_length_HVG <- merge(ZT8,Gene_size,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,length) %>%
  filter(TypeV=="Variable")

ZT8_length_LVG <- merge(ZT8,Gene_size,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,length) %>%
  filter(TypeV=="Stable")



HVG_LVG_distribution_length_ZT8=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_length_ZT8=data.frame(HVG_LVG_distribution_length_ZT8)

d_background <-  density(Gene_size$length, from=153, to=26435, na.rm=TRUE)
HVG_LVG_distribution_length_ZT8[,1] <- d_background$x

d_HVG <- density(ZT8_length_HVG$length, from=153, to=26435, na.rm=TRUE)
d_LVG <- density(ZT8_length_LVG$length, from=153, to=26435, na.rm=TRUE)

HVG_LVG_distribution_length_ZT8[,2] <- d_HVG$y
HVG_LVG_distribution_length_ZT8[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_length_ZT8) <- c("length", "HVG", "LVG")
HVG_LVG_distribution_length_ZT8_long <- gather(HVG_LVG_distribution_length_ZT8, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_length_ZT8 <- mutate(HVG_LVG_distribution_length_ZT8, mean_bootstrat=rowMeans(random_distribution_length_ZT8_1000[,2:1001]),
                                          bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_length_ZT8_1000[,2:1001])),
                                          top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_length_ZT8_1000[,2:1001])),
                                          ZT= "ZT8")


# wilcoxon test

Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(ZT8_length_HVG$length,random_length_ZT8_1000[,i])$p.value
}

# All values between 3.849197e-19 and 1.601820e-34 #


# ZT10 ----
random_distribution_length_ZT10_1000=matrix(nrow=512,ncol=1000)
random_distribution_length_ZT10_1000=data.frame(random_distribution_length_ZT10_1000)

random_length_ZT10_1000=matrix(nrow=310,ncol=1000)
random_length_ZT10_1000=data.frame(random_length_ZT10_1000)


for (i in 1:1000) {
  random_set <- Gene_size[sample(nrow(Gene_size), 310),]
  random_length_ZT10_1000[,i] <- random_set[,6]
  d <- density(random_set[,6], from=153, to=26435, na.rm=TRUE)
  random_distribution_length_ZT10_1000[,i] <- d$y
}

colnames(random_length_ZT10_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(Gene_size$length, from=153, to=26435, na.rm=TRUE)
random_distribution_length_ZT10_1000 <- cbind(d_background$x, random_distribution_length_ZT10_1000)
colnames(random_distribution_length_ZT10_1000) <- c("length",paste("random_set",1:100, sep="_"))


ZT10 <- read.table("Gene-length_intron-number/ZT10_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT10)

ZT10_length_HVG <- merge(ZT10,Gene_size,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,length) %>%
  filter(TypeV=="Variable")

ZT10_length_LVG <- merge(ZT10,Gene_size,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,length) %>%
  filter(TypeV=="Stable")



HVG_LVG_distribution_length_ZT10=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_length_ZT10=data.frame(HVG_LVG_distribution_length_ZT10)

d_background <-  density(Gene_size$length, from=153, to=26435, na.rm=TRUE)
HVG_LVG_distribution_length_ZT10[,1] <- d_background$x

d_HVG <- density(ZT10_length_HVG$length, from=153, to=26435, na.rm=TRUE)
d_LVG <- density(ZT10_length_LVG$length, from=153, to=26435, na.rm=TRUE)

HVG_LVG_distribution_length_ZT10[,2] <- d_HVG$y
HVG_LVG_distribution_length_ZT10[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_length_ZT10) <- c("length", "HVG", "LVG")
HVG_LVG_distribution_length_ZT10_long <- gather(HVG_LVG_distribution_length_ZT10, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_length_ZT10 <- mutate(HVG_LVG_distribution_length_ZT10, mean_bootstrat=rowMeans(random_distribution_length_ZT10_1000[,2:1001]),
                                           bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_length_ZT10_1000[,2:1001])),
                                           top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_length_ZT10_1000[,2:1001])),
                                           ZT= "ZT10")


# wilcoxon test

Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(ZT10_length_HVG$length,random_length_ZT10_1000[,i])$p.value
}

# All values between 9.707502e-16 and 2.177487e-32 #



# ZT12 ----
random_distribution_length_ZT12_1000=matrix(nrow=512,ncol=1000)
random_distribution_length_ZT12_1000=data.frame(random_distribution_length_ZT12_1000)

random_length_ZT12_1000=matrix(nrow=636,ncol=1000)
random_length_ZT12_1000=data.frame(random_length_ZT12_1000)


for (i in 1:1000) {
  random_set <- Gene_size[sample(nrow(Gene_size), 636),]
  random_length_ZT12_1000[,i] <- random_set[,6]
  d <- density(random_set[,6], from=153, to=26435, na.rm=TRUE)
  random_distribution_length_ZT12_1000[,i] <- d$y
}

colnames(random_length_ZT12_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(Gene_size$length, from=153, to=26435, na.rm=TRUE)
random_distribution_length_ZT12_1000 <- cbind(d_background$x, random_distribution_length_ZT12_1000)
colnames(random_distribution_length_ZT12_1000) <- c("length",paste("random_set",1:100, sep="_"))


ZT12 <- read.table("Gene-length_intron-number/ZT12_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT12)

ZT12_length_HVG <- merge(ZT12,Gene_size,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,length) %>%
  filter(TypeV=="Variable")

ZT12_length_LVG <- merge(ZT12,Gene_size,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,length) %>%
  filter(TypeV=="Stable")



HVG_LVG_distribution_length_ZT12=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_length_ZT12=data.frame(HVG_LVG_distribution_length_ZT12)

d_background <-  density(Gene_size$length, from=153, to=26435, na.rm=TRUE)
HVG_LVG_distribution_length_ZT12[,1] <- d_background$x

d_HVG <- density(ZT12_length_HVG$length, from=153, to=26435, na.rm=TRUE)
d_LVG <- density(ZT12_length_LVG$length, from=153, to=26435, na.rm=TRUE)

HVG_LVG_distribution_length_ZT12[,2] <- d_HVG$y
HVG_LVG_distribution_length_ZT12[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_length_ZT12) <- c("length", "HVG", "LVG")
HVG_LVG_distribution_length_ZT12_long <- gather(HVG_LVG_distribution_length_ZT12, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_length_ZT12 <- mutate(HVG_LVG_distribution_length_ZT12, mean_bootstrat=rowMeans(random_distribution_length_ZT12_1000[,2:1001]),
                                           bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_length_ZT12_1000[,2:1001])),
                                           top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_length_ZT12_1000[,2:1001])),
                                           ZT= "ZT12")


# wilcoxon test

Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(ZT12_length_HVG$length,random_length_ZT12_1000[,i])$p.value
}

# All values between 2.257185e-30 and 4.065248e-53 #



# ZT14 ----
random_distribution_length_ZT14_1000=matrix(nrow=512,ncol=1000)
random_distribution_length_ZT14_1000=data.frame(random_distribution_length_ZT14_1000)

random_length_ZT14_1000=matrix(nrow=418,ncol=1000)
random_length_ZT14_1000=data.frame(random_length_ZT14_1000)


for (i in 1:1000) {
  random_set <- Gene_size[sample(nrow(Gene_size), 418),]
  random_length_ZT14_1000[,i] <- random_set[,6]
  d <- density(random_set[,6], from=153, to=26435, na.rm=TRUE)
  random_distribution_length_ZT14_1000[,i] <- d$y
}

colnames(random_length_ZT14_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(Gene_size$length, from=153, to=26435, na.rm=TRUE)
random_distribution_length_ZT14_1000 <- cbind(d_background$x, random_distribution_length_ZT14_1000)
colnames(random_distribution_length_ZT14_1000) <- c("length",paste("random_set",1:100, sep="_"))


ZT14 <- read.table("ZT14_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT14)

ZT14_length_HVG <- merge(ZT14,Gene_size,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,length) %>%
  filter(TypeV=="Variable")

ZT14_length_LVG <- merge(ZT14,Gene_size,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,length) %>%
  filter(TypeV=="Stable")



HVG_LVG_distribution_length_ZT14=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_length_ZT14=data.frame(HVG_LVG_distribution_length_ZT14)

d_background <-  density(Gene_size$length, from=153, to=26435, na.rm=TRUE)
HVG_LVG_distribution_length_ZT14[,1] <- d_background$x

d_HVG <- density(ZT14_length_HVG$length, from=153, to=26435, na.rm=TRUE)
d_LVG <- density(ZT14_length_LVG$length, from=153, to=26435, na.rm=TRUE)

HVG_LVG_distribution_length_ZT14[,2] <- d_HVG$y
HVG_LVG_distribution_length_ZT14[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_length_ZT14) <- c("length", "HVG", "LVG")
HVG_LVG_distribution_length_ZT14_long <- gather(HVG_LVG_distribution_length_ZT14, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_length_ZT14 <- mutate(HVG_LVG_distribution_length_ZT14, mean_bootstrat=rowMeans(random_distribution_length_ZT14_1000[,2:1001]),
                                           bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_length_ZT14_1000[,2:1001])),
                                           top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_length_ZT14_1000[,2:1001])),
                                           ZT= "ZT14")


# wilcoxon test

Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(ZT14_length_HVG$length,random_length_ZT14_1000[,i])$p.value
}

# All values between 2.893868e-19 and 1.006388e-36 #



# ZT16 ----
random_distribution_length_ZT16_1000=matrix(nrow=512,ncol=1000)
random_distribution_length_ZT16_1000=data.frame(random_distribution_length_ZT16_1000)

random_length_ZT16_1000=matrix(nrow=630,ncol=1000)
random_length_ZT16_1000=data.frame(random_length_ZT16_1000)


for (i in 1:1000) {
  random_set <- Gene_size[sample(nrow(Gene_size), 630),]
  random_length_ZT16_1000[,i] <- random_set[,6]
  d <- density(random_set[,6], from=153, to=26435, na.rm=TRUE)
  random_distribution_length_ZT16_1000[,i] <- d$y
}

colnames(random_length_ZT16_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(Gene_size$length, from=153, to=26435, na.rm=TRUE)
random_distribution_length_ZT16_1000 <- cbind(d_background$x, random_distribution_length_ZT16_1000)
colnames(random_distribution_length_ZT16_1000) <- c("length",paste("random_set",1:100, sep="_"))


ZT16 <- read.table("Gene-length_intron-number/ZT16_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT16)

ZT16_length_HVG <- merge(ZT16,Gene_size,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,length) %>%
  filter(TypeV=="Variable")

ZT16_length_LVG <- merge(ZT16,Gene_size,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,length) %>%
  filter(TypeV=="Stable")



HVG_LVG_distribution_length_ZT16=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_length_ZT16=data.frame(HVG_LVG_distribution_length_ZT16)

d_background <-  density(Gene_size$length, from=153, to=26435, na.rm=TRUE)
HVG_LVG_distribution_length_ZT16[,1] <- d_background$x

d_HVG <- density(ZT16_length_HVG$length, from=153, to=26435, na.rm=TRUE)
d_LVG <- density(ZT16_length_LVG$length, from=153, to=26435, na.rm=TRUE)

HVG_LVG_distribution_length_ZT16[,2] <- d_HVG$y
HVG_LVG_distribution_length_ZT16[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_length_ZT16) <- c("length", "HVG", "LVG")
HVG_LVG_distribution_length_ZT16_long <- gather(HVG_LVG_distribution_length_ZT16, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_length_ZT16 <- mutate(HVG_LVG_distribution_length_ZT16, mean_bootstrat=rowMeans(random_distribution_length_ZT16_1000[,2:1001]),
                                           bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_length_ZT16_1000[,2:1001])),
                                           top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_length_ZT16_1000[,2:1001])),
                                           ZT= "ZT16")


# wilcoxon test

Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(ZT16_length_HVG$length,random_length_ZT16_1000[,i])$p.value
}

# All values between 1.788187e-26 and 5.163659e-44 #



# ZT18 ----
random_distribution_length_ZT18_1000=matrix(nrow=512,ncol=1000)
random_distribution_length_ZT18_1000=data.frame(random_distribution_length_ZT18_1000)

random_length_ZT18_1000=matrix(nrow=408,ncol=1000)
random_length_ZT18_1000=data.frame(random_length_ZT18_1000)


for (i in 1:1000) {
  random_set <- Gene_size[sample(nrow(Gene_size), 408),]
  random_length_ZT18_1000[,i] <- random_set[,6]
  d <- density(random_set[,6], from=153, to=26435, na.rm=TRUE)
  random_distribution_length_ZT18_1000[,i] <- d$y
}

colnames(random_length_ZT18_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(Gene_size$length, from=153, to=26435, na.rm=TRUE)
random_distribution_length_ZT18_1000 <- cbind(d_background$x, random_distribution_length_ZT18_1000)
colnames(random_distribution_length_ZT18_1000) <- c("length",paste("random_set",1:100, sep="_"))


ZT18 <- read.table("Gene-length_intron-number/ZT18_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT18)

ZT18_length_HVG <- merge(ZT18,Gene_size,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,length) %>%
  filter(TypeV=="Variable")

ZT18_length_LVG <- merge(ZT18,Gene_size,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,length) %>%
  filter(TypeV=="Stable")



HVG_LVG_distribution_length_ZT18=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_length_ZT18=data.frame(HVG_LVG_distribution_length_ZT18)

d_background <-  density(Gene_size$length, from=153, to=26435, na.rm=TRUE)
HVG_LVG_distribution_length_ZT18[,1] <- d_background$x

d_HVG <- density(ZT18_length_HVG$length, from=153, to=26435, na.rm=TRUE)
d_LVG <- density(ZT18_length_LVG$length, from=153, to=26435, na.rm=TRUE)

HVG_LVG_distribution_length_ZT18[,2] <- d_HVG$y
HVG_LVG_distribution_length_ZT18[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_length_ZT18) <- c("length", "HVG", "LVG")
HVG_LVG_distribution_length_ZT18_long <- gather(HVG_LVG_distribution_length_ZT18, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_length_ZT18 <- mutate(HVG_LVG_distribution_length_ZT18, mean_bootstrat=rowMeans(random_distribution_length_ZT18_1000[,2:1001]),
                                           bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_length_ZT18_1000[,2:1001])),
                                           top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_length_ZT18_1000[,2:1001])),
                                           ZT= "ZT18")


# wilcoxon test

Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(ZT18_length_HVG$length,random_length_ZT18_1000[,i])$p.value
}

# All values between 3.562388e-19 and 2.728502e-37 #


# ZT20 ----
random_distribution_length_ZT20_1000=matrix(nrow=512,ncol=1000)
random_distribution_length_ZT20_1000=data.frame(random_distribution_length_ZT20_1000)

random_length_ZT20_1000=matrix(nrow=716,ncol=1000)
random_length_ZT20_1000=data.frame(random_length_ZT20_1000)


for (i in 1:1000) {
  random_set <- Gene_size[sample(nrow(Gene_size), 716),]
  random_length_ZT20_1000[,i] <- random_set[,6]
  d <- density(random_set[,6], from=153, to=26435, na.rm=TRUE)
  random_distribution_length_ZT20_1000[,i] <- d$y
}

colnames(random_length_ZT20_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(Gene_size$length, from=153, to=26435, na.rm=TRUE)
random_distribution_length_ZT20_1000 <- cbind(d_background$x, random_distribution_length_ZT20_1000)
colnames(random_distribution_length_ZT20_1000) <- c("length",paste("random_set",1:100, sep="_"))


ZT20 <- read.table("Gene-length_intron-number/ZT20_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT20)

ZT20_length_HVG <- merge(ZT20,Gene_size,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,length) %>%
  filter(TypeV=="Variable")

ZT20_length_LVG <- merge(ZT20,Gene_size,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,length) %>%
  filter(TypeV=="Stable")



HVG_LVG_distribution_length_ZT20=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_length_ZT20=data.frame(HVG_LVG_distribution_length_ZT20)

d_background <-  density(Gene_size$length, from=153, to=26435, na.rm=TRUE)
HVG_LVG_distribution_length_ZT20[,1] <- d_background$x

d_HVG <- density(ZT20_length_HVG$length, from=153, to=26435, na.rm=TRUE)
d_LVG <- density(ZT20_length_LVG$length, from=153, to=26435, na.rm=TRUE)

HVG_LVG_distribution_length_ZT20[,2] <- d_HVG$y
HVG_LVG_distribution_length_ZT20[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_length_ZT20) <- c("length", "HVG", "LVG")
HVG_LVG_distribution_length_ZT20_long <- gather(HVG_LVG_distribution_length_ZT20, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_length_ZT20 <- mutate(HVG_LVG_distribution_length_ZT20, mean_bootstrat=rowMeans(random_distribution_length_ZT20_1000[,2:1001]),
                                           bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_length_ZT20_1000[,2:1001])),
                                           top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_length_ZT20_1000[,2:1001])),
                                           ZT= "ZT20")


# wilcoxon test

Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(ZT20_length_HVG$length,random_length_ZT20_1000[,i])$p.value
}

# All values between 7.763010e-34 and 1.714399e-56 #


# ZT22 ----
random_distribution_length_ZT22_1000=matrix(nrow=512,ncol=1000)
random_distribution_length_ZT22_1000=data.frame(random_distribution_length_ZT22_1000)

random_length_ZT22_1000=matrix(nrow=566,ncol=1000)
random_length_ZT22_1000=data.frame(random_length_ZT22_1000)


for (i in 1:1000) {
  random_set <- Gene_size[sample(nrow(Gene_size), 566),]
  random_length_ZT22_1000[,i] <- random_set[,6]
  d <- density(random_set[,6], from=153, to=26435, na.rm=TRUE)
  random_distribution_length_ZT22_1000[,i] <- d$y
}

colnames(random_length_ZT22_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(Gene_size$length, from=153, to=26435, na.rm=TRUE)
random_distribution_length_ZT22_1000 <- cbind(d_background$x, random_distribution_length_ZT22_1000)
colnames(random_distribution_length_ZT22_1000) <- c("length",paste("random_set",1:100, sep="_"))


ZT22 <- read.table("Gene-length_intron-number/ZT22_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT22)

ZT22_length_HVG <- merge(ZT22,Gene_size,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,length) %>%
  filter(TypeV=="Variable")

ZT22_length_LVG <- merge(ZT22,Gene_size,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,length) %>%
  filter(TypeV=="Stable")



HVG_LVG_distribution_length_ZT22=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_length_ZT22=data.frame(HVG_LVG_distribution_length_ZT22)

d_background <-  density(Gene_size$length, from=153, to=26435, na.rm=TRUE)
HVG_LVG_distribution_length_ZT22[,1] <- d_background$x

d_HVG <- density(ZT22_length_HVG$length, from=153, to=26435, na.rm=TRUE)
d_LVG <- density(ZT22_length_LVG$length, from=153, to=26435, na.rm=TRUE)

HVG_LVG_distribution_length_ZT22[,2] <- d_HVG$y
HVG_LVG_distribution_length_ZT22[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_length_ZT22) <- c("length", "HVG", "LVG")
HVG_LVG_distribution_length_ZT22_long <- gather(HVG_LVG_distribution_length_ZT22, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_length_ZT22 <- mutate(HVG_LVG_distribution_length_ZT22, mean_bootstrat=rowMeans(random_distribution_length_ZT22_1000[,2:1001]),
                                           bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_length_ZT22_1000[,2:1001])),
                                           top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_length_ZT22_1000[,2:1001])),
                                           ZT= "ZT22")


# wilcoxon test

Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(ZT22_length_HVG$length,random_length_ZT22_1000[,i])$p.value
}

# All values between 6.991256e-38 and 1.498199e-59 #




# ZT24 ----
random_distribution_length_ZT24_1000=matrix(nrow=512,ncol=1000)
random_distribution_length_ZT24_1000=data.frame(random_distribution_length_ZT24_1000)

random_length_ZT24_1000=matrix(nrow=377,ncol=1000)
random_length_ZT24_1000=data.frame(random_length_ZT24_1000)


for (i in 1:1000) {
  random_set <- Gene_size[sample(nrow(Gene_size), 377),]
  random_length_ZT24_1000[,i] <- random_set[,6]
  d <- density(random_set[,6], from=153, to=26435, na.rm=TRUE)
  random_distribution_length_ZT24_1000[,i] <- d$y
}

colnames(random_length_ZT24_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(Gene_size$length, from=153, to=26435, na.rm=TRUE)
random_distribution_length_ZT24_1000 <- cbind(d_background$x, random_distribution_length_ZT24_1000)
colnames(random_distribution_length_ZT24_1000) <- c("length",paste("random_set",1:100, sep="_"))


ZT24 <- read.table("Gene-length_intron-number/ZT0_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT24)

ZT24_length_HVG <- merge(ZT24,Gene_size,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,length) %>%
  filter(TypeV=="Variable")

ZT24_length_LVG <- merge(ZT24,Gene_size,by.x="Gene",by.y="gene_name") %>%
  select(Gene,TypeV,ZT,length) %>%
  filter(TypeV=="Stable")



HVG_LVG_distribution_length_ZT24=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_length_ZT24=data.frame(HVG_LVG_distribution_length_ZT24)

d_background <-  density(Gene_size$length, from=153, to=26435, na.rm=TRUE)
HVG_LVG_distribution_length_ZT24[,1] <- d_background$x

d_HVG <- density(ZT24_length_HVG$length, from=153, to=26435, na.rm=TRUE)
d_LVG <- density(ZT24_length_LVG$length, from=153, to=26435, na.rm=TRUE)

HVG_LVG_distribution_length_ZT24[,2] <- d_HVG$y
HVG_LVG_distribution_length_ZT24[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_length_ZT24) <- c("length", "HVG", "LVG")
HVG_LVG_distribution_length_ZT24_long <- gather(HVG_LVG_distribution_length_ZT24, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_length_ZT24 <- mutate(HVG_LVG_distribution_length_ZT24, mean_bootstrat=rowMeans(random_distribution_length_ZT24_1000[,2:1001]),
                                           bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_length_ZT24_1000[,2:1001])),
                                           top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_length_ZT24_1000[,2:1001])),
                                           ZT= "ZT24")


# wilcoxon test

Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(ZT24_length_HVG$length,random_length_ZT24_1000[,i])$p.value
}

# All values between 1.590812e-29 and 6.277734e-48 #




# all ZT combined and graph


allZT_boostrap_length <- bind_rows(HVG_LVG_distribution_length_ZT2,
                                   HVG_LVG_distribution_length_ZT4,
                                   HVG_LVG_distribution_length_ZT6,
                                   HVG_LVG_distribution_length_ZT8,
                                   HVG_LVG_distribution_length_ZT10,
                                   HVG_LVG_distribution_length_ZT12,
                                   HVG_LVG_distribution_length_ZT14,
                                   HVG_LVG_distribution_length_ZT16,
                                   HVG_LVG_distribution_length_ZT18,
                                   HVG_LVG_distribution_length_ZT20,
                                   HVG_LVG_distribution_length_ZT22,
                                   HVG_LVG_distribution_length_ZT24)



gather(allZT_boostrap_length, type, distr, HVG:top_CI, factor_key=TRUE) %>%
  ggplot(aes(x=length, y=distr, color=type, linetype=type)) +
  geom_line(size=1) +
  scale_x_continuous(limits = c(0, 8000)) +
  facet_wrap(~ZT) +
  xlab("gene length (bp)") +
  ylab("Density") +
  scale_color_manual(values=c("blue","green4","black","grey30","grey30")) +
  scale_linetype_manual(values=c("solid","solid","solid", "dashed", "dashed")) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line( size=.1, color="black" ),
        panel.background = element_rect(fill="white",colour="grey50"), 
        text=element_text(size=22),
        axis.text.x = element_text(angle=90))



### Quantitative comparison of gene length and Corrrected CV2 ----

Gene_size <- read.table("Gene-length_intron-number/TAIR10_Gene_sizes.txt", header=TRUE, sep='\t')

CorCV <- read.table("Gene-length_intron-number/BrenneckeMod_allZT_RUV2.txt",header=T,sep='\t') %>%
  select(Gene, contains("Mean"), contains("BioVar"))

CorCV_gene_size <- merge( CorCV, Gene_size,
                          by.x="Gene",by.y="gene_name")

CorCV_gene_size %>%
  gather("sample", "BioVar", c(BioVar_ZT2:BioVar_ZT0),
         factor_key=TRUE) %>%
  ggplot(aes(x=BioVar, y=length)) + 
  geom_point() + 
  facet_wrap(~ sample) +
  stat_bin2d(bins=50) +
  scale_fill_gradient(limits=c(0,140)) +
  xlim(-5,6) +
  ylim(0,10000)






