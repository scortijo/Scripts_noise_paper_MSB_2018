
library(tidyverse)

# Distribution of CV2 at each ZT for HVGs, LVGs and the 1000 sets of random genes

# data preparation ----

CorCV2 <- read.table("Distribution_CV2_HVG_LVG_random/CorCV2_allZT_allGenes.txt",header=TRUE,sep='\t')

CorCV2_allZT <- gather(CorCV2, ZT, CorCV2, -tracking_id)

# bootstrap for each ZT ----
# ZT2 ----
random_distribution_CorCV2_ZT2_1000=matrix(nrow=512,ncol=1000)
random_distribution_CorCV2_ZT2_1000=data.frame(random_distribution_CorCV2_ZT2_1000)

random_CorCV2_ZT2_1000=matrix(nrow=313,ncol=1000)
random_CorCV2_ZT2_1000=data.frame(random_CorCV2_ZT2_1000)


for (i in 1:1000) {
  random_set <- CorCV2[sample(nrow(CorCV2), 313),]
  random_CorCV2_ZT2_1000[,i] <- random_set[,2]
  d <- density(random_set[,2], from=-6.5, to=8.1, na.rm=TRUE)
  random_distribution_CorCV2_ZT2_1000[,i] <- d$y
}

colnames(random_CorCV2_ZT2_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(CorCV2_allZT$CorCV2, from=-6.5, to=8.1, na.rm=TRUE)
random_distribution_CorCV2_ZT2_1000 <- cbind(d_background$x, random_distribution_CorCV2_ZT2_1000)
colnames(random_distribution_CorCV2_ZT2_1000) <- c("CorCV2",paste("random_set",1:100, sep="_"))


ZT2 <- read.table("Distribution_CV2_HVG_LVG_random/ZT2_HVG_LVG_random.txt",header=TRUE,sep='\t')

ZT2_CorCV2_HVG <- merge(ZT2,CorCV2,by.x="Gene",by.y="tracking_id") %>%
  select(Gene,TypeV,ZT,BioVar_ZT2) %>%
  filter(TypeV=="Variable")
names(ZT2_CorCV2_HVG) <-c("Gene","TypeV","ZT","BioVar_correspondingZT")

ZT2_CorCV2_LVG <- merge(ZT2,CorCV2,by.x="Gene",by.y="tracking_id") %>%
  select(Gene,TypeV,ZT,BioVar_ZT2) %>%
  filter(TypeV=="Stable")
names(ZT2_CorCV2_LVG) <-c("Gene","TypeV","ZT","BioVar_correspondingZT")


HVG_LVG_distribution_CorCV2_ZT2=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_CorCV2_ZT2=data.frame(HVG_LVG_distribution_CorCV2_ZT2)

d_background <-  density(CorCV2_allZT$CorCV2, from=-6.5, to=8.1, na.rm=TRUE)
HVG_LVG_distribution_CorCV2_ZT2[,1] <- d_background$x

d_HVG <- density(ZT2_CorCV2_HVG$BioVar_correspondingZT, from=-6.5, to=8.1, na.rm=TRUE)
d_LVG <- density(ZT2_CorCV2_LVG$BioVar_correspondingZT, from=-6.5, to=8.1, na.rm=TRUE)

HVG_LVG_distribution_CorCV2_ZT2[,2] <- d_HVG$y
HVG_LVG_distribution_CorCV2_ZT2[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_CorCV2_ZT2) <- c("CorCV2", "HVG", "LVG")
HVG_LVG_distribution_CorCV2_ZT2_long <- gather(HVG_LVG_distribution_CorCV2_ZT2, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_CorCV2_ZT2 <- mutate(HVG_LVG_distribution_CorCV2_ZT2, mean_bootstrat=rowMeans(random_distribution_CorCV2_ZT2_1000[,2:1001]),
                                          bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_CorCV2_ZT2_1000[,2:1001])),
                                          top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_CorCV2_ZT2_1000[,2:1001])),
                                          ZT= "ZT2")


# ZT4 ----
random_distribution_CorCV2_ZT4_1000=matrix(nrow=512,ncol=1000)
random_distribution_CorCV2_ZT4_1000=data.frame(random_distribution_CorCV2_ZT4_1000)

random_CorCV2_ZT4_1000=matrix(nrow=257,ncol=1000)
random_CorCV2_ZT4_1000=data.frame(random_CorCV2_ZT4_1000)


for (i in 1:1000) {
  random_set <- CorCV2[sample(nrow(CorCV2), 257),]
  random_CorCV2_ZT4_1000[,i] <- random_set[,3]
  d <- density(random_set[,3], from=-6.5, to=8.1, na.rm=TRUE)
  random_distribution_CorCV2_ZT4_1000[,i] <- d$y
}

colnames(random_CorCV2_ZT4_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(CorCV2_allZT$CorCV2, from=-6.5, to=8.1, na.rm=TRUE)
random_distribution_CorCV2_ZT4_1000 <- cbind(d_background$x, random_distribution_CorCV2_ZT4_1000)
colnames(random_distribution_CorCV2_ZT4_1000) <- c("CorCV2",paste("random_set",1:100, sep="_"))


ZT4 <- read.table("Distribution_CV2_HVG_LVG_random/ZT4_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT4)

ZT4_CorCV2_HVG <- merge(ZT4,CorCV2,by.x="Gene",by.y="tracking_id") %>%
  select(Gene,TypeV,ZT,BioVar_ZT4) %>%
  filter(TypeV=="Variable")
names(ZT4_CorCV2_HVG) <-c("Gene","TypeV","ZT","BioVar_correspondingZT")

ZT4_CorCV2_LVG <- merge(ZT4,CorCV2,by.x="Gene",by.y="tracking_id") %>%
  select(Gene,TypeV,ZT,BioVar_ZT4) %>%
  filter(TypeV=="Stable")
names(ZT4_CorCV2_LVG) <-c("Gene","TypeV","ZT","BioVar_correspondingZT")


HVG_LVG_distribution_CorCV2_ZT4=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_CorCV2_ZT4=data.frame(HVG_LVG_distribution_CorCV2_ZT4)

d_background <-  density(CorCV2_allZT$CorCV2, from=-6.5, to=8.1, na.rm=TRUE)
HVG_LVG_distribution_CorCV2_ZT4[,1] <- d_background$x

d_HVG <- density(ZT4_CorCV2_HVG$BioVar_correspondingZT, from=-6.5, to=8.1, na.rm=TRUE)
d_LVG <- density(ZT4_CorCV2_LVG$BioVar_correspondingZT, from=-6.5, to=8.1, na.rm=TRUE)

HVG_LVG_distribution_CorCV2_ZT4[,2] <- d_HVG$y
HVG_LVG_distribution_CorCV2_ZT4[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_CorCV2_ZT4) <- c("CorCV2", "HVG", "LVG")
HVG_LVG_distribution_CorCV2_ZT4_long <- gather(HVG_LVG_distribution_CorCV2_ZT4, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_CorCV2_ZT4 <- mutate(HVG_LVG_distribution_CorCV2_ZT4, mean_bootstrat=rowMeans(random_distribution_CorCV2_ZT4_1000[,2:1001]),
                                          bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_CorCV2_ZT4_1000[,2:1001])),
                                          top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_CorCV2_ZT4_1000[,2:1001])),
                                          ZT= "ZT4")


# ZT6 ----
random_distribution_CorCV2_ZT6_1000=matrix(nrow=512,ncol=1000)
random_distribution_CorCV2_ZT6_1000=data.frame(random_distribution_CorCV2_ZT6_1000)

random_CorCV2_ZT6_1000=matrix(nrow=374,ncol=1000)
random_CorCV2_ZT6_1000=data.frame(random_CorCV2_ZT6_1000)


for (i in 1:1000) {
  random_set <- CorCV2[sample(nrow(CorCV2), 374),]
  random_CorCV2_ZT6_1000[,i] <- random_set[,4]
  d <- density(random_set[,4], from=-6.5, to=8.1, na.rm=TRUE)
  random_distribution_CorCV2_ZT6_1000[,i] <- d$y
}

colnames(random_CorCV2_ZT6_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(CorCV2_allZT$CorCV2, from=-6.5, to=8.1, na.rm=TRUE)
random_distribution_CorCV2_ZT6_1000 <- cbind(d_background$x, random_distribution_CorCV2_ZT6_1000)
colnames(random_distribution_CorCV2_ZT6_1000) <- c("CorCV2",paste("random_set",1:100, sep="_"))


ZT6 <- read.table("Distribution_CV2_HVG_LVG_random/ZT6_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT6)

ZT6_CorCV2_HVG <- merge(ZT6,CorCV2,by.x="Gene",by.y="tracking_id") %>%
  select(Gene,TypeV,ZT,BioVar_ZT6) %>%
  filter(TypeV=="Variable")
names(ZT6_CorCV2_HVG) <-c("Gene","TypeV","ZT","BioVar_correspondingZT")

ZT6_CorCV2_LVG <- merge(ZT6,CorCV2,by.x="Gene",by.y="tracking_id") %>%
  select(Gene,TypeV,ZT,BioVar_ZT6) %>%
  filter(TypeV=="Stable")
names(ZT6_CorCV2_LVG) <-c("Gene","TypeV","ZT","BioVar_correspondingZT")


HVG_LVG_distribution_CorCV2_ZT6=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_CorCV2_ZT6=data.frame(HVG_LVG_distribution_CorCV2_ZT6)

d_background <-  density(CorCV2_allZT$CorCV2, from=-6.5, to=8.1, na.rm=TRUE)
HVG_LVG_distribution_CorCV2_ZT6[,1] <- d_background$x

d_HVG <- density(ZT6_CorCV2_HVG$BioVar_correspondingZT, from=-6.5, to=8.1, na.rm=TRUE)
d_LVG <- density(ZT6_CorCV2_LVG$BioVar_correspondingZT, from=-6.5, to=8.1, na.rm=TRUE)

HVG_LVG_distribution_CorCV2_ZT6[,2] <- d_HVG$y
HVG_LVG_distribution_CorCV2_ZT6[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_CorCV2_ZT6) <- c("CorCV2", "HVG", "LVG")
HVG_LVG_distribution_CorCV2_ZT6_long <- gather(HVG_LVG_distribution_CorCV2_ZT6, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_CorCV2_ZT6 <- mutate(HVG_LVG_distribution_CorCV2_ZT6, mean_bootstrat=rowMeans(random_distribution_CorCV2_ZT6_1000[,2:1001]),
                                          bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_CorCV2_ZT6_1000[,2:1001])),
                                          top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_CorCV2_ZT6_1000[,2:1001])),
                                          ZT= "ZT6")



# ZT8 ----
random_distribution_CorCV2_ZT8_1000=matrix(nrow=512,ncol=1000)
random_distribution_CorCV2_ZT8_1000=data.frame(random_distribution_CorCV2_ZT8_1000)

random_CorCV2_ZT8_1000=matrix(nrow=327,ncol=1000)
random_CorCV2_ZT8_1000=data.frame(random_CorCV2_ZT8_1000)


for (i in 1:1000) {
  random_set <- CorCV2[sample(nrow(CorCV2), 327),]
  random_CorCV2_ZT8_1000[,i] <- random_set[,5]
  d <- density(random_set[,5], from=-6.5, to=8.1, na.rm=TRUE)
  random_distribution_CorCV2_ZT8_1000[,i] <- d$y
}

colnames(random_CorCV2_ZT8_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(CorCV2_allZT$CorCV2, from=-6.5, to=8.1, na.rm=TRUE)
random_distribution_CorCV2_ZT8_1000 <- cbind(d_background$x, random_distribution_CorCV2_ZT8_1000)
colnames(random_distribution_CorCV2_ZT8_1000) <- c("CorCV2",paste("random_set",1:100, sep="_"))


ZT8 <- read.table("Distribution_CV2_HVG_LVG_random/ZT8_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT8)

ZT8_CorCV2_HVG <- merge(ZT8,CorCV2,by.x="Gene",by.y="tracking_id") %>%
  select(Gene,TypeV,ZT,BioVar_ZT8) %>%
  filter(TypeV=="Variable")
names(ZT8_CorCV2_HVG) <-c("Gene","TypeV","ZT","BioVar_correspondingZT")

ZT8_CorCV2_LVG <- merge(ZT8,CorCV2,by.x="Gene",by.y="tracking_id") %>%
  select(Gene,TypeV,ZT,BioVar_ZT8) %>%
  filter(TypeV=="Stable")
names(ZT8_CorCV2_LVG) <-c("Gene","TypeV","ZT","BioVar_correspondingZT")


HVG_LVG_distribution_CorCV2_ZT8=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_CorCV2_ZT8=data.frame(HVG_LVG_distribution_CorCV2_ZT8)

d_background <-  density(CorCV2_allZT$CorCV2, from=-6.5, to=8.1, na.rm=TRUE)
HVG_LVG_distribution_CorCV2_ZT8[,1] <- d_background$x

d_HVG <- density(ZT8_CorCV2_HVG$BioVar_correspondingZT, from=-6.5, to=8.1, na.rm=TRUE)
d_LVG <- density(ZT8_CorCV2_LVG$BioVar_correspondingZT, from=-6.5, to=8.1, na.rm=TRUE)

HVG_LVG_distribution_CorCV2_ZT8[,2] <- d_HVG$y
HVG_LVG_distribution_CorCV2_ZT8[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_CorCV2_ZT8) <- c("CorCV2", "HVG", "LVG")
HVG_LVG_distribution_CorCV2_ZT8_long <- gather(HVG_LVG_distribution_CorCV2_ZT8, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_CorCV2_ZT8 <- mutate(HVG_LVG_distribution_CorCV2_ZT8, mean_bootstrat=rowMeans(random_distribution_CorCV2_ZT8_1000[,2:1001]),
                                          bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_CorCV2_ZT8_1000[,2:1001])),
                                          top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_CorCV2_ZT8_1000[,2:1001])),
                                          ZT= "ZT8")



# ZT10 ----
random_distribution_CorCV2_ZT10_1000=matrix(nrow=512,ncol=1000)
random_distribution_CorCV2_ZT10_1000=data.frame(random_distribution_CorCV2_ZT10_1000)

random_CorCV2_ZT10_1000=matrix(nrow=310,ncol=1000)
random_CorCV2_ZT10_1000=data.frame(random_CorCV2_ZT10_1000)


for (i in 1:1000) {
  random_set <- CorCV2[sample(nrow(CorCV2), 310),]
  random_CorCV2_ZT10_1000[,i] <- random_set[,6]
  d <- density(random_set[,6], from=-6.5, to=8.1, na.rm=TRUE)
  random_distribution_CorCV2_ZT10_1000[,i] <- d$y
}

colnames(random_CorCV2_ZT10_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(CorCV2_allZT$CorCV2, from=-6.5, to=8.1, na.rm=TRUE)
random_distribution_CorCV2_ZT10_1000 <- cbind(d_background$x, random_distribution_CorCV2_ZT10_1000)
colnames(random_distribution_CorCV2_ZT10_1000) <- c("CorCV2",paste("random_set",1:100, sep="_"))


ZT10 <- read.table("Distribution_CV2_HVG_LVG_random/ZT10_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT10)

ZT10_CorCV2_HVG <- merge(ZT10,CorCV2,by.x="Gene",by.y="tracking_id") %>%
  select(Gene,TypeV,ZT,BioVar_ZT10) %>%
  filter(TypeV=="Variable")
names(ZT10_CorCV2_HVG) <-c("Gene","TypeV","ZT","BioVar_correspondingZT")

ZT10_CorCV2_LVG <- merge(ZT10,CorCV2,by.x="Gene",by.y="tracking_id") %>%
  select(Gene,TypeV,ZT,BioVar_ZT10) %>%
  filter(TypeV=="Stable")
names(ZT10_CorCV2_LVG) <-c("Gene","TypeV","ZT","BioVar_correspondingZT")


HVG_LVG_distribution_CorCV2_ZT10=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_CorCV2_ZT10=data.frame(HVG_LVG_distribution_CorCV2_ZT10)

d_background <-  density(CorCV2_allZT$CorCV2, from=-6.5, to=8.1, na.rm=TRUE)
HVG_LVG_distribution_CorCV2_ZT10[,1] <- d_background$x

d_HVG <- density(ZT10_CorCV2_HVG$BioVar_correspondingZT, from=-6.5, to=8.1, na.rm=TRUE)
d_LVG <- density(ZT10_CorCV2_LVG$BioVar_correspondingZT, from=-6.5, to=8.1, na.rm=TRUE)

HVG_LVG_distribution_CorCV2_ZT10[,2] <- d_HVG$y
HVG_LVG_distribution_CorCV2_ZT10[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_CorCV2_ZT10) <- c("CorCV2", "HVG", "LVG")
HVG_LVG_distribution_CorCV2_ZT10_long <- gather(HVG_LVG_distribution_CorCV2_ZT10, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_CorCV2_ZT10 <- mutate(HVG_LVG_distribution_CorCV2_ZT10, mean_bootstrat=rowMeans(random_distribution_CorCV2_ZT10_1000[,2:1001]),
                                           bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_CorCV2_ZT10_1000[,2:1001])),
                                           top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_CorCV2_ZT10_1000[,2:1001])),
                                           ZT= "ZT10")



# ZT12 ----
random_distribution_CorCV2_ZT12_1000=matrix(nrow=512,ncol=1000)
random_distribution_CorCV2_ZT12_1000=data.frame(random_distribution_CorCV2_ZT12_1000)

random_CorCV2_ZT12_1000=matrix(nrow=636,ncol=1000)
random_CorCV2_ZT12_1000=data.frame(random_CorCV2_ZT12_1000)


for (i in 1:1000) {
  random_set <- CorCV2[sample(nrow(CorCV2), 636),]
  random_CorCV2_ZT12_1000[,i] <- random_set[,7]
  d <- density(random_set[,7], from=-6.5, to=8.1, na.rm=TRUE)
  random_distribution_CorCV2_ZT12_1000[,i] <- d$y
}

colnames(random_CorCV2_ZT12_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(CorCV2_allZT$CorCV2, from=-6.5, to=8.1, na.rm=TRUE)
random_distribution_CorCV2_ZT12_1000 <- cbind(d_background$x, random_distribution_CorCV2_ZT12_1000)
colnames(random_distribution_CorCV2_ZT12_1000) <- c("CorCV2",paste("random_set",1:100, sep="_"))


ZT12 <- read.table("Distribution_CV2_HVG_LVG_random/ZT12_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT12)

ZT12_CorCV2_HVG <- merge(ZT12,CorCV2,by.x="Gene",by.y="tracking_id") %>%
  select(Gene,TypeV,ZT,BioVar_ZT12) %>%
  filter(TypeV=="Variable")
names(ZT12_CorCV2_HVG) <-c("Gene","TypeV","ZT","BioVar_correspondingZT")

ZT12_CorCV2_LVG <- merge(ZT12,CorCV2,by.x="Gene",by.y="tracking_id") %>%
  select(Gene,TypeV,ZT,BioVar_ZT12) %>%
  filter(TypeV=="Stable")
names(ZT12_CorCV2_LVG) <-c("Gene","TypeV","ZT","BioVar_correspondingZT")


HVG_LVG_distribution_CorCV2_ZT12=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_CorCV2_ZT12=data.frame(HVG_LVG_distribution_CorCV2_ZT12)

d_background <-  density(CorCV2_allZT$CorCV2, from=-6.5, to=8.1, na.rm=TRUE)
HVG_LVG_distribution_CorCV2_ZT12[,1] <- d_background$x

d_HVG <- density(ZT12_CorCV2_HVG$BioVar_correspondingZT, from=-6.5, to=8.1, na.rm=TRUE)
d_LVG <- density(ZT12_CorCV2_LVG$BioVar_correspondingZT, from=-6.5, to=8.1, na.rm=TRUE)

HVG_LVG_distribution_CorCV2_ZT12[,2] <- d_HVG$y
HVG_LVG_distribution_CorCV2_ZT12[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_CorCV2_ZT12) <- c("CorCV2", "HVG", "LVG")
HVG_LVG_distribution_CorCV2_ZT12_long <- gather(HVG_LVG_distribution_CorCV2_ZT12, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_CorCV2_ZT12 <- mutate(HVG_LVG_distribution_CorCV2_ZT12, mean_bootstrat=rowMeans(random_distribution_CorCV2_ZT12_1000[,2:1001]),
                                           bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_CorCV2_ZT12_1000[,2:1001])),
                                           top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_CorCV2_ZT12_1000[,2:1001])),
                                           ZT= "ZT12")


# ZT14 ----
random_distribution_CorCV2_ZT14_1000=matrix(nrow=512,ncol=1000)
random_distribution_CorCV2_ZT14_1000=data.frame(random_distribution_CorCV2_ZT14_1000)

random_CorCV2_ZT14_1000=matrix(nrow=418,ncol=1000)
random_CorCV2_ZT14_1000=data.frame(random_CorCV2_ZT14_1000)


for (i in 1:1000) {
  random_set <- CorCV2[sample(nrow(CorCV2), 418),]
  random_CorCV2_ZT14_1000[,i] <- random_set[,8]
  d <- density(random_set[,8], from=-6.5, to=8.1, na.rm=TRUE)
  random_distribution_CorCV2_ZT14_1000[,i] <- d$y
}

colnames(random_CorCV2_ZT14_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(CorCV2_allZT$CorCV2, from=-6.5, to=8.1, na.rm=TRUE)
random_distribution_CorCV2_ZT14_1000 <- cbind(d_background$x, random_distribution_CorCV2_ZT14_1000)
colnames(random_distribution_CorCV2_ZT14_1000) <- c("CorCV2",paste("random_set",1:100, sep="_"))


ZT14 <- read.table("Distribution_CV2_HVG_LVG_random/ZT14_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT14)

ZT14_CorCV2_HVG <- merge(ZT14,CorCV2,by.x="Gene",by.y="tracking_id") %>%
  select(Gene,TypeV,ZT,BioVar_ZT14) %>%
  filter(TypeV=="Variable")
names(ZT14_CorCV2_HVG) <-c("Gene","TypeV","ZT","BioVar_correspondingZT")

ZT14_CorCV2_LVG <- merge(ZT14,CorCV2,by.x="Gene",by.y="tracking_id") %>%
  select(Gene,TypeV,ZT,BioVar_ZT14) %>%
  filter(TypeV=="Stable")
names(ZT14_CorCV2_LVG) <-c("Gene","TypeV","ZT","BioVar_correspondingZT")


HVG_LVG_distribution_CorCV2_ZT14=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_CorCV2_ZT14=data.frame(HVG_LVG_distribution_CorCV2_ZT14)

d_background <-  density(CorCV2_allZT$CorCV2, from=-6.5, to=8.1, na.rm=TRUE)
HVG_LVG_distribution_CorCV2_ZT14[,1] <- d_background$x

d_HVG <- density(ZT14_CorCV2_HVG$BioVar_correspondingZT, from=-6.5, to=8.1, na.rm=TRUE)
d_LVG <- density(ZT14_CorCV2_LVG$BioVar_correspondingZT, from=-6.5, to=8.1, na.rm=TRUE)

HVG_LVG_distribution_CorCV2_ZT14[,2] <- d_HVG$y
HVG_LVG_distribution_CorCV2_ZT14[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_CorCV2_ZT14) <- c("CorCV2", "HVG", "LVG")
HVG_LVG_distribution_CorCV2_ZT14_long <- gather(HVG_LVG_distribution_CorCV2_ZT14, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_CorCV2_ZT14 <- mutate(HVG_LVG_distribution_CorCV2_ZT14, mean_bootstrat=rowMeans(random_distribution_CorCV2_ZT14_1000[,2:1001]),
                                           bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_CorCV2_ZT14_1000[,2:1001])),
                                           top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_CorCV2_ZT14_1000[,2:1001])),
                                           ZT= "ZT14")


# ZT16 ----
random_distribution_CorCV2_ZT16_1000=matrix(nrow=512,ncol=1000)
random_distribution_CorCV2_ZT16_1000=data.frame(random_distribution_CorCV2_ZT16_1000)

random_CorCV2_ZT16_1000=matrix(nrow=630,ncol=1000)
random_CorCV2_ZT16_1000=data.frame(random_CorCV2_ZT16_1000)


for (i in 1:1000) {
  random_set <- CorCV2[sample(nrow(CorCV2), 630),]
  random_CorCV2_ZT16_1000[,i] <- random_set[,9]
  d <- density(random_set[,9], from=-6.5, to=8.1, na.rm=TRUE)
  random_distribution_CorCV2_ZT16_1000[,i] <- d$y
}

colnames(random_CorCV2_ZT16_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(CorCV2_allZT$CorCV2, from=-6.5, to=8.1, na.rm=TRUE)
random_distribution_CorCV2_ZT16_1000 <- cbind(d_background$x, random_distribution_CorCV2_ZT16_1000)
colnames(random_distribution_CorCV2_ZT16_1000) <- c("CorCV2",paste("random_set",1:100, sep="_"))


ZT16 <- read.table("Distribution_CV2_HVG_LVG_random/ZT16_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT16)

ZT16_CorCV2_HVG <- merge(ZT16,CorCV2,by.x="Gene",by.y="tracking_id") %>%
  select(Gene,TypeV,ZT,BioVar_ZT16) %>%
  filter(TypeV=="Variable")
names(ZT16_CorCV2_HVG) <-c("Gene","TypeV","ZT","BioVar_correspondingZT")

ZT16_CorCV2_LVG <- merge(ZT16,CorCV2,by.x="Gene",by.y="tracking_id") %>%
  select(Gene,TypeV,ZT,BioVar_ZT16) %>%
  filter(TypeV=="Stable")
names(ZT16_CorCV2_LVG) <-c("Gene","TypeV","ZT","BioVar_correspondingZT")


HVG_LVG_distribution_CorCV2_ZT16=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_CorCV2_ZT16=data.frame(HVG_LVG_distribution_CorCV2_ZT16)

d_background <-  density(CorCV2_allZT$CorCV2, from=-6.5, to=8.1, na.rm=TRUE)
HVG_LVG_distribution_CorCV2_ZT16[,1] <- d_background$x

d_HVG <- density(ZT16_CorCV2_HVG$BioVar_correspondingZT, from=-6.5, to=8.1, na.rm=TRUE)
d_LVG <- density(ZT16_CorCV2_LVG$BioVar_correspondingZT, from=-6.5, to=8.1, na.rm=TRUE)

HVG_LVG_distribution_CorCV2_ZT16[,2] <- d_HVG$y
HVG_LVG_distribution_CorCV2_ZT16[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_CorCV2_ZT16) <- c("CorCV2", "HVG", "LVG")
HVG_LVG_distribution_CorCV2_ZT16_long <- gather(HVG_LVG_distribution_CorCV2_ZT16, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_CorCV2_ZT16 <- mutate(HVG_LVG_distribution_CorCV2_ZT16, mean_bootstrat=rowMeans(random_distribution_CorCV2_ZT16_1000[,2:1001]),
                                           bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_CorCV2_ZT16_1000[,2:1001])),
                                           top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_CorCV2_ZT16_1000[,2:1001])),
                                           ZT= "ZT16")


# ZT18 ----
random_distribution_CorCV2_ZT18_1000=matrix(nrow=512,ncol=1000)
random_distribution_CorCV2_ZT18_1000=data.frame(random_distribution_CorCV2_ZT18_1000)

random_CorCV2_ZT18_1000=matrix(nrow=408,ncol=1000)
random_CorCV2_ZT18_1000=data.frame(random_CorCV2_ZT18_1000)


for (i in 1:1000) {
  random_set <- CorCV2[sample(nrow(CorCV2), 408),]
  random_CorCV2_ZT18_1000[,i] <- random_set[,10]
  d <- density(random_set[,10], from=-6.5, to=8.1, na.rm=TRUE)
  random_distribution_CorCV2_ZT18_1000[,i] <- d$y
}

colnames(random_CorCV2_ZT18_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(CorCV2_allZT$CorCV2, from=-6.5, to=8.1, na.rm=TRUE)
random_distribution_CorCV2_ZT18_1000 <- cbind(d_background$x, random_distribution_CorCV2_ZT18_1000)
colnames(random_distribution_CorCV2_ZT18_1000) <- c("CorCV2",paste("random_set",1:100, sep="_"))


ZT18 <- read.table("Distribution_CV2_HVG_LVG_random/ZT18_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT18)

ZT18_CorCV2_HVG <- merge(ZT18,CorCV2,by.x="Gene",by.y="tracking_id") %>%
  select(Gene,TypeV,ZT,BioVar_ZT18) %>%
  filter(TypeV=="Variable")
names(ZT18_CorCV2_HVG) <-c("Gene","TypeV","ZT","BioVar_correspondingZT")

ZT18_CorCV2_LVG <- merge(ZT18,CorCV2,by.x="Gene",by.y="tracking_id") %>%
  select(Gene,TypeV,ZT,BioVar_ZT18) %>%
  filter(TypeV=="Stable")
names(ZT18_CorCV2_LVG) <-c("Gene","TypeV","ZT","BioVar_correspondingZT")


HVG_LVG_distribution_CorCV2_ZT18=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_CorCV2_ZT18=data.frame(HVG_LVG_distribution_CorCV2_ZT18)

d_background <-  density(CorCV2_allZT$CorCV2, from=-6.5, to=8.1, na.rm=TRUE)
HVG_LVG_distribution_CorCV2_ZT18[,1] <- d_background$x

d_HVG <- density(ZT18_CorCV2_HVG$BioVar_correspondingZT, from=-6.5, to=8.1, na.rm=TRUE)
d_LVG <- density(ZT18_CorCV2_LVG$BioVar_correspondingZT, from=-6.5, to=8.1, na.rm=TRUE)

HVG_LVG_distribution_CorCV2_ZT18[,2] <- d_HVG$y
HVG_LVG_distribution_CorCV2_ZT18[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_CorCV2_ZT18) <- c("CorCV2", "HVG", "LVG")
HVG_LVG_distribution_CorCV2_ZT18_long <- gather(HVG_LVG_distribution_CorCV2_ZT18, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_CorCV2_ZT18 <- mutate(HVG_LVG_distribution_CorCV2_ZT18, mean_bootstrat=rowMeans(random_distribution_CorCV2_ZT18_1000[,2:1001]),
                                           bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_CorCV2_ZT18_1000[,2:1001])),
                                           top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_CorCV2_ZT18_1000[,2:1001])),
                                           ZT= "ZT18")



# ZT20 ----
random_distribution_CorCV2_ZT20_1000=matrix(nrow=512,ncol=1000)
random_distribution_CorCV2_ZT20_1000=data.frame(random_distribution_CorCV2_ZT20_1000)

random_CorCV2_ZT20_1000=matrix(nrow=716,ncol=1000)
random_CorCV2_ZT20_1000=data.frame(random_CorCV2_ZT20_1000)


for (i in 1:1000) {
  random_set <- CorCV2[sample(nrow(CorCV2), 716),]
  random_CorCV2_ZT20_1000[,i] <- random_set[,11]
  d <- density(random_set[,11], from=-6.5, to=8.1, na.rm=TRUE)
  random_distribution_CorCV2_ZT20_1000[,i] <- d$y
}

colnames(random_CorCV2_ZT20_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(CorCV2_allZT$CorCV2, from=-6.5, to=8.1, na.rm=TRUE)
random_distribution_CorCV2_ZT20_1000 <- cbind(d_background$x, random_distribution_CorCV2_ZT20_1000)
colnames(random_distribution_CorCV2_ZT20_1000) <- c("CorCV2",paste("random_set",1:100, sep="_"))


ZT20 <- read.table("Distribution_CV2_HVG_LVG_random/ZT20_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT20)

ZT20_CorCV2_HVG <- merge(ZT20,CorCV2,by.x="Gene",by.y="tracking_id") %>%
  select(Gene,TypeV,ZT,BioVar_ZT20) %>%
  filter(TypeV=="Variable")
names(ZT20_CorCV2_HVG) <-c("Gene","TypeV","ZT","BioVar_correspondingZT")

ZT20_CorCV2_LVG <- merge(ZT20,CorCV2,by.x="Gene",by.y="tracking_id") %>%
  select(Gene,TypeV,ZT,BioVar_ZT20) %>%
  filter(TypeV=="Stable")
names(ZT20_CorCV2_LVG) <-c("Gene","TypeV","ZT","BioVar_correspondingZT")


HVG_LVG_distribution_CorCV2_ZT20=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_CorCV2_ZT20=data.frame(HVG_LVG_distribution_CorCV2_ZT20)

d_background <-  density(CorCV2_allZT$CorCV2, from=-6.5, to=8.1, na.rm=TRUE)
HVG_LVG_distribution_CorCV2_ZT20[,1] <- d_background$x

d_HVG <- density(ZT20_CorCV2_HVG$BioVar_correspondingZT, from=-6.5, to=8.1, na.rm=TRUE)
d_LVG <- density(ZT20_CorCV2_LVG$BioVar_correspondingZT, from=-6.5, to=8.1, na.rm=TRUE)

HVG_LVG_distribution_CorCV2_ZT20[,2] <- d_HVG$y
HVG_LVG_distribution_CorCV2_ZT20[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_CorCV2_ZT20) <- c("CorCV2", "HVG", "LVG")
HVG_LVG_distribution_CorCV2_ZT20_long <- gather(HVG_LVG_distribution_CorCV2_ZT20, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_CorCV2_ZT20 <- mutate(HVG_LVG_distribution_CorCV2_ZT20, mean_bootstrat=rowMeans(random_distribution_CorCV2_ZT20_1000[,2:1001]),
                                           bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_CorCV2_ZT20_1000[,2:1001])),
                                           top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_CorCV2_ZT20_1000[,2:1001])),
                                           ZT= "ZT20")



# ZT22 ----
random_distribution_CorCV2_ZT22_1000=matrix(nrow=512,ncol=1000)
random_distribution_CorCV2_ZT22_1000=data.frame(random_distribution_CorCV2_ZT22_1000)

random_CorCV2_ZT22_1000=matrix(nrow=566,ncol=1000)
random_CorCV2_ZT22_1000=data.frame(random_CorCV2_ZT22_1000)


for (i in 1:1000) {
  random_set <- CorCV2[sample(nrow(CorCV2), 566),]
  random_CorCV2_ZT22_1000[,i] <- random_set[,12]
  d <- density(random_set[,12], from=-6.5, to=8.1, na.rm=TRUE)
  random_distribution_CorCV2_ZT22_1000[,i] <- d$y
}

colnames(random_CorCV2_ZT22_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(CorCV2_allZT$CorCV2, from=-6.5, to=8.1, na.rm=TRUE)
random_distribution_CorCV2_ZT22_1000 <- cbind(d_background$x, random_distribution_CorCV2_ZT22_1000)
colnames(random_distribution_CorCV2_ZT22_1000) <- c("CorCV2",paste("random_set",1:100, sep="_"))


ZT22 <- read.table("Distribution_CV2_HVG_LVG_random/ZT22_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT22)

ZT22_CorCV2_HVG <- merge(ZT22,CorCV2,by.x="Gene",by.y="tracking_id") %>%
  select(Gene,TypeV,ZT,BioVar_ZT22) %>%
  filter(TypeV=="Variable")
names(ZT22_CorCV2_HVG) <-c("Gene","TypeV","ZT","BioVar_correspondingZT")

ZT22_CorCV2_LVG <- merge(ZT22,CorCV2,by.x="Gene",by.y="tracking_id") %>%
  select(Gene,TypeV,ZT,BioVar_ZT22) %>%
  filter(TypeV=="Stable")
names(ZT22_CorCV2_LVG) <-c("Gene","TypeV","ZT","BioVar_correspondingZT")


HVG_LVG_distribution_CorCV2_ZT22=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_CorCV2_ZT22=data.frame(HVG_LVG_distribution_CorCV2_ZT22)

d_background <-  density(CorCV2_allZT$CorCV2, from=-6.5, to=8.1, na.rm=TRUE)
HVG_LVG_distribution_CorCV2_ZT22[,1] <- d_background$x

d_HVG <- density(ZT22_CorCV2_HVG$BioVar_correspondingZT, from=-6.5, to=8.1, na.rm=TRUE)
d_LVG <- density(ZT22_CorCV2_LVG$BioVar_correspondingZT, from=-6.5, to=8.1, na.rm=TRUE)

HVG_LVG_distribution_CorCV2_ZT22[,2] <- d_HVG$y
HVG_LVG_distribution_CorCV2_ZT22[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_CorCV2_ZT22) <- c("CorCV2", "HVG", "LVG")
HVG_LVG_distribution_CorCV2_ZT22_long <- gather(HVG_LVG_distribution_CorCV2_ZT22, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_CorCV2_ZT22 <- mutate(HVG_LVG_distribution_CorCV2_ZT22, mean_bootstrat=rowMeans(random_distribution_CorCV2_ZT22_1000[,2:1001]),
                                           bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_CorCV2_ZT22_1000[,2:1001])),
                                           top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_CorCV2_ZT22_1000[,2:1001])),
                                           ZT= "ZT22")


# ZT24 ----
random_distribution_CorCV2_ZT24_1000=matrix(nrow=512,ncol=1000)
random_distribution_CorCV2_ZT24_1000=data.frame(random_distribution_CorCV2_ZT24_1000)

random_CorCV2_ZT24_1000=matrix(nrow=377,ncol=1000)
random_CorCV2_ZT24_1000=data.frame(random_CorCV2_ZT24_1000)


for (i in 1:1000) {
  random_set <- CorCV2[sample(nrow(CorCV2), 377),]
  random_CorCV2_ZT24_1000[,i] <- random_set[,13]
  d <- density(random_set[,13], from=-6.5, to=8.1, na.rm=TRUE)
  random_distribution_CorCV2_ZT24_1000[,i] <- d$y
}

colnames(random_CorCV2_ZT24_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(CorCV2_allZT$CorCV2, from=-6.5, to=8.1, na.rm=TRUE)
random_distribution_CorCV2_ZT24_1000 <- cbind(d_background$x, random_distribution_CorCV2_ZT24_1000)
colnames(random_distribution_CorCV2_ZT24_1000) <- c("CorCV2",paste("random_set",1:100, sep="_"))


ZT24 <- read.table("Distribution_CV2_HVG_LVG_random/ZT0_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT24)

ZT24_CorCV2_HVG <- merge(ZT24,CorCV2,by.x="Gene",by.y="tracking_id") %>%
  select(Gene,TypeV,ZT,BioVar_ZT0) %>%
  filter(TypeV=="Variable")
names(ZT24_CorCV2_HVG) <-c("Gene","TypeV","ZT","BioVar_correspondingZT")

ZT24_CorCV2_LVG <- merge(ZT24,CorCV2,by.x="Gene",by.y="tracking_id") %>%
  select(Gene,TypeV,ZT,BioVar_ZT0) %>%
  filter(TypeV=="Stable")
names(ZT24_CorCV2_LVG) <-c("Gene","TypeV","ZT","BioVar_correspondingZT")


HVG_LVG_distribution_CorCV2_ZT24=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_CorCV2_ZT24=data.frame(HVG_LVG_distribution_CorCV2_ZT24)

d_background <-  density(CorCV2_allZT$CorCV2, from=-6.5, to=8.1, na.rm=TRUE)
HVG_LVG_distribution_CorCV2_ZT24[,1] <- d_background$x

d_HVG <- density(ZT24_CorCV2_HVG$BioVar_correspondingZT, from=-6.5, to=8.1, na.rm=TRUE)
d_LVG <- density(ZT24_CorCV2_LVG$BioVar_correspondingZT, from=-6.5, to=8.1, na.rm=TRUE)

HVG_LVG_distribution_CorCV2_ZT24[,2] <- d_HVG$y
HVG_LVG_distribution_CorCV2_ZT24[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_CorCV2_ZT24) <- c("CorCV2", "HVG", "LVG")
HVG_LVG_distribution_CorCV2_ZT24_long <- gather(HVG_LVG_distribution_CorCV2_ZT24, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_CorCV2_ZT24 <- mutate(HVG_LVG_distribution_CorCV2_ZT24, mean_bootstrat=rowMeans(random_distribution_CorCV2_ZT24_1000[,2:1001]),
                                           bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_CorCV2_ZT24_1000[,2:1001])),
                                           top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_CorCV2_ZT24_1000[,2:1001])),
                                           ZT= "ZT24")


# all ZT combined and graph


allZT_boostrap_CorCV2 <- bind_rows(HVG_LVG_distribution_CorCV2_ZT2, HVG_LVG_distribution_CorCV2_ZT4,
                                   HVG_LVG_distribution_CorCV2_ZT6, HVG_LVG_distribution_CorCV2_ZT8,
                                   HVG_LVG_distribution_CorCV2_ZT10, HVG_LVG_distribution_CorCV2_ZT12,
                                   HVG_LVG_distribution_CorCV2_ZT14, HVG_LVG_distribution_CorCV2_ZT16,
                                   HVG_LVG_distribution_CorCV2_ZT18, HVG_LVG_distribution_CorCV2_ZT20,
                                   HVG_LVG_distribution_CorCV2_ZT22, HVG_LVG_distribution_CorCV2_ZT24)



gather(allZT_boostrap_CorCV2, type, distr, HVG:top_CI, factor_key=TRUE) %>%
  ggplot(aes(x=CorCV2, y=distr, color=type, linetype=type)) +
  geom_line(size=1) +
  facet_wrap(~ZT) +
  scale_x_continuous(limits = c(-6, 8)) +
  facet_wrap(~ZT) +
  xlab("Corrected CV2 ZT24") +
  ylab("Density") +
  scale_color_manual(values=c("blue","green4","black","grey30","grey30")) +
  scale_linetype_manual(values=c("solid","solid","solid", "dashed", "dashed")) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line( size=.1, color="black" ),
        panel.background = element_rect(fill="white",colour="grey50"), 
        text=element_text(size=22))





