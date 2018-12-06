
library(tidyverse)

## Expression level at each ZT ----

# data preparation

Expression <- read.table("Expression_level/TPM2_14Seedlings_average_eachZT.txt",header=TRUE,sep='\t')

Expression_allZT <- gather(Expression, ZT, TPM, -Gene)

# bootstrap for each ZT
# ZT2 ----
random_distribution_Expression_ZT2_1000=matrix(nrow=512,ncol=1000)
random_distribution_Expression_ZT2_1000=data.frame(random_distribution_Expression_ZT2_1000)

random_Expression_ZT2_1000=matrix(nrow=313,ncol=1000)
random_Expression_ZT2_1000=data.frame(random_Expression_ZT2_1000)


for (i in 1:1000) {
  random_set <- Expression[sample(nrow(Expression), 313),]
  random_Expression_ZT2_1000[,i] <- random_set[,2]
  d <- density(log2(random_set[,2]), from=1.5, to=10, na.rm=TRUE)
  random_distribution_Expression_ZT2_1000[,i] <- d$y
}

colnames(random_Expression_ZT2_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(log2(Expression_allZT$TPM), from=1.5, to=10, na.rm=TRUE)
random_distribution_Expression_ZT2_1000 <- cbind(d_background$x, random_distribution_Expression_ZT2_1000)
colnames(random_distribution_Expression_ZT2_1000) <- c("TPM",paste("random_set",1:100, sep="_"))


ZT2 <- read.table("Expression_level/ZT2_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT2)

ZT2_Expression_HVG <- merge(ZT2,Expression,by.x="Gene",by.y="Gene") %>%
  select(Gene,TypeV,ZT,TPM_average_ZT2) %>%
  filter(TypeV=="Variable")
names(ZT2_Expression_HVG) <-c("Gene","TypeV","ZT","TPM_correspondingZT")

ZT2_Expression_LVG <- merge(ZT2,Expression,by.x="Gene",by.y="Gene") %>%
  select(Gene,TypeV,ZT,TPM_average_ZT2) %>%
  filter(TypeV=="Stable")
names(ZT2_Expression_LVG) <-c("Gene","TypeV","ZT","TPM_correspondingZT")


HVG_LVG_distribution_Expression_ZT2=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_Expression_ZT2=data.frame(HVG_LVG_distribution_Expression_ZT2)

d_background <-  density(log2(Expression_allZT$TPM), from=1.5, to=10, na.rm=TRUE)
HVG_LVG_distribution_Expression_ZT2[,1] <- d_background$x

d_HVG <- density(log2(ZT2_Expression_HVG$TPM_correspondingZT), from=1.5, to=10, na.rm=TRUE)
d_LVG <- density(log2(ZT2_Expression_LVG$TPM_correspondingZT), from=1.5, to=10, na.rm=TRUE)

HVG_LVG_distribution_Expression_ZT2[,2] <- d_HVG$y
HVG_LVG_distribution_Expression_ZT2[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_Expression_ZT2) <- c("TPM", "HVG", "LVG")
HVG_LVG_distribution_Expression_ZT2_long <- gather(HVG_LVG_distribution_Expression_ZT2, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_Expression_ZT2 <- mutate(HVG_LVG_distribution_Expression_ZT2, mean_bootstrat=rowMeans(random_distribution_Expression_ZT2_1000[,2:1001]),
                                              bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_Expression_ZT2_1000[,2:1001])),
                                              top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_Expression_ZT2_1000[,2:1001])),
                                              ZT= "ZT2")


# wilcoxon test

Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(ZT2_Expression_HVG$TPM_correspondingZT,random_Expression_ZT2_1000[,i])$p.value
}

names(Wilcox_pvalues) <- "result"

ggplot(Wilcox_pvalues, aes(x=result)) +
  geom_density()

# All values between 0.67793785 and 1.043956e-06 #



# ZT4 ----
random_distribution_Expression_ZT4_1000=matrix(nrow=512,ncol=1000)
random_distribution_Expression_ZT4_1000=data.frame(random_distribution_Expression_ZT4_1000)

random_Expression_ZT4_1000=matrix(nrow=257,ncol=1000)
random_Expression_ZT4_1000=data.frame(random_Expression_ZT4_1000)


for (i in 1:1000) {
  random_set <- Expression[sample(nrow(Expression), 257),]
  random_Expression_ZT4_1000[,i] <- random_set[,3]
  d <- density(log2(random_set[,3]), from=1.5, to=10, na.rm=TRUE)
  random_distribution_Expression_ZT4_1000[,i] <- d$y
}

colnames(random_Expression_ZT4_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(log2(Expression_allZT$TPM), from=1.5, to=10, na.rm=TRUE)
random_distribution_Expression_ZT4_1000 <- cbind(d_background$x, random_distribution_Expression_ZT4_1000)
colnames(random_distribution_Expression_ZT4_1000) <- c("TPM",paste("random_set",1:100, sep="_"))


ZT4 <- read.table("Expression_level/ZT4_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT4)

ZT4_Expression_HVG <- merge(ZT4,Expression,by.x="Gene",by.y="Gene") %>%
  select(Gene,TypeV,ZT,TPM_average_ZT4) %>%
  filter(TypeV=="Variable")
names(ZT4_Expression_HVG) <-c("Gene","TypeV","ZT","TPM_correspondingZT")

ZT4_Expression_LVG <- merge(ZT4,Expression,by.x="Gene",by.y="Gene") %>%
  select(Gene,TypeV,ZT,TPM_average_ZT4) %>%
  filter(TypeV=="Stable")
names(ZT4_Expression_LVG) <-c("Gene","TypeV","ZT","TPM_correspondingZT")


HVG_LVG_distribution_Expression_ZT4=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_Expression_ZT4=data.frame(HVG_LVG_distribution_Expression_ZT4)

d_background <-  density(log2(Expression_allZT$TPM), from=1.5, to=10, na.rm=TRUE)
HVG_LVG_distribution_Expression_ZT4[,1] <- d_background$x

d_HVG <- density(log2(ZT4_Expression_HVG$TPM_correspondingZT), from=1.5, to=10, na.rm=TRUE)
d_LVG <- density(log2(ZT4_Expression_LVG$TPM_correspondingZT), from=1.5, to=10, na.rm=TRUE)

HVG_LVG_distribution_Expression_ZT4[,2] <- d_HVG$y
HVG_LVG_distribution_Expression_ZT4[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_Expression_ZT4) <- c("TPM", "HVG", "LVG")
HVG_LVG_distribution_Expression_ZT4_long <- gather(HVG_LVG_distribution_Expression_ZT4, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_Expression_ZT4 <- mutate(HVG_LVG_distribution_Expression_ZT4, mean_bootstrat=rowMeans(random_distribution_Expression_ZT4_1000[,2:1001]),
                                              bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_Expression_ZT4_1000[,2:1001])),
                                              top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_Expression_ZT4_1000[,2:1001])),
                                              ZT= "ZT4")


# wilcoxon test

Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(ZT4_Expression_HVG$TPM_correspondingZT,random_Expression_ZT4_1000[,i])$p.value
}

names(Wilcox_pvalues) <- "result"

ggplot(Wilcox_pvalues, aes(x=result)) +
  geom_density()

# All values between 0.32291453 and 4.129703e-07 #



# ZT6 ----
random_distribution_Expression_ZT6_1000=matrix(nrow=512,ncol=1000)
random_distribution_Expression_ZT6_1000=data.frame(random_distribution_Expression_ZT6_1000)

random_Expression_ZT6_1000=matrix(nrow=374,ncol=1000)
random_Expression_ZT6_1000=data.frame(random_Expression_ZT6_1000)


for (i in 1:1000) {
  random_set <- Expression[sample(nrow(Expression), 374),]
  random_Expression_ZT6_1000[,i] <- random_set[,4]
  d <- density(log2(random_set[,4]), from=1.5, to=10, na.rm=TRUE)
  random_distribution_Expression_ZT6_1000[,i] <- d$y
}

colnames(random_Expression_ZT6_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(log2(Expression_allZT$TPM), from=1.5, to=10, na.rm=TRUE)
random_distribution_Expression_ZT6_1000 <- cbind(d_background$x, random_distribution_Expression_ZT6_1000)
colnames(random_distribution_Expression_ZT6_1000) <- c("TPM",paste("random_set",1:100, sep="_"))


ZT6 <- read.table("Expression_level/ZT6_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT6)

ZT6_Expression_HVG <- merge(ZT6,Expression,by.x="Gene",by.y="Gene") %>%
  select(Gene,TypeV,ZT,TPM_average_ZT6) %>%
  filter(TypeV=="Variable")
names(ZT6_Expression_HVG) <-c("Gene","TypeV","ZT","TPM_correspondingZT")

ZT6_Expression_LVG <- merge(ZT6,Expression,by.x="Gene",by.y="Gene") %>%
  select(Gene,TypeV,ZT,TPM_average_ZT6) %>%
  filter(TypeV=="Stable")
names(ZT6_Expression_LVG) <-c("Gene","TypeV","ZT","TPM_correspondingZT")


HVG_LVG_distribution_Expression_ZT6=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_Expression_ZT6=data.frame(HVG_LVG_distribution_Expression_ZT6)

d_background <-  density(log2(Expression_allZT$TPM), from=1.5, to=10, na.rm=TRUE)
HVG_LVG_distribution_Expression_ZT6[,1] <- d_background$x

d_HVG <- density(log2(ZT6_Expression_HVG$TPM_correspondingZT), from=1.5, to=10, na.rm=TRUE)
d_LVG <- density(log2(ZT6_Expression_LVG$TPM_correspondingZT), from=1.5, to=10, na.rm=TRUE)

HVG_LVG_distribution_Expression_ZT6[,2] <- d_HVG$y
HVG_LVG_distribution_Expression_ZT6[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_Expression_ZT6) <- c("TPM", "HVG", "LVG")
HVG_LVG_distribution_Expression_ZT6_long <- gather(HVG_LVG_distribution_Expression_ZT6, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_Expression_ZT6 <- mutate(HVG_LVG_distribution_Expression_ZT6, mean_bootstrat=rowMeans(random_distribution_Expression_ZT6_1000[,2:1001]),
                                              bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_Expression_ZT6_1000[,2:1001])),
                                              top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_Expression_ZT6_1000[,2:1001])),
                                              ZT= "ZT6")



# wilcoxon test

Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(ZT6_Expression_HVG$TPM_correspondingZT,random_Expression_ZT6_1000[,i])$p.value
}

names(Wilcox_pvalues) <- "result"

ggplot(Wilcox_pvalues, aes(x=result)) +
  geom_density()

# All values between 0.0108541657 and 3.187834e-12 #



# ZT8 ----
random_distribution_Expression_ZT8_1000=matrix(nrow=512,ncol=1000)
random_distribution_Expression_ZT8_1000=data.frame(random_distribution_Expression_ZT8_1000)

random_Expression_ZT8_1000=matrix(nrow=327,ncol=1000)
random_Expression_ZT8_1000=data.frame(random_Expression_ZT8_1000)


for (i in 1:1000) {
  random_set <- Expression[sample(nrow(Expression), 327),]
  random_Expression_ZT8_1000[,i] <- random_set[,5]
  d <- density(log2(random_set[,5]), from=1.5, to=10, na.rm=TRUE)
  random_distribution_Expression_ZT8_1000[,i] <- d$y
}

colnames(random_Expression_ZT8_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(log2(Expression_allZT$TPM), from=1.5, to=10, na.rm=TRUE)
random_distribution_Expression_ZT8_1000 <- cbind(d_background$x, random_distribution_Expression_ZT8_1000)
colnames(random_distribution_Expression_ZT8_1000) <- c("TPM",paste("random_set",1:100, sep="_"))


ZT8 <- read.table("Expression_level/ZT8_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT8)

ZT8_Expression_HVG <- merge(ZT8,Expression,by.x="Gene",by.y="Gene") %>%
  select(Gene,TypeV,ZT,TPM_average_ZT8) %>%
  filter(TypeV=="Variable")
names(ZT8_Expression_HVG) <-c("Gene","TypeV","ZT","TPM_correspondingZT")

ZT8_Expression_LVG <- merge(ZT8,Expression,by.x="Gene",by.y="Gene") %>%
  select(Gene,TypeV,ZT,TPM_average_ZT8) %>%
  filter(TypeV=="Stable")
names(ZT8_Expression_LVG) <-c("Gene","TypeV","ZT","TPM_correspondingZT")


HVG_LVG_distribution_Expression_ZT8=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_Expression_ZT8=data.frame(HVG_LVG_distribution_Expression_ZT8)

d_background <-  density(log2(Expression_allZT$TPM), from=1.5, to=10, na.rm=TRUE)
HVG_LVG_distribution_Expression_ZT8[,1] <- d_background$x

d_HVG <- density(log2(ZT8_Expression_HVG$TPM_correspondingZT), from=1.5, to=10, na.rm=TRUE)
d_LVG <- density(log2(ZT8_Expression_LVG$TPM_correspondingZT), from=1.5, to=10, na.rm=TRUE)

HVG_LVG_distribution_Expression_ZT8[,2] <- d_HVG$y
HVG_LVG_distribution_Expression_ZT8[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_Expression_ZT8) <- c("TPM", "HVG", "LVG")
HVG_LVG_distribution_Expression_ZT8_long <- gather(HVG_LVG_distribution_Expression_ZT8, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_Expression_ZT8 <- mutate(HVG_LVG_distribution_Expression_ZT8, mean_bootstrat=rowMeans(random_distribution_Expression_ZT8_1000[,2:1001]),
                                              bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_Expression_ZT8_1000[,2:1001])),
                                              top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_Expression_ZT8_1000[,2:1001])),
                                              ZT= "ZT8")


# wilcoxon test

Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(ZT8_Expression_HVG$TPM_correspondingZT,random_Expression_ZT8_1000[,i])$p.value
}

names(Wilcox_pvalues) <- "result"

ggplot(Wilcox_pvalues, aes(x=result)) +
  geom_density()

# All values between 0.053156629 and 2.124110e-09 #


# ZT10 ----
random_distribution_Expression_ZT10_1000=matrix(nrow=512,ncol=1000)
random_distribution_Expression_ZT10_1000=data.frame(random_distribution_Expression_ZT10_1000)

random_Expression_ZT10_1000=matrix(nrow=310,ncol=1000)
random_Expression_ZT10_1000=data.frame(random_Expression_ZT10_1000)


for (i in 1:1000) {
  random_set <- Expression[sample(nrow(Expression), 310),]
  random_Expression_ZT10_1000[,i] <- random_set[,6]
  d <- density(log2(random_set[,6]), from=1.5, to=10, na.rm=TRUE)
  random_distribution_Expression_ZT10_1000[,i] <- d$y
}

colnames(random_Expression_ZT10_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(log2(Expression_allZT$TPM), from=1.5, to=10, na.rm=TRUE)
random_distribution_Expression_ZT10_1000 <- cbind(d_background$x, random_distribution_Expression_ZT10_1000)
colnames(random_distribution_Expression_ZT10_1000) <- c("TPM",paste("random_set",1:100, sep="_"))


ZT10 <- read.table("Expression_level/ZT10_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT10)

ZT10_Expression_HVG <- merge(ZT10,Expression,by.x="Gene",by.y="Gene") %>%
  select(Gene,TypeV,ZT,TPM_average_ZT10) %>%
  filter(TypeV=="Variable")
names(ZT10_Expression_HVG) <-c("Gene","TypeV","ZT","TPM_correspondingZT")

ZT10_Expression_LVG <- merge(ZT10,Expression,by.x="Gene",by.y="Gene") %>%
  select(Gene,TypeV,ZT,TPM_average_ZT10) %>%
  filter(TypeV=="Stable")
names(ZT10_Expression_LVG) <-c("Gene","TypeV","ZT","TPM_correspondingZT")


HVG_LVG_distribution_Expression_ZT10=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_Expression_ZT10=data.frame(HVG_LVG_distribution_Expression_ZT10)

d_background <-  density(log2(Expression_allZT$TPM), from=1.5, to=10, na.rm=TRUE)
HVG_LVG_distribution_Expression_ZT10[,1] <- d_background$x

d_HVG <- density(log2(ZT10_Expression_HVG$TPM_correspondingZT), from=1.5, to=10, na.rm=TRUE)
d_LVG <- density(log2(ZT10_Expression_LVG$TPM_correspondingZT), from=1.5, to=10, na.rm=TRUE)

HVG_LVG_distribution_Expression_ZT10[,2] <- d_HVG$y
HVG_LVG_distribution_Expression_ZT10[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_Expression_ZT10) <- c("TPM", "HVG", "LVG")
HVG_LVG_distribution_Expression_ZT10_long <- gather(HVG_LVG_distribution_Expression_ZT10, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_Expression_ZT10 <- mutate(HVG_LVG_distribution_Expression_ZT10, mean_bootstrat=rowMeans(random_distribution_Expression_ZT10_1000[,2:1001]),
                                               bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_Expression_ZT10_1000[,2:1001])),
                                               top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_Expression_ZT10_1000[,2:1001])),
                                               ZT= "ZT10")


# wilcoxon test

Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(ZT10_Expression_HVG$TPM_correspondingZT,random_Expression_ZT10_1000[,i])$p.value
}

names(Wilcox_pvalues) <- "result"

ggplot(Wilcox_pvalues, aes(x=result)) +
  geom_density()

# All values between 2.852049e-04 and 1.027666e-14 #



# ZT12 ----
random_distribution_Expression_ZT12_1000=matrix(nrow=512,ncol=1000)
random_distribution_Expression_ZT12_1000=data.frame(random_distribution_Expression_ZT12_1000)

random_Expression_ZT12_1000=matrix(nrow=636,ncol=1000)
random_Expression_ZT12_1000=data.frame(random_Expression_ZT12_1000)


for (i in 1:1000) {
  random_set <- Expression[sample(nrow(Expression), 636),]
  random_Expression_ZT12_1000[,i] <- random_set[,7]
  d <- density(log2(random_set[,7]), from=1.5, to=10, na.rm=TRUE)
  random_distribution_Expression_ZT12_1000[,i] <- d$y
}

colnames(random_Expression_ZT12_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(log2(Expression_allZT$TPM), from=1.5, to=10, na.rm=TRUE)
random_distribution_Expression_ZT12_1000 <- cbind(d_background$x, random_distribution_Expression_ZT12_1000)
colnames(random_distribution_Expression_ZT12_1000) <- c("TPM",paste("random_set",1:100, sep="_"))


ZT12 <- read.table("Expression_level/ZT12_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT12)

ZT12_Expression_HVG <- merge(ZT12,Expression,by.x="Gene",by.y="Gene") %>%
  select(Gene,TypeV,ZT,TPM_average_ZT12) %>%
  filter(TypeV=="Variable")
names(ZT12_Expression_HVG) <-c("Gene","TypeV","ZT","TPM_correspondingZT")

ZT12_Expression_LVG <- merge(ZT12,Expression,by.x="Gene",by.y="Gene") %>%
  select(Gene,TypeV,ZT,TPM_average_ZT12) %>%
  filter(TypeV=="Stable")
names(ZT12_Expression_LVG) <-c("Gene","TypeV","ZT","TPM_correspondingZT")


HVG_LVG_distribution_Expression_ZT12=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_Expression_ZT12=data.frame(HVG_LVG_distribution_Expression_ZT12)

d_background <-  density(log2(Expression_allZT$TPM), from=1.5, to=10, na.rm=TRUE)
HVG_LVG_distribution_Expression_ZT12[,1] <- d_background$x

d_HVG <- density(log2(ZT12_Expression_HVG$TPM_correspondingZT), from=1.5, to=10, na.rm=TRUE)
d_LVG <- density(log2(ZT12_Expression_LVG$TPM_correspondingZT), from=1.5, to=10, na.rm=TRUE)

HVG_LVG_distribution_Expression_ZT12[,2] <- d_HVG$y
HVG_LVG_distribution_Expression_ZT12[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_Expression_ZT12) <- c("TPM", "HVG", "LVG")
HVG_LVG_distribution_Expression_ZT12_long <- gather(HVG_LVG_distribution_Expression_ZT12, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_Expression_ZT12 <- mutate(HVG_LVG_distribution_Expression_ZT12, mean_bootstrat=rowMeans(random_distribution_Expression_ZT12_1000[,2:1001]),
                                               bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_Expression_ZT12_1000[,2:1001])),
                                               top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_Expression_ZT12_1000[,2:1001])),
                                               ZT= "ZT12")

# wilcoxon test

Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(ZT12_Expression_HVG$TPM_correspondingZT,random_Expression_ZT12_1000[,i])$p.value
}

names(Wilcox_pvalues) <- "result"

ggplot(Wilcox_pvalues, aes(x=result)) +
  geom_density()

# All values between 1.686746e-03 and 2.774496e-14 #


# ZT14 ----
random_distribution_Expression_ZT14_1000=matrix(nrow=512,ncol=1000)
random_distribution_Expression_ZT14_1000=data.frame(random_distribution_Expression_ZT14_1000)

random_Expression_ZT14_1000=matrix(nrow=418,ncol=1000)
random_Expression_ZT14_1000=data.frame(random_Expression_ZT14_1000)


for (i in 1:1000) {
  random_set <- Expression[sample(nrow(Expression), 418),]
  random_Expression_ZT14_1000[,i] <- random_set[,8]
  d <- density(log2(random_set[,8]), from=1.5, to=10, na.rm=TRUE)
  random_distribution_Expression_ZT14_1000[,i] <- d$y
}

colnames(random_Expression_ZT14_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(log2(Expression_allZT$TPM), from=1.5, to=10, na.rm=TRUE)
random_distribution_Expression_ZT14_1000 <- cbind(d_background$x, random_distribution_Expression_ZT14_1000)
colnames(random_distribution_Expression_ZT14_1000) <- c("TPM",paste("random_set",1:100, sep="_"))


ZT14 <- read.table("Expression_level/ZT14_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT14)

ZT14_Expression_HVG <- merge(ZT14,Expression,by.x="Gene",by.y="Gene") %>%
  select(Gene,TypeV,ZT,TPM_average_ZT14) %>%
  filter(TypeV=="Variable")
names(ZT14_Expression_HVG) <-c("Gene","TypeV","ZT","TPM_correspondingZT")

ZT14_Expression_LVG <- merge(ZT14,Expression,by.x="Gene",by.y="Gene") %>%
  select(Gene,TypeV,ZT,TPM_average_ZT14) %>%
  filter(TypeV=="Stable")
names(ZT14_Expression_LVG) <-c("Gene","TypeV","ZT","TPM_correspondingZT")


HVG_LVG_distribution_Expression_ZT14=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_Expression_ZT14=data.frame(HVG_LVG_distribution_Expression_ZT14)

d_background <-  density(log2(Expression_allZT$TPM), from=1.5, to=10, na.rm=TRUE)
HVG_LVG_distribution_Expression_ZT14[,1] <- d_background$x

d_HVG <- density(log2(ZT14_Expression_HVG$TPM_correspondingZT), from=1.5, to=10, na.rm=TRUE)
d_LVG <- density(log2(ZT14_Expression_LVG$TPM_correspondingZT), from=1.5, to=10, na.rm=TRUE)

HVG_LVG_distribution_Expression_ZT14[,2] <- d_HVG$y
HVG_LVG_distribution_Expression_ZT14[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_Expression_ZT14) <- c("TPM", "HVG", "LVG")
HVG_LVG_distribution_Expression_ZT14_long <- gather(HVG_LVG_distribution_Expression_ZT14, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_Expression_ZT14 <- mutate(HVG_LVG_distribution_Expression_ZT14, mean_bootstrat=rowMeans(random_distribution_Expression_ZT14_1000[,2:1001]),
                                               bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_Expression_ZT14_1000[,2:1001])),
                                               top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_Expression_ZT14_1000[,2:1001])),
                                               ZT= "ZT14")


# wilcoxon test

Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(ZT14_Expression_HVG$TPM_correspondingZT,random_Expression_ZT14_1000[,i])$p.value
}

names(Wilcox_pvalues) <- "result"

ggplot(Wilcox_pvalues, aes(x=result)) +
  geom_density()

# All values between 2.923199e-08 and 3.008820e-23 #




# ZT16 ----
random_distribution_Expression_ZT16_1000=matrix(nrow=512,ncol=1000)
random_distribution_Expression_ZT16_1000=data.frame(random_distribution_Expression_ZT16_1000)

random_Expression_ZT16_1000=matrix(nrow=630,ncol=1000)
random_Expression_ZT16_1000=data.frame(random_Expression_ZT16_1000)


for (i in 1:1000) {
  random_set <- Expression[sample(nrow(Expression), 630),]
  random_Expression_ZT16_1000[,i] <- random_set[,9]
  d <- density(log2(random_set[,9]), from=1.5, to=10, na.rm=TRUE)
  random_distribution_Expression_ZT16_1000[,i] <- d$y
}

colnames(random_Expression_ZT16_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(log2(Expression_allZT$TPM), from=1.5, to=10, na.rm=TRUE)
random_distribution_Expression_ZT16_1000 <- cbind(d_background$x, random_distribution_Expression_ZT16_1000)
colnames(random_distribution_Expression_ZT16_1000) <- c("TPM",paste("random_set",1:100, sep="_"))


ZT16 <- read.table("Expression_level/ZT16_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT16)

ZT16_Expression_HVG <- merge(ZT16,Expression,by.x="Gene",by.y="Gene") %>%
  select(Gene,TypeV,ZT,TPM_average_ZT16) %>%
  filter(TypeV=="Variable")
names(ZT16_Expression_HVG) <-c("Gene","TypeV","ZT","TPM_correspondingZT")

ZT16_Expression_LVG <- merge(ZT16,Expression,by.x="Gene",by.y="Gene") %>%
  select(Gene,TypeV,ZT,TPM_average_ZT16) %>%
  filter(TypeV=="Stable")
names(ZT16_Expression_LVG) <-c("Gene","TypeV","ZT","TPM_correspondingZT")


HVG_LVG_distribution_Expression_ZT16=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_Expression_ZT16=data.frame(HVG_LVG_distribution_Expression_ZT16)

d_background <-  density(log2(Expression_allZT$TPM), from=1.5, to=10, na.rm=TRUE)
HVG_LVG_distribution_Expression_ZT16[,1] <- d_background$x

d_HVG <- density(log2(ZT16_Expression_HVG$TPM_correspondingZT), from=1.5, to=10, na.rm=TRUE)
d_LVG <- density(log2(ZT16_Expression_LVG$TPM_correspondingZT), from=1.5, to=10, na.rm=TRUE)

HVG_LVG_distribution_Expression_ZT16[,2] <- d_HVG$y
HVG_LVG_distribution_Expression_ZT16[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_Expression_ZT16) <- c("TPM", "HVG", "LVG")
HVG_LVG_distribution_Expression_ZT16_long <- gather(HVG_LVG_distribution_Expression_ZT16, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_Expression_ZT16 <- mutate(HVG_LVG_distribution_Expression_ZT16, mean_bootstrat=rowMeans(random_distribution_Expression_ZT16_1000[,2:1001]),
                                               bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_Expression_ZT16_1000[,2:1001])),
                                               top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_Expression_ZT16_1000[,2:1001])),
                                               ZT= "ZT16")


# wilcoxon test

Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(ZT16_Expression_HVG$TPM_correspondingZT,random_Expression_ZT16_1000[,i])$p.value
}

names(Wilcox_pvalues) <- "result"

ggplot(Wilcox_pvalues, aes(x=result)) +
  geom_density()

# All values between 7.904203e-13 and 2.151787e-28 #



# ZT18 ----
random_distribution_Expression_ZT18_1000=matrix(nrow=512,ncol=1000)
random_distribution_Expression_ZT18_1000=data.frame(random_distribution_Expression_ZT18_1000)

random_Expression_ZT18_1000=matrix(nrow=408,ncol=1000)
random_Expression_ZT18_1000=data.frame(random_Expression_ZT18_1000)


for (i in 1:1000) {
  random_set <- Expression[sample(nrow(Expression), 408),]
  random_Expression_ZT18_1000[,i] <- random_set[,10]
  d <- density(log2(random_set[,10]), from=1.5, to=10, na.rm=TRUE)
  random_distribution_Expression_ZT18_1000[,i] <- d$y
}

colnames(random_Expression_ZT18_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(log2(Expression_allZT$TPM), from=1.5, to=10, na.rm=TRUE)
random_distribution_Expression_ZT18_1000 <- cbind(d_background$x, random_distribution_Expression_ZT18_1000)
colnames(random_distribution_Expression_ZT18_1000) <- c("TPM",paste("random_set",1:100, sep="_"))


ZT18 <- read.table("Expression_level/ZT18_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT18)

ZT18_Expression_HVG <- merge(ZT18,Expression,by.x="Gene",by.y="Gene") %>%
  select(Gene,TypeV,ZT,TPM_average_ZT18) %>%
  filter(TypeV=="Variable")
names(ZT18_Expression_HVG) <-c("Gene","TypeV","ZT","TPM_correspondingZT")

ZT18_Expression_LVG <- merge(ZT18,Expression,by.x="Gene",by.y="Gene") %>%
  select(Gene,TypeV,ZT,TPM_average_ZT18) %>%
  filter(TypeV=="Stable")
names(ZT18_Expression_LVG) <-c("Gene","TypeV","ZT","TPM_correspondingZT")


HVG_LVG_distribution_Expression_ZT18=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_Expression_ZT18=data.frame(HVG_LVG_distribution_Expression_ZT18)

d_background <-  density(log2(Expression_allZT$TPM), from=1.5, to=10, na.rm=TRUE)
HVG_LVG_distribution_Expression_ZT18[,1] <- d_background$x

d_HVG <- density(log2(ZT18_Expression_HVG$TPM_correspondingZT), from=1.5, to=10, na.rm=TRUE)
d_LVG <- density(log2(ZT18_Expression_LVG$TPM_correspondingZT), from=1.5, to=10, na.rm=TRUE)

HVG_LVG_distribution_Expression_ZT18[,2] <- d_HVG$y
HVG_LVG_distribution_Expression_ZT18[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_Expression_ZT18) <- c("TPM", "HVG", "LVG")
HVG_LVG_distribution_Expression_ZT18_long <- gather(HVG_LVG_distribution_Expression_ZT18, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_Expression_ZT18 <- mutate(HVG_LVG_distribution_Expression_ZT18, mean_bootstrat=rowMeans(random_distribution_Expression_ZT18_1000[,2:1001]),
                                               bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_Expression_ZT18_1000[,2:1001])),
                                               top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_Expression_ZT18_1000[,2:1001])),
                                               ZT= "ZT18")


# wilcoxon test

Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(ZT18_Expression_HVG$TPM_correspondingZT,random_Expression_ZT18_1000[,i])$p.value
}

names(Wilcox_pvalues) <- "result"

ggplot(Wilcox_pvalues, aes(x=result)) +
  geom_density()

# All values between 1.077009e-08 and 1.073335e-24 #


# ZT20 ----
random_distribution_Expression_ZT20_1000=matrix(nrow=512,ncol=1000)
random_distribution_Expression_ZT20_1000=data.frame(random_distribution_Expression_ZT20_1000)

random_Expression_ZT20_1000=matrix(nrow=716,ncol=1000)
random_Expression_ZT20_1000=data.frame(random_Expression_ZT20_1000)


for (i in 1:1000) {
  random_set <- Expression[sample(nrow(Expression), 716),]
  random_Expression_ZT20_1000[,i] <- random_set[,11]
  d <- density(log2(random_set[,11]), from=1.5, to=10, na.rm=TRUE)
  random_distribution_Expression_ZT20_1000[,i] <- d$y
}

colnames(random_Expression_ZT20_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(log2(Expression_allZT$TPM), from=1.5, to=10, na.rm=TRUE)
random_distribution_Expression_ZT20_1000 <- cbind(d_background$x, random_distribution_Expression_ZT20_1000)
colnames(random_distribution_Expression_ZT20_1000) <- c("TPM",paste("random_set",1:100, sep="_"))


ZT20 <- read.table("Expression_level/ZT20_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT20)

ZT20_Expression_HVG <- merge(ZT20,Expression,by.x="Gene",by.y="Gene") %>%
  select(Gene,TypeV,ZT,TPM_average_ZT20) %>%
  filter(TypeV=="Variable")
names(ZT20_Expression_HVG) <-c("Gene","TypeV","ZT","TPM_correspondingZT")

ZT20_Expression_LVG <- merge(ZT20,Expression,by.x="Gene",by.y="Gene") %>%
  select(Gene,TypeV,ZT,TPM_average_ZT20) %>%
  filter(TypeV=="Stable")
names(ZT20_Expression_LVG) <-c("Gene","TypeV","ZT","TPM_correspondingZT")


HVG_LVG_distribution_Expression_ZT20=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_Expression_ZT20=data.frame(HVG_LVG_distribution_Expression_ZT20)

d_background <-  density(log2(Expression_allZT$TPM), from=1.5, to=10, na.rm=TRUE)
HVG_LVG_distribution_Expression_ZT20[,1] <- d_background$x

d_HVG <- density(log2(ZT20_Expression_HVG$TPM_correspondingZT), from=1.5, to=10, na.rm=TRUE)
d_LVG <- density(log2(ZT20_Expression_LVG$TPM_correspondingZT), from=1.5, to=10, na.rm=TRUE)

HVG_LVG_distribution_Expression_ZT20[,2] <- d_HVG$y
HVG_LVG_distribution_Expression_ZT20[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_Expression_ZT20) <- c("TPM", "HVG", "LVG")
HVG_LVG_distribution_Expression_ZT20_long <- gather(HVG_LVG_distribution_Expression_ZT20, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_Expression_ZT20 <- mutate(HVG_LVG_distribution_Expression_ZT20, mean_bootstrat=rowMeans(random_distribution_Expression_ZT20_1000[,2:1001]),
                                               bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_Expression_ZT20_1000[,2:1001])),
                                               top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_Expression_ZT20_1000[,2:1001])),
                                               ZT= "ZT20")



# wilcoxon test

Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(ZT20_Expression_HVG$TPM_correspondingZT,random_Expression_ZT20_1000[,i])$p.value
}

names(Wilcox_pvalues) <- "result"

ggplot(Wilcox_pvalues, aes(x=result)) +
  geom_density()

# All values between 3.247710e-07 and 2.011718e-21 #



# ZT22 ----
random_distribution_Expression_ZT22_1000=matrix(nrow=512,ncol=1000)
random_distribution_Expression_ZT22_1000=data.frame(random_distribution_Expression_ZT22_1000)

random_Expression_ZT22_1000=matrix(nrow=566,ncol=1000)
random_Expression_ZT22_1000=data.frame(random_Expression_ZT22_1000)


for (i in 1:1000) {
  random_set <- Expression[sample(nrow(Expression), 566),]
  random_Expression_ZT22_1000[,i] <- random_set[,12]
  d <- density(log2(random_set[,12]), from=1.5, to=10, na.rm=TRUE)
  random_distribution_Expression_ZT22_1000[,i] <- d$y
}

colnames(random_Expression_ZT22_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(log2(Expression_allZT$TPM), from=1.5, to=10, na.rm=TRUE)
random_distribution_Expression_ZT22_1000 <- cbind(d_background$x, random_distribution_Expression_ZT22_1000)
colnames(random_distribution_Expression_ZT22_1000) <- c("TPM",paste("random_set",1:100, sep="_"))


ZT22 <- read.table("Expression_level/ZT22_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT22)

ZT22_Expression_HVG <- merge(ZT22,Expression,by.x="Gene",by.y="Gene") %>%
  select(Gene,TypeV,ZT,TPM_average_ZT22) %>%
  filter(TypeV=="Variable")
names(ZT22_Expression_HVG) <-c("Gene","TypeV","ZT","TPM_correspondingZT")

ZT22_Expression_LVG <- merge(ZT22,Expression,by.x="Gene",by.y="Gene") %>%
  select(Gene,TypeV,ZT,TPM_average_ZT22) %>%
  filter(TypeV=="Stable")
names(ZT22_Expression_LVG) <-c("Gene","TypeV","ZT","TPM_correspondingZT")


HVG_LVG_distribution_Expression_ZT22=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_Expression_ZT22=data.frame(HVG_LVG_distribution_Expression_ZT22)

d_background <-  density(log2(Expression_allZT$TPM), from=1.5, to=10, na.rm=TRUE)
HVG_LVG_distribution_Expression_ZT22[,1] <- d_background$x

d_HVG <- density(log2(ZT22_Expression_HVG$TPM_correspondingZT), from=1.5, to=10, na.rm=TRUE)
d_LVG <- density(log2(ZT22_Expression_LVG$TPM_correspondingZT), from=1.5, to=10, na.rm=TRUE)

HVG_LVG_distribution_Expression_ZT22[,2] <- d_HVG$y
HVG_LVG_distribution_Expression_ZT22[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_Expression_ZT22) <- c("TPM", "HVG", "LVG")
HVG_LVG_distribution_Expression_ZT22_long <- gather(HVG_LVG_distribution_Expression_ZT22, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_Expression_ZT22 <- mutate(HVG_LVG_distribution_Expression_ZT22, mean_bootstrat=rowMeans(random_distribution_Expression_ZT22_1000[,2:1001]),
                                               bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_Expression_ZT22_1000[,2:1001])),
                                               top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_Expression_ZT22_1000[,2:1001])),
                                               ZT= "ZT22")


# wilcoxon test

Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(ZT22_Expression_HVG$TPM_correspondingZT,random_Expression_ZT22_1000[,i])$p.value
}

names(Wilcox_pvalues) <- "result"

ggplot(Wilcox_pvalues, aes(x=result)) +
  geom_density()

# All values between 3.431072e-04 and 2.462566e-15 #




# ZT24 ----
random_distribution_Expression_ZT24_1000=matrix(nrow=512,ncol=1000)
random_distribution_Expression_ZT24_1000=data.frame(random_distribution_Expression_ZT24_1000)

random_Expression_ZT24_1000=matrix(nrow=377,ncol=1000)
random_Expression_ZT24_1000=data.frame(random_Expression_ZT24_1000)


for (i in 1:1000) {
  random_set <- Expression[sample(nrow(Expression), 377),]
  random_Expression_ZT24_1000[,i] <- random_set[,13]
  d <- density(log2(random_set[,13]), from=1.5, to=10, na.rm=TRUE)
  random_distribution_Expression_ZT24_1000[,i] <- d$y
}

colnames(random_Expression_ZT24_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(log2(Expression_allZT$TPM), from=1.5, to=10, na.rm=TRUE)
random_distribution_Expression_ZT24_1000 <- cbind(d_background$x, random_distribution_Expression_ZT24_1000)
colnames(random_distribution_Expression_ZT24_1000) <- c("TPM",paste("random_set",1:100, sep="_"))


ZT24 <- read.table("Expression_level/ZT0_HVG_LVG_random.txt",header=TRUE,sep='\t')
dim(ZT24)

ZT24_Expression_HVG <- merge(ZT24,Expression,by.x="Gene",by.y="Gene") %>%
  select(Gene,TypeV,ZT,TPM_average_ZT24) %>%
  filter(TypeV=="Variable")
names(ZT24_Expression_HVG) <-c("Gene","TypeV","ZT","TPM_correspondingZT")

ZT24_Expression_LVG <- merge(ZT24,Expression,by.x="Gene",by.y="Gene") %>%
  select(Gene,TypeV,ZT,TPM_average_ZT24) %>%
  filter(TypeV=="Stable")
names(ZT24_Expression_LVG) <-c("Gene","TypeV","ZT","TPM_correspondingZT")


HVG_LVG_distribution_Expression_ZT24=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_Expression_ZT24=data.frame(HVG_LVG_distribution_Expression_ZT24)

d_background <-  density(log2(Expression_allZT$TPM), from=1.5, to=10, na.rm=TRUE)
HVG_LVG_distribution_Expression_ZT24[,1] <- d_background$x

d_HVG <- density(log2(ZT24_Expression_HVG$TPM_correspondingZT), from=1.5, to=10, na.rm=TRUE)
d_LVG <- density(log2(ZT24_Expression_LVG$TPM_correspondingZT), from=1.5, to=10, na.rm=TRUE)

HVG_LVG_distribution_Expression_ZT24[,2] <- d_HVG$y
HVG_LVG_distribution_Expression_ZT24[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_Expression_ZT24) <- c("TPM", "HVG", "LVG")
HVG_LVG_distribution_Expression_ZT24_long <- gather(HVG_LVG_distribution_Expression_ZT24, type, distr, HVG:LVG)

# distribution

HVG_LVG_distribution_Expression_ZT24 <- mutate(HVG_LVG_distribution_Expression_ZT24, mean_bootstrat=rowMeans(random_distribution_Expression_ZT24_1000[,2:1001]),
                                               bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_Expression_ZT24_1000[,2:1001])),
                                               top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_Expression_ZT24_1000[,2:1001])),
                                               ZT= "ZT24")


# wilcoxon test

Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(ZT24_Expression_HVG$TPM_correspondingZT,random_Expression_ZT24_1000[,i])$p.value
}

names(Wilcox_pvalues) <- "result"

ggplot(Wilcox_pvalues, aes(x=result)) +
  geom_density()

# All values between 2.556529e-02 and 6.257754e-12 #



# all ZT combined and graph


allZT_boostrap_Expression <- bind_rows(HVG_LVG_distribution_Expression_ZT2,
                                       HVG_LVG_distribution_Expression_ZT4,
                                       HVG_LVG_distribution_Expression_ZT6,
                                       HVG_LVG_distribution_Expression_ZT8,
                                       HVG_LVG_distribution_Expression_ZT10,
                                       HVG_LVG_distribution_Expression_ZT12,
                                       HVG_LVG_distribution_Expression_ZT14,
                                       HVG_LVG_distribution_Expression_ZT16,
                                       HVG_LVG_distribution_Expression_ZT18,
                                       HVG_LVG_distribution_Expression_ZT20,
                                       HVG_LVG_distribution_Expression_ZT22,
                                       HVG_LVG_distribution_Expression_ZT24)



gather(allZT_boostrap_Expression, type, distr, HVG:top_CI, factor_key=TRUE) %>%
  ggplot(aes(x=TPM, y=distr, color=type, linetype=type)) +
  geom_line(size=1) +
  scale_x_continuous(limits = c(1.5, 10)) +
  facet_wrap(~ZT) +
  xlab("Corrected CV2 ZT24") +
  ylab("Density") +
  scale_color_manual(values=c("blue","green4","black","grey30","grey30")) +
  scale_linetype_manual(values=c("solid","solid","solid", "dashed", "dashed")) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line( size=.1, color="black" ),
        panel.background = element_rect(fill="white",colour="grey50"), 
        text=element_text(size=22))




### Plotting of the correlation between gene expression levels and variability profiles during the time course  ----

## Correlation expression vs variability. Spearman correlation ----

#data preparation


All_data <- read.table("Expression_level/BrenneckeMod_allZT_RUV2.txt", header=TRUE, sep='\t')
dim(All_data)

CorCV2 <- select(All_data, Gene, contains("BioVar"))
Ave_expr <- select(All_data, Gene, contains("Mean"))

Ave_expr <- dplyr::mutate(Ave_expr, 
                          MeanNorm_ZT2 = Mean_ZT2/rowMeans(Ave_expr[,2:13],na.rm = TRUE),
                          MeanNorm_ZT4 = Mean_ZT4/rowMeans(Ave_expr[,2:13],na.rm = TRUE),
                          MeanNorm_ZT6 = Mean_ZT6/rowMeans(Ave_expr[,2:13],na.rm = TRUE),
                          MeanNorm_ZT8 = Mean_ZT8/rowMeans(Ave_expr[,2:13],na.rm = TRUE),
                          MeanNorm_ZT10 = Mean_ZT10/rowMeans(Ave_expr[,2:13],na.rm = TRUE),
                          MeanNorm_ZT12 = Mean_ZT12/rowMeans(Ave_expr[,2:13],na.rm = TRUE),
                          MeanNorm_ZT14 = Mean_ZT14/rowMeans(Ave_expr[,2:13],na.rm = TRUE),
                          MeanNorm_ZT16 = Mean_ZT16/rowMeans(Ave_expr[,2:13],na.rm = TRUE),
                          MeanNorm_ZT18 = Mean_ZT18/rowMeans(Ave_expr[,2:13],na.rm = TRUE),
                          MeanNorm_ZT20 = Mean_ZT20/rowMeans(Ave_expr[,2:13],na.rm = TRUE),
                          MeanNorm_ZT22 = Mean_ZT22/rowMeans(Ave_expr[,2:13],na.rm = TRUE),
                          MeanNorm_ZT0 = Mean_ZT0/rowMeans(Ave_expr[,2:13],na.rm = TRUE))


Ave_expr_norm <- select(Ave_expr, Gene, contains("MeanNorm"))

CorCV2_bis <- column_to_rownames(CorCV2, var="Gene")
Ave_expr_bis <- column_to_rownames(Ave_expr, var="Gene")
Ave_expr_norm_bis <- column_to_rownames(Ave_expr_norm, var="Gene")


Correlation_BioVar_Average<- diag(cor(t(Ave_expr_norm_bis), t(CorCV2_bis), use="pairwise.complete.obs", method="spearman"))

Correlation_BioVar_Average2=as.data.frame(Correlation_BioVar_Average)

Correlation_BioVar_Average2<-rownames_to_column(Correlation_BioVar_Average2,var="tracking_id")



# extract 1000 sets of random genes ----


random_distribution_Correlation_1000=matrix(nrow=512,ncol=1000)
random_distribution_Correlation_1000=data.frame(random_distribution_Correlation_1000)

random_Correlation_1000=matrix(nrow=1358,ncol=1000)
random_Correlation_1000=data.frame(random_Correlation_1000)


for (i in 1:1000) {
  random_set <- Correlation_BioVar_Average2[sample(nrow(Correlation_BioVar_Average2), 1358),]
  random_Correlation_1000[,i] <- random_set[,2]
  d <- density(random_set[,2], from=-1, to=1)
  random_distribution_Correlation_1000[,i] <- d$y
}

colnames(random_Correlation_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(Correlation_BioVar_Average2$Correlation_BioVar_Average, from=-1, to=1)
random_distribution_Correlation_1000 <- cbind(d_background$x, random_distribution_Correlation_1000)
colnames(random_distribution_Correlation_1000) <- c("Correlation_BioVar_Average",paste("random_set",1:100, sep="_"))


HVG <- read.table("Expression_level/HVG_allZT_background.txt", header=TRUE, sep='\t')
dim(HVG)
[1] 1358    2

Variable_Correlation <- inner_join(HVG, Correlation_BioVar_Average2, by=c("allZT"="tracking_id")) 


LVG <- read.table("Expression_level/LVG1000_allZT_background.txt", header=TRUE, sep='\t')
dim(LVG)
[1] 5727    2

Stable_Correlation <- inner_join(LVG, Correlation_BioVar_Average2, by=c("allZT"="tracking_id")) 



HVG_LVG_distribution_Correlation=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_Correlation=data.frame(HVG_LVG_distribution_Correlation)

d_background <-  density(Correlation_BioVar_Average2$Correlation_BioVar_Average, from=-1, to=1)
HVG_LVG_distribution_Correlation[,1] <- d_background$x

d_HVG <- density(Variable_Correlation$Correlation_BioVar_Average, from=-1, to=1)
d_LVG <- density(Stable_Correlation$Correlation_BioVar_Average,from=-1, to=1)

HVG_LVG_distribution_Correlation[,2] <- d_HVG$y
HVG_LVG_distribution_Correlation[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_Correlation) <- c("Correlation", "HVG", "LVG")
HVG_LVG_distribution_Correlation_long <- gather(HVG_LVG_distribution_Correlation, type, distr, HVG:LVG)


# graph

HVG_LVG_distribution_Correlation <- mutate(HVG_LVG_distribution_Correlation, mean_bootstrat=rowMeans(random_distribution_Correlation_1000[,2:1001]),
                                           bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_Correlation_1000[,2:1001])),
                                           top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_Correlation_1000[,2:1001])))



gather(HVG_LVG_distribution_Correlation, type, distr, HVG:top_CI, factor_key=TRUE) %>%
  ggplot(aes(x=Correlation, y=distr, color=type, linetype=type)) +
  geom_line(size=1) +
  scale_x_continuous(limits = c(-1, 1)) +
  xlab("Spearman correlation") +
  ylab("Density") +
  scale_color_manual(values=c("blue","green4","black","grey30","grey30")) +
  scale_linetype_manual(values=c("solid","solid","solid", "dashed", "dashed")) +
  theme(panel.grid.major.x = element_blank(),panel.grid.major.y = element_line( size=.1, color="black" ) 
        , panel.background = element_rect(fill="white",colour="grey50"), text=element_text(size=22))

median(Variable_Correlation$Correlation_BioVar_Average)
[1] 0.1825758

median(Stable_Correlation$Correlation_BioVar_Average)
[1] -0.06315828


# median correlation


Median_correlation=matrix(nrow=1000,ncol=1)
Median_correlation=data.frame(Median_correlation)


for (i in 1:1000) {
  Median_correlation[i,] <- median(random_Correlation_1000[,i])
}

mean(Median_correlation$Median_correlation, na.rm=TRUE)
[1] -0.005166325

mean(Median_correlation$Median_correlation, na.rm=TRUE) - sd(Median_correlation$Median_correlation, na.rm=TRUE)
[1] -0.01863709

mean(Median_correlation$Median_correlation, na.rm=TRUE) + sd(Median_correlation$Median_correlation, na.rm=TRUE)
[1] 0.008304439

# wilcoxon test

Wilcox_stat=matrix(nrow=1000,ncol=1)
Wilcox_stat=data.frame(Wilcox_stat)


for (i in 1:1000) {
  Wilcox_stat[i,] <- wilcox.test(Variable_Correlation$Correlation_BioVar_Average,random_Correlation_1000[,i])$statistic
}

names(Wilcox_stat) <- "result"

ggplot(Wilcox_stat, aes(x=result)) +
  geom_density()


Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(Variable_Correlation$Correlation_BioVar_Average,random_Correlation_1000[,i])$p.value
}

names(Wilcox_pvalues) <- "result"

ggplot(Wilcox_pvalues, aes(x=result)) +
  geom_density()

#All values between 6.017999e-13 and 7.973866e-32, but for whatever reason doesn't plot it correclty#



# Analysis of genes with Cor<0.1 or Cor >0.1 ----

# Data preparation


All_data <- read.table("Expression_level/BrenneckeMod_allZT_RUV2.txt", header=TRUE, sep='\t')
dim(All_data)

HVG <- read.table("Expression_level/HVG_allZT_background.txt", header=TRUE, sep='\t')

HVG_expression <- inner_join(All_data, HVG, by=c("Gene"="allZT"))

CorCV2 <- select(HVG_expression, Gene, contains("BioVar"))
Ave_expr <- select(HVG_expression, Gene, contains("Mean"))

Ave_expr <- dplyr::mutate(Ave_expr, 
                          MeanNorm_ZT2 = Mean_ZT2/rowMeans(Ave_expr[,2:13],na.rm = TRUE),
                          MeanNorm_ZT4 = Mean_ZT4/rowMeans(Ave_expr[,2:13],na.rm = TRUE),
                          MeanNorm_ZT6 = Mean_ZT6/rowMeans(Ave_expr[,2:13],na.rm = TRUE),
                          MeanNorm_ZT8 = Mean_ZT8/rowMeans(Ave_expr[,2:13],na.rm = TRUE),
                          MeanNorm_ZT10 = Mean_ZT10/rowMeans(Ave_expr[,2:13],na.rm = TRUE),
                          MeanNorm_ZT12 = Mean_ZT12/rowMeans(Ave_expr[,2:13],na.rm = TRUE),
                          MeanNorm_ZT14 = Mean_ZT14/rowMeans(Ave_expr[,2:13],na.rm = TRUE),
                          MeanNorm_ZT16 = Mean_ZT16/rowMeans(Ave_expr[,2:13],na.rm = TRUE),
                          MeanNorm_ZT18 = Mean_ZT18/rowMeans(Ave_expr[,2:13],na.rm = TRUE),
                          MeanNorm_ZT20 = Mean_ZT20/rowMeans(Ave_expr[,2:13],na.rm = TRUE),
                          MeanNorm_ZT22 = Mean_ZT22/rowMeans(Ave_expr[,2:13],na.rm = TRUE),
                          MeanNorm_ZT0 = Mean_ZT0/rowMeans(Ave_expr[,2:13],na.rm = TRUE))


Ave_expr_norm <- select(Ave_expr, Gene, contains("MeanNorm"))

CorCV2_bis <- column_to_rownames(CorCV2, var="Gene")
Ave_expr_bis <- column_to_rownames(Ave_expr, var="Gene")
Ave_expr_norm_bis <- column_to_rownames(Ave_expr_norm, var="Gene")


Correlation_BioVar_Average<- diag(cor(t(Ave_expr_norm_bis), t(CorCV2_bis), use="pairwise.complete.obs", method="spearman"))

Correlation_BioVar_Average2=as.data.frame(Correlation_BioVar_Average)

Correlation_BioVar_Average2<-rownames_to_column(Correlation_BioVar_Average2,var="tracking_id")



HVG_left_peak <-  filter(Correlation_BioVar_Average2,  Correlation_BioVar_Average <= 0.1)

HVG_right_peak <-  filter(Correlation_BioVar_Average2, Correlation_BioVar_Average > 0.1)

# heatmap expression profiles
purpleblackyellow <- colorRampPalette(c("purple","black", "yellow"))(n = 100)

ColDay=c("orange","orange","orange","orange","orange","red","black","black","black","black","black","blue")

ZT24 <- read.table("Expression_level/BrenneckeMod_ZT0_RUV2.txt",header=TRUE,sep='\t')
ZT2 <- read.table("Expression_level/BrenneckeMod_ZT2_RUV2.txt",header=TRUE,sep='\t')
ZT4 <- read.table("Expression_level/BrenneckeMod_ZT4_RUV2.txt",header=TRUE,sep='\t')
ZT6 <- read.table("Expression_level/BrenneckeMod_ZT6_RUV2.txt",header=TRUE,sep='\t')
ZT8 <- read.table("Expression_level/BrenneckeMod_ZT8_RUV2.txt",header=TRUE,sep='\t')
ZT10 <- read.table("Expression_level/BrenneckeMod_ZT10_RUV2.txt",header=TRUE,sep='\t')
ZT12 <- read.table("Expression_level/BrenneckeMod_ZT12_RUV2.txt",header=TRUE,sep='\t')
ZT14 <- read.table("Expression_level/BrenneckeMod_ZT14_RUV2.txt",header=TRUE,sep='\t')
ZT16 <- read.table("Expression_level/BrenneckeMod_ZT16_RUV2.txt",header=TRUE,sep='\t')
ZT18 <- read.table("Expression_level/BrenneckeMod_ZT18_RUV2.txt",header=TRUE,sep='\t')
ZT20 <- read.table("Expression_level/BrenneckeMod_ZT20_RUV2.txt",header=TRUE,sep='\t')
ZT22 <- read.table("Expression_level/BrenneckeMod_ZT22_RUV2.txt",header=TRUE,sep='\t')

AllZT <- Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2, by="Gene"),
                list(ZT2,ZT4,ZT6,ZT8,ZT10,ZT12,
                     ZT14,ZT16,ZT18,ZT20,ZT22,ZT24))

AllZT_Mean <- select(AllZT, Gene, contains("Mean"))
dim(AllZT_Mean)
[1] 15646    13

AllZT_Mean_meanNorm <- mutate(AllZT_Mean,
                              Mean_ZT2_meanNorm=Mean_ZT2/rowMeans(AllZT_Mean[,2:13],na.rm=TRUE),
                              Mean_ZT4_meanNorm=Mean_ZT4/rowMeans(AllZT_Mean[,2:13],na.rm=TRUE),
                              Mean_ZT6_meanNorm=Mean_ZT6/rowMeans(AllZT_Mean[,2:13],na.rm=TRUE),
                              Mean_ZT8_meanNorm=Mean_ZT8/rowMeans(AllZT_Mean[,2:13],na.rm=TRUE),
                              Mean_ZT10_meanNorm=Mean_ZT10/rowMeans(AllZT_Mean[,2:13],na.rm=TRUE),
                              Mean_ZT12_meanNorm=Mean_ZT12/rowMeans(AllZT_Mean[,2:13],na.rm=TRUE),
                              Mean_ZT14_meanNorm=Mean_ZT14/rowMeans(AllZT_Mean[,2:13],na.rm=TRUE),
                              Mean_ZT16_meanNorm=Mean_ZT16/rowMeans(AllZT_Mean[,2:13],na.rm=TRUE),
                              Mean_ZT18_meanNorm=Mean_ZT18/rowMeans(AllZT_Mean[,2:13],na.rm=TRUE),
                              Mean_ZT20_meanNorm=Mean_ZT20/rowMeans(AllZT_Mean[,2:13],na.rm=TRUE),
                              Mean_ZT22_meanNorm=Mean_ZT22/rowMeans(AllZT_Mean[,2:13],na.rm=TRUE),
                              Mean_ZT0_meanNorm=Mean_ZT0/rowMeans(AllZT_Mean[,2:13],na.rm=TRUE))


HVG_left_peak_exp_meanNorm <- inner_join(HVG_left_peak, AllZT_Mean_meanNorm,
                                         b=c("tracking_id"="Gene"))


HVG_left_peak_BioVar.pearson <- cor(t(HVG_left_peak_exp_meanNorm[,c(15:26)]), method = "pearson")
HVG_left_peak_BioVar.pearson[is.na(HVG_left_peak_BioVar.pearson)]<-0
HVG_left_peak_BioVar.pearson.dist<-as.dist((1-HVG_left_peak_BioVar.pearson))
HVG_left_peak_BioVar.pearson.dist[is.na(HVG_left_peak_BioVar.pearson.dist)]<-0
HVG_left_peak_BioVar.pearson.dist.clust.complete <- hclust(HVG_left_peak_BioVar.pearson.dist, method="complete")
HVG_left_peak_BioVar.pearson.dist.clust.complete.den <- as.dendrogram(HVG_left_peak_BioVar.pearson.dist.clust.complete)
HVG_left_peak_BioVar.pearson.dist.clust.complete.col=brewer.pal(12, 'Set3')[cutree(HVG_left_peak_BioVar.pearson.dist.clust.complete, k = 5)]
HVG_left_peak_exp_meanNorm$cluster_5_pearson_color <-brewer.pal(12, 'Set3')[cutree(HVG_left_peak_BioVar.pearson.dist.clust.complete, k = 5)]


heatmap.2(as.matrix(HVG_left_peak_exp_meanNorm[,15:26]), Colv=NA, dendrogram='row', 
          Rowv=HVG_left_peak_BioVar.pearson.dist.clust.complete.den, 
          col= purpleblackyellow,breaks=seq(0.8,1.2,length.out=101), trace = 'none' ,
          labRow=NA, labCol=names(HVG_left_peak_exp_meanNorm[,16:27]), 
          main="Expression, HVG pearson corr <0.1",ColSideColors=ColDay,
          RowSideColors = HVG_left_peak_BioVar.pearson.dist.clust.complete.col,
          na.rm=F, na.color="black")





HVG_right_peak_exp_meanNorm <- inner_join(HVG_right_peak, AllZT_Mean_meanNorm,
                                          b=c("tracking_id"="Gene"))


HVG_right_peak_BioVar.pearson <- cor(t(HVG_right_peak_exp_meanNorm[,c(15:26)]), method = "pearson")
HVG_right_peak_BioVar.pearson[is.na(HVG_right_peak_BioVar.pearson)]<-0
HVG_right_peak_BioVar.pearson.dist<-as.dist((1-HVG_right_peak_BioVar.pearson))
HVG_right_peak_BioVar.pearson.dist[is.na(HVG_right_peak_BioVar.pearson.dist)]<-0
HVG_right_peak_BioVar.pearson.dist.clust.complete <- hclust(HVG_right_peak_BioVar.pearson.dist, method="complete")
HVG_right_peak_BioVar.pearson.dist.clust.complete.den <- as.dendrogram(HVG_right_peak_BioVar.pearson.dist.clust.complete)
HVG_right_peak_BioVar.pearson.dist.clust.complete.col=brewer.pal(12, 'Set3')[cutree(HVG_right_peak_BioVar.pearson.dist.clust.complete, k = 5)]
HVG_right_peak_exp_meanNorm$cluster_5_pearson_color <-brewer.pal(12, 'Set3')[cutree(HVG_right_peak_BioVar.pearson.dist.clust.complete, k = 5)]


heatmap.2(as.matrix(HVG_right_peak_exp_meanNorm[,15:26]), Colv=NA, dendrogram='row', 
          Rowv=HVG_right_peak_BioVar.pearson.dist.clust.complete.den, 
          col= purpleblackyellow,breaks=seq(0.8,1.2,length.out=101), trace = 'none' ,
          labRow=NA, labCol=names(HVG_right_peak_exp_meanNorm[,16:27]), 
          main="Expression, HVG pearson corr >0.1",ColSideColors=ColDay,
          RowSideColors = HVG_right_peak_BioVar.pearson.dist.clust.complete.col,
          na.rm=F, na.color="black")



# boxplot expression level for each ZT


Correlation_BioVar_Average2 <-  mutate(Correlation_BioVar_Average2, left_peak = Correlation_BioVar_Average <= 0.1)

HVG_exp_meanNorm <- inner_join(Correlation_BioVar_Average2, AllZT_Mean_meanNorm,
                               b=c("tracking_id"="Gene"))

gather(HVG_exp_meanNorm, ZT, expression, Mean_ZT2:Mean_ZT0, factor_key = T) %>%
  ggplot(aes(x=ZT, y=expression, fill=left_peak)) +
  geom_boxplot(outlier.alpha = 0, notch = TRUE) +
  scale_y_log10(limits=c(5,2000)) +
  scale_fill_manual("", values=c("gray30", "gray80"),
                    labels=c("FALSE"="Cor > 0.1",
                             "TRUE"="Cor < 0.1")) +
  scale_x_discrete(labels=c("ZT2","ZT4","ZT6","ZT8","ZT10","ZT12",
                            "ZT14","ZT16","ZT18","ZT20","ZT22","ZT24")) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line( size=.1, color="black" ) , 
        panel.background = element_rect(fill="white",colour="grey50"), 
        text=element_text(size=21),
        axis.text.x = element_text(angle = 90))


# Distribution nb ZT variable for each group


Correlation_BioVar_Average2 <-  mutate(Correlation_BioVar_Average2, left_peak = Correlation_BioVar_Average <= 0.1)

nbZT_variable <- read.table("Expression_level/HVG_allZT_nbZT.txt",header=TRUE, sep='\t')

HVG_nbZT_variable <- inner_join(Correlation_BioVar_Average2, nbZT_variable,
                                by=c("tracking_id"="allZT"))

ggplot(HVG_nbZT_variable, aes(x=nbZT, color=left_peak)) +
  geom_freqpoly(bins=12,size=1.5) +
  xlab("Number of ZT highly variable") +
  xlim(0,11.5) +
  ylab("Proportion of genes") +
  ggtitle("") +
  scale_color_manual("", values=c("gray30", "gray80"),
                     labels=c("FALSE"="Cor > 0.1",
                              "TRUE"="Cor < 0.1")) +
  theme_bw() +
  theme(legend.text=element_text(size=22),legend.title = element_text(size=22, face="bold"),
        axis.text=element_text(size=22), axis.title=element_text(size=22),
        title=element_text(size=22)) 


