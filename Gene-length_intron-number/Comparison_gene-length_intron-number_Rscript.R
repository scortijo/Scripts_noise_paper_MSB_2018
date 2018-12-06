
library(tidyverse)

### Direct comparison of gene length and intron number ----

# Data import and formatting

Gene_nbIntrons <- read.table("Gene-length_intron-number/TAIR10_Gene_nbIntrons.txt", header=TRUE, sep='\t')

Gene_size <- read.table("Gene-length_intron-number/TAIR10_Gene_sizes.txt", header=TRUE, sep='\t')

Gene_size_nbIntrons <- merge(Gene_nbIntrons, Gene_size,
                                    by.x="gene_name",by.y="gene_name")

# Plotting

ggplot(Gene_size_nbIntrons, aes(x=Number_of_introns,y=length)) +
  geom_point()+
  xlim(0,150)+
  ylim(0,20000) +
  xlab("Number of Introns") +
  ylab("Gene length")+
  theme(text=element_text(size=24)) 

cor(Gene_size_nbIntrons$Number_of_introns,Gene_size_nbIntrons$length)


### Plotting of the gene length for HVGs, LVGs and a 1000 sets of random genes with a fixed number of introns ----

#data preparation

Gene_size <- read.table("Gene-length_intron-number/TAIR10_Gene_sizes.txt", header=TRUE, sep='\t')
dim(Gene_size)

Background <- read.table("Gene-length_intron-number/All_genes_background.txt", header=TRUE, sep='\t')
dim(Background)

Background_gene_size <- inner_join(Background, Gene_size, by=c("Gene"="gene_name")) %>%
  select(Gene, length)


Gene_nbIntrons <- read.table("Gene-length_intron-number/TAIR10_Gene_nbIntrons.txt", header=TRUE, sep='\t')
dim(Gene_nbIntrons)


Background_gene_size_nbIntron <- inner_join(Background_gene_size, Gene_nbIntrons,
                                            by=c("Gene"="gene_name"))

# no intron ----

subset_0intron <- filter(Background_gene_size_nbIntron,
                         Number_of_introns==0)
dim(subset_0intron)

# extract 1000 sets of random genes


random_distribution_0intron_1000=matrix(nrow=512,ncol=1000)
random_distribution_0intron_1000=data.frame(random_distribution_0intron_1000)

random_0intron_1000=matrix(nrow=32,ncol=1000)
random_0intron_1000=data.frame(random_0intron_1000)


for (i in 1:1000) {
  random_set <- subset_0intron[sample(nrow(subset_0intron), 32),]
  random_0intron_1000[,i] <- random_set[,2]
  d <- density(random_set[,2], from=153, to=26435)
  random_distribution_0intron_1000[,i] <- d$y
}

colnames(random_0intron_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(subset_0intron$length, from=153, to=26435)
random_distribution_0intron_1000 <- cbind(d_background$x, random_distribution_0intron_1000)
colnames(random_distribution_0intron_1000) <- c("length",paste("random_set",1:100, sep="_"))


HVG <- read.table("Gene-length_intron-number/HVG_allZT_background.txt", header=TRUE, sep='\t')
dim(HVG)
[1] 1358    2

Variable_0intron <- inner_join(HVG, subset_0intron, by=c("allZT"="Gene")) %>%
  select(allZT, length) 
dim(Variable_0intron)
[1] 32  2

LVG <- read.table("Gene-length_intron-number/LVG1000_allZT_background.txt", header=TRUE, sep='\t')
dim(LVG)
[1] 5727    2

Stable_0intron <- inner_join(LVG, subset_0intron, by=c("allZT"="Gene")) %>%
  select(allZT, length) 
dim(Stable_0intron)
[1] 23  2


HVG_LVG_distribution_0intron=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_0intron=data.frame(HVG_LVG_distribution_0intron)

d_background <-  density(subset_0intron$length, from=153, to=26435)
HVG_LVG_distribution_0intron[,1] <- d_background$x

d_HVG <- density(Variable_0intron$length, from=153, to=26435)
d_LVG <- density(Stable_0intron$length, from=153, to=26435)

HVG_LVG_distribution_0intron[,2] <- d_HVG$y
HVG_LVG_distribution_0intron[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_0intron) <- c("length", "HVG", "LVG")
HVG_LVG_distribution_0intron_long <- gather(HVG_LVG_distribution_0intron, type, distr, HVG:LVG)

# graph

HVG_LVG_distribution_0intron <- mutate(HVG_LVG_distribution_0intron, mean_bootstrat=rowMeans(random_distribution_0intron_1000[,2:1001]),
                                       bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_0intron_1000[,2:1001])),
                                       top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_0intron_1000[,2:1001])))



gather(HVG_LVG_distribution_0intron, type, distr, HVG:top_CI, factor_key=TRUE) %>%
  ggplot(aes(x=length, y=distr, color=type, linetype=type)) +
  geom_line(size=1) +
  scale_x_continuous(limits = c(0, 8000)) +
  xlab("Gene length (bp)") +
  ylab("Density") +
  ggtitle(" no intron") +
  scale_color_manual(values=c("blue","green4","black","grey30","grey30")) +
  scale_linetype_manual(values=c("solid","solid","solid", "dashed", "dashed")) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line( size=.1, color="black" ), 
        panel.background = element_rect(fill="white",colour="grey50"), 
        text=element_text(size=38))

# wilcoxon test


Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(Variable_0intron$length,random_0intron_1000[,i])$p.value
}

wilcox.test(Variable_0intron$length,Stable_0intron$length)$p.value
[1] 5.032923e-07

# All values between 0.088115329 and 2.692207e-06 #




# 1 intron ----

subset_1intron <- filter(Background_gene_size_nbIntron,
                         Number_of_introns==1)
dim(subset_1intron)

# extract 1000 sets of random genes


random_distribution_1intron_1000=matrix(nrow=512,ncol=1000)
random_distribution_1intron_1000=data.frame(random_distribution_1intron_1000)

random_1intron_1000=matrix(nrow=362,ncol=1000)
random_1intron_1000=data.frame(random_1intron_1000)


for (i in 1:1000) {
  random_set <- subset_1intron[sample(nrow(subset_1intron), 362),]
  random_1intron_1000[,i] <- random_set[,2]
  d <- density(random_set[,2], from=153, to=26435)
  random_distribution_1intron_1000[,i] <- d$y
}

colnames(random_1intron_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(subset_1intron$length, from=153, to=26435)
random_distribution_1intron_1000 <- cbind(d_background$x, random_distribution_1intron_1000)
colnames(random_distribution_1intron_1000) <- c("length",paste("random_set",1:100, sep="_"))


HVG <- read.table("Gene-length_intron-number/HVG_allZT_background.txt", header=TRUE, sep='\t')
dim(HVG)
[1] 1358    2

Variable_1intron <- inner_join(HVG, subset_1intron, by=c("allZT"="Gene")) %>%
  select(allZT, length) 
dim(Variable_1intron)
[1] 362  2

LVG <- read.table("Gene-length_intron-number/LVG1000_allZT_background.txt", header=TRUE, sep='\t')
dim(LVG)
[1] 5727    2

Stable_1intron <- inner_join(LVG, subset_1intron, by=c("allZT"="Gene")) %>%
  select(allZT, length) 
dim(Stable_1intron)
[1] 594  2


HVG_LVG_distribution_1intron=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_1intron=data.frame(HVG_LVG_distribution_1intron)

d_background <-  density(subset_1intron$length, from=153, to=26435)
HVG_LVG_distribution_1intron[,1] <- d_background$x

d_HVG <- density(Variable_1intron$length, from=153, to=26435)
d_LVG <- density(Stable_1intron$length, from=153, to=26435)

HVG_LVG_distribution_1intron[,2] <- d_HVG$y
HVG_LVG_distribution_1intron[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_1intron) <- c("length", "HVG", "LVG")
HVG_LVG_distribution_1intron_long <- gather(HVG_LVG_distribution_1intron, type, distr, HVG:LVG)

# graph

HVG_LVG_distribution_1intron <- mutate(HVG_LVG_distribution_1intron, mean_bootstrat=rowMeans(random_distribution_1intron_1000[,2:1001]),
                                       bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_1intron_1000[,2:1001])),
                                       top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_1intron_1000[,2:1001])))



gather(HVG_LVG_distribution_1intron, type, distr, HVG:top_CI, factor_key=TRUE) %>%
  ggplot(aes(x=length, y=distr, color=type, linetype=type)) +
  geom_line(size=1) +
  scale_x_continuous(limits = c(0, 8000)) +
  xlab("Gene length (bp)") +
  ylab("Density") +
  ggtitle("1 intron") +
  scale_color_manual(values=c("blue","green4","black","grey30","grey30")) +
  scale_linetype_manual(values=c("solid","solid","solid", "dashed", "dashed")) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line( size=.1, color="black" ), 
        panel.background = element_rect(fill="white",colour="grey50"), 
        text=element_text(size=38))


# wilcoxon test


Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(Variable_1intron$length,random_1intron_1000[,i])$p.value
}

wilcox.test(Variable_1intron$length,Stable_1intron$length)$p.value
[1] 2.735142e-76

# All values between 2.775528e-11 and 9.102639e-26 #


# 2 introns ----

subset_2intron <- filter(Background_gene_size_nbIntron,
                         Number_of_introns==2)
dim(subset_2intron)

# extract 1000 sets of random genes


random_distribution_2intron_1000=matrix(nrow=512,ncol=1000)
random_distribution_2intron_1000=data.frame(random_distribution_2intron_1000)

random_2intron_1000=matrix(nrow=255,ncol=1000)
random_2intron_1000=data.frame(random_2intron_1000)


for (i in 1:1000) {
  random_set <- subset_2intron[sample(nrow(subset_2intron), 255),]
  random_2intron_1000[,i] <- random_set[,2]
  d <- density(random_set[,2], from=153, to=26435)
  random_distribution_2intron_1000[,i] <- d$y
}

colnames(random_2intron_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(subset_2intron$length, from=153, to=26435)
random_distribution_2intron_1000 <- cbind(d_background$x, random_distribution_2intron_1000)
colnames(random_distribution_2intron_1000) <- c("length",paste("random_set",1:100, sep="_"))


HVG <- read.table("Gene-length_intron-number/HVG_allZT_background.txt", header=TRUE, sep='\t')
dim(HVG)
[1] 1358    2

Variable_2intron <- inner_join(HVG, subset_2intron, by=c("allZT"="Gene")) %>%
  select(allZT, length) 
dim(Variable_2intron)
[1] 255  2

LVG <- read.table("Gene-length_intron-number/LVG1000_allZT_background.txt", header=TRUE, sep='\t')
dim(LVG)
[1] 5727    2

Stable_2intron <- inner_join(LVG, subset_2intron, by=c("allZT"="Gene")) %>%
  select(allZT, length) 
dim(Stable_2intron)
[1] 428  2


HVG_LVG_distribution_2intron=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_2intron=data.frame(HVG_LVG_distribution_2intron)

d_background <-  density(subset_2intron$length, from=153, to=26435)
HVG_LVG_distribution_2intron[,1] <- d_background$x

d_HVG <- density(Variable_2intron$length, from=153, to=26435)
d_LVG <- density(Stable_2intron$length, from=153, to=26435)

HVG_LVG_distribution_2intron[,2] <- d_HVG$y
HVG_LVG_distribution_2intron[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_2intron) <- c("length", "HVG", "LVG")
HVG_LVG_distribution_2intron_long <- gather(HVG_LVG_distribution_2intron, type, distr, HVG:LVG)

# graph

HVG_LVG_distribution_2intron <- mutate(HVG_LVG_distribution_2intron, mean_bootstrat=rowMeans(random_distribution_2intron_1000[,2:1001]),
                                       bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_2intron_1000[,2:1001])),
                                       top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_2intron_1000[,2:1001])))



gather(HVG_LVG_distribution_2intron, type, distr, HVG:top_CI, factor_key=TRUE) %>%
  ggplot(aes(x=length, y=distr, color=type, linetype=type)) +
  geom_line(size=1) +
  scale_x_continuous(limits = c(0, 8000)) +
  xlab("Gene length (bp)") +
  ylab("Density") +
  ggtitle("2 introns") +
  scale_color_manual(values=c("blue","green4","black","grey30","grey30")) +
  scale_linetype_manual(values=c("solid","solid","solid", "dashed", "dashed")) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line( size=.1, color="black" ), 
        panel.background = element_rect(fill="white",colour="grey50"), 
        text=element_text(size=38))


# wilcoxon test

wilcox.test(Variable_2intron$length,Stable_2intron$length)$p.value
[1] 4.087408e-50


Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(Variable_2intron$length,random_2intron_1000[,i])$p.value
}

# All values between 1.440280e-06 and 1.671117e-18 #



# 3 introns ----

subset_3intron <- filter(Background_gene_size_nbIntron,
                         Number_of_introns==3)
dim(subset_3intron)

# extract 1000 sets of random genes


random_distribution_3intron_1000=matrix(nrow=512,ncol=1000)
random_distribution_3intron_1000=data.frame(random_distribution_3intron_1000)

random_3intron_1000=matrix(nrow=172,ncol=1000)
random_3intron_1000=data.frame(random_3intron_1000)


for (i in 1:1000) {
  random_set <- subset_3intron[sample(nrow(subset_3intron), 172),]
  random_3intron_1000[,i] <- random_set[,2]
  d <- density(random_set[,2], from=153, to=26435)
  random_distribution_3intron_1000[,i] <- d$y
}

colnames(random_3intron_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(subset_3intron$length, from=153, to=26435)
random_distribution_3intron_1000 <- cbind(d_background$x, random_distribution_3intron_1000)
colnames(random_distribution_3intron_1000) <- c("length",paste("random_set",1:100, sep="_"))


HVG <- read.table("Gene-length_intron-number/HVG_allZT_background.txt", header=TRUE, sep='\t')
dim(HVG)
[1] 1358    2

Variable_3intron <- inner_join(HVG, subset_3intron, by=c("allZT"="Gene")) %>%
  select(allZT, length) 
dim(Variable_3intron)
[1] 172  2

LVG <- read.table("Gene-length_intron-number/LVG1000_allZT_background.txt", header=TRUE, sep='\t')
dim(LVG)
[1] 5727    2

Stable_3intron <- inner_join(LVG, subset_3intron, by=c("allZT"="Gene")) %>%
  select(allZT, length) 
dim(Stable_3intron)
[1] 366  2


HVG_LVG_distribution_3intron=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_3intron=data.frame(HVG_LVG_distribution_3intron)

d_background <-  density(subset_3intron$length, from=153, to=26435)
HVG_LVG_distribution_3intron[,1] <- d_background$x

d_HVG <- density(Variable_3intron$length, from=153, to=26435)
d_LVG <- density(Stable_3intron$length, from=153, to=26435)

HVG_LVG_distribution_3intron[,2] <- d_HVG$y
HVG_LVG_distribution_3intron[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_3intron) <- c("length", "HVG", "LVG")
HVG_LVG_distribution_3intron_long <- gather(HVG_LVG_distribution_3intron, type, distr, HVG:LVG)

# graph

HVG_LVG_distribution_3intron <- mutate(HVG_LVG_distribution_3intron, mean_bootstrat=rowMeans(random_distribution_3intron_1000[,2:1001]),
                                       bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_3intron_1000[,2:1001])),
                                       top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_3intron_1000[,2:1001])))



gather(HVG_LVG_distribution_3intron, type, distr, HVG:top_CI, factor_key=TRUE) %>%
  ggplot(aes(x=length, y=distr, color=type, linetype=type)) +
  geom_line(size=1) +
  scale_x_continuous(limits = c(0, 8000)) +
  xlab("Gene length (bp)") +
  ylab("Density") +
  ggtitle("3 introns") +
  scale_color_manual(values=c("blue","green4","black","grey30","grey30")) +
  scale_linetype_manual(values=c("solid","solid","solid", "dashed", "dashed")) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line( size=.1, color="black" ), 
        panel.background = element_rect(fill="white",colour="grey50"), 
        text=element_text(size=38))


# wilcoxon test

wilcox.test(Variable_3intron$length,Stable_3intron$length)$p.value
[1] 1.506918e-24


Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(Variable_3intron$length,random_3intron_1000[,i])$p.value
}

# All values between 0.0371065926 and 1.992590e-09 #

# 4 introns ----

subset_4intron <- filter(Background_gene_size_nbIntron,
                         Number_of_introns==4)
dim(subset_4intron)

# extract 1000 sets of random genes


random_distribution_4intron_1000=matrix(nrow=512,ncol=1000)
random_distribution_4intron_1000=data.frame(random_distribution_4intron_1000)

random_4intron_1000=matrix(nrow=123,ncol=1000)
random_4intron_1000=data.frame(random_4intron_1000)


for (i in 1:1000) {
  random_set <- subset_4intron[sample(nrow(subset_4intron), 123),]
  random_4intron_1000[,i] <- random_set[,2]
  d <- density(random_set[,2], from=153, to=26435)
  random_distribution_4intron_1000[,i] <- d$y
}

colnames(random_4intron_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(subset_4intron$length, from=153, to=26435)
random_distribution_4intron_1000 <- cbind(d_background$x, random_distribution_4intron_1000)
colnames(random_distribution_4intron_1000) <- c("length",paste("random_set",1:100, sep="_"))


HVG <- read.table("Gene-length_intron-number/HVG_allZT_background.txt", header=TRUE, sep='\t')
dim(HVG)
[1] 1358    2

Variable_4intron <- inner_join(HVG, subset_4intron, by=c("allZT"="Gene")) %>%
  select(allZT, length) 
dim(Variable_4intron)
[1] 123  2

LVG <- read.table("Gene-length_intron-number/LVG1000_allZT_background.txt", header=TRUE, sep='\t')
dim(LVG)
[1] 5727    2

Stable_4intron <- inner_join(LVG, subset_4intron, by=c("allZT"="Gene")) %>%
  select(allZT, length) 
dim(Stable_4intron)
[1] 384  2


HVG_LVG_distribution_4intron=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_4intron=data.frame(HVG_LVG_distribution_4intron)

d_background <-  density(subset_4intron$length, from=153, to=26435)
HVG_LVG_distribution_4intron[,1] <- d_background$x

d_HVG <- density(Variable_4intron$length, from=153, to=26435)
d_LVG <- density(Stable_4intron$length, from=153, to=26435)

HVG_LVG_distribution_4intron[,2] <- d_HVG$y
HVG_LVG_distribution_4intron[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_4intron) <- c("length", "HVG", "LVG")
HVG_LVG_distribution_4intron_long <- gather(HVG_LVG_distribution_4intron, type, distr, HVG:LVG)

# graph

HVG_LVG_distribution_4intron <- mutate(HVG_LVG_distribution_4intron, mean_bootstrat=rowMeans(random_distribution_4intron_1000[,2:1001]),
                                       bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_4intron_1000[,2:1001])),
                                       top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_4intron_1000[,2:1001])))



gather(HVG_LVG_distribution_4intron, type, distr, HVG:top_CI, factor_key=TRUE) %>%
  ggplot(aes(x=length, y=distr, color=type, linetype=type)) +
  geom_line(size=1) +
  scale_x_continuous(limits = c(0, 8000)) +
  xlab("Gene length (bp)") +
  ylab("Density") +
  ggtitle("4 introns") +
  scale_color_manual(values=c("blue","green4","black","grey30","grey30")) +
  scale_linetype_manual(values=c("solid","solid","solid", "dashed", "dashed")) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line( size=.1, color="black" ), 
        panel.background = element_rect(fill="white",colour="grey50"), 
        text=element_text(size=38))


# wilcoxon test

wilcox.test(Variable_4intron$length,Stable_4intron$length)$p.value
[1] 1.223043e-12


Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(Variable_4intron$length,random_4intron_1000[,i])$p.value
}

# All values between 0.57731267 and 4.278270e-07 #



# 5 introns ----

subset_5intron <- filter(Background_gene_size_nbIntron,
                         Number_of_introns==5)
dim(subset_5intron)

# extract 1000 sets of random genes


random_distribution_5intron_1000=matrix(nrow=512,ncol=1000)
random_distribution_5intron_1000=data.frame(random_distribution_5intron_1000)

random_5intron_1000=matrix(nrow=67,ncol=1000)
random_5intron_1000=data.frame(random_5intron_1000)


for (i in 1:1000) {
  random_set <- subset_5intron[sample(nrow(subset_5intron), 67),]
  random_5intron_1000[,i] <- random_set[,2]
  d <- density(random_set[,2], from=153, to=26435)
  random_distribution_5intron_1000[,i] <- d$y
}

colnames(random_5intron_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(subset_5intron$length, from=153, to=26435)
random_distribution_5intron_1000 <- cbind(d_background$x, random_distribution_5intron_1000)
colnames(random_distribution_5intron_1000) <- c("length",paste("random_set",1:100, sep="_"))


HVG <- read.table("Gene-length_intron-number/HVG_allZT_background.txt", header=TRUE, sep='\t')
dim(HVG)
[1] 1358    2

Variable_5intron <- inner_join(HVG, subset_5intron, by=c("allZT"="Gene")) %>%
  select(allZT, length) 
dim(Variable_5intron)
[1] 67  2

LVG <- read.table("Gene-length_intron-number/LVG1000_allZT_background.txt", header=TRUE, sep='\t')
dim(LVG)
[1] 5727    2

Stable_5intron <- inner_join(LVG, subset_5intron, by=c("allZT"="Gene")) %>%
  select(allZT, length) 
dim(Stable_5intron)
[1] 354  2


HVG_LVG_distribution_5intron=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_5intron=data.frame(HVG_LVG_distribution_5intron)

d_background <-  density(subset_5intron$length, from=153, to=26435)
HVG_LVG_distribution_5intron[,1] <- d_background$x

d_HVG <- density(Variable_5intron$length, from=153, to=26435)
d_LVG <- density(Stable_5intron$length, from=153, to=26435)

HVG_LVG_distribution_5intron[,2] <- d_HVG$y
HVG_LVG_distribution_5intron[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_5intron) <- c("length", "HVG", "LVG")
HVG_LVG_distribution_5intron_long <- gather(HVG_LVG_distribution_5intron, type, distr, HVG:LVG)

# graph

HVG_LVG_distribution_5intron <- mutate(HVG_LVG_distribution_5intron, mean_bootstrat=rowMeans(random_distribution_5intron_1000[,2:1001]),
                                       bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_5intron_1000[,2:1001])),
                                       top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_5intron_1000[,2:1001])))



gather(HVG_LVG_distribution_5intron, type, distr, HVG:top_CI, factor_key=TRUE) %>%
  ggplot(aes(x=length, y=distr, color=type, linetype=type)) +
  geom_line(size=1) +
  scale_x_continuous(limits = c(0, 8000)) +
  xlab("Gene length (bp)") +
  ylab("Density") +
  ggtitle("5 introns") +
  scale_color_manual(values=c("blue","green4","black","grey30","grey30")) +
  scale_linetype_manual(values=c("solid","solid","solid", "dashed", "dashed")) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line( size=.1, color="black" ), 
        panel.background = element_rect(fill="white",colour="grey50"), 
        text=element_text(size=38))


# wilcoxon test

wilcox.test(Variable_5intron$length,Stable_5intron$length)$p.value
[1] 1.223043e-12


Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(Variable_5intron$length,random_5intron_1000[,i])$p.value
}

# All values between 0.57731267 and 4.278270e-07 #



# 6 introns ----

subset_6intron <- filter(Background_gene_size_nbIntron,
                         Number_of_introns==6)
dim(subset_6intron)

# extract 1000 sets of random genes


random_distribution_6intron_1000=matrix(nrow=512,ncol=1000)
random_distribution_6intron_1000=data.frame(random_distribution_6intron_1000)

random_6intron_1000=matrix(nrow=62,ncol=1000)
random_6intron_1000=data.frame(random_6intron_1000)


for (i in 1:1000) {
  random_set <- subset_6intron[sample(nrow(subset_6intron), 62),]
  random_6intron_1000[,i] <- random_set[,2]
  d <- density(random_set[,2], from=153, to=26435)
  random_distribution_6intron_1000[,i] <- d$y
}

colnames(random_6intron_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(subset_6intron$length, from=153, to=26435)
random_distribution_6intron_1000 <- cbind(d_background$x, random_distribution_6intron_1000)
colnames(random_distribution_6intron_1000) <- c("length",paste("random_set",1:100, sep="_"))


HVG <- read.table("Gene-length_intron-number/HVG_allZT_background.txt", header=TRUE, sep='\t')
dim(HVG)
[1] 1358    2

Variable_6intron <- inner_join(HVG, subset_6intron, by=c("allZT"="Gene")) %>%
  select(allZT, length) 
dim(Variable_6intron)
[1] 62  2

LVG <- read.table("Gene-length_intron-number/LVG1000_allZT_background.txt", header=TRUE, sep='\t')
dim(LVG)
[1] 5727    2

Stable_6intron <- inner_join(LVG, subset_6intron, by=c("allZT"="Gene")) %>%
  select(allZT, length) 
dim(Stable_6intron)
[1] 338  2


HVG_LVG_distribution_6intron=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_6intron=data.frame(HVG_LVG_distribution_6intron)

d_background <-  density(subset_6intron$length, from=153, to=26435)
HVG_LVG_distribution_6intron[,1] <- d_background$x

d_HVG <- density(Variable_6intron$length, from=153, to=26435)
d_LVG <- density(Stable_6intron$length, from=153, to=26435)

HVG_LVG_distribution_6intron[,2] <- d_HVG$y
HVG_LVG_distribution_6intron[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_6intron) <- c("length", "HVG", "LVG")
HVG_LVG_distribution_6intron_long <- gather(HVG_LVG_distribution_6intron, type, distr, HVG:LVG)

# graph

HVG_LVG_distribution_6intron <- mutate(HVG_LVG_distribution_6intron, mean_bootstrat=rowMeans(random_distribution_6intron_1000[,2:1001]),
                                       bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_6intron_1000[,2:1001])),
                                       top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_6intron_1000[,2:1001])))



gather(HVG_LVG_distribution_6intron, type, distr, HVG:top_CI, factor_key=TRUE) %>%
  ggplot(aes(x=length, y=distr, color=type, linetype=type)) +
  geom_line(size=1) +
  scale_x_continuous(limits = c(0, 8000)) +
  xlab("Gene length (bp)") +
  ylab("Density") +
  ggtitle("6 introns") +
  scale_color_manual(values=c("blue","green4","black","grey30","grey30")) +
  scale_linetype_manual(values=c("solid","solid","solid", "dashed", "dashed")) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line( size=.1, color="black" ), 
        panel.background = element_rect(fill="white",colour="grey50"), 
        text=element_text(size=38))


# wilcoxon test

wilcox.test(Variable_6intron$length,Stable_6intron$length)$p.value
[1] 3.276595e-07


Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(Variable_6intron$length,random_6intron_1000[,i])$p.value
}

# All values between 0.73398898 and 2.694564e-05 #





# 7 introns ----

subset_7intron <- filter(Background_gene_size_nbIntron,
                         Number_of_introns==7)
dim(subset_7intron)

# extract 1000 sets of random genes


random_distribution_7intron_1000=matrix(nrow=512,ncol=1000)
random_distribution_7intron_1000=data.frame(random_distribution_7intron_1000)

random_7intron_1000=matrix(nrow=35,ncol=1000)
random_7intron_1000=data.frame(random_7intron_1000)


for (i in 1:1000) {
  random_set <- subset_7intron[sample(nrow(subset_7intron), 35),]
  random_7intron_1000[,i] <- random_set[,2]
  d <- density(random_set[,2], from=153, to=26435)
  random_distribution_7intron_1000[,i] <- d$y
}

colnames(random_7intron_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(subset_7intron$length, from=153, to=26435)
random_distribution_7intron_1000 <- cbind(d_background$x, random_distribution_7intron_1000)
colnames(random_distribution_7intron_1000) <- c("length",paste("random_set",1:100, sep="_"))


HVG <- read.table("Gene-length_intron-number/HVG_allZT_background.txt", header=TRUE, sep='\t')
dim(HVG)
[1] 1358    2

Variable_7intron <- inner_join(HVG, subset_7intron, by=c("allZT"="Gene")) %>%
  select(allZT, length) 
dim(Variable_7intron)
[1] 35  2

LVG <- read.table("Gene-length_intron-number/LVG1000_allZT_background.txt", header=TRUE, sep='\t')
dim(LVG)
[1] 5727    2

Stable_7intron <- inner_join(LVG, subset_7intron, by=c("allZT"="Gene")) %>%
  select(allZT, length) 
dim(Stable_7intron)
[1] 313  2


HVG_LVG_distribution_7intron=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_7intron=data.frame(HVG_LVG_distribution_7intron)

d_background <-  density(subset_7intron$length, from=153, to=26435)
HVG_LVG_distribution_7intron[,1] <- d_background$x

d_HVG <- density(Variable_7intron$length, from=153, to=26435)
d_LVG <- density(Stable_7intron$length, from=153, to=26435)

HVG_LVG_distribution_7intron[,2] <- d_HVG$y
HVG_LVG_distribution_7intron[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_7intron) <- c("length", "HVG", "LVG")
HVG_LVG_distribution_7intron_long <- gather(HVG_LVG_distribution_7intron, type, distr, HVG:LVG)

# graph

HVG_LVG_distribution_7intron <- mutate(HVG_LVG_distribution_7intron, mean_bootstrat=rowMeans(random_distribution_7intron_1000[,2:1001]),
                                       bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_7intron_1000[,2:1001])),
                                       top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_7intron_1000[,2:1001])))



gather(HVG_LVG_distribution_7intron, type, distr, HVG:top_CI, factor_key=TRUE) %>%
  ggplot(aes(x=length, y=distr, color=type, linetype=type)) +
  geom_line(size=1) +
  scale_x_continuous(limits = c(0, 8000)) +
  xlab("Gene length (bp)") +
  ylab("Density") +
  ggtitle("7 introns") +
  scale_color_manual(values=c("blue","green4","black","grey30","grey30")) +
  scale_linetype_manual(values=c("solid","solid","solid", "dashed", "dashed")) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line( size=.1, color="black" ), 
        panel.background = element_rect(fill="white",colour="grey50"), 
        text=element_text(size=38))


# wilcoxon test

wilcox.test(Variable_7intron$length,Stable_7intron$length)$p.value
[1] 0.008406366


Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(Variable_7intron$length,random_7intron_1000[,i])$p.value
}

# All values between 1.0000000 and 0.003857455 #




# 8 to 10 introns ----

subset_8_10intron <- filter(Background_gene_size_nbIntron,
                            Number_of_introns>=8 & Number_of_introns<=10)
dim(subset_8_10intron)

# extract 1000 sets of random genes


random_distribution_8_10intron_1000=matrix(nrow=512,ncol=1000)
random_distribution_8_10intron_1000=data.frame(random_distribution_8_10intron_1000)

random_8_10intron_1000=matrix(nrow=80,ncol=1000)
random_8_10intron_1000=data.frame(random_8_10intron_1000)


for (i in 1:1000) {
  random_set <- subset_8_10intron[sample(nrow(subset_8_10intron), 80),]
  random_8_10intron_1000[,i] <- random_set[,2]
  d <- density(random_set[,2], from=153, to=26435)
  random_distribution_8_10intron_1000[,i] <- d$y
}

colnames(random_8_10intron_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(subset_8_10intron$length, from=153, to=26435)
random_distribution_8_10intron_1000 <- cbind(d_background$x, random_distribution_8_10intron_1000)
colnames(random_distribution_8_10intron_1000) <- c("length",paste("random_set",1:100, sep="_"))


HVG <- read.table("Gene-length_intron-number/HVG_allZT_background.txt", header=TRUE, sep='\t')
dim(HVG)
[1] 1358    2

Variable_8_10intron <- inner_join(HVG, subset_8_10intron, by=c("allZT"="Gene")) %>%
  select(allZT, length) 
dim(Variable_8_10intron)
[1] 80  2

LVG <- read.table("Gene-length_intron-number/LVG1000_allZT_background.txt", header=TRUE, sep='\t')
dim(LVG)
[1] 5727    2

Stable_8_10intron <- inner_join(LVG, subset_8_10intron, by=c("allZT"="Gene")) %>%
  select(allZT, length) 
dim(Stable_8_10intron)
[1] 814  2


HVG_LVG_distribution_8_10intron=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_8_10intron=data.frame(HVG_LVG_distribution_8_10intron)

d_background <-  density(subset_8_10intron$length, from=153, to=26435)
HVG_LVG_distribution_8_10intron[,1] <- d_background$x

d_HVG <- density(Variable_8_10intron$length, from=153, to=26435)
d_LVG <- density(Stable_8_10intron$length, from=153, to=26435)

HVG_LVG_distribution_8_10intron[,2] <- d_HVG$y
HVG_LVG_distribution_8_10intron[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_8_10intron) <- c("length", "HVG", "LVG")
HVG_LVG_distribution_8_10intron_long <- gather(HVG_LVG_distribution_8_10intron, type, distr, HVG:LVG)

# graph

HVG_LVG_distribution_8_10intron <- mutate(HVG_LVG_distribution_8_10intron, mean_bootstrat=rowMeans(random_distribution_8_10intron_1000[,2:1001]),
                                          bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_8_10intron_1000[,2:1001])),
                                          top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_8_10intron_1000[,2:1001])))



gather(HVG_LVG_distribution_8_10intron, type, distr, HVG:top_CI, factor_key=TRUE) %>%
  ggplot(aes(x=length, y=distr, color=type, linetype=type)) +
  geom_line(size=1) +
  scale_x_continuous(limits = c(0, 8000)) +
  xlab("Gene length (bp)") +
  ylab("Density") +
  ggtitle("8 to 10 introns") +
  scale_color_manual(values=c("blue","green4","black","grey30","grey30")) +
  scale_linetype_manual(values=c("solid","solid","solid", "dashed", "dashed")) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line( size=.1, color="black" ), 
        panel.background = element_rect(fill="white",colour="grey50"), 
        text=element_text(size=38))


# wilcoxon test

wilcox.test(Variable_8_10intron$length,Stable_8_10intron$length)$p.value
[1] 6.855292e-05


Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(Variable_8_10intron$length,random_8_10intron_1000[,i])$p.value
}

# All values between 0.9931929 and 0.0002312065 #


# 11 to 15 introns ----

subset_11_15intron <- filter(Background_gene_size_nbIntron,
                             Number_of_introns>=11 & Number_of_introns<=15)
dim(subset_11_15intron)

# extract 1000 sets of random genes


random_distribution_11_15intron_1000=matrix(nrow=512,ncol=1000)
random_distribution_11_15intron_1000=data.frame(random_distribution_11_15intron_1000)

random_11_15intron_1000=matrix(nrow=79,ncol=1000)
random_11_15intron_1000=data.frame(random_11_15intron_1000)


for (i in 1:1000) {
  random_set <- subset_11_15intron[sample(nrow(subset_11_15intron), 79),]
  random_11_15intron_1000[,i] <- random_set[,2]
  d <- density(random_set[,2], from=153, to=26435)
  random_distribution_11_15intron_1000[,i] <- d$y
}

colnames(random_11_15intron_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(subset_11_15intron$length, from=153, to=26435)
random_distribution_11_15intron_1000 <- cbind(d_background$x, random_distribution_11_15intron_1000)
colnames(random_distribution_11_15intron_1000) <- c("length",paste("random_set",1:100, sep="_"))


HVG <- read.table("Gene-length_intron-number/HVG_allZT_background.txt", header=TRUE, sep='\t')
dim(HVG)
[1] 1358    2

Variable_11_15intron <- inner_join(HVG, subset_11_15intron, by=c("allZT"="Gene")) %>%
  select(allZT, length) 
dim(Variable_11_15intron)
[1] 79  2

LVG <- read.table("Gene-length_intron-number/LVG1000_allZT_background.txt", header=TRUE, sep='\t')
dim(LVG)
[1] 5727    2

Stable_11_15intron <- inner_join(LVG, subset_11_15intron, by=c("allZT"="Gene")) %>%
  select(allZT, length) 
dim(Stable_11_15intron)
[1] 849  2


HVG_LVG_distribution_11_15intron=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_11_15intron=data.frame(HVG_LVG_distribution_11_15intron)

d_background <-  density(subset_11_15intron$length, from=153, to=26435)
HVG_LVG_distribution_11_15intron[,1] <- d_background$x

d_HVG <- density(Variable_11_15intron$length, from=153, to=26435)
d_LVG <- density(Stable_11_15intron$length, from=153, to=26435)

HVG_LVG_distribution_11_15intron[,2] <- d_HVG$y
HVG_LVG_distribution_11_15intron[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_11_15intron) <- c("length", "HVG", "LVG")
HVG_LVG_distribution_11_15intron_long <- gather(HVG_LVG_distribution_11_15intron, type, distr, HVG:LVG)

# graph

HVG_LVG_distribution_11_15intron <- mutate(HVG_LVG_distribution_11_15intron, mean_bootstrat=rowMeans(random_distribution_11_15intron_1000[,2:1001]),
                                           bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_11_15intron_1000[,2:1001])),
                                           top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_11_15intron_1000[,2:1001])))



gather(HVG_LVG_distribution_11_15intron, type, distr, HVG:top_CI, factor_key=TRUE) %>%
  ggplot(aes(x=length, y=distr, color=type, linetype=type)) +
  geom_line(size=1) +
  scale_x_continuous(limits = c(0, 8000)) +
  xlab("Gene length (bp)") +
  ylab("Density") +
  ggtitle("11 to 15 introns") +
  scale_color_manual(values=c("blue","green4","black","grey30","grey30")) +
  scale_linetype_manual(values=c("solid","solid","solid", "dashed", "dashed")) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line( size=.1, color="black" ), 
        panel.background = element_rect(fill="white",colour="grey50"), 
        text=element_text(size=38))

# wilcoxon test

wilcox.test(Variable_11_15intron$length,Stable_11_15intron$length)$p.value
[1] 9.457579e-08


Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(Variable_11_15intron$length,random_11_15intron_1000[,i])$p.value
}

# All values between 0.54165998 and 2.751205e-07 #



# 16 to 20 introns ----

subset_16_20intron <- filter(Background_gene_size_nbIntron,
                             Number_of_introns>=16 & Number_of_introns<=20)
dim(subset_16_20intron)

# extract 1000 sets of random genes


random_distribution_16_20intron_1000=matrix(nrow=512,ncol=1000)
random_distribution_16_20intron_1000=data.frame(random_distribution_16_20intron_1000)

random_16_20intron_1000=matrix(nrow=25,ncol=1000)
random_16_20intron_1000=data.frame(random_16_20intron_1000)


for (i in 1:1000) {
  random_set <- subset_11_15intron[sample(nrow(subset_16_20intron), 25),]
  random_16_20intron_1000[,i] <- random_set[,2]
  d <- density(random_set[,2], from=153, to=26435)
  random_distribution_16_20intron_1000[,i] <- d$y
}

colnames(random_16_20intron_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(subset_16_20intron$length, from=153, to=26435)
random_distribution_16_20intron_1000 <- cbind(d_background$x, random_distribution_16_20intron_1000)
colnames(random_distribution_16_20intron_1000) <- c("length",paste("random_set",1:100, sep="_"))


HVG <- read.table("Gene-length_intron-number/HVG_allZT_background.txt", header=TRUE, sep='\t')
dim(HVG)
[1] 1358    2

Variable_16_20intron <- inner_join(HVG, subset_16_20intron, by=c("allZT"="Gene")) %>%
  select(allZT, length) 
dim(Variable_16_20intron)
[1] 25  2

LVG <- read.table("Gene-length_intron-number/LVG1000_allZT_background.txt", header=TRUE, sep='\t')
dim(LVG)
[1] 5727    2

Stable_16_20intron <- inner_join(LVG, subset_16_20intron, by=c("allZT"="Gene")) %>%
  select(allZT, length) 
dim(Stable_16_20intron)
[1] 501  2


HVG_LVG_distribution_16_20intron=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_16_20intron=data.frame(HVG_LVG_distribution_16_20intron)

d_background <-  density(subset_16_20intron$length, from=153, to=26435)
HVG_LVG_distribution_16_20intron[,1] <- d_background$x

d_HVG <- density(Variable_16_20intron$length, from=153, to=26435)
d_LVG <- density(Stable_16_20intron$length, from=153, to=26435)

HVG_LVG_distribution_16_20intron[,2] <- d_HVG$y
HVG_LVG_distribution_16_20intron[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_16_20intron) <- c("length", "HVG", "LVG")
HVG_LVG_distribution_16_20intron_long <- gather(HVG_LVG_distribution_16_20intron, type, distr, HVG:LVG)

# graph

HVG_LVG_distribution_16_20intron <- mutate(HVG_LVG_distribution_16_20intron, mean_bootstrat=rowMeans(random_distribution_16_20intron_1000[,2:1001]),
                                           bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_16_20intron_1000[,2:1001])),
                                           top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_16_20intron_1000[,2:1001])))



gather(HVG_LVG_distribution_16_20intron, type, distr, HVG:top_CI, factor_key=TRUE) %>%
  ggplot(aes(x=length, y=distr, color=type, linetype=type)) +
  geom_line(size=1) +
  scale_x_continuous(limits = c(0, 8000)) +
  xlab("Gene length (bp)") +
  ylab("Density") +
  ggtitle("16 to 20 introns") +
  scale_color_manual(values=c("blue","green4","black","grey30","grey30")) +
  scale_linetype_manual(values=c("solid","solid","solid", "dashed", "dashed")) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line( size=.1, color="black" ), 
        panel.background = element_rect(fill="white",colour="grey50"), 
        text=element_text(size=38))


# wilcoxon test

wilcox.test(Variable_16_20intron$length,Stable_16_20intron$length)$p.value
[1] 3.451906e-05


Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(Variable_16_20intron$length,random_16_20intron_1000[,i])$p.value
}

# All values between 1.0000000 and 0.003477174 #



# 21 introns and more ----

subset_21intron <- filter(Background_gene_size_nbIntron,
                          Number_of_introns>=21)
dim(subset_21intron)

# extract 1000 sets of random genes


random_distribution_21intron_1000=matrix(nrow=512,ncol=1000)
random_distribution_21intron_1000=data.frame(random_distribution_21intron_1000)

random_21intron_1000=matrix(nrow=66,ncol=1000)
random_21intron_1000=data.frame(random_21intron_1000)


for (i in 1:1000) {
  random_set <- subset_11_15intron[sample(nrow(subset_21intron), 66),]
  random_21intron_1000[,i] <- random_set[,2]
  d <- density(random_set[,2], from=153, to=26435)
  random_distribution_21intron_1000[,i] <- d$y
}

colnames(random_21intron_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(subset_21intron$length, from=153, to=26435)
random_distribution_21intron_1000 <- cbind(d_background$x, random_distribution_21intron_1000)
colnames(random_distribution_21intron_1000) <- c("length",paste("random_set",1:100, sep="_"))


HVG <- read.table("Gene-length_intron-number/HVG_allZT_background.txt", header=TRUE, sep='\t')
dim(HVG)
[1] 1358    2

Variable_21intron <- inner_join(HVG, subset_21intron, by=c("allZT"="Gene")) %>%
  select(allZT, length) 
dim(Variable_21intron)
[1] 66  2

LVG <- read.table("Gene-length_intron-number/LVG1000_allZT_background.txt", header=TRUE, sep='\t')
dim(LVG)
[1] 5727    2

Stable_21intron <- inner_join(LVG, subset_21intron, by=c("allZT"="Gene")) %>%
  select(allZT, length) 
dim(Stable_21intron)
[1] 763  2


HVG_LVG_distribution_21intron=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_21intron=data.frame(HVG_LVG_distribution_21intron)

d_background <-  density(subset_21intron$length, from=153, to=26435)
HVG_LVG_distribution_21intron[,1] <- d_background$x

d_HVG <- density(Variable_21intron$length, from=153, to=26435)
d_LVG <- density(Stable_21intron$length, from=153, to=26435)

HVG_LVG_distribution_21intron[,2] <- d_HVG$y
HVG_LVG_distribution_21intron[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_21intron) <- c("length", "HVG", "LVG")
HVG_LVG_distribution_21intron_long <- gather(HVG_LVG_distribution_21intron, type, distr, HVG:LVG)

# graph

HVG_LVG_distribution_21intron <- mutate(HVG_LVG_distribution_21intron, mean_bootstrat=rowMeans(random_distribution_21intron_1000[,2:1001]),
                                        bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_21intron_1000[,2:1001])),
                                        top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_21intron_1000[,2:1001])))



gather(HVG_LVG_distribution_21intron, type, distr, HVG:top_CI, factor_key=TRUE) %>%
  ggplot(aes(x=length, y=distr, color=type, linetype=type)) +
  geom_line(size=1) +
  scale_x_continuous(limits = c(0, 8000)) +
  xlab("Gene length (bp)") +
  ylab("Density") +
  ggtitle("21 introns and more") +
  scale_color_manual(values=c("blue","green4","black","grey30","grey30")) +
  scale_linetype_manual(values=c("solid","solid","solid", "dashed", "dashed")) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line( size=.1, color="black" ), 
        panel.background = element_rect(fill="white",colour="grey50"), 
        text=element_text(size=38))


# wilcoxon test

wilcox.test(Variable_21intron$length,Stable_21intron$length)$p.value
[1] 0.005837275


Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)


for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(Variable_21intron$length,random_21intron_1000[,i])$p.value
}

# All values between 0.62789073 and 1.075565e-07 #





### Plotting of the number of introns for HVGs, LVGs and random genes with a fixed gene lenght ----

# Data preparation

Gene_size <- read.table("Gene-length_intron-number/TAIR10_Gene_sizes.txt", header=TRUE, sep='\t')
dim(Gene_size)

Background <- read.table("Gene-length_intron-number/All_genes_background.txt", header=TRUE, sep='\t')
dim(Background)

Background_gene_size <- inner_join(Background, Gene_size, by=c("Gene"="gene_name")) %>%
  select(Gene, length)


Gene_nbIntrons <- read.table("Gene-length_intron-number/TAIR10_Gene_nbIntrons.txt", header=TRUE, sep='\t')
dim(Gene_nbIntrons)


Background_gene_size_nbIntron <- inner_join(Background_gene_size, Gene_nbIntrons,
                                            by=c("Gene"="gene_name"))

# Genes less than 1000bp ----

subset_0to1000bp <- filter(Background_gene_size_nbIntron,
                           length<1000)
dim(subset_0to1000bp)


# extract 1000 sets of random genes 


random_distribution_0to1000bp_1000=matrix(nrow=512,ncol=1000)
random_distribution_0to1000bp_1000=data.frame(random_distribution_0to1000bp_1000)

random_0to1000bp_1000=matrix(nrow=425,ncol=1000)
random_0to1000bp_1000=data.frame(random_0to1000bp_1000)


for (i in 1:1000) {
  random_set <- subset_0to1000bp[sample(nrow(subset_0to1000bp), 425),]
  random_0to1000bp_1000[,i] <- random_set[,3]
  d <- density(random_set[,3], from=0, to=196)
  random_distribution_0to1000bp_1000[,i] <- d$y
}

colnames(random_0to1000bp_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(subset_0to1000bp$Number_of_introns, from=0, to=196)
random_distribution_0to1000bp_1000 <- cbind(d_background$x, random_distribution_0to1000bp_1000)
colnames(random_distribution_0to1000bp_1000) <- c("Number_of_introns",paste("random_set",1:100, sep="_"))


HVG <- read.table("Gene-length_intron-number/HVG_allZT_background.txt", header=TRUE, sep='\t')
dim(HVG)
[1] 1358    2

Variable_0to1000bp <- inner_join(HVG, subset_0to1000bp, by=c("allZT"="Gene")) %>%
  select(allZT, Number_of_introns) 
dim(Variable_0to1000bp)
[1] 425   2

LVG <- read.table("Gene-length_intron-number/LVG1000_allZT_background.txt", header=TRUE, sep='\t')
dim(LVG)
[1] 5727    2

Stable_0to1000bp <- inner_join(LVG, subset_0to1000bp, by=c("allZT"="Gene")) %>%
  select(allZT, Number_of_introns) 
dim(Stable_0to1000bp)
[1] 73  2

HVG_LVG_distribution_0to1000bp=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_0to1000bp=data.frame(HVG_LVG_distribution_0to1000bp)

d_background <-  density(subset_0to1000bp$Number_of_introns, from=0, to=196)
HVG_LVG_distribution_0to1000bp[,1] <- d_background$x

d_HVG <- density(Variable_0to1000bp$Number_of_introns, from=0, to=196)
d_LVG <- density(Stable_0to1000bp$Number_of_introns, from=0, to=196)

HVG_LVG_distribution_0to1000bp[,2] <- d_HVG$y
HVG_LVG_distribution_0to1000bp[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_0to1000bp) <- c("Number_of_introns", "HVG", "LVG")
HVG_LVG_distribution_0to1000bp_long <- gather(HVG_LVG_distribution_0to1000bp, type, distr, HVG:LVG)

# graph

HVG_LVG_distribution_0to1000bp <- mutate(HVG_LVG_distribution_0to1000bp, mean_bootstrat=rowMeans(random_distribution_0to1000bp_1000[,2:1001]),
                                         bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_0to1000bp_1000[,2:1001])),
                                         top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_0to1000bp_1000[,2:1001])))



gather(HVG_LVG_distribution_0to1000bp, type, distr, HVG:top_CI, factor_key=TRUE) %>%
  ggplot(aes(x=Number_of_introns, y=distr, color=type, linetype=type)) +
  geom_line(size=1) +
  scale_x_continuous(limits = c(0, 15)) +
  xlab("Number of introns") +
  ylab("Density") +
  ggtitle("Genes less than 1000bp") +
  scale_color_manual(values=c("blue","green4","black","grey30","grey30")) +
  scale_linetype_manual(values=c("solid","solid","solid", "dashed", "dashed")) +
  theme(panel.grid.major.x = element_blank(),panel.grid.major.y = element_line( size=.1, color="black" ) 
        , panel.background = element_rect(fill="white",colour="grey50"), text=element_text(size=38))


# wilcoxon test

wilcox.test(Variable_0to1000bp$Number_of_introns,Stable_0to1000bp$Number_of_introns)$p.value
[1] 0.006280032


Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)

for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(Variable_0to1000bp$Number_of_introns,random_0to1000bp_1000[,i])$p.value
}


#All values between 0.9985146 and 0.007649847 #




# Genes 1000bp to 2000bp ----

subset_1000to2000bp <- filter(Background_gene_size_nbIntron,
                              length>=1000 & length<2000)
dim(subset_1000to2000bp)


# extract 1000 sets of random genes 


random_distribution_1000to2000bp_1000=matrix(nrow=512,ncol=1000)
random_distribution_1000to2000bp_1000=data.frame(random_distribution_1000to2000bp_1000)

random_1000to2000bp_1000=matrix(nrow=527,ncol=1000)
random_1000to2000bp_1000=data.frame(random_1000to2000bp_1000)


for (i in 1:1000) {
  random_set <- subset_1000to2000bp[sample(nrow(subset_1000to2000bp), 527),]
  random_1000to2000bp_1000[,i] <- random_set[,3]
  d <- density(random_set[,3], from=0, to=196)
  random_distribution_1000to2000bp_1000[,i] <- d$y
}

colnames(random_1000to2000bp_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(subset_1000to2000bp$Number_of_introns, from=0, to=196)
random_distribution_1000to2000bp_1000 <- cbind(d_background$x, random_distribution_1000to2000bp_1000)
colnames(random_distribution_1000to2000bp_1000) <- c("Number_of_introns",paste("random_set",1:100, sep="_"))


HVG <- read.table("Gene-length_intron-number/HVG_allZT_background.txt", header=TRUE, sep='\t')
dim(HVG)
[1] 1358    2

Variable_1000to2000bp <- inner_join(HVG, subset_1000to2000bp, by=c("allZT"="Gene")) %>%
  select(allZT, Number_of_introns) 
dim(Variable_1000to2000bp)
[1] 527   2

LVG <- read.table("Gene-length_intron-number/LVG1000_allZT_background.txt", header=TRUE, sep='\t')
dim(LVG)
[1] 5727    2

Stable_1000to2000bp <- inner_join(LVG, subset_1000to2000bp, by=c("allZT"="Gene")) %>%
  select(allZT, Number_of_introns) 
dim(Stable_1000to2000bp)
[1] 1287  2

HVG_LVG_distribution_1000to2000bp=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_1000to2000bp=data.frame(HVG_LVG_distribution_1000to2000bp)

d_background <-  density(subset_1000to2000bp$Number_of_introns, from=0, to=196)
HVG_LVG_distribution_1000to2000bp[,1] <- d_background$x

d_HVG <- density(Variable_1000to2000bp$Number_of_introns, from=0, to=196)
d_LVG <- density(Stable_1000to2000bp$Number_of_introns, from=0, to=196)

HVG_LVG_distribution_1000to2000bp[,2] <- d_HVG$y
HVG_LVG_distribution_1000to2000bp[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_1000to2000bp) <- c("Number_of_introns", "HVG", "LVG")
HVG_LVG_distribution_1000to2000bp_long <- gather(HVG_LVG_distribution_1000to2000bp, type, distr, HVG:LVG)

# graph

HVG_LVG_distribution_1000to2000bp <- mutate(HVG_LVG_distribution_1000to2000bp, mean_bootstrat=rowMeans(random_distribution_1000to2000bp_1000[,2:1001]),
                                            bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_1000to2000bp_1000[,2:1001])),
                                            top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_1000to2000bp_1000[,2:1001])))



gather(HVG_LVG_distribution_1000to2000bp, type, distr, HVG:top_CI, factor_key=TRUE) %>%
  ggplot(aes(x=Number_of_introns, y=distr, color=type, linetype=type)) +
  geom_line(size=1) +
  scale_x_continuous(limits = c(0, 40)) +
  xlab("Number of introns") +
  ylab("Density") +
  ggtitle("Genes between 1000 and 2000bp") +
  scale_color_manual(values=c("blue","green4","black","grey30","grey30")) +
  scale_linetype_manual(values=c("solid","solid","solid", "dashed", "dashed")) +
  theme(panel.grid.major.x = element_blank(),panel.grid.major.y = element_line( size=.1, color="black" ) 
        , panel.background = element_rect(fill="white",colour="grey50"), text=element_text(size=38))


# wilcoxon test

wilcox.test(Variable_1000to2000bp$Number_of_introns,Stable_1000to2000bp$Number_of_introns)$p.value
[1] 0.3548196


Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)

for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(Variable_1000to2000bp$Number_of_introns,random_1000to2000bp_1000[,i])$p.value
}


#All values between 0.9985146 and 0.007649847 #



# Genes 2000bp to 3000bp ----

subset_2000to3000bp <- filter(Background_gene_size_nbIntron,
                              length>=2000 & length<3000)
dim(subset_2000to3000bp)


# extract 1000 sets of random genes 


random_distribution_2000to3000bp_1000=matrix(nrow=512,ncol=1000)
random_distribution_2000to3000bp_1000=data.frame(random_distribution_2000to3000bp_1000)

random_2000to3000bp_1000=matrix(nrow=250,ncol=1000)
random_2000to3000bp_1000=data.frame(random_2000to3000bp_1000)


for (i in 1:1000) {
  random_set <- subset_2000to3000bp[sample(nrow(subset_2000to3000bp), 250),]
  random_2000to3000bp_1000[,i] <- random_set[,3]
  d <- density(random_set[,3], from=0, to=196)
  random_distribution_2000to3000bp_1000[,i] <- d$y
}

colnames(random_2000to3000bp_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(subset_2000to3000bp$Number_of_introns, from=0, to=196)
random_distribution_2000to3000bp_1000 <- cbind(d_background$x, random_distribution_2000to3000bp_1000)
colnames(random_distribution_2000to3000bp_1000) <- c("Number_of_introns",paste("random_set",1:100, sep="_"))


HVG <- read.table("Gene-length_intron-number/HVG_allZT_background.txt", header=TRUE, sep='\t')
dim(HVG)
[1] 1358    2

Variable_2000to3000bp <- inner_join(HVG, subset_2000to3000bp, by=c("allZT"="Gene")) %>%
  select(allZT, Number_of_introns) 
dim(Variable_2000to3000bp)
[1] 250   2

LVG <- read.table("Gene-length_intron-number/LVG1000_allZT_background.txt", header=TRUE, sep='\t')
dim(LVG)
[1] 5727    2

Stable_2000to3000bp <- inner_join(LVG, subset_2000to3000bp, by=c("allZT"="Gene")) %>%
  select(allZT, Number_of_introns) 
dim(Stable_2000to3000bp)
[1] 2150  2

HVG_LVG_distribution_2000to3000bp=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_2000to3000bp=data.frame(HVG_LVG_distribution_2000to3000bp)

d_background <-  density(subset_2000to3000bp$Number_of_introns, from=0, to=196)
HVG_LVG_distribution_2000to3000bp[,1] <- d_background$x

d_HVG <- density(Variable_2000to3000bp$Number_of_introns, from=0, to=196)
d_LVG <- density(Stable_2000to3000bp$Number_of_introns, from=0, to=196)

HVG_LVG_distribution_2000to3000bp[,2] <- d_HVG$y
HVG_LVG_distribution_2000to3000bp[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_2000to3000bp) <- c("Number_of_introns", "HVG", "LVG")
HVG_LVG_distribution_2000to3000bp_long <- gather(HVG_LVG_distribution_2000to3000bp, type, distr, HVG:LVG)

# graph

HVG_LVG_distribution_2000to3000bp <- mutate(HVG_LVG_distribution_2000to3000bp, mean_bootstrat=rowMeans(random_distribution_2000to3000bp_1000[,2:1001]),
                                            bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_2000to3000bp_1000[,2:1001])),
                                            top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_2000to3000bp_1000[,2:1001])))



gather(HVG_LVG_distribution_2000to3000bp, type, distr, HVG:top_CI, factor_key=TRUE) %>%
  ggplot(aes(x=Number_of_introns, y=distr, color=type, linetype=type)) +
  geom_line(size=1) +
  scale_x_continuous(limits = c(0, 70)) +
  xlab("Number of introns") +
  ylab("Density") +
  ggtitle("Genes between 2000 and 3000bp") +
  scale_color_manual(values=c("blue","green4","black","grey30","grey30")) +
  scale_linetype_manual(values=c("solid","solid","solid", "dashed", "dashed")) +
  theme(panel.grid.major.x = element_blank(),panel.grid.major.y = element_line( size=.1, color="black" ) 
        , panel.background = element_rect(fill="white",colour="grey50"), text=element_text(size=38))

# wilcoxon test

wilcox.test(Variable_2000to3000bp$Number_of_introns,Stable_2000to3000bp$Number_of_introns)$p.value
[1] 0.0001901308


Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)

for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(Variable_2000to3000bp$Number_of_introns,random_2000to3000bp_1000[,i])$p.value
}


#All values between 0.323600817 and 2.437231e-07 #



# Genes 3000bp to 4000bp ----

subset_3000to4000bp <- filter(Background_gene_size_nbIntron,
                              length>=3000 & length<4000)
dim(subset_3000to4000bp)


# extract 1000 sets of random genes 


random_distribution_3000to4000bp_1000=matrix(nrow=512,ncol=1000)
random_distribution_3000to4000bp_1000=data.frame(random_distribution_3000to4000bp_1000)

random_3000to4000bp_1000=matrix(nrow=103,ncol=1000)
random_3000to4000bp_1000=data.frame(random_3000to4000bp_1000)


for (i in 1:1000) {
  random_set <- subset_3000to4000bp[sample(nrow(subset_3000to4000bp), 103),]
  random_3000to4000bp_1000[,i] <- random_set[,3]
  d <- density(random_set[,3], from=0, to=196)
  random_distribution_3000to4000bp_1000[,i] <- d$y
}

colnames(random_3000to4000bp_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(subset_3000to4000bp$Number_of_introns, from=0, to=196)
random_distribution_3000to4000bp_1000 <- cbind(d_background$x, random_distribution_3000to4000bp_1000)
colnames(random_distribution_3000to4000bp_1000) <- c("Number_of_introns",paste("random_set",1:100, sep="_"))


HVG <- read.table("Gene-length_intron-number/HVG_allZT_background.txt", header=TRUE, sep='\t')
dim(HVG)
[1] 1358    2

Variable_3000to4000bp <- inner_join(HVG, subset_3000to4000bp, by=c("allZT"="Gene")) %>%
  select(allZT, Number_of_introns) 
dim(Variable_3000to4000bp)
[1] 103   2

LVG <- read.table("Gene-length_intron-number/LVG1000_allZT_background.txt", header=TRUE, sep='\t')
dim(LVG)
[1] 5727    2

Stable_3000to4000bp <- inner_join(LVG, subset_3000to4000bp, by=c("allZT"="Gene")) %>%
  select(allZT, Number_of_introns) 
dim(Stable_3000to4000bp)
[1] 1139  2

HVG_LVG_distribution_3000to4000bp=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_3000to4000bp=data.frame(HVG_LVG_distribution_3000to4000bp)

d_background <-  density(subset_3000to4000bp$Number_of_introns, from=0, to=196)
HVG_LVG_distribution_3000to4000bp[,1] <- d_background$x

d_HVG <- density(Variable_3000to4000bp$Number_of_introns, from=0, to=196)
d_LVG <- density(Stable_3000to4000bp$Number_of_introns, from=0, to=196)

HVG_LVG_distribution_3000to4000bp[,2] <- d_HVG$y
HVG_LVG_distribution_3000to4000bp[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_3000to4000bp) <- c("Number_of_introns", "HVG", "LVG")
HVG_LVG_distribution_3000to4000bp_long <- gather(HVG_LVG_distribution_3000to4000bp, type, distr, HVG:LVG)

# graph

HVG_LVG_distribution_3000to4000bp <- mutate(HVG_LVG_distribution_3000to4000bp, mean_bootstrat=rowMeans(random_distribution_3000to4000bp_1000[,2:1001]),
                                            bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_3000to4000bp_1000[,2:1001])),
                                            top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_3000to4000bp_1000[,2:1001])))



gather(HVG_LVG_distribution_3000to4000bp, type, distr, HVG:top_CI, factor_key=TRUE) %>%
  ggplot(aes(x=Number_of_introns, y=distr, color=type, linetype=type)) +
  geom_line(size=1) +
  scale_x_continuous(limits = c(0, 80)) +
  xlab("Number of introns") +
  ylab("Density") +
  ggtitle("Genes between 3000 and 4000bp") +
  scale_color_manual(values=c("blue","green4","black","grey30","grey30")) +
  scale_linetype_manual(values=c("solid","solid","solid", "dashed", "dashed")) +
  theme(panel.grid.major.x = element_blank(),panel.grid.major.y = element_line( size=.1, color="black" ) 
        , panel.background = element_rect(fill="white",colour="grey50"), text=element_text(size=38))


# wilcoxon test

wilcox.test(Variable_3000to4000bp$Number_of_introns,Stable_3000to4000bp$Number_of_introns)$p.value
[1] 0.5843233


Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)

for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(Variable_3000to4000bp$Number_of_introns,random_3000to4000bp_1000[,i])$p.value
}


#All values between 1.0000000 and 0.03158386 #


# Genes 4000bp to 5000bp ----

subset_4000to5000bp <- filter(Background_gene_size_nbIntron,
                              length>=4000 & length<5000)
dim(subset_4000to5000bp)


# extract 1000 sets of random genes 


random_distribution_4000to5000bp_1000=matrix(nrow=512,ncol=1000)
random_distribution_4000to5000bp_1000=data.frame(random_distribution_4000to5000bp_1000)

random_4000to5000bp_1000=matrix(nrow=34,ncol=1000)
random_4000to5000bp_1000=data.frame(random_4000to5000bp_1000)


for (i in 1:1000) {
  random_set <- subset_4000to5000bp[sample(nrow(subset_4000to5000bp), 34),]
  random_4000to5000bp_1000[,i] <- random_set[,3]
  d <- density(random_set[,3], from=0, to=196)
  random_distribution_4000to5000bp_1000[,i] <- d$y
}

colnames(random_4000to5000bp_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(subset_4000to5000bp$Number_of_introns, from=0, to=196)
random_distribution_4000to5000bp_1000 <- cbind(d_background$x, random_distribution_4000to5000bp_1000)
colnames(random_distribution_4000to5000bp_1000) <- c("Number_of_introns",paste("random_set",1:100, sep="_"))


HVG <- read.table("Gene-length_intron-number/HVG_allZT_background.txt", header=TRUE, sep='\t')
dim(HVG)
[1] 1358    2

Variable_4000to5000bp <- inner_join(HVG, subset_4000to5000bp, by=c("allZT"="Gene")) %>%
  select(allZT, Number_of_introns) 
dim(Variable_4000to5000bp)
[1] 34   2

LVG <- read.table("Gene-length_intron-number/LVG1000_allZT_background.txt", header=TRUE, sep='\t')
dim(LVG)
[1] 5727    2

Stable_4000to5000bp <- inner_join(LVG, subset_4000to5000bp, by=c("allZT"="Gene")) %>%
  select(allZT, Number_of_introns) 
dim(Stable_4000to5000bp)
[1] 532  2

HVG_LVG_distribution_4000to5000bp=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_4000to5000bp=data.frame(HVG_LVG_distribution_4000to5000bp)

d_background <-  density(subset_4000to5000bp$Number_of_introns, from=0, to=196)
HVG_LVG_distribution_4000to5000bp[,1] <- d_background$x

d_HVG <- density(Variable_4000to5000bp$Number_of_introns, from=0, to=196)
d_LVG <- density(Stable_4000to5000bp$Number_of_introns, from=0, to=196)

HVG_LVG_distribution_4000to5000bp[,2] <- d_HVG$y
HVG_LVG_distribution_4000to5000bp[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_4000to5000bp) <- c("Number_of_introns", "HVG", "LVG")
HVG_LVG_distribution_4000to5000bp_long <- gather(HVG_LVG_distribution_4000to5000bp, type, distr, HVG:LVG)

# graph

HVG_LVG_distribution_4000to5000bp <- mutate(HVG_LVG_distribution_4000to5000bp, mean_bootstrat=rowMeans(random_distribution_4000to5000bp_1000[,2:1001]),
                                            bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_4000to5000bp_1000[,2:1001])),
                                            top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_4000to5000bp_1000[,2:1001])))



gather(HVG_LVG_distribution_4000to5000bp, type, distr, HVG:top_CI, factor_key=TRUE) %>%
  ggplot(aes(x=Number_of_introns, y=distr, color=type, linetype=type)) +
  geom_line(size=1) +
  scale_x_continuous(limits = c(0, 80)) +
  xlab("Number of introns") +
  ylab("Density") +
  ggtitle("Genes between 4000 and 5000bp") +
  scale_color_manual(values=c("blue","green4","black","grey30","grey30")) +
  scale_linetype_manual(values=c("solid","solid","solid", "dashed", "dashed")) +
  theme(panel.grid.major.x = element_blank(),panel.grid.major.y = element_line( size=.1, color="black" ) 
        , panel.background = element_rect(fill="white",colour="grey50"), text=element_text(size=38))

# wilcoxon test

wilcox.test(Variable_4000to5000bp$Number_of_introns,Stable_4000to5000bp$Number_of_introns)$p.value
[1] 0.9533834


Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)

for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(Variable_4000to5000bp$Number_of_introns,random_4000to5000bp_1000[,i])$p.value
}


#All values between 1.0000000 and 0.06028156 #



# Genes more than 5000bp ----

subset_5000bp <- filter(Background_gene_size_nbIntron,
                        length>5000)
dim(subset_5000bp)


# extract 1000 sets of random genes 


random_distribution_5000bp_1000=matrix(nrow=512,ncol=1000)
random_distribution_5000bp_1000=data.frame(random_distribution_5000bp_1000)

random_5000bp_1000=matrix(nrow=19,ncol=1000)
random_5000bp_1000=data.frame(random_5000bp_1000)


for (i in 1:1000) {
  random_set <- subset_5000bp[sample(nrow(subset_5000bp), 19),]
  random_5000bp_1000[,i] <- random_set[,3]
  d <- density(random_set[,3], from=0, to=196)
  random_distribution_5000bp_1000[,i] <- d$y
}

colnames(random_5000bp_1000) <- paste("random_set",1:1000, sep="_")

d_background <-  density(subset_5000bp$Number_of_introns, from=0, to=196)
random_distribution_5000bp_1000 <- cbind(d_background$x, random_distribution_5000bp_1000)
colnames(random_distribution_5000bp_1000) <- c("Number_of_introns",paste("random_set",1:100, sep="_"))


HVG <- read.table("Gene-length_intron-number/HVG_allZT_background.txt", header=TRUE, sep='\t')
dim(HVG)
[1] 1358    2

Variable_5000bp <- inner_join(HVG, subset_5000bp, by=c("allZT"="Gene")) %>%
  select(allZT, Number_of_introns) 
dim(Variable_5000bp)
[1] 19   2

LVG <- read.table("Gene-length_intron-number/LVG1000_allZT_background.txt", header=TRUE, sep='\t')
dim(LVG)
[1] 5727    2

Stable_5000bp <- inner_join(LVG, subset_5000bp, by=c("allZT"="Gene")) %>%
  select(allZT, Number_of_introns) 
dim(Stable_5000bp)
[1] 546  2

HVG_LVG_distribution_5000bp=matrix(nrow=512,ncol=3)
HVG_LVG_distribution_5000bp=data.frame(HVG_LVG_distribution_5000bp)

d_background <-  density(subset_5000bp$Number_of_introns, from=0, to=196)
HVG_LVG_distribution_5000bp[,1] <- d_background$x

d_HVG <- density(Variable_5000bp$Number_of_introns, from=0, to=196)
d_LVG <- density(Stable_5000bp$Number_of_introns, from=0, to=196)

HVG_LVG_distribution_5000bp[,2] <- d_HVG$y
HVG_LVG_distribution_5000bp[,3] <- d_LVG$y

colnames(HVG_LVG_distribution_5000bp) <- c("Number_of_introns", "HVG", "LVG")
HVG_LVG_distribution_5000bp_long <- gather(HVG_LVG_distribution_5000bp, type, distr, HVG:LVG)

# graph

HVG_LVG_distribution_5000bp <- mutate(HVG_LVG_distribution_5000bp, mean_bootstrat=rowMeans(random_distribution_5000bp_1000[,2:1001]),
                                      bottom_CI=mean_bootstrat-rowSds(as.matrix(random_distribution_5000bp_1000[,2:1001])),
                                      top_CI=mean_bootstrat+rowSds(as.matrix(random_distribution_5000bp_1000[,2:1001])))



gather(HVG_LVG_distribution_5000bp, type, distr, HVG:top_CI, factor_key=TRUE) %>%
  ggplot(aes(x=Number_of_introns, y=distr, color=type, linetype=type)) +
  geom_line(size=1) +
  scale_x_continuous(limits = c(0, 130)) +
  xlab("Number of introns") +
  ylab("Density") +
  ggtitle("Genes more than 5000bp") +
  scale_color_manual(values=c("blue","green4","black","grey30","grey30")) +
  scale_linetype_manual(values=c("solid","solid","solid", "dashed", "dashed")) +
  theme(panel.grid.major.x = element_blank(),panel.grid.major.y = element_line( size=.1, color="black" ) 
        , panel.background = element_rect(fill="white",colour="grey50"), text=element_text(size=38))


# wilcoxon test

wilcox.test(Variable_5000bp$Number_of_introns,Stable_5000bp$Number_of_introns)$p.value
[1] 0.933325


Wilcox_pvalues=matrix(nrow=1000,ncol=1)
Wilcox_pvalues=data.frame(Wilcox_pvalues)

for (i in 1:1000) {
  Wilcox_pvalues[i,] <- wilcox.test(Variable_5000bp$Number_of_introns,random_5000bp_1000[,i])$p.value
}


#All values between 1.0000000 and 0.03543054 #



