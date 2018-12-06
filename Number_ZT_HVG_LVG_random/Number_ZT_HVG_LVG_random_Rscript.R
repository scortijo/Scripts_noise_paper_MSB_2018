
library(tidyverse)



## Number of ZT variable for HVGs, LVGs and 1000 sets of random genes----


Background <- read.table("Number_ZT_HVG_LVG_random/All_genes_background.txt", header=TRUE, sep='\t')
dim(Background)
Background$type <- "background"



# extract 1000 sets of random genes 


random_nbZT_1000=matrix(nrow=5332,ncol=1000)
random_nbZT_1000=data.frame(random_nbZT_1000)

random_set <- Background[sample(nrow(Background), 313),]
random_set <- bind_rows(random_set,Background[sample(nrow(Background), 257),])
random_set <- bind_rows(random_set,Background[sample(nrow(Background), 374),])
random_set <- bind_rows(random_set,Background[sample(nrow(Background), 327),])
random_set <- bind_rows(random_set,Background[sample(nrow(Background), 310),])
random_set <- bind_rows(random_set,Background[sample(nrow(Background), 636),])
random_set <- bind_rows(random_set,Background[sample(nrow(Background), 418),])
random_set <- bind_rows(random_set,Background[sample(nrow(Background), 630),])
random_set <- bind_rows(random_set,Background[sample(nrow(Background), 408),])
random_set <- bind_rows(random_set,Background[sample(nrow(Background), 716),])
random_set <- bind_rows(random_set,Background[sample(nrow(Background), 566),])
random_set <- bind_rows(random_set,Background[sample(nrow(Background), 377),])
Nb_ZT <- dplyr::count(random_set, Gene)
Nb_ZT_compact <- dplyr::count(Nb_ZT,n)

for (i in 1:999) {
  random_set <- Background[sample(nrow(Background), 313),]
  random_set <- bind_rows(random_set,Background[sample(nrow(Background), 257),])
  random_set <- bind_rows(random_set,Background[sample(nrow(Background), 374),])
  random_set <- bind_rows(random_set,Background[sample(nrow(Background), 327),])
  random_set <- bind_rows(random_set,Background[sample(nrow(Background), 310),])
  random_set <- bind_rows(random_set,Background[sample(nrow(Background), 636),])
  random_set <- bind_rows(random_set,Background[sample(nrow(Background), 418),])
  random_set <- bind_rows(random_set,Background[sample(nrow(Background), 630),])
  random_set <- bind_rows(random_set,Background[sample(nrow(Background), 408),])
  random_set <- bind_rows(random_set,Background[sample(nrow(Background), 716),])
  random_set <- bind_rows(random_set,Background[sample(nrow(Background), 566),])
  random_set <- bind_rows(random_set,Background[sample(nrow(Background), 377),])
  Nb_ZT <- dplyr::count(random_set, Gene)
  Nb_ZT_compact <- dplyr::count(Nb_ZT,n) %>%
    full_join(Nb_ZT_compact, ., by="n")
}


names(Nb_ZT_compact) = c("nbZT",paste("nb_genes_set",1:1000, sep=""))

Nb_ZT_compact[is.na(Nb_ZT_compact)] <- 0




nbZT_summary_1000sets=matrix(nrow=6,ncol=1)
nbZT_summary_1000sets=data.frame(nbZT_summary_1000sets)

nbZT_summary_1000sets[,1] <- Nb_ZT_compact$nbZT

nbZT_summary_1000sets <- mutate(nbZT_summary_1000sets, mean_bootstrat=rowMeans(Nb_ZT_compact[,2:1001]),
                                bottom_CI=mean_bootstrat-rowSds(as.matrix(Nb_ZT_compact[,2:1001])),
                                top_CI=mean_bootstrat+rowSds(as.matrix(Nb_ZT_compact[,2:1001])))


AllVariable <- read.table("Number_ZT_HVG_LVG_random/HVG_allZT_nbZT.txt",header=TRUE, sep='\t')
dim(AllVariable)
1358    2
names(AllVariable)
[1] "allZT" "nbZT" 
AllVariable <- mutate(AllVariable, TypeV = "HVG")

AllVariable_compact <- dplyr::count(AllVariable, nbZT)

nbZT_summary_1000sets <- full_join(nbZT_summary_1000sets, AllVariable_compact,
                                   by=c("nbZT_summary_1000sets"="nbZT"))

names(nbZT_summary_1000sets) <- c("nbZT_summary_1000sets", "mean_bootstrat",
                                  "bottom_CI", "top_CI", "HVG")

AllStable <- read.table("Number_ZT_HVG_LVG_random/LVG1000_allZT_nbZT.txt",header=TRUE, sep='\t')
dim(AllStable)
[1] 5727     2
names(AllStable)
[1] "allZT" "nbZT"  
AllStable <- mutate(AllStable, TypeV = "LVG")


AllStable_compact <- dplyr::count(AllStable, nbZT)

nbZT_summary_1000sets <- full_join(nbZT_summary_1000sets, AllStable_compact,
                                   by=c("nbZT_summary_1000sets"="nbZT"))

names(nbZT_summary_1000sets) <- c("nbZT_summary_1000sets", "mean_bootstrat",
                                  "bottom_CI", "top_CI", "HVG", "LVG")


nbZT_summary_1000sets[is.na(nbZT_summary_1000sets)] <- 0





gather(nbZT_summary_1000sets, type, distr, mean_bootstrat:LVG, factor_key=TRUE) %>%
  ggplot(aes(x=nbZT_summary_1000sets, y=distr, color=type, linetype=type)) +
  geom_line(size=1) +
  scale_x_continuous(limits = c(0, 12),
                     breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12)) +
  xlab("Number of ZT selected") +
  ylab("Number of genes") +
  ggtitle("") +
  scale_color_manual(values=c("black","grey30","grey30","blue","green4")) +
  scale_linetype_manual(values=c("solid", "dashed", "dashed","solid","solid")) +
  theme(panel.grid.major.x = element_blank(),panel.grid.major.y = element_line( size=.1, color="black" ) 
        , panel.background = element_rect(fill="white",colour="grey50"), text=element_text(size=22))


