



## Percentage of chromatin marks----


Background <- read.table("Chromatin_proportion/All_genes_background.txt", header=TRUE, sep='\t')
dim(Background)
Background$type <- "background"

Chromatin_marks <- read.table("Chromatin_proportion/Gene_list_chromatin_marks.txt", 
                              header=TRUE, sep='\t')
names(Chromatin_marks)
[1] "H3K4me2"       "H3K4me3"       "H3K27me1"      "H3K27me3"      "H3K36me3"     
[6] "H2Bub"         "X5mC"          "Gene_highH2AZ" "Gene_lowH2AZ" 

# extract 1000 sets of random genes 


random_percentage_chromatin_1000=matrix(nrow=1000,ncol=9)
random_percentage_chromatin_1000=data.frame(random_percentage_chromatin_1000)
names(random_percentage_chromatin_1000) <- names(Chromatin_marks)

for (i in 1:1000) {
  random_set <- Background[sample(nrow(Background), 1358),]
  random_percentage_chromatin_1000[i,1] <- nrow(inner_join(random_set, Chromatin_marks,
                                                           by=c("Gene"="H3K4me2")))
  random_percentage_chromatin_1000[i,2] <- nrow(inner_join(random_set, Chromatin_marks,
                                                           by=c("Gene"="H3K4me3")))
  random_percentage_chromatin_1000[i,3] <- nrow(inner_join(random_set, Chromatin_marks,
                                                           by=c("Gene"="H3K27me1")))
  random_percentage_chromatin_1000[i,4] <- nrow(inner_join(random_set, Chromatin_marks,
                                                           by=c("Gene"="H3K27me3")))
  random_percentage_chromatin_1000[i,5] <- nrow(inner_join(random_set, Chromatin_marks,
                                                           by=c("Gene"="H3K36me3")))
  random_percentage_chromatin_1000[i,6] <- nrow(inner_join(random_set, Chromatin_marks,
                                                           by=c("Gene"="H2Bub")))
  random_percentage_chromatin_1000[i,7] <- nrow(inner_join(random_set, Chromatin_marks,
                                                           by=c("Gene"="X5mC")))
  random_percentage_chromatin_1000[i,8] <- nrow(inner_join(random_set, Chromatin_marks,
                                                           by=c("Gene"="Gene_highH2AZ")))
  random_percentage_chromatin_1000[i,9] <- nrow(inner_join(random_set, Chromatin_marks,
                                                           by=c("Gene"="Gene_lowH2AZ")))
}


random_percentage_chromatin_mean_CI=matrix(nrow=3,ncol=9)
random_percentage_chromatin_mean_CI=data.frame(random_percentage_chromatin_mean_CI)
names(random_percentage_chromatin_mean_CI) <- names(Chromatin_marks)
row.names(random_percentage_chromatin_mean_CI) <- c("mean", "bottom_CI", "top_CI")

for (i in 1:9) {
  random_percentage_chromatin_mean_CI[1,i] <- mean(random_percentage_chromatin_1000[,i])
  random_percentage_chromatin_mean_CI[2,i] <- mean(random_percentage_chromatin_1000[,i])-sd(random_percentage_chromatin_1000[,i])
  random_percentage_chromatin_mean_CI[3,i] <- mean(random_percentage_chromatin_1000[,i])+sd(random_percentage_chromatin_1000[,i])
}

write.table(random_percentage_chromatin_mean_CI, "Chromatin_proportion/random_percentage_chromatin_mean_CI.txt", sep='\t')


