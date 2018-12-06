

library(tidyverse)

# Combining all files in raw_data folder, which contains TPM, FPKM and raw reads for each of the 14 seedslings for all 12 time points (168 files in total)

filenames <-  list.files(path="1.raw_data", full.names=TRUE)

datalist  <-  lapply(filenames, function(x){read.table(file=x,header=T,sep='\t')})

TPM_FPKM_Raw_14seedling <- Reduce(function(x,y) {merge(x,y)}, datalist)


# Extracting TPM and saving it as TPM_14seedling.txt, which is the basis for all analyses

TPM_14seedling <- select(TPM_FPKM_Raw_14seedling, tracking_id, contains("TPM_"))

write.table(TPM_14seedling, file="2.combining_data/TPM_14seedling.txt", row.names=FALSE, sep='\t')

