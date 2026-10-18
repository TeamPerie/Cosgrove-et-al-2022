
df <- read.csv("/Users/jasoncosgrove/Dropbox (Team_Perie)/Jason/Experiments/Dry_Lab/MetaFate_final_pipeline/data/Step4_RNA_BC_processing/split_VDJs/Split_VDJs_expts_merged.csv")


colnames(df)

hist(nchar(df[,2]) + nchar(df[,4]), main = "number of insertions", xlab = "insertions per barcode")
hist(nchar(df[,6]) + nchar(df[,7]) + nchar(df[,8]) + nchar(df[,9]), main = "number of deletions", xlab = "deletions per barcode")


bc.pairs <- read.csv('/Users/jasoncosgrove/Dropbox (Team_Perie)/Jason/Experiments/Dry_Lab/MetaFate_final_pipeline/data/Step6_10X_Data_integration/DNA_RNA_BC_PAIRS/all_dna_rna_bc_pairs.csv')

hist(nchar(bc.pairs$Con.seq), main = "length of DRAG barcodes", xlab = "barcode length (bases)")



#now lets look at DNA sequences.



#samples 5-24 are the ones we keep in our final analyses

setwd("/Users/jasoncosgrove/Dropbox (Team_Perie)/Jason/Experiments/Dry_Lab/MetaFate_final_pipeline/data/Step3_DNA_BC_processing/outputs/threshold_1/")


files <- list.files()[grepl("split" ,list.files())]
files.subset <- files[5:24]


test <- files.subset %>% lapply(read.table,fill = TRUE,header = 1)

df <- rbind(test[1][[1]],
test[2][[1]],
test[3][[1]],
test[4][[1]],
test[5][[1]],
test[6][[1]],
test[7][[1]],
test[8][[1]],
test[9][[1]],
test[10][[1]]

)



    


hist(nchar(df[,2]) + nchar(df[,4]), main = "number of insertions", xlab = "insertions per DNA barcode")
hist(nchar(df[,6]) + nchar(df[,7]) + nchar(df[,8]) + nchar(df[,9]), main = "number of deletions", xlab = "deletions per DNA barcode")
