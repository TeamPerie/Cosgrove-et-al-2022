
rm(list =ls())
setwd("/Users/jasoncosgrove/Dropbox (Team_Perie)/Jason/Experiments/Dry_Lab/MetaFate_final_pipeline/pipeline_files/Step4_RNA_BC_processing/RNA_barcodes/3_045_1_045")


file.names <- list.files()

for(i in 1:length(file.names)){
  f <- read.csv(file.names[i])
  
  print("sample: ")
  print(file.names[i])
  
  print("number of cells: ")
  print(length(unique(f$Cell.bc)))
  
  print("number of barcodes: ")
  print(length(unique(f$Con.seq)))
  
}