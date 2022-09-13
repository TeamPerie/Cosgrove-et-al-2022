#!/usr/bin/env Rscript

# if in the cluster
#.libPaths(setdiff(.libPaths(), normalizePath(Sys.getenv("R_LIBS_USER"))))



library(dplyr)
library(plyr)
library(readr)

args=commandArgs(trailingOnly=TRUE)
analysis=as.character(args[1])
lib=as.character(args[2])

file_list <- list.files("./outputs/tmp/")
i=0
for(file in file_list){
  data<-read.csv(paste("./outputs/tmp/", file, sep = ''), header = T, sep = "\t")
  if(lib=="new"){
    # sum columns in each file
    filename=sub('.csv', '', file)
    data$tot<-rowSums(data[,-1])
    # if 2 indices: replace data[,c(1,6)] by data[,c(1,4)]
    data<-data[,c(1,6)]
    names(data)[names(data) == "tot"] <-filename
    #combine matrices
    if(i==0){ 
      previous_file <- data                         
    }else{ 
      previous_file <- merge(previous_file, data, by = "name", all = TRUE)
    }
    i=i+1
  }else{
    previous_file<-select(data,-contains("mismatch"))
  }
}

previous_file<-replace(previous_file,is.na(previous_file),0)
previous_file<-previous_file[rowSums(previous_file[, -1])>0, ]
previous_file<-ddply(previous_file, "name", unique)
previous_file<-previous_file[-which(previous_file$name=="nohit_rows"),]

write_csv(previous_file, paste("./outputs/final_matrices/CountMatrix_", analysis, ".csv", sep = ""))



