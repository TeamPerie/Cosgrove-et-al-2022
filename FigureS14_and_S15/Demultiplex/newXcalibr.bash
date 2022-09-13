#!/bin/bash
################# PROCESS THE OUTPUT OF A CELLULAR BARCODING EXPERIMENT USING XCALIBR  #####################
# author: Jason Cosgrove (jason.cosgrove@curie.fr, UMR168 Perie team)
# date: 06/07/2018
############################################################################################################

# if you have issues running the file it might be because the line endings are not compatible with mac
# tr -d '\r' <XcalibrProcessing.sh >newXcalibr.bash

#Name of the barcode reference library
LIB="LG22"

#path to the plate index library and the barcode reference library
LIB_INDEX_PATH="/Users/jasoncosgrove/Desktop/L400_STB11_STB13/Scripts"

#path to xcalibr
XCALIBR_PATH="/Users/jasoncosgrove/Desktop/L400_STB11_STB13/Scripts/xcalibr-master"

#path to sequencing data, there should be a separate fasta for each sample
DATA_PATH="/Users/jasoncosgrove/Desktop/L400_STB11_STB13/KDI_2013335_2021-05-03_16-05-22"

#Add the prefix for the experiment
PREFIX="L400R"

#add the min and max values for the samples
SAMPLE_MIN=10
SAMPLE_MAX=64



#Add location where the app is to PATH
PATH=$PATH:$XCALIBR_PATH

for ((i=$SAMPLE_MIN;i<=$SAMPLE_MAX;i++))
do

  echo $i

  f1=L400R$i

  echo $f1

  #change working directory to folder where data is
  cd $DATA_PATH/$f1
  
  #unzip the fasta file
  gunzip $f1.R1.fastq.gz

  #if you get a permission denied error make sure to change permissions for xcalibr first using chmod +x 
  cat $f1.R1.fastq | xcalibr hash $f1.bin
  
  #check the statistics for the plate indices, to check the quality of our data
  xcalibr analyze --number 50 $f1.bin 1-8 > $f1-plate_index_stats.txt
  
  #check the statistics for the barcodes, to check the quality of our data
  xcalibr analyze --number 50 $f1.bin 28-47 > $f1-barcode_stats.txt

  #generate the barcode matrix for this sample  
  xcalibr extract --template X8N19Y20 --matchX $LIB_INDEX_PATH/plate_index_library.txt --matchY $LIB_INDEX_PATH/$LIB"_filtered_updated.txt" --printemptyY $f1.bin > $f1-res-all.txt

  #zip the file when you are done to keep the folder size as small as possible
  gzip $f1.R1.fastq

done

cd $LIB_INDEX_PATH

#TODO we need to output some kind of summary file, how many reads do not map to a barcode, or how many times do we get a barcode but no sample index and vice versa. 
Rscript processXcalibrOutput.R $DATA_PATH $LIB_INDEX_PATH $PREFIX SAMPLE_MIN SAMPLE_MAX
