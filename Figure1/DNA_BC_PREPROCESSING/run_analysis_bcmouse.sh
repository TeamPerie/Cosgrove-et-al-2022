#!/bin/bash
#
# Adapted A Lyne 31/10/2019
# script to run various parts of barcode pipeline. Extracting barcodes


# working directory: TPerie/barcode_mouse_DNA_pipeline_AML/VDJ_V616 (symbolic link to /data/kdi_prod/project_result/1056/08.00)


# directories where raw data is stored, each sequencing cohort separated by ";"
DATA_DIRS="/data/kdi_prod/project_result/1056/08.00/data_long_reads"
# avant "/data/tmp/example_fastqs" mais à disparu comme toutes les data de ce projet ...

# sample name stem for each sequencing cohort (assumes files are named: "$SAMP_STEM$i".R1.fastq.gz)
# each sequencing cohort separated by ";"
SAMP_STEMS="V618"
NSEQ=18 #number of sequencing runs

#output directories to be made by user, one per sequencing cohort, each separated by ";"
#preferably somewhere permanent and backed up 
#structure should be as below, one main directory e.g. Emilie_Jacques_bc_mouse, then one subdirectory per sequencing cohort named with cohort stem
PERM_DIRS="/data/kdi_prod/project_result/1056/08.00/outputs"

i=1 #run code for this sequencing cohort
#extract relevant directories for this cohort
DATA_DIR=$(echo $DATA_DIRS | awk -v l="$i" 'BEGIN {FS=";"} {print $l}') #/data/kdi_prod/project_result/1056/08.00/data/2013758"
SAMP_STEM=$(echo $SAMP_STEMS | awk -v l="$i" 'BEGIN {FS=";"} {print $l}') #V616
PERM_DIR=$(echo $PERM_DIRS | awk -v l="$i" 'BEGIN {FS=";"} {print $l}') #/data/kdi_prod/project_result/1056/08.00/outputs

# directory where barcode mouse perl scripts are saved
BC_MOUSE_DIR="/data/kdi_prod/project_result/1056/08.00/bc_mouse"

#directory on cluster where this script is saved
SCRIPT_DIR="/data/kdi_prod/project_result/1056/08.00/Emilie_Jacques_bc_mouse/"

#command to execute an R script
R_EXEC="/bioinfo/local/build/Centos/R/R-3.4.0/bin/Rscript"

#find how many samples were sequenced
NSAMP_SPLIT=$(wc -l "$PERM_DIR"/"$SAMP_STEM".sampleDescription.txt | awk '{print $1}')
#half to get number biological samples
NSAMP=$((NSAMP_SPLIT/2))

THRESH=2 #read threshold to use for barcode filtering analysis
TRIMMED=1 #set to 1 if you want to use trimmed reads, 0 if not
FULL_LENGTH=0 #set to 1 if you only want to use full length reads, 0 if not

UMI_L=16 #length of UMI
CONST_ST=17 #start of constant part in read
CONST_EN=40 #end of constant part in read
READ_L=150 #length of sequencing read, best to check this in one of the fastq files

####°°°°°°°°°°°°°°°°°°°°°°°°°°°  INSTALLATION  °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°##
#run code to install IGoR for barcode probability analysis
#check paths and version match in install_igor_cluster.sh
#only need to do this once
# qsub -N igor_inst -l nodes=1:ppn=1 -l walltime=5:00:00 -l mem=5gb -k oe install_igor_cluster.sh

#run code to install perl on cluster
#only need to do this once
# qsub -N igor_inst -l nodes=1:ppn=1 -l walltime=5:00:00 -l mem=5gb -k oe install_perl_cluster.sh
####°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°##


####°°°°°°°°°°°°°°°°°°°°°°°°°°°  PIPELINE STEPS  °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°##

##################_______________________________Part 1_________________________________#####################
##----------------- STEP 1 ----------------##
## run xcalibr to extract barcodes from reads, one task per sequenced sample (so if sample=R14 use -t 14 or if all: -t 1-"$NSAMP_SPLIT")

RUNID=$(qsub -N bcm_xcalibr -l nodes=1:ppn=6 -l walltime=1:00:00 -l mem=45gb -t 1-"$NSAMP_SPLIT" -v UMI_L="$UMI_L",CONST_ST="$CONST_ST",CONST_EN="$CONST_EN",READ_L="$READ_L",DATA_DIR="$DATA_DIR",BC_MOUSE_DIR="$BC_MOUSE_DIR",PERM_DIR="$PERM_DIR",SAMP_STEM="$SAMP_STEM",TRIMMED="$TRIMMED",FULL_LENGTH="$FULL_LENGTH" bcmouse_xcalibr.sh)

##################_______________________________Part 2_________________________________######################
#run Joost's pipeline to get list of observed barcodes, counts and alignment of VDJ
RUNID=$(qsub -N bcm_clean_tags -l nodes=1:ppn=4 -l walltime=1:00:00 -l mem=40gb -t 1-"$NSAMP_SPLIT" -k oe -v PERM_DIR="$PERM_DIR",R_EXEC="$R_EXEC",SCRIPT_DIR="$SCRIPT_DIR",SAMP_STEM="$SAMP_STEM",THRESH="$THRESH" get_cleaned_tags.sh)
