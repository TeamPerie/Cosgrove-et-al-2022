#!/bin/bash
# run_10x_analysis_vdj_Jason.sh
#
#  Created by A Lyne on 18/04/2018.
#  

#second version of grep, using fewer bases to try and get more barcodes

#make id for analysis, used to name cluster jobs
ID="Jason"
jid="_${ID}"

#MAKE DIRECTORY FOR FINAL OUTPUT
OUT_DIR="/data/users/alyne/lin_trace/svn/analyse/data/Jason_bc_mouse3"
mkdir $OUT_DIR

#MAKE DIRECTORY IN /data/tmp FOR TEMPORARY OUTPUT
TARG_TMP="/data/tmp/Jason_targeted_CR"

#MAKE DIRECTORY IN $HOME FOR SCRIPTS
SCRIPT_DIR="/data/users/alyne/scrnaseq/Jason_10x_barcode_analysis"

#software that is called e.g. samtools in /bioinfo/local/build/Centos, needs to be added to PATH variable in .profile in your home directory on the cluster
# e.g. export PATH=$PATH:/bioinfo/local/build/Centos/samtools/samtools-1.9/bin

#TRANSFER ALL SCRIPTS TO $SCRIPT_DIR

#CHANGE PERMISSIONS OF THIS FILE TO EXECUTABLE chmod +x run_10x_analysis_vdj.sh

#UNCOMMENT EACH QSUB COMMAND AS REQUIRED

#cd to "$SCRIPT_DIR" and run using ./run_10x_analysis_bc_mouse.sh

#location of cellranger count output directory for each sample separated by ;
COUNT_DIRS="/data/kdi_prod/dataset_all/2015177/export/user/count_M1;/data/kdi_prod/dataset_all/2015177/export/user/count_M2;/data/kdi_prod/dataset_all/2015322/export/user/count_M3;/data/kdi_prod/dataset_all/2015322/export/user/count_M4;/data/kdi_prod/dataset_all/2011098/export/user/count_mouse_454;/data/kdi_prod/dataset_all/2011098/export/user/count_mouse_534;/data/kdi_prod/dataset_all/2012664/export/user/count_ET27_Dec2020_LSKs;/data/kdi_prod/dataset_all/2007899/export/user/count_LSK_GFPplus__Mus_musculus;/data/kdi_prod/dataset_all/2012299/export/user/count_ET27_LSKs;/data/kdi_prod/dataset_all/2012299/export/user/count_ET27_LSKs"

#name of each sample
SAMPLES="M1;M2;M3;M4;m454;m534;ET27_Dec2020_LSKs;GFP_plus;ET27_LSKs1;ET27_LSKs2"

#barcode list from cellranger output
BC_LISTS="/data/kdi_prod/dataset_all/2015177/export/user/count_M1/filtered_feature_bc_matrix/barcodes.tsv.gz;/data/kdi_prod/dataset_all/2015177/export/user/count_M2/filtered_feature_bc_matrix/barcodes.tsv.gz;/data/kdi_prod/dataset_all/2015322/export/user/count_M3/filtered_feature_bc_matrix/barcodes.tsv.gz;/data/kdi_prod/dataset_all/2015322/export/user/count_M4/filtered_feature_bc_matrix/barcodes.tsv.gz;/data/kdi_prod/dataset_all/2011098/export/user/count_mouse_454/filtered_feature_bc_matrix/barcodes.tsv.gz;/data/kdi_prod/dataset_all/2011098/export/user/count_mouse_534/filtered_feature_bc_matrix/barcodes.tsv.gz;/data/kdi_prod/dataset_all/2012664/export/user/count_ET27_Dec2020_LSKs/filtered_feature_bc_matrix/barcodes.tsv.gz;/data/kdi_prod/dataset_all/2007899/export/user/count_LSK_GFPplus__Mus_musculus/filtered_gene_bc_matrices/mm10/barcodes.tsv;/data/kdi_prod/dataset_all/2012299/export/user/count_ET27_LSKs/filtered_feature_bc_matrix/barcodes.tsv.gz;/data/kdi_prod/dataset_all/2012299/export/user/count_ET27_LSKs/filtered_feature_bc_matrix/barcodes.tsv.gz"

#sequencing run info
KDI_NOS="2015730;2015327" #KDI number for targeted sequencing
SEQ_NOS="L438;V638"
TARGS="L438T01;L438T02;L438T03;L438T04;L438T05;L438T06;L438T07;L438T08;V638T1;V638T2"
NTARGS=10

#path to Rscript executable
R_EXEC="/bioinfo/local/build/Centos/R/R-4.1.0/bin/Rscript"

#flanking sequences and primer sequences at start/end of V/J
FL_V="ATTGTAATACGACTCACTATAGGGAGACGCGTGTT"
PR_V="ACCTCCTCGAGGTCATCGAAGTATCAAG"
FL_J="CTACTGGAATCAGACCGCCACCATGGTGAGCAAGG"
PR_J="TAGCAAGCTCGAGAGTAGAC"
#number of bases to match outside of flanking region 
#we choose 52 for V in order to match primer, and 10 for J
NMATCH_V=52
NMATCH_J=10

#barcode filtering parameters
MIN_READ=3 # minimum number of reads to retain a barcode
MIN_PROP="0.45" # minimum proportion of read for each UMI assigned to top barcode
MIN_PROP_UMI="0.75"
MIN_UMI=1 # minimum number of UMIs required to retain a cell
MIN_PROP_1UMI="0.45" # minimum proportion of cell reads assigned to top barcode if cell has only one UMI after filtering



##################
## run cell ranger on targeted sequencing data to add CB and UB
##################

for ((i=1;i<="$NTARGS";i++)); do

TARG=$(echo "$TARGS" | awk -v l="$i" 'BEGIN {FS=";";OFS=";"} {print $l}')

#first eight samples are from one sequencing batch
if [[ "$i" -le 8 ]]; then
KDI_NO=$(echo "$KDI_NOS" | awk 'print $1')
else
KDI_NO=$(echo "$KDI_NOS" | awk 'print $2')
fi

RUNID1=$(qsub -l nodes=1:ppn=8 -l walltime=12:00:00 -l mem=48gb -N 10x_count"$jid" -v TARG_TMP="$TARG_TMP",TARG="$TARG",KDI_NO="$KDI_NO" -k oe run_CR_count_targeted_seq.sh)

done





##################
#run barcode extraction on targeted sequencing
##################

for ((i=1;i<="$NTARGS";i++)); do

if [[ "$i" -le 8 ]]; then
KDI_NO=$(echo "$KDI_NOS" | awk 'print $1')
SEQ_NO=$(echo "$SEQ_NOS" | awk 'print $1')
else
KDI_NO=$(echo "$KDI_NOS" | awk 'print $2')
SEQ_NO=$(echo "$SEQ_NOS" | awk 'print $1')
fi

#get sample name according to 10x experiment
SAMPLE=$(echo "$SAMPLES" | awk -v l="$i" 'BEGIN {FS=";";OFS=";"} {print $l}')

#use conversion file to get sample name according to targeted sequencing
SAMPLE2=$(awk -v samp="$SAMPLE" '$2==samp {print $1}' "$SCRIPT_DIR"/sample_conversion.txt)

#use sequencing platform file to get TARG
TARG=$(awk -v samp="$SAMPLE2" 'BEGIN {FS="|"} $2==samp {print $1}' "/data/kdi_prod/dataset_all/${KDI_NO}/export/user/${SEQ_NO}.sampleDescription.txt")

COUNT_DIR=$(echo "$COUNT_DIRS" | awk -v l="$i" 'BEGIN {FS=";";OFS=";"} {print $l}')

BC_LIST=$(echo "$BC_LISTS" | awk -v l="$i" 'BEGIN {FS=";";OFS=";"} {print $l}')

AGREP_DIR="$OUT_DIR"/agrep_reads_"$SAMPLE"_"$TARG"/match_"$NMATCH_V"_"$NMATCH_J"
mkdir "$OUT_DIR"/agrep_reads_"$SAMPLE"_"$TARG"
mkdir "$AGREP_DIR"

# #extract unmapped reads with match to VDJ flanking region
RUNID=$(qsub -l nodes=1:ppn=1 -l walltime=6:00:00 -l mem=10gb -N 10x_grep"$jid" -k oe -v SAMPLE="$SAMPLE",COUNT_DIR="$COUNT_DIR",FL_V="$FL_V",PR_V="$PR_V",FL_J="$FL_J",PR_J="$PR_J",NMATCH_V="$NMATCH_V",NMATCH_J="$NMATCH_J",R_EXEC="$R_EXEC",AGREP_DIR="$AGREP_DIR",BC_LIST="$BC_LIST",TARG_TMP="$TARG_TMP",TARG="$TARG",KDI_NO="$KDI_NO",SEQ_NO="$SEQ_NO" grep_flanking_both_10x_vdj2_plus_targeted.sh)

# #filter to keep high quality UMIs and output one VDJ barcode for each cell
qsub -W depend=afterany:"$RUNID" -l nodes=1:ppn=1 -l walltime=2:00:00 -l mem=5gb -N 10x_con_VDJ"$jid" -k oe -v AGREP_DIR="$AGREP_DIR",NMATCH_V="$NMATCH_V",NMATCH_J="$NMATCH_J",R_EXEC="$R_EXEC",SAMPLE="$SAMPLE",MIN_READ="$MIN_READ",MIN_PROP="$MIN_PROP",MIN_PROP_UMI="$MIN_PROP_UMI",MIN_UMI="$MIN_UMI",MIN_PROP_1UMI="$MIN_PROP_1UMI",SCRIPT_DIR="$SCRIPT_DIR" run_get_consensus_VDJ_reads.sh

done





##################
### Igor pgen computation
##################

#before running this, need to make list of all barcodes using make_bc_list_igor.sh
#from file DRAG_DNA_barcodes_normalised_and_ab_filtered2.csv

#run igor to infer model parameters and compute generation probability
qsub -l nodes=1:ppn=1 -l walltime=5:00:00 -l mem=50gb -N 10x_igor"$jid" -k oe -v OUT_DIR="$OUT_DIR",R_EXEC="$R_EXEC",SCRIPT_DIR="$SCRIPT_DIR" infer_igor_tags.sh



