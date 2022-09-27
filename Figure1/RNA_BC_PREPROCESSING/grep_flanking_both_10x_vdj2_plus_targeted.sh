#!/bin/bash
#  grep_flanking_10x_vdj.sh
#
#  Created by A Lyne on 26/04/2018.
#

LINE=$(awk -v l="$TARG" 'BEGIN {FS="|"} $1==l {print $0}' "/data/kdi_prod/dataset_all/${KDI_NO}/export/user/${SEQ_NO}.sampleDescription.txt")

echo $LINE
echo $TARG
echo $SAMPLE

cd "$AGREP_DIR"

#extract unmapped reads
#this step is quite slow and only needs doing once, you can comment this line out
#when trying different parameters in the parts below
#extract unmapped reads from original 10x sequencing
samtools view -f 4 "$COUNT_DIR"/possorted_genome_bam.bam > ../all_reads_"${SAMPLE}".txt
#add unmapped reads from targeted sequencing
samtools view -f 4 "$TARG_TMP"/"$TARG"/CR_output/CR_"$TARG"/outs/possorted_genome_bam.bam >> ../all_reads_"${SAMPLE}".txt


#make sequences to grep bases on start/end of V/J plus some flanking sequence
#add flanking and primer together
TOT_V="${FL_V}${PR_V}"
NCHAR=${#TOT_V}
START=$((NCHAR-NMATCH_V+1))
GREP_L=$(echo $TOT_V | cut -c "$START"-"$NCHAR")

TOT_J="${PR_J}${FL_J}"
NCHAR=${#TOT_J}
GREP_R=$(echo $TOT_J | cut -c 1-"$NMATCH_J")

echo $GREP_L
echo $GREP_R


#make regular expression with 5' constant part, 10-40 bases variable and then 3' constant part
regex="${GREP_L}[ACGT]{10,40}${GREP_R}"

#extract reads with regex match
awk -v regex="$regex" '$0 ~ regex' ../all_reads_"${SAMPLE}".txt > reads_pattern_match_"${SAMPLE}".txt

wc -l reads_pattern_match_"${SAMPLE}".txt


# take reads which have10x BC and UMI, and print sequence, read name, CB and UB
awk 'BEGIN {FS="\t";OFS="\t"} { for (i=10; i<=NF; ++i) { for (j=10; j<=NF; ++j) { if ($i ~ "CB:" && $j ~ "UB:") print substr($i,6),$10,$1,substr($j,6) } } }' reads_pattern_match_"${SAMPLE}".txt > unaligned_10xbc_seq_"$SAMPLE".txt

wc -l unaligned_10xbc_seq_"$SAMPLE".txt


#keep only reads which are in cells (as defined by cellranger)
#check whether bc list file is zipped or not
NCHAR=${#BC_LIST}
ST=$((NCHAR-1))
FIN2=$(echo $BC_LIST | cut -c "$ST"-"$NCHAR")

if [[ "$FIN2" == "gz" ]]; then
awk 'NR == FNR { x[$1] = 1; next; } ($1 in x)' <(gzip -dc "$BC_LIST") FS=" " unaligned_10xbc_seq_"$SAMPLE".txt > unaligned_10xbc_seq_"$SAMPLE"_cells.txt
else
awk 'NR == FNR { x[$1] = 1; next; } ($1 in x)' "$BC_LIST" FS=" " unaligned_10xbc_seq_"$SAMPLE".txt > unaligned_10xbc_seq_"$SAMPLE"_cells.txt
fi

wc -l unaligned_10xbc_seq_"$SAMPLE"_cells.txt


#extract regex match (VDJ barcodes plus flanking regions) from reads
grep -Eo "$regex" unaligned_10xbc_seq_"$SAMPLE"_cells.txt > regexs_10xbc_"$SAMPLE"_cells.txt


#add VDJ barcodes/flanking as a column, remove flanking, rearrange columns
NMATCH_TOT=$((NMATCH_V+NMATCH_J))

paste unaligned_10xbc_seq_"$SAMPLE"_cells.txt regexs_10xbc_"$SAMPLE"_cells.txt | awk -v NMATCH_V="$NMATCH_V" -v NMATCH_TOT="$NMATCH_TOT" 'BEGIN {FS="\t";OFS="\t"} {NCHAR=length($5); BCLEN=NCHAR-NMATCH_TOT; $5=substr($5,NMATCH_V+1,BCLEN)} {print $1,$4,$5,$3,$2}' > agrep_10xbc_and_vbc_"$SAMPLE"_both.txt

rm regexs_10xbc_"$SAMPLE"_cells.txt

gzip ../all_reads_"${SAMPLE}".txt









