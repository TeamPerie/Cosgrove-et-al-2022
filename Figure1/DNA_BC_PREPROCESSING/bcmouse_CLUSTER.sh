#!/bin/bash
#  bcmouse.sh
#
#  Created by Arno Velds
# Adapted A Lyne 31/10/2019
# script to extract VJD barcodes from sequencing reads

#PATH=$PATH:/bioinfo/local/build/xcalibr/xcalibr-20151125/xcalibr
#PATH=$PATH:/bioinfo/local/build/xcalibr/xcalibr-20151125/bin/


export PATH=$PATH:/data/kdi_prod/project_result/1056/08.00/blast2/bin

export LD_LIBRARY_PATH=/data/kdi_prod/project_result/1056/08.00/blast2/lib:$LD_LIBRARY_PATH

echo "UMI has length $UMI_L"
echo "Constant part is bases $CONST_ST to $CONST_EN"

#get sequencing platform ID
k="$PBS_ARRAYID"
if [[ ${#k} -eq 1 ]]; then k="0${k}"; fi
SAMPLE="${SAMP_STEM}R${k}"

#move to temporary directory
cd "$TMPDIR"

#copy trimmed (or not) file from data directory to temp dir
if [[ "$TRIMMED" == 0 ]]; then

cp "$DATA_DIR"/"$SAMPLE"/"$SAMPLE".R1.fastq.gz .

else

cp "$DATA_DIR"/"$SAMPLE"/"$SAMPLE"_trimmed.R1.fastq.gz "$SAMPLE".R1.fastq.gz

fi


#remove short reads if FULL_LENGTH==1
if [[ "$FULL_LENGTH" == 1 ]]; then

#unzip fastq
gzip -d "$SAMPLE".R1.fastq.gz

paste - - - - < "$SAMPLE".R1.fastq | awk -v l="$READ_L" 'BEGIN {FS="\t"} length($2)==l' | tr '\t' '\n' | gzip > "$SAMPLE".R1.fastq.gz
rm "$SAMPLE".R1.fastq

fi


#find index sequence from fastq file
#NB this assumes lines in fastq file with index are in format:
# @M04391:341:000000000-G3W64:1:1101:19012:1881 1:N:0:TTGGTTG
#and that the most abundant one is the correct one
INDEX=$(zcat "$SAMPLE".R1.fastq.gz | head -n 10000 | grep '1:N:0:' | awk '{print $2}' | awk 'BEGIN {FS=":"} {print $4}'| sort | uniq -c | sort -nr -k1,1 | awk 'NR==1 {print $2}')
echo "index is ${INDEX}"

zcat "$SAMPLE".R1.fastq.gz | head -n 10000 | grep '1:N:0:' | awk '{print $2}' | awk 'BEGIN {FS=":"} {print $4}' | sort | uniq -c | sort -nr -k1,1 | head

# Normal demultiplexing allows 1bp distance,
# but we collect only the perfect matches with the index sequences.
# these reads are piped through 'xcalibr hash' to create the analysis data
zcat "$SAMPLE".R1.fastq.gz | perl "$BC_MOUSE_DIR"/perfect.pl "$INDEX" | /bioinfo/local/build/xcalibr/xcalibr-20151125/xcalibr hash "$SAMPLE"-perfect.bin
ls -l "$SAMPLE"-perfect.bin

# we test the hash output for the expected constant part. The illumina sequencer tends to make mistakes here. We
# collect the top hits to create a constant template sequence that discards as little data as possible
/bioinfo/local/build/xcalibr/xcalibr-20151125/xcalibr analyze "$SAMPLE"-perfect.bin "$CONST_ST"-"$CONST_EN" > "$SAMPLE"-constant.txt
cat "$SAMPLE"-constant.txt
cat *constant.txt | grep "^C.CG" | cut -f1 | sort | uniq -c | sort -n | tail -n20
#67		CTCGAGGCCATCGAAGTATCAAGTCC
#112	CTCGAGGTCATCGAAGAATCAAGTCC
#131	CNCGAGGNCATCGAAGTATCAAGTNC
#131	CNCGAGGNCATCGAAGTATCAAGTCC
#132	CTCGAGGTCATCGAAGTATCAAGTCC
#132	CTCGAGGNCATCGAAGTATCAAGTCC

# use constant: CTCGAGGNCATCGAAGTATCAAG
# note the N's in the variable positions

# extract a count table from the hash using the template. We use the sparse output to avoid creating a huge matrix with al lot of 0's
/bioinfo/local/build/xcalibr/xcalibr-20151125/xcalibr extract --template N1X"$UMI_L"CTCGAGGNCATCGAAGTATCAAGYn --rows Y --cols X --requireseq --sparse --output "$SAMPLE"-perfect-long.txt "$SAMPLE"-perfect.bin &> "$SAMPLE".log
cat "$SAMPLE".log
ls
# apply post-processing to the sparse count table.
# we expect a constant sequence (template2.fa) to start somewhere in every sequence
# we use ncbi blast2 to find these hits
# we hope to find only one hit per sequence, but this is not always the case.
# a perl script uses the blast results to process the long output

F="$SAMPLE"-perfect-long.txt
BASE=${F%.txt}
#Blast against template clip tags
cut -f1,2 --output-delimiter=Y $F | perl "$BC_MOUSE_DIR"/count2fasta.pl > $BASE.fa
blast2 -p blastn -i $BASE.fa -j "$BC_MOUSE_DIR"/template2.fa  -W7 -K1  -e 1e-6 -m8 >$BASE-blast-template2-e6.txt
#check counts
echo "Tailing $BASE, should be 1 count:"
cut -f1 $BASE-blast-template2-e6.txt | sort | uniq -c | sort -n | tail -5
#cut tags
perl "$BC_MOUSE_DIR"/cuttagslong.pl $F $BASE-blast-template2-e6.txt > $BASE-clipped-template2-e6.txt

rm "$SAMPLE".R1.fastq.gz "$SAMPLE"-perfect.bin "$SAMPLE".R1.fastq "$SAMPLE".R1.fastq.gz "$SAMPLE"-perfect-long.fa "$SAMPLE"-perfect-long.txt "$SAMPLE".log


#make directory for permanent output, then copy all output there
if [[ "$TRIMMED" == 0 ]]; then

if [[ "$FULL_LENGTH" == 0 ]]; then

mkdir "$PERM_DIR"/"$SAMPLE"
cp * "$PERM_DIR"/"$SAMPLE"/

else

mkdir "$PERM_DIR"/"$SAMPLE"_full_length
cp * "$PERM_DIR"/"$SAMPLE"_full_length/

fi

else

if [[ "$FULL_LENGTH" == 0 ]]; then

mkdir "$PERM_DIR"/"$SAMPLE"_trimmed
cp * "$PERM_DIR"/"$SAMPLE"_trimmed/

else

mkdir "$PERM_DIR"/"$SAMPLE"_trimmed_full_length
cp * "$PERM_DIR"/"$SAMPLE"_trimmed_full_length/

fi

fi





