#!/bin/bash
#  run_get_consensus_VDJ.sh
#
#  Created by A Lyne on 20/08/2018.
#

#make file with all VDJ barcodes

cd "/Users/alyne/Dropbox (Personal)/Curie/Curie_Paris/anne-marie/Jason_barcode_mouse"

#remove quotation marks
awk 'BEGIN {FS=","} NR>1 {print $1}' DNA_barcodes/DRAG_DNA_barcodes_normalised_and_ab_filtered2.csv | awk '{print substr($0,2,length()-2);}' > tags_complete/jason_vdj_barcodes_temp.txt




cd "/Users/alyne/Dropbox (Personal)/Curie/Curie_Paris/anne-marie/Jason_barcode_mouse/output_Jason_032022/Jason_bc_mouse4_no_mm"

NSAMP=$(ls | grep -e '^agrep_reads' | wc -l)

#iterate through each sample
for (( i=1; i<=$NSAMP; i++ )); do

DIR=$(ls | grep -e '^agrep_reads' | awk -v l=$i 'NR==l')

#find file top_bc_per_cell_VDJ_barcodes_*_52_10_both_splitVDJ.csv
if [ "$DIR" != "agrep_reads_ET27_LSKs_V638T1" ]; then

FILE=$(ls "$DIR"/output_3_045_1_045 | grep "both_splitVDJ")

#print second column which stores VDJ barcode
awk 'BEGIN {FS=","} NR>1 && $8==1 {print $2}' "$DIR"/output_3_045_1_045/"$FILE" >> ../../tags_complete/jason_vdj_barcodes_temp.txt

fi

done





cd "/Users/alyne/Dropbox (Personal)/Curie/Curie_Paris/anne-marie/Jason_barcode_mouse/tags_complete"

#sort and take unique
sort -u jason_vdj_barcodes_temp.txt > jason_vdj_barcodes.txt





