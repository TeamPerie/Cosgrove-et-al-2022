#!/bin/bash
#  run_igor_tags.sh
#
#  script to run IGoR to infer probabilistic model for VDJ barcode generation

FASTA_DIR="/data/users/alyne/igor_install/VDJ_fasta_files"

cd "$OUT_DIR"

NAME="Jason_igor"
mkdir "$NAME"
mkdir "${NAME}"/logs

cd "$NAME"

#implement pipeline in IGoR as written by Yuval

#align full reads to V and J
igor -set_wd "$OUT_DIR"/"${NAME}" -stdout_f logs/"${NAME}".out -read_seqs jason_vdj_barcodes.txt -batch "${NAME}" -set_genomic --V "$FASTA_DIR"/V.fasta --D "$FASTA_DIR"/D.fasta --J "$FASTA_DIR"/J.fasta -set_CDR3_anchors --V "$FASTA_DIR"/V_gene_anchor.csv --J "$FASTA_DIR"/J_gene_anchor.csv -align --ntCDR3 --V --J

echo "aligning step done"


#create empty D alignments file
>aligns/"${NAME}"_D_alignments.csv


#infer generation model
igor -set_wd "$OUT_DIR"/"${NAME}" -stdout_f logs/"${NAME}"_infer.out -batch "${NAME}" -chain beta -species human -set_genomic --V "$FASTA_DIR"/V.fasta --D "$FASTA_DIR"/D.fasta --J "$FASTA_DIR"/J.fasta -infer --N_iter 100 --P_ratio_thresh 0.0001

echo "inference step done"


#calculate pgen of aligned sequences using previously learned model
igor -set_wd "$OUT_DIR"/"${NAME}" -stdout_f logs/"${NAME}"_pgen.out -batch "${NAME}" -set_genomic --V "$FASTA_DIR"/V.fasta --D "$FASTA_DIR"/D.fasta --J "$FASTA_DIR"/J.fasta -set_custom_model "${NAME}"_inference/final_parms.txt "${NAME}"_inference/final_marginals.txt -evaluate -output --Pgen

echo "pgen computation done"


R_FILE="$SCRIPT_DIR"/make_igor_output_file.R
$R_EXEC $R_FILE "$OUT_DIR"/Jason_igor/jason_vdj_barcodes.txt "$OUT_DIR"/Jason_igor/Jason_igor_output/Pgen_counts.csv "$OUT_DIR"/Jason_igor/jason_bcs_pgens.csv




