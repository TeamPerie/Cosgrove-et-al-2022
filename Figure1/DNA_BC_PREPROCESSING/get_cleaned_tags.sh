#!/bin/bash
#
#  
#script to run Joost's pipeline to get cleaned "tags" (barcodes) from UMI-collapsed output

k="$PBS_ARRAYID"
if [[ ${#k} -eq 1 ]]; then k="0${k}"; fi
SAMP_DIR="${SAMP_STEM}R${k}_trimmed"
#A1202R{number}_trimmed

cd "$PERM_DIR"
#/data/tmp/lhadjabe/barcode_mouse_DNA_pipeline_AML/outputs/A1202

#make directory for output specific to a given threshold
mkdir output_thresh"$THRESH"

#run R file
R_FILE="$SCRIPT_DIR"/tagcleanup.R
$R_EXEC $R_FILE "$PERM_DIR"/"$SAMP_DIR" "$THRESH" "$SCRIPT_DIR"


#"$PERM_DIR"/"$SAMP_DIR" = /data/tmp/lhadjabe/barcode_mouse_DNA_pipeline_AML/outputs/A1202/A1202R{number}_trimmed


