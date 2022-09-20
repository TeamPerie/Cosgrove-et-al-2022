#!/bin/bash
#  run_CR_count_targeted_seq.sh
#
#  Created by A Lyne on 01/04/2022.
# rename files for CR and then run count


mkdir "$TARG_TMP"
cd "$TARG_TMP"

mkdir "$TARG"
cd "$TARG"

mkdir fastqs

#copy fastq files from permanent location
cp /data/kdi_prod/dataset_all/"$KDI_NO"/export/user/"$TARG"/"$TARG".R1.fastq.gz fastqs/"$TARG"_S1_L001_R1_001.fastq.gz

cp /data/kdi_prod/dataset_all/"$KDI_NO"/export/user/"$TARG"/"$TARG".R2.fastq.gz fastqs/"$TARG"_S1_L001_R2_001.fastq.gz

mkdir CR_output
cd CR_output

#run cellranger
cellranger count --id="CR_${TARG}" \
--fastqs=/data/tmp/Jason_targeted_CR/"$TARG"/fastqs \
--sample="$TARG" \
--transcriptome=/bioinfo/local/build/Centos/cellranger/refdata/refdata-cellranger-mm10-3.0.0 \
--localmem=45 \
--localcores=8


