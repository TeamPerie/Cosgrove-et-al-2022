# DRAG barcode extraction from 10x scRNA-seq data and targeted sequencing
--------------

This is the code to extract and analyse VDJ DRAG barcodes from 10x single cell RNA-seq data as used in the paper:

Metabolic Heterogeneity in Hematopoietic Multipotent Progenitors Fuels Innate Immune Cell Production (Cosgrove et al., 2022).

Full details of the approach are given in the methods section of the paper.


### Running the code

Originally the code was run on a cluster using the main script: run_10x_analysis_vdj.sh

If the code is to be run this way, there are several variables to be set at the beginning of the script (this is also annotated in the script):

* OUT_DIR – writable directory in which all output will go
* TARG_TMP - temporary directory for cellranger output for targeted sequencing
* SCRIPT_DIR – directory where all the bash and R scripts are stored on the cluster
* COUNT_DIRS – list of full paths to count directories of cellranger output for 10x experiments, separated by ;
* SAMPLES – sample names, separated by ; (should match those used for cellranger)
* BC_LISTS – list of full paths to barcode lists for each 10x experiment, separated by ;
* Other parameters related to the output from our sequencing platform
* NTARGS – number of samples
* DATA_DIR – directory where separate 10x bam files are stored on cluster (each in a subdirectory SAMPLE_NAME/possorted_genome_bam.bam)
* R_EXEC – path to Rscript
* Parameters which can be changed, such as the number of bases to match in the flanking regions, which are currently set to the values used in the paper.
* Parameters related to barcode filtering


At least part of the code probably needs to be run on a cluster as running cellranger for multiple samples takes some time. On the other hand, some parts of the code can also be run locally e.g. the R code to filter barcodes, and then the input directories and paths can be changed appropriately.


### Requirements

* Samtools
* R (libraries seqinr, Biostrings)
* [IGoR](https://github.com/qmarcou/IGoR)


### How does it work?

Experimentally, a standard 10x experiment is carried out on cells which are expressing DRAG barcodes, and then cellranger is run to obtain count matrices and bam files. In addition, for each sample, targeted sequencing is used to amplify reads containing the DRAG barcodes and then these reads are sequenced such that the 10x cell and UMI barcodes are retained.

The code

* Runs cellranger on the targeted sequencing so as to get a bam file with aligned (and unaligned) reads with CB and UB tags for the 10x cell and UMI barcodes (run_CR_count_targeted_seq.sh).

* Extract DRAG barcodes from bam files using script grep_flanking_both_10x_vdj2_plus_targeted.sh. First extract unaligned reads from both bams using samtools, then use a grep match to the flanking regions to extract DRAG barcodes from reads, as well as 10x CB and UB. Keep only DRAG barcodes in cells listed in cellranger barcode file.

* Filter DRAG barcodes to get one per cell using consensus_per_UMI_and_cell_reads.R. Briefly, this is done by filtering UMIs with < MIN_READ reads, and then those with < MIN_PROP reads assigned to top barcode. The barcode with most reads is then assigned to each UMI. Cells are then filtered to keep only those with MIN_UMI umis (although we use MIN_UMI=1 in the paper). For cells with more than one UMI, we then filter to keep cells with at least MIN_PROP_UMI umis with dominant barcode i.e. cells with good consensus across UMIs. Finally, for cells with only one UMI remaining, we keep those where the final barcode is observed in at least MIN_PROP_1UMI of the original reads for this cell.

* Filter DRAG barcodes according to whether they have expected VDJ structure using code splitVDJ.R.

* Collate all barcodes and run IGoR analysis using script infer_igor_tags.sh. Fit statistical model of VDJ barcode generation and compute generation probability of each barcode.



### Input and example output

10x scRNA-seq data from Cosgrove et al can be found [here](******LINK*****) and targeted sequencing data [here](******LINK*****). Input and output for the barcode extraction and filtering parts of the code are provided for sample m534 (files have been zipped to save space but input to code should be unzipped). 


