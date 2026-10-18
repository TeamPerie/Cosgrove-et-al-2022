# Cosgrove et al. – The pentose phosphate pathway modulates myelo-lymphoid lineage specification

Analysis code for:

> Cosgrove J., Donada A., Lyne A.-M., *et al.* **The pentose phosphate pathway modulates myelo-lymphoid lineage specification.** *Nature Communications* (2026). DOI: [to be added on publication]

An earlier version of this work was posted as a preprint under the title *"Metabolic priming of multipotent hematopoietic progenitors fuels innate immunity"* (bioRxiv, 2023). This repository was created for that preprint and has been updated for the published article.

> **Figure numbering.** The top-level `Figure*/` folders were named after the **preprint** figures and do **not** match the published article. Use the [figure-to-code table](#where-to-find-the-code-for-each-figure) below to find the code for each panel of the published paper. New analyses added during revision live in descriptively named folders (`scCHICseq_H3K4me1/`, `inflammation_reanalysis/`, `PAGA_Dahlin2018/`, `Figure1/DRAG_barcode_QC/`).

## Repository layout

| Folder | Analysis | Main scripts |
|---|---|---|
| `Figure1/DNA_BC_PREPROCESSING/` | DRAG barcode extraction, clean-up and QC from bulk DNA sequencing of mature myeloid/erythroid cells | `run_analysis_bcmouse.sh`, `bcmouse*.sh`, `tagcleanup.R`, `splitVDJ.R`, `DNA_BC_QC.Rmd` |
| `Figure1/RNA_BC_PREPROCESSING/` | Recovery of DRAG barcodes from 10x scRNA-seq libraries (targeted amplification), consensus calling, IGoR generation probabilities | `run_10x_analysis_vdj.sh`, `consensus_per_UMI_and_cell_reads.R`, `splitVDJ.R`, `infer_igor_tags.sh` (see folder README) |
| `Figure1/INTEGRATION/` | MetaFate: scRNA-seq QC and integration, lineage-bias scoring, DEGs, pathway enrichment, signatures, permutation tests, 4-fold cross-validation, validation on public datasets | `IntegrateDataset_QCandPreprocessing.Rmd`, `IntegrateDataset_Analysis.Rmd`, `IntegrateDataset_Validation.Rmd`, `helper_methods_refactored.R` |
| `Figure1/DRAG_barcode_QC/` *(new)* | Barcode length / insertion / deletion summaries; RNA barcode recovery summary | `splitVDJ_analysis.R`, `summarise-RNA-BC-recovery.R` |
| `Figure5/` (identical copies in `Figure6/`, `FigureS11_S12_S13/`) | Lentiviral barcoding of CD62L<sup>neg/high</sup> MPPs (3 weeks post-transplant): demultiplexing, QC, clone size and lineage bias | `Demultiplex/`, `QC/lentiviralbarcoding_QC_WT.Rmd`, `Analysis/CD62L_MPP_analysis.Rmd` |
| `Figure3/` (identical copies in `Figure4/`, `FigureS10/`, `FigureS9/`) | Lentiviral barcoding of WT vs G6PD-Tg MPPs; metabolomics | `QC/lentiviralbarcoding_QC_{WT,TG}.Rmd`, `Analysis/Analysis.Rmd`, `metabolomics_anlaysis.R` |
| `FigureS14_and_S15/` | Lentiviral barcoding of SLAM HSCs, 12 months post-transplant | `Demultiplex/newXcalibr.bash`, `QC/step2_QCandPreprocessing.Rmd`, `Analysis/step3_Analysis.Rmd` |
| `FigureS6/` | Metabolic pathway signature dynamics along haematopoietic trajectories (Dahlin et al. 2018) | `S6.Rmd` |
| `FigureS7/` | Haemopedia bulk RNA-seq: PCA, lineage DEGs, MIIC input and network heatmap | `FigureS7.Rmd`, `MIIC_heatmap.Rmd`, `S7_helpermethods.R` |
| `PAGA_Dahlin2018/` *(new)* | PAGA trajectory model of the Dahlin et al. 2018 cKit<sup>+</sup> dataset | `process_paga_dataset_dahlin.py`, `dahlin18.ipynb` |
| `inflammation_reanalysis/` *(new)* | Re-analysis of public scRNA-seq datasets: G-CSF (Fast et al. 2022), EPO (Tusi et al. 2018), PyMT breast cancer (Gerber-Ferder et al. 2023), *P. berghei* infection (Haltalli et al. 2020); HSC/MPP signature scoring | one sub-folder per dataset; `HSPC_signatures/calculate_HSPC_sigs.Rmd` |
| `scCHICseq_H3K4me1/` *(new)* | Single-cell ChIC-seq (H3K4me1) of CD62L<sup>neg/high</sup> MPPs after 6-AN treatment and of WT vs G6PD-Tg MPPs: QC, clustering, differential enrichment, motif analysis | `ChIC_Analysis_1.Rmd`, `global_var.R` |
| repository root | MetaFate per-cell metadata and normalised expression matrix | `S5a_metafate_seurat_metadata.csv`, `S5b_normalised_expression_matrix.zip` |

## Where to find the code for each figure

"Flow / Prism" means the panel was quantified in FlowJo and plotted and tested in GraphPad Prism, with no custom code. The numbers behind every plot are in the article's Source Data file.

| Published figure | Panels | Code |
|---|---|---|
| **Fig. 1** | b, c | `Figure1/DNA_BC_PREPROCESSING/DNA_BC_QC.Rmd`; `Figure1/INTEGRATION/IntegrateDataset_Analysis.Rmd` |
| | d | `Figure1/RNA_BC_PREPROCESSING/`; `Figure1/INTEGRATION/IntegrateDataset_QCandPreprocessing.Rmd`, `IntegrateDataset_Analysis.Rmd` (UMAP, density maps) |
| | e, f, g | `Figure1/INTEGRATION/IntegrateDataset_Analysis.Rmd` (DEGs/volcano, enrichR KEGG, signature permutation tests) |
| **Fig. 2** | a | MetaFate Seurat object from `Figure1/INTEGRATION/`; plotted in Prism |
| | b, g | Flow / Prism |
| | d, e | `Figure5/Analysis/CD62L_MPP_analysis.Rmd` (QC: `Figure5/QC/`) |
| | h, i | `inflammation_reanalysis/GCSF_Fast2022/`, `EPO_Tusi2018/`, `PyMT_GerberFerder2023/` |
| **Fig. 3** | b–e | Flow (SCENITH, SPICE-Met, TMRE) / Prism |
| | h–j | `Figure3/metabolomics_anlaysis.R` |
| **Fig. 4** | a (left) | `PAGA_Dahlin2018/` |
| | a (right) | `FigureS6/S6.Rmd` |
| | b | `FigureS7/FigureS7.Rmd` |
| | c, d | `FigureS7/MIIC_heatmap.Rmd`, `S7_helpermethods.R` (input matrix = Supplementary Table 10; network inferred with MIIC, https://miic.curie.fr) |
| **Fig. 5** | a | Flow / Prism |
| | c–e, h | `scCHICseq_H3K4me1/ChIC_Analysis_1.Rmd` |
| **Fig. 6** | b | Flow / Prism |
| | c–e | `Figure3/Analysis/Analysis.Rmd` (QC: `Figure3/QC/`) |
| **Supplementary Fig. 1** | b–g | `Figure1/DNA_BC_PREPROCESSING/splitVDJ.R`, `Figure1/RNA_BC_PREPROCESSING/splitVDJ.R`, `infer_igor_tags.sh` (generation probabilities), `Figure1/DRAG_barcode_QC/splitVDJ_analysis.R` |
| **Supplementary Fig. 2** | – | Flow (gating) |
| **Supplementary Fig. 3** | a | `Figure1/DNA_BC_PREPROCESSING/DNA_BC_QC.Rmd` |
| | b | `Figure1/INTEGRATION/IntegrateDataset_Analysis.Rmd` |
| **Supplementary Fig. 4** | a–d | `Figure1/INTEGRATION/IntegrateDataset_QCandPreprocessing.Rmd`, `IntegrateDataset_Analysis.Rmd` |
| | e–g | `IntegrateDataset_Analysis.Rmd` (bias-threshold sensitivity analysis) |
| | h | `IntegrateDataset_Analysis.Rmd` (4-fold cross-validation) |
| **Supplementary Fig. 5** | a, b | `Figure1/INTEGRATION/IntegrateDataset_Analysis.Rmd` (a); `IntegrateDataset_Validation.Rmd` (b: Dahlin, Tusi, Weinreb) |
| **Supplementary Fig. 6** | a–c | Flow / Prism |
| **Supplementary Fig. 7** | b, c | `Figure5/QC/lentiviralbarcoding_QC_WT.Rmd` |
| | d–f | `Figure5/Analysis/CD62L_MPP_analysis.Rmd` |
| **Supplementary Fig. 8** | c–f | `Figure5/Analysis/CD62L_MPP_analysis.Rmd` (a, b: flow) |
| **Supplementary Figs. 9, 10** | all | `FigureS14_and_S15/` |
| **Supplementary Fig. 11** | a, b | Flow |
| **Supplementary Fig. 12** | a–d | `Figure1/INTEGRATION/IntegrateDataset_Validation.Rmd` (Haltalli et al. section); `inflammation_reanalysis/Plasmodium_Haltalli2020/` |
| **Supplementary Fig. 13** | a–d | `inflammation_reanalysis/HSPC_signatures/calculate_HSPC_sigs.Rmd` (with the dataset-specific scripts in `inflammation_reanalysis/`) |
| **Supplementary Figs. 14, 15, 17, 21** | all | Flow / Prism |
| **Supplementary Fig. 16** | a–d | `Figure3/metabolomics_anlaysis.R` |
| **Supplementary Figs. 18, 23** | all except b (FRiP) | `scCHICseq_H3K4me1/ChIC_Analysis_1.Rmd` |
| **Supplementary Fig. 19** | b–e | `Figure3/QC/lentiviralbarcoding_QC_{WT,TG}.Rmd` |
| **Supplementary Fig. 20** | all | `Figure3/Analysis/Analysis.Rmd` |
| **Supplementary Fig. 22** | a, b | `Figure3/Analysis/Analysis.Rmd` |
| | c–e | Flow / Prism |

Fig. 7 is a schematic.

## Data

| Data | Location |
|---|---|
| MetaFate scRNA-seq and DRAG barcode sequencing | GEO / SRA: [accessions to be added] |
| Lentiviral barcode sequencing | SRA: [accession to be added]; processed count matrices are Supplementary Tables 8, 9, 11 and 12 of the article |
| scCHIC-seq | GEO: [accession to be added] |
| Metabolomics (LC-MS/MS) | [repository / accession to be added]; per-sample peak heights in the article's Source Data |
| MetaFate per-cell metadata and normalised expression | this repository (`S5a_*`, `S5b_*`) |
| Public datasets re-analysed | Dahlin et al. 2018; Tusi et al. 2018; Weinreb et al. 2020; Fast et al. 2022; Haltalli et al. 2020; Gerber-Ferder et al. 2023; Haemopedia (Choi et al. 2019). See the article's Data Availability statement for accession codes. |

## Software

The analyses were run in R (≥ 4.0) and Python 3. Key packages:

- **scRNA-seq**: Seurat, SeuratDisk, scran, scater, Rmagic, destiny, monocle3, clustree, enrichR, msigdbr, fgsea, decoupleR
- **Bulk RNA-seq / networks**: limma, edgeR, DESeq2, Glimma, genefilter; MIIC (web server)
- **Barcoding**: stringdist, seqinr, Biostrings; xcalibr and NCBI BLAST 2 (bundled under `*/Demultiplex/xcalibr-master` and `Figure1/DNA_BC_PREPROCESSING/blast2`); IGoR
- **scCHIC-seq**: ChromSCape, Signac, chromVAR, motifmatchr, TFBSTools, JASPAR2020, GenomicRanges, BSgenome.Mmusculus.UCSC.mm10, DropletUtils, irlba
- **Plotting and utilities**: ggplot2, ggpubr, ggrepel, cowplot, gridExtra, pheatmap, gplots, plotgardener, dplyr, tidyr, plyr, reshape2, data.table
- **Python**: scanpy, anndata, numpy, pandas, scipy, matplotlib (PAGA analysis)

## Reproducing the analyses

- Scripts were run interactively on the authors' machines and the Institut Curie cluster. Most contain **absolute paths** (`setwd()`, `read.csv("/Users/...")`). Change these to your local copies of the input data before running.
- Inputs that are not in this repository (raw FASTQ, Cell Ranger outputs, Seurat objects of public datasets, gene-set `.Rda` files) can be obtained from the repositories listed above or regenerated from the raw data.
- `Figure4/`, `Figure6/`, `FigureS9/`, `FigureS10/` and `FigureS11_S12_S13/` contain byte-identical copies of scripts in `Figure3/` and `Figure5/`. They are kept so that links to the preprint version remain valid.

## Credits

- **DRAG barcode processing**: adapted from pipelines by Arno Velds and Anne-Marie Lyne.
- **scCHIC-seq analysis**: Grégoire Jouault.
- **All other analyses**: Jason Cosgrove (jason.cosgrove@curie.fr), Team Perié, Institut Curie.

## License

MIT – see `LICENSE`.
