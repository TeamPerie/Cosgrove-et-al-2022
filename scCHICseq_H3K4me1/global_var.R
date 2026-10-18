
###### LIBRARIES ######

# Packages
library(ChromSCape)          
library(devtools)            
library(DropletUtils)        
library(irlba)               
library(corrplot)           
library(R.utils)            
library(scater)             
library(Rtsne)               
library(ccRemover)           
library(viridis)             
library(colorRamps)          
library(RColorBrewer)       
library(edgeR)                
library(gplots)              
library(ggplot2)             
library(RColorBrewer)       
library(genefilter)          
library(xtable)              
library(WriteXLS)             
library(data.table)           
library(stringr)             
library(limma)                
library(edgeR)               
library(monocle3)           
library(Seurat)              
library(dendextend)          
library(ape)                
library(ConsensusClusterPlus) 
library(Matrix)              
library(genefilter)          
library(ggpubr)              
library(eulerr)              
library(GenomicRanges)      
library(kableExtra)           
library(colorspace)          
library(forcats)             
library(dplyr)
library(tidyr)
library(scTools)
library(EnsDb.Hsapiens.v86)
library(Signac)
library(decoupleR)
library(factoextra)
library(ggdendro)
library(msigdbr)
library(fgsea)
library(Nebulosa)
library(memes)
library(universalmotif)


###### VARIABLES ######

set.seed(1)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

###### FUNCTIONS ######

Associate_Annot = function(masterPeak, Annotations){
  Peak$peak_loca = paste0(Peak@seqnames, "_",Peak@ranges)
  genes.coords.all = read.table(file.path(Annotations, "EnsDB.All.bed.df"), header = T)
  genes.coords.all$seqnames = paste0("chr", genes.coords.all$seqnames)
  genes.coords.all = plyranges::as_granges(genes.coords.all)
  
  Peak = plyranges::join_overlap_inner(Peak, genes.coords.all)
  Peak = as.data.frame(Peak)
  Order_loca = c("Promoter_1K","Promoter_2K","Promoter_3K","Promoter_5K","TSS_1K","TSS_2K","TSS_5K","TSS_10K","Enhancers","Genes","Intergenic")
  Peak$seq_coord_system = factor(Peak$seq_coord_system, levels = Order_loca)
  Peak = Peak %>% dplyr::group_by(peak_loca) %>% dplyr::arrange(seq_coord_system, .by_group = T) %>% dplyr::slice_head(n = 1)
  Peak = plyranges::as_granges(Peak)
  
  genes.coords.all_onlygenes = genes.coords.all %>% dplyr::filter(seq_coord_system == "Genes")
  Intergenic_matdiff = Peak %>% dplyr::filter(seq_coord_system == "Intergenic") %>% plyranges::as_granges()
  Intergenic_matdiff$Nearest_gene = genes.coords.all_onlygenes$gene_name[nearest(Intergenic_matdiff, genes.coords.all_onlygenes)]
  Intergenic_matdiff = as.data.frame(Intergenic_matdiff)
  gene_MatDiff_fullannot = Peak %>% dplyr::filter(seq_coord_system != "Intergenic")
  gene_MatDiff_fullannot$Nearest_gene = gene_MatDiff_fullannot$gene_name
  gene_MatDiff_fullannot = as.data.frame(gene_MatDiff_fullannot)
  gene_MatDiff_fullannot = rbind(gene_MatDiff_fullannot, Intergenic_matdiff)
  gene_MatDiff_fullannot = plyranges::as_granges(gene_MatDiff_fullannot)
  return(gene_MatDiff_fullannot)
}

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


###### CHROMSCAPE FUNCTIONS ######
## Use macs3 instead of macs2 
call_macs2_merge_peaks <- function(affectation, odir, p.value,
                                   format = c("scBED","BAM")[1],
                                   ref, peak_distance_to_merge,min_peak_size, 
                                   conda_env = "Env_macs3"){
  suffix = ".bam"
  if(format != "BAM") suffix = ".bed"
  if(format != "BAM") format = "BED"
  merged_peaks=list()
  for (class in levels(factor(affectation$cell_cluster))){
    
    macs2_options = paste0(" -f ",format," --nomodel --extsize 300 ")
    
    command <- paste0(
      "conda run -n ", conda_env, " macs3 callpeak ", p.value,
      macs2_options, " --keep-dup all --broad -t '",
      file.path(odir, paste0(class, suffix)),
      "' --outdir '", odir, "' --name ", class
    )
    
    message("Running: ", command)
    system(command)
    merged_peaks[[class]] <-
      merge_MACS2_peaks(odir, class, peak_distance_to_merge,min_peak_size, ref)
  }
  unlink(file.path(odir, "*.xls"))
  unlink(file.path(odir, "*.gappedPeak"))
  unlink(file.path(odir, "*_model.r"))
  unlink(file.path(odir, "header.sam"))
  unlink(file.path(odir, 'merged.bam'))
  return(merged_peaks)
}

# Put merged instead of merged_bam in Rsamtools::mergeBam
subset_bam_call_peaks <- function(scExp, odir, input, format = "BAM", p.value = 0.005,
                                  ref = "hg38", peak_distance_to_merge = 10000,min_peak_size = 200,
                                  geneTSS_annotation = NULL, run_coverage = FALSE,
                                  progress = NULL)
{
  stopifnot(is(scExp, "SingleCellExperiment"), dir.exists(odir),
            is.numeric(p.value), is.character(ref),
            is.numeric(peak_distance_to_merge))
  if (!is.null(progress)) progress$set(detail = "Checking for annotation...", value = 0.15)
  if (!ref %in% c("hg38", "mm10", "ce11")) 
    stop("ChromSCape::subset_bam_call_peaks - 
            Reference genome (ref) must be 'hg38' or 'mm10'.")
  if (!"barcode" %in% colnames(SingleCellExperiment::colData(scExp))) 
    stop("ChromSCape::subset_bam_call_peaks - 
            Colnames of cell annotation must contain 'barcode'.")
  if (!"cell_cluster" %in% colnames(SingleCellExperiment::colData(scExp))) 
    stop("ChromSCape::subset_bam_call_peaks -
            Colnames of cell annotation must contain 'cell_cluster'.")
  
  if (is.null(geneTSS_annotation))
  {
    message("ChromSCape::subset_bam_call_peaks - Selecting ", ref,
            " genes from Gencode.")
    eval(parse(text = paste0("data(", ref, ".GeneTSS)")))
    geneTSS_annotation = 
      eval(parse(text = paste0("", ref, ".GeneTSS")))
    start = geneTSS_annotation$start
    geneTSS_annotation$start = ifelse(geneTSS_annotation$strand == "+",
                                      geneTSS_annotation$start,
                                      geneTSS_annotation$end)
    geneTSS_annotation$end = ifelse(geneTSS_annotation$strand == "+",
                                    start +1,
                                    geneTSS_annotation$end +1)
    geneTSS_annotation$strand = NULL
    geneTSS_annotation = as(geneTSS_annotation,"GRanges")
  } else geneTSS_annotation = as(geneTSS_annotation, "GRanges")
  if(format == "BAM"){
    if (!is.null(progress)) progress$set(detail = "Merging BAM files...", value = 0.25)
    if (length(input) > 1)
    {
      merged = file.path(odir, "merged.bam")
      message( "Merging BAM files...")
      Rsamtools::mergeBam(input, destination = merged,
                          overwrite = TRUE, indexDestination = TRUE)
      
    } else merged = input[1]
    if(!file.exists(paste0(merged,".bai"))) 
      Rsamtools::indexBam(merged, overwrite= TRUE)
  } 
  if (!is.null(progress)) progress$set(detail = "Generating pseudo-bulk of each cluster...", value = 0.35)
  affectation = SingleCellExperiment::colData(scExp)
  affectation$cell_cluster = as.factor(affectation$cell_cluster)
  message("Generating pseudo-bulk of each cluster...")
  p.value = paste0(" -p ", p.value, " ")  # format for system call to macs2
  
  if(format == "BAM") separate_BAM_into_clusters(affectation, odir, merged)
  else concatenate_scBed_into_clusters(affectation, input, odir)
  
  if (!is.null(progress)) progress$set(detail = "Calling peaks with MACS3...", value = 0.45)
  message("Calling peaks with MACS3...")
  merged_peaks <- call_macs2_merge_peaks(affectation, odir, p.value, 
                                         format, ref, peak_distance_to_merge, min_peak_size)
  if (!is.null(progress)) progress$set(detail = "Refining annotation based on merged peaks...", value = 0.65)
  message("Refining annotation based on merged peaks...")
  scExp@metadata$refined_annotation <- 
    annotation_from_merged_peaks(scExp, odir,
                                 merged_peaks, geneTSS_annotation)
  
  return(scExp)
}

# Not present in the functions of the package but inside ChromSCape/R/peak_calling.R
separate_BAM_into_clusters <- function(affectation, odir, merged_bam){
  
  for (class in levels(factor(affectation$cell_cluster)))
  {
    message("ChromSCape:::separate_BAM_into_clusters - generating ",class,
            " BAM file...")
    class_barcodes <- affectation$barcode[which(affectation$cell_cluster == 
                                                  as.character(class))]
    filt <- Rsamtools::ScanBamParam(tag=c("XB"),
                                    tagFilter=list(XB=class_barcodes))
    
    Rsamtools::filterBam(merged_bam,file.path(odir,paste0(class,".bam")),
                         indexDestination = TRUE, overwrite=TRUE, param=filt)
  }
}


# Not present in the functions of the package but inside ChromSCape/R/peak_calling.R
# Modified the Input
# Remove perged_peaks list initiation
merge_MACS2_peaks <- function(odir, class, peak_distance_to_merge = 1000, 
                              min_peak_size = 200, ref = "hg38"){
  # merged_peaks = list()
  suffix = "_peaks.broadpeak"
  peak_file = file.path(odir, paste0(class, suffix))
  peaks = read.table(
    file = peak_file,
    colClasses = c("character", "integer", "integer", rep("NULL", 6)))
  colnames(peaks) = c("chr", "start", "end")
  peaks = GenomicRanges::GRanges(peaks)
  peaks = peaks[which(width(GenomicRanges::ranges(peaks)) >= min_peak_size), ]
  peaks = GenomicRanges::reduce(peaks, min.gapwidth = peak_distance_to_merge, 
                                ignore.strand = TRUE)
  
  eval(parse(text = paste0("data(", ref, ".chromosomes)")))
  ref_chromosomes <- eval(parse(text = paste0("", ref, ".chromosomes")))
  ref_chromosomes = GenomicRanges::GRanges(ref_chromosomes)
  peaks = IRanges::subsetByOverlaps(peaks, ref_chromosomes, 
                                    ignore.strand = TRUE)
  return(peaks)
}


# Pour gérer les empty cells, avec le threshold à 100 on ne peut pas mettre une limite en dessous de 100
filter_scExp =  function (
    scExp, min_cov_cell = 1600, quant_removal = 95, min_count_per_feature = 10,
    verbose = TRUE){
  stopifnot(is(scExp, "SingleCellExperiment"), is.numeric(min_cov_cell),
            is.numeric(quant_removal),
            is.numeric(min_count_per_feature), verbose %in% c(FALSE, TRUE))
  if (is.null(scExp)) warning(
    "ChromSCape::filter_scExp - Please specify a  SingleCellExperiment")
  
  cellCounts <- Matrix::colSums(SingleCellExperiment::counts(scExp))
  thresh <- stats::quantile(cellCounts, probs = seq(0, 1, 0.01))
  sel <- (cellCounts > min_cov_cell & cellCounts <= thresh[quant_removal + 1])
  
  scExp <- scExp[, sel]
  counts <- SingleCellExperiment::counts(scExp)
  sel_feature <- (Matrix::rowSums(counts) >= min_count_per_feature)
  
  scExp <- scExp[sel_feature,]
  empty_cells = (Matrix::colSums(SingleCellExperiment::counts(scExp)) < min_cov_cell)
  if(any(empty_cells)) scExp <- scExp[,!empty_cells]
  empty_features = (Matrix::rowSums(counts(scExp)) < min_count_per_feature)
  if(any(empty_features)) scExp <- scExp[!empty_features,]
  
  SummarizedExperiment::colData(scExp)$total_counts = 
    Matrix::colSums(SingleCellExperiment::counts(scExp))
  
  if (verbose) message(
    "ChromSCape::filter_scExp - ", ncol(scExp), " cells pass the ",
    " threshold of ", min_cov_cell, " minimum reads and are lower than ",
    "the ", quant_removal, "th centile of library size ~= ",
    round(thresh[quant_removal + 1]), " reads.")
  
  if (verbose) message(
    "ChromSCape::filter_scExp - ",nrow(scExp)," features pass ",
    "the threshold of ", min_count_per_feature," count per feature.")
  
  return(scExp)
}

reduce_dims_scExp <-
  function(scExp, dimension_reductions = c("PCA", "UMAP"), n = 10,
           batch_correction = FALSE, batch_list = NULL,
           remove_PC = NULL, verbose = TRUE){
    stopifnot(is(scExp, "SingleCellExperiment"), is.numeric(n),
              dimension_reductions[1] %in% c("PCA", "TSNE", "UMAP"))
    
    if(!is.null(remove_PC)) stopifnot(is.character(remove_PC), 
                                      all(remove_PC %in% paste0("Component_",1:100)))
    if (!"normcounts" %in% names(SummarizedExperiment::assays(scExp))){
      warning("ChromSCape::reduce_dims_scExp - The raw counts are not
                normalized, running dimensionality reduction on raw counts.")
      mat <- SingleCellExperiment::counts(scExp)
    } else{
      mat <- SingleCellExperiment::normcounts(scExp)
    }
    if (batch_correction && !is.list(batch_list)){
      stop("ChromSCape::reduce_dims_scExp - If doing batch correction,
                batch_list must be a list of the samples IDs of each batch.")}
    if (batch_correction){
      out <- reduce_dim_batch_correction(scExp, mat, batch_list, n)
      scExp <- out$scExp
      pca <- out$pca
    } else{
      scExp$batch_id <- factor(1)
      if ( is(mat, "dgCMatrix") | is(mat, "dgTMatrix")) {
        pca <- pca_irlba_for_sparseMatrix(Matrix::t(mat), n)
      } else{
        pca <- stats::prcomp(Matrix::t(mat),center = TRUE,
                             scale. = FALSE)
        pca <- pca$x[, seq_len(n)]}}
    pca <- as.data.frame(as.matrix(pca))
    colnames(pca) <- paste0("Component_", seq_len(ncol(pca)))
    if(!is.null(remove_PC)) {
      pca <- pca[,-which(colnames(pca) %in% remove_PC)]
      if(verbose) message("ChromSCape::reduce_dims_scExp - removing ",
                          paste(remove_PC, collapse = ", "), " as probably correlated to library size.")
    }
    rownames(pca) <- colnames(scExp)
    if ("TSNE" %in% dimension_reductions){
      tsne <- Rtsne::Rtsne(pca, dims = 2, pca = FALSE, theta = 0,
                           perplexity = choose_perplexity(pca), verbose = verbose,
                           max_iter = 1000)
      tsne <- as.data.frame(tsne$Y)
      rownames(tsne) <- rownames(pca)
      colnames(tsne) <- c("Component_1", "Component_2")
    }
    if ("UMAP" %in% dimension_reductions){
      config <- umap::umap.defaults
      # config$metric <- "cosine"
      umap <- umap::umap(pca, config = config, method = "naive")
      umap <- as.data.frame(umap$layout)
      colnames(umap) <- c("Component_1", "Component_2")
    }
    listReducedDim <- list(PCA = pca)
    if ("TSNE" %in% dimension_reductions) listReducedDim$TSNE <- tsne
    if ("UMAP" %in% dimension_reductions) listReducedDim$UMAP <- umap
    SingleCellExperiment::reducedDims(scExp) <- listReducedDim
    return(scExp)
  }


pca_irlba_for_sparseMatrix <- function(x, n_comp, work = 3 * n_comp)
{
  x.means <- Matrix::colMeans(x)
  svd.0 <- irlba::irlba(x, center = x.means, nv = n_comp, work = work)
  pca <- svd.0$u %*% diag(svd.0$d)
  return(pca)
}



choose_perplexity <- function(dataset)
{
  stopifnot(!is.null(dataset), !is.null(dim(dataset)))
  perplexity <- 30
  if (nrow(dataset) <= 200)
  {
    perplexity <- 20
  }
  if (nrow(dataset) <= 250)
  {
    perplexity <- 25
  }
  if (nrow(dataset) <= 150)
  {
    perplexity <- 15
  }
  if (nrow(dataset) <= 100)
  {
    perplexity <- 10
  }
  if (nrow(dataset) <= 50)
  {
    perplexity <- 5
  }
  perplexity
}


define_feature <- function(ref = c("hg38","mm10", "ce11")[1], peak_file = NULL,
                           bin_width  = NULL, genebody = FALSE,
                           extendPromoter = 2500){
  eval(parse(text = paste0("data(", ref, ".chromosomes)")))
  chr <- eval(parse(text = paste0("", ref, ".chromosomes")))
  chr <- GenomicRanges::GRanges(chr)
  if (!is.null(peak_file)) {
    message('ChromSCape::define_feature - ',
            'Reading in peaks file...')
    features <- rtracklayer::import(peak_file, format="bed")
    which <- IRanges::subsetByOverlaps(features,chr)
    if(length(which(GenomicRanges::duplicated(which)))>0) 
      message("ChromSCape::define_feature - Removing ",
              length(which(GenomicRanges::duplicated(which))), 
              " duplicated regions from peak file.")
    which = which[!GenomicRanges::duplicated(which)]
  } else if (!is.null(bin_width)) {
    message('ChromSCape::define_feature - ',
            'Counting on genomic bins...')
    which <- GenomicRanges::tileGenome(
      setNames(
        GenomicRanges::width(chr),
        GenomicRanges::seqnames(chr)
      ),
      tilewidth = bin_width, cut.last.tile.in.chrom = T
    )
  } else if (genebody == TRUE) {
    # Retrieve gene TSS from ref and create GRanges
    eval(parse(text = paste0("data(", ref, ".GeneTSS)")))
    geneTSS_df <- eval(parse(text = paste0("", ref, ".GeneTSS")))
    
    geneTSS_df$start = ifelse(geneTSS_df$strand == "+",
                              geneTSS_df$start - extendPromoter,
                              geneTSS_df$start)
    
    geneTSS_df$end = ifelse(geneTSS_df$strand == "-",
                            geneTSS_df$end + extendPromoter,
                            geneTSS_df$end)
    
    geneTSS_df$strand = NULL
    geneTSS_df$distanceToTSS = 0
    which = GenomicRanges::GRanges(geneTSS_df)
    which = which[!GenomicRanges::duplicated(which)]
    names(which@ranges) = NULL
  }
  message("Using custom function adapted from ChromSCape")
  return(which)
}
