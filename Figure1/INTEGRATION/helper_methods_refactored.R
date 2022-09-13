##############################################################
# Helper methods used in Metafate analysis
#
# @Author: Jason Cosgrove (jason.cosgrove@curie.fr)
# @Date: 30/08/2022
###############################################################


# is X is in the vector Y
# inputs : "X" is an object 
#          "Y" is a vector
# outputs: A boolean stating TRUE if X is in Y
'%ni%' <- function(x,y)!('%in%'(x,y))


#calculate the average expression score of a list of genes
geneSignatureScoring <- function(sobj,  gene.sets,gene.set.names,assay){
  
  for(i in 1:length(gene.set.names)){
    genes.of.interest <- paste(gene.sets[,gene.set.names[i]])
    genes.of.interest <- genes.of.interest[genes.of.interest != ""]
    genes.of.interest <- intersect(genes.of.interest,rownames(sobj@assays$RNA))
    sobj <- AddModuleScore(sobj , features = list(genes.of.interest),
                           name = gene.set.names[i],replace = TRUE,assay = assay)
  }
  
  return(sobj)
  
}




removeUninformativeGenes <- function(dataset){
  genes.to.remove <- c(rownames(dataset)[grepl("Gm",rownames(dataset))],
                       rownames(dataset)[grepl("Rik",rownames(dataset))],
                       rownames(dataset)[grepl("Rp",rownames(dataset))])
  
  
  genes  <- rownames(dataset)[rownames(dataset) %ni% genes.to.remove]
  
  dataset <- subset(dataset, features = genes)
  return(dataset)
}


#load in the metadata for the experiment
makeMetaData <- function(mat,expt,bc.metadata){
  
  bc.metadata.subset <- bc.metadata[bc.metadata$expt == expt,]
  
  metadata <- data.frame(matrix(0,nrow = ncol(mat), ncol = 7))
  colnames(metadata) <- colnames(bc.metadata)[c(2,9:14)]
  rownames(metadata) <- colnames(mat)
  
  for(i in 1:nrow(bc.metadata.subset)){
    cell <- bc.metadata.subset$Cell.bc[i]
    metadata[cell,] <- bc.metadata.subset[i,c(2,9:14)]
  }
  
  #convert NAs to zeros
  metadata[is.na(metadata)] = 0
  
  metadata$Con.seq[metadata$Con.seq == 0] = "nb"
  
  return(metadata)
  
}



updatePGEN <- function(bc.metadata){
  
  #get the pgen for each barcode
  pgen <- read.csv('/Users/jasoncosgrove/Dropbox (Team_Perie)/Jason/Experiments/Dry_Lab/MetaFate_final_pipeline/data/Step4_RNA_BC_processing/IGOR/Jason_igor/jason_bcs_pgens.csv')
  
  bc.metadata$pgen <- NA
  
  for(i in 1:nrow(bc.metadata)){
    
    bc <- bc.metadata$Con.seq[i]
    
    if(nrow(pgen[pgen$VDJ == bc,]) > 0){
      bc.metadata$pgen[i] <- pgen[pgen$VDJ == bc,]$igor_pgen
    }
  }
  
  return(bc.metadata)
}



printBarcodeSummaryStatistics <- function(bc.metadata){
  #number of cells with barcodes
  print("number of cells with barcodes: ")
  print(nrow(bc.metadata))
  
  #unique number of barcodes
  print("number of unique barcodes: ")
  print(length(unique(bc.metadata$Con.seq_perexpt)))
  
  #distribution of clone sizes
  hist(table(bc.metadata$Con.seq_perexpt), main = "distribution of clone sizes for RNA barcodes")
  
  #number of cells with a DNA-RNA pairing
  print("number of cells with a DNA-RNA barcode pairing: ")
  bc.metadata[is.na(bc.metadata)] = 0
  table(bc.metadata$Myeloid > 0 | bc.metadata$Erythroid > 0)
  
}



CreateSeuratObjects <- function(m1,m2,m3,m4,m5,bc.metadata){
  m1.sobj <- CreateSeuratObject(counts = m1, min.cells = 50, min.features = 300,project = "m1",meta.data = makeMetaData(m1,"ET27a_M1",bc.metadata))
  
  m2.sobj <- CreateSeuratObject(counts = m2, min.cells = 50, min.features = 300,project = "m2",meta.data = makeMetaData(m2,"ET27b_M1",bc.metadata))
  
  m3.sobj <- CreateSeuratObject(counts = m3, min.cells = 50, min.features = 300,project = "m3",meta.data = makeMetaData(m3,"JCW26_M2",bc.metadata))
  
  m4.sobj <- CreateSeuratObject(counts = m4, min.cells = 50, min.features = 300,project = "m4",meta.data = makeMetaData(m4,"JCW26_M3",bc.metadata))
  
  m5.sobj <- CreateSeuratObject(counts = m5, min.cells = 50, min.features = 300,project = "m5",meta.data = makeMetaData(m5,"JCW26_M4",bc.metadata))
  
  dataset.list <- c(m1.sobj, m2.sobj, m3.sobj, m4.sobj, m5.sobj)
  return(dataset.list)
  
}

#add the percentage of reads that map to mitochondrial reads, this will be used as a QC metric
update_QC_metadata <- function(dataset.list){
  for (i in 1:length(dataset.list)) {
    sobj <- dataset.list[[i]]
    
    #update mitochondrial content as a QC metric
    sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^mt-")
    #also update whether the cell is barcoded which can help when assessing QC filters
    sobj@meta.data$barcoded <- (sobj@meta.data$Con.seq != "nb")
    dataset.list[[i]] <- sobj
  }
  return(dataset.list)
}


cellCycleAnalysis <- function(dataset.integrated){
  
  mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
  ensembl <- mapIds(org.Mm.eg.db, keys=rownames(dataset.integrated), keytype="SYMBOL", column="ENSEMBL")
  assignments <- cyclone(dataset.integrated@assays$integrated@data, mm.pairs, gene.names=ensembl)
  
  dataset.integrated@meta.data$phases <- assignments$phases
  dataset.integrated@meta.data$G1_score <- assignments$normalized.scores$G1
  dataset.integrated@meta.data$S_score <- assignments$normalized.scores$S
  dataset.integrated@meta.data$G2M_score <- assignments$normalized.scores$G2M
  
  return(dataset.integrated)
  
}



#this should be done mouse by mouse

setLineageBias <- function(diff.active,E_threshold = 0.25,M_threshold = 0.75){
  
  
  bias <- rep("unbiased", ncol(diff.active))
  bias[diff.active$bias.score >= quantile(diff.active$bias.score,M_threshold)[[1]]] <- "M"
  bias[diff.active$bias.score <= quantile(diff.active$bias.score,E_threshold)[[1]]] <- "E"
  diff.active@meta.data$bias <- bias
  diff.active@meta.data$bias.score <- diff.active$bias.score
  
  return(diff.active)
  
}



setDiffBias <- function(bc){
  
  bc@meta.data$bias <- rep("diff_inactive", length(all.barcoded.cells))
  
  
  for(i in 1:length(all.barcoded.cells)){
    
    cell <- all.barcoded.cells[i]
    #if the cell is differentiation active, update the bias variable using information from the 
    # diff.active seurat object
    active <- sum(bc@meta.data[cell,c("Myeloid","Erythroid")]) > 0
    if(active){
      bc@meta.data[cell,]$bias <- diff.active@meta.data[cell,]$bias
    }
  }
  
  return(bc)
}


createDensityPlot <- function(sobj,bins = 10){
  
  tmp.all<-as.data.frame(Embeddings(object = sobj, reduction = "umap"))
  p <- ggplot(tmp.all, aes(x = UMAP_1, y = UMAP_2)) + geom_point(colour="#00000000") + 
    stat_density_2d(aes(fill = stat(level)), geom = "polygon", bins=bins) + 
    scale_fill_gradientn(colors = c("#4169E100","royalblue", "darkolivegreen3","goldenrod1","red")) +
    theme_classic() + 
    theme(legend.position="bottom")
  
  return(p)
}

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

number_of_DEGs_per_threshold <- function(fp){
  
  threshold_0_55 <- read.csv(paste(fp,'cloneDEGs_bias_threshold_0_55.csv',sep = ""))
  
  threshold_0_66 <- read.csv(paste(fp,'cloneDEGs_bias_threshold_0_66.csv',sep = ""))
  
  threshold_0_75 <- read.csv(paste(fp,'cloneDEGs_bias_threshold_0_75.csv',sep = ""))
  
  threshold_0_9 <- read.csv(paste(fp,'cloneDEGs_bias_threshold_0_9.csv',sep = ""))
  
  #how many genes in each case
  print('55% threshold')
  print(nrow(threshold_0_55))
  
  print('66% threshold')
  print(nrow(threshold_0_66))
  
  print('75% threshold')
  print(nrow(threshold_0_75))
  
  print('90% threshold')
  print(nrow(threshold_0_9))
  
  barplot(c(nrow(threshold_0_55),
            nrow(threshold_0_66),
            nrow(threshold_0_75),
            nrow(threshold_0_9)),col =c("grey10","grey30","grey50","grey70"),
          main = "number of DEGs per threshold",
          xlab = "threshold value",ylab = "number of DEGs",
          names = c("55%","66%","75%","90%"))
  
}


DEGs_overlap <- function(fp){
  
  threshold_0_55 <- read.csv(paste(fp,'cloneDEGs_bias_threshold_0_55.csv',sep = ""))
  
  threshold_0_66 <- read.csv(paste(fp,'cloneDEGs_bias_threshold_0_66.csv',sep = ""))
  
  threshold_0_75 <- read.csv(paste(fp,'cloneDEGs_bias_threshold_0_75.csv',sep = ""))
  
  threshold_0_9 <- read.csv(paste(fp,'cloneDEGs_bias_threshold_0_9.csv',sep = ""))
  
  #how many genes in each case
  #overlap between the different thresholds
  
  print('DEG overlap 55% and 66%')
  print(jaccard(threshold_0_55$X,threshold_0_66$X))
  
  print('DEG overlap 66% and 75%')
  print(jaccard(threshold_0_66$X,threshold_0_75$X))
  
  print('DEG overlap 75% and 90%')
  print(jaccard(threshold_0_75$X,threshold_0_9$X))
  
  
}



log2FC_per_threshold <- function(fp){
  
  threshold_0_55 <- read.csv(paste(fp,'cloneDEGs_bias_threshold_0_55.csv',sep = ""))
  
  threshold_0_66 <- read.csv(paste(fp,'cloneDEGs_bias_threshold_0_66.csv',sep = ""))
  
  threshold_0_75 <- read.csv(paste(fp,'cloneDEGs_bias_threshold_0_75.csv',sep = ""))
  
  threshold_0_9 <- read.csv(paste(fp,'cloneDEGs_bias_threshold_0_9.csv',sep = ""))
  #plot the logfc in expression markers take that for our top markers
  threshold_0_55$avg_log2FC <-  rowMeans(threshold_0_55[,grepl('log2FC',colnames(threshold_0_55))])
  
  threshold_0_66$avg_log2FC <-  rowMeans(threshold_0_66[,grepl('log2FC',colnames(threshold_0_66))])
  
  threshold_0_75$avg_log2FC <-  rowMeans(threshold_0_75[,grepl('log2FC',colnames(threshold_0_75))])
  
  threshold_0_9$avg_log2FC <-  rowMeans(threshold_0_9[,grepl('log2FC',colnames(threshold_0_9))])
  
  
  threshold_0_55$threshold <- "55%"
  threshold_0_66$threshold <- "66%"
  threshold_0_75$threshold <- "75%"
  threshold_0_9$threshold <- "90%"
  
  combined <- rbind(threshold_0_55[,c("avg_log2FC","threshold")], threshold_0_66[,c("avg_log2FC","threshold")], 
                    threshold_0_75[,c("avg_log2FC","threshold")], threshold_0_9[,c("avg_log2FC","threshold")])
  
  
  p <- ggplot(combined, aes(x=threshold, y=avg_log2FC,color=threshold)) + 
    geom_boxplot(lwd = 1.5,show.legend = FALSE)
  p + coord_flip() + scale_color_manual(values=c("grey10","grey30","grey50","grey70")) + theme_classic() +theme(axis.text=element_text(size=12,face = 'bold'),
                                                                                                                axis.title=element_text(size=14,face="bold"),
                                                                                                                axis.line = element_line(size = 1.0, colour = "black")) + labs(y = "distribution of log2FC values for DEGs" , x= "bias classifier threshold value")
  
}





number_of_cells_per_threshold <- function(diff.active){
  
  bias <- rep("unbiased", ncol(diff.active))
  bias.score <- diff.active$bias.score
  bias[bias.score >= quantile(bias.score,0.55)[[1]]] <- "M"
  bias[bias.score <= quantile(bias.score,0.45)[[1]]] <- "E"
  a <- table(bias)
  
  bias <- rep("unbiased", length(diff.active.barcoded.cells))
  bias[bias.score >= quantile(bias.score,0.66)[[1]]] <- "M"
  bias[bias.score <= quantile(bias.score,0.33)[[1]]] <- "E"
  b <- table(bias)
  
  bias <- rep("unbiased", length(diff.active.barcoded.cells))
  bias[bias.score >= quantile(bias.score,0.75)[[1]]] <- "M"
  bias[bias.score <= quantile(bias.score,0.25)[[1]]] <- "E"
  c <- table(bias)
  
  bias <- rep("unbiased", length(diff.active.barcoded.cells))
  bias[bias.score >= quantile(bias.score,0.9)[[1]]] <- "M"
  bias[bias.score <= quantile(bias.score,0.1)[[1]]] <- "E"
  d <- table(bias)
  
  combined <- data.frame(rbind(a,b,c,d))
  
  
  par(mfrow=c(1,2))
  barplot(combined$E,col =c("grey10","grey30","grey50","grey70"),xlab = "threshold value",ylab = "number of erythroid-biased cells",
          names = c("55%","66%","75%","90%"))
  
  
  barplot(combined$M,col =c("grey10","grey30","grey50","grey70"),xlab = "threshold value",ylab = "number of myeloid-biased cells",
          names = c("55%","66%","75%","90%"))
  
  
}



#generate a volcano plot to visualise DEGs
volcanoPlotHSPC <- function(res,gene.list){
  res$gene <- rownames(res)
  
  res$sign <- 0
  
  
  res$sign[which(res$minimump_p_val < 0.05 & res$avg_log2FC > 0.1)] <- 2
  res$sign[which(res$minimump_p_val < 0.05 & res$avg_log2FC < -0.1)] <- 1
  res$minimump_p_val[res$minimump_p_val == 0] <- 1e-312
  
  p <- ggplot(data=res, aes(x=avg_log2FC, y=-log10(minimump_p_val), colour=as.factor(sign))) + geom_point( size=1) +
    scale_color_manual(name="", values=c("1" = "red","2" = rgb(67/255,138/255,201/255,0.4),"3"=rgb(128/255,128/255, 128/255,0.4), "0"=rgb(220/255,220/255, 220/255,0.1))) +  
    theme(legend.position = "none") + xlim(-1.4,1.4) + 
    xlab("log2 fold change") + ylab("-log10 pvalue") + 
    geom_vline(xintercept=c(-0.1, 0.1), linetype=2,colour = "grey70") + 
    geom_hline(yintercept=-log10(0.05), linetype=2,colour = "grey70") 
  
  p <- p  + theme(
    text = element_text(size = 20),
    # Remove panel border
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "grey")
  ) +geom_text_repel(data=res[gene.list,],aes(label=gene),fontface = "bold",color = "black",min.segment.length = unit(0, 'lines'),size = 5)
  
  
  return(p)
  
}

permutationTest <- function(gene.set,sobj,metric,simRuns = 100){
  
  output.bias.rho <- matrix(0,ncol = 1, nrow = simRuns)
  gene.set.size <- length(gene.set)
  
  for(i in 1:simRuns){
    
    Sys.sleep(0.1)
    set.seed(1245 * i)
    
    gene.index <- sample(1:nrow(sobj), gene.set.size, replace=FALSE)
    temp.gene.set <- rownames(dataset.integrated)[gene.index]
    
    sobj <- AddModuleScore(sobj,features = list(temp.gene.set), name = "temp")
    cor <- cor.test( sobj$temp1,sobj@meta.data[,metric], method = "spearman")
    
    output.bias.rho[i] <- cor$estimate[[1]]
    
  }
  
  return(output.bias.rho)
  
}




plotPermutationResults <- function(){
  df <- data.frame(x = rep("Fate-Myeloid",200),y = m.perm)
  df <- rbind(df, data.frame(x = rep("MetaFate-Myeloid",200),y = m.metab.perm))
  df <- rbind(df, data.frame(x = rep("Transcription Factor-Myeloid",200),y = m.tf.perm))
  df <- rbind(df, data.frame(x = rep("MPP3",200),y =  mpp3.perm))
  
  
  pointdata <- data.frame(x = c("Fate-Myeloid","MetaFate-Myeloid","Transcription Factor-Myeloid","MPP3"),y = c(0.34,0.29,0.27,0.2))
  
  df$x <- factor(df$x,levels = c("Fate-Myeloid","MetaFate-Myeloid","Transcription Factor-Myeloid","MPP3"))
  
  ggplot(df,aes(x,y,color = x )) + geom_jitter(width = 0.25,color = rgb(0,0,1,0.3),size = 2)  
  + coord_flip() + theme_classic() +theme(axis.text=element_text(size=15,face = 'bold',colour = "black"),
    axis.title=element_text(size=0,face="bold"), axis.line = element_line(size = 1.0, colour = "black"),
    axis.ticks.length=unit(.25, "cm")) + geom_point(data = pointdata,  mapping = aes(x = x, y = y,size = 3) , 
                                                    color = "red") + NoLegend()
  
  
}
