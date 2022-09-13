############################ Helper Methods ################################
#
# @author: Jason Cosgrove (jason.cosgrove@curie.fr)
# @date:   23/12/2019
#
############################################################################


# input:  two vectors
# output: a vector of booleans stating whether the element of x occurs in y
'%ni%' <- function(x,y)!('%in%'(x,y))


# given a lineage, return the appropriate weighting for the classifier
# input:   a string "Lymphoid", "Myeloid", "Erythroid"
# outputs: a number representing the weighting of the lineage 
#          for the machine learning classifier
calculateLineageWeights <- function(input){
	
	lineages <- c("Lymphoid", "Myeloid", "Erythroid")
	
	if(input %ni% lineages){ stop('input is not a recognised lineage')}
	
	else{
       if(input == "Lymphoid"){ return(0.65) } 
       if(input == "Myeloid"){  return(0.28) } 
       if(input == "Erythroid"){ return(1.0) } 
  }
}


# plot a kegg pathway overlaying gene expression data
# inputs:  the kegg pathway ID, the output suffix, and foldchanges in expression between two conditions
# outputs: sends kegg pathway images to file, in the current working directory
plot_pathway = function(pid, out_suffix,foldchanges) 
			   pathview(gene.data=foldchanges, pathway.id=pid, species="mmu", 
                new.signature=FALSE,out.suffix = out_suffix,low = list(gene = "blue", cpd = "blue"), 
                mid = list(gene = "gray", cpd = "gray"), high = list(gene = "red", cpd ="yellow"))


# this method processes the scRNAseq dataset from Tusi et al Nature (2018)
# returns a processed dataset in the form of a seurat object
createTusiSeuratObject <- function(){
  
  # load in data from the Tusi paper
  dat.ctrl <- read.csv("datasets/Tusi/basal_bone_marrow/GSM2388072_basal_bone_marrow.filtered_normalized_counts.csv")
  dat.epo <- read.csv("datasets/Tusi/epo_bone_marrow/GSM2388073_epo_bone_marrow.filtered_normalized_counts.csv")
  
  #Keep only the columns that we need for the analysis. 
  dat.ctrl.filtered <- t(dat.ctrl[,6:28210])
  colnames(dat.ctrl.filtered) <- dat.ctrl$cell_id
  
  dat.epo.filtered <- t(dat.epo[,6:28210])
  colnames(dat.epo.filtered) <- dat.epo$cell_id
  
  #now combine the WT and EPO datasets into an individual datafrane
  dataset <- cbind(dat.epo.filtered,dat.ctrl.filtered)
  
  # generate the metadata for the dataset saying whether a cell came from a control BM or one treated with EPO
  Batch <- c(rep("EPO_treated",ncol(dat.epo.filtered)),rep("Control",ncol(dat.ctrl.filtered)))
  cell_anns <- data.frame(batch = Batch)
  rownames(cell_anns) <- colnames(dataset)
  
  # convert the dataset into a seurat dataset
  TusiEPO <- CreateSeuratObject(counts= dataset, min.cells = 20, min.features = 750,meta.data = cell_anns, 
                                project = "EPO")
  
  #the data are already normalsied so we just need to scale for visualisation
  TusiEPO <- ScaleData(TusiEPO)
  
  #add the spring visualisation as metadata to the seurat object
  ctrl.metadata <- read.csv("datasets/Tusi/basal_bone_marrow/basal_bone_marrow.metadata.csv")
  epo.metadata <- read.csv("datasets/Tusi/epo_bone_marrow/epo_bone_marrow.metadata.csv")
  metadata.subset <- rbind(epo.metadata[,6:7],ctrl.metadata[,6:7])
  colnames(metadata.subset) <- c("spring1","spring2")
  rownames(metadata.subset) <- c(paste(epo.metadata$cell_id),paste(ctrl.metadata$cell_id))
  metadata.subset <- metadata.subset[colnames(TusiEPO),]
  TusiEPO[["spring"]] <- CreateDimReducObject(embeddings = as.matrix(metadata.subset), key = "spring_", assay = DefaultAssay(TusiEPO))
  
  return(TusiEPO)
  
}



# given a list of genes see if the expression is higher in 
# both positive samples than in negative samples
filterOnBulkProgenitors <- function(sig,pos,neg,dataset){
  
  toKeep <- c()
  
  #iterate over the signature setting a boolean to TRUE if expression 
  # is higher in the pos group than the negative fraction. 
  for(i in 1:length(sig)){  
    #expression must be higher in both samples
    if(dataset[sig[i],pos[1]] < 
       max(dataset[sig[i],neg])
       & dataset[sig[i],pos[2]] < 
       max(dataset[sig[i],neg]))
    {
      toKeep <- c(toKeep, FALSE)
    }
    else{
      toKeep <- c(toKeep, TRUE)
    }
  }
  
  #keep only the genes that have higher expression in the pos group
  return(sig[toKeep])
}


# convert a gene name to an entrez gene name
# inputs: A single gene name and an entrezgene object generated from the biomart getBM object
# outputs: A vector of entrezgene names with the conventional gene names as headers
convertToEntrezGenes <- function(gene_name,entrezgenes){
  if(length(entrezgenes[entrezgenes$external_gene_name == gene_name,2]) > 0){
    return(entrezgenes[entrezgenes$external_gene_name == gene_name,2][[1]])
  }
  else{return(NA)}
}



#calculate the shapley score for all variables in a classifier model
# inputs: all_features = all genes used in the classifier model
#		  dict   = a dictionary that links lineage to its numeric entry
#				   in our classifier model
#         nclass = the number of output classes in your classifier model
#		  data   = the input data used in training the classifier
#
# outputs: A dataframe with the shapley score for each instance (sample in the training data) 
#          and each variable and each class
calculateSHAPScores <- function(all_features,dict,nclass_input,data,nfolds,nrounds){
  
  topn   = length(all_features)
  top_features <- all_features
  nclass = nclass_input - 1 #numbers of classes to predict thqt you have, we do minus 1 because we start from zero
  trees0 = seq(from=0, by=nclass, length.out=nrounds)
  all_res_SHAP = data.frame(value=numeric(), SHAP=numeric(), index=numeric(), fold=numeric(), class=factor(levels = c(0:nclass)), feature=character())
  keys = all_features[1:50]  #Arbitrary keys for finding the original index of permuted data
  
  
  #for each sample get a shap score for each feature for each predictive class. 
  for(fold in 1:nfolds){
    print(fold)
    for(class in 0:nclass){
      
      res_SHAP = xgb.plot.shap(data, model = cv$models[[fold]],
                               features = top_features[1:topn],
                               trees = trees0 + class,
                               target_class = class, plot = F)
      
      
      
      # The data is randomly permuted, so we need to find a consistent index. 
      #is that the cells that are permited?
      index = apply(res_SHAP$data[,keys], 1, function(row){
        which(apply(data[,keys], 1, function(x){all(x==row)}))
      })
      
      
      d <- lapply(top_features[1:topn], function(feature) {
        data.frame(value   = res_SHAP$data[,feature], 
                   SHAP    = res_SHAP$shap_contrib[,feature],
                   index   = index,
                   fold    = rep(fold, nrow(res_SHAP$data)), 
                   class   = factor(rep(dict[[class + 1]], nrow(res_SHAP$data))),
                   feature = rep(feature, nrow(res_SHAP$data)))
      })
      
      all_res_SHAP <- rbind(all_res_SHAP,do.call(rbind,d))
      
    }
  }
  
  #format the results 
  SHAP_formatted<- all_res_SHAP %>% group_by(class,feature,index) %>% summarise(value = mean(value), SHAP = mean(SHAP)) %>%
    dplyr::mutate(medianAbs = mean(abs(SHAP))) %>% dplyr::arrange(medianAbs) 
  
  return(data.frame(SHAP_formatted))
}


# format the ShapleyScore analysis to make the outputs more interpretable. 
# input: the output of calculateSHAPScores()
# output: for a given variable for a given class take the mean 
# shapley score across all instances
createSHAPSummary <- function(SHAP_results){
  genes <- paste(unique(SHAP_results$feature))
  
  e.shapscores <- matrix(0, nrow = length(genes), ncol = 1)
  colnames(e.shapscores) <- c( "SHAP_E")
  rownames(e.shapscores) <- (genes)
  
  m.shapscores <- matrix(0, nrow = length(genes), ncol = 1)
  colnames(m.shapscores) <- c( "SHAP_M")
  rownames(m.shapscores) <- (genes)
  
  l.shapscores <- matrix(0, nrow = length(genes), ncol = 1)
  colnames(l.shapscores) <- c( "SHAP_L")
  rownames(l.shapscores) <- (genes)
  
  
  #for each gene
  for(i in 1:length(genes)){
    e.shapscores[i,1] <- mean(SHAP_results[SHAP_results$feature ==  genes[i] & paste(SHAP_results$class) == "Erythroid",6])
    m.shapscores[i,1] <- mean(SHAP_results[SHAP_results$feature ==  genes[i] & paste(SHAP_results$class) == "Myeloid",6])
    l.shapscores[i,1] <- mean(SHAP_results[SHAP_results$feature ==  genes[i] & paste(SHAP_results$class) == "Lymphoid",6])
  }
  
  shapscores <- cbind(e.shapscores, m.shapscores, l.shapscores)
  
  return(shapscores)
  
}



# see how our classifier model performs in the RPP compartment
# inputs: 
#			haem: the processed haemosphere dataset in the form of a seurat object
#			cv: the outputs of xgb.cv containing the trained classifier models
#
# outputs:
#       the classifier predictions for the RPP compartment
#
assessPerformanceOnPrecursors <- function(haem,cv){
  
  #create a dictionary that links lineage to its numeric entry
  labels.ordered <- colnames(table(xgb.preds, labels))
  dict <- vector(mode="list", length=length(labels.ordered))
  names(dict) <- 1:(length(labels.ordered))
  for(i in 1:length(labels.ordered)){dict[[i]] <- labels.ordered[i]}
  
  #create an input matrix with gene expression profiles just for lineage biased precursors
  mat <- as.matrix(haem@assays$RNA@scale.data)
  test.dat <- mat[interesting.genes,haem@meta.data$lineage == "RPP" | 
                    haem@meta.data$lineage == "MPP"]
  test.data <- as.data.frame(t(test.dat))
  
  # for each cell sample, lets generate a lineage prediction from each of our 10 models
  # and store it in a matrix
  # we have 10 models because we did 10-fold cross validation
  output.mat <- matrix(0, nrow = length(rownames(test.data)), ncol = nfolds)
  output.mat.numeric <- matrix(0, nrow = length(rownames(test.data)), ncol = nfolds)
  rownames(output.mat) <- rownames(test.data)
  for(i in 1:length(rownames(test.data))){
    for(j in 1:nfolds){
      index <- which.max(predict(cv$models[[j]],as.matrix(test.data[i,])))
      output.mat[i,j] <- dict[[index]]
      output.mat.numeric[i,j] <- index
    }
  }
  
  #for each cell type find out what the most popular lineage classification is and store it in a vector
  Consensus <- function(x){return(names(which.max(table(x))))}
  consensus.prediction <- data.frame(apply(output.mat, 1,Consensus))
  x <- consensus.prediction[c("CLP.1","CLP.2", "GMP.1", "GMP.2", "MEP.1","MEP.2"),] # we dont have CD9s because they're used in the training
  
  #CD9s are: Lin:- cKit:+ Sca1:- CD150:+ FcgRII/III (CD16/32):- Endoglin:lo CD9:hi
  #see table S1 at https://www.pnas.org/content/pnas/suppl/2012/01/25/1121385109.DCSupplemental/pnas.201121385SI.pdf?targetid=nameddest%3DST1
  names(x) <- c("CLP.1","CLP.2", "GMP.1", "GMP.2", "MEP.1","MEP.2")
  
  return(x)
  
}




# preps the data matrix so it can be used as an input for the classifier
# Inputs: 
#			haem: the seurat object for the haemopedia dataset
#			interesting.genes: the genes that you want to use in training the classifier
# 			MPPONLY: a boolean stating whether you want to use RPPs in the training
#					if FALSE the RPPs are not used in training
#
# Outputs: Returns a dataset that can be used as input into our classifer training algorithm
prepData <- function(interesting.genes,haem,MPPONLY){
  
  
  if(MPPONLY){
    mat <- as.matrix(haem@assays$RNA@scale.data)
    dat <- mat[interesting.genes, haem@meta.data$lineage != "MPP" & colnames(haem) %ni% c("CMP.1","CMP.2")]
    data.withlineage <- as.data.frame(t(dat))
    
    #add lineage as an extra variable
    data.withlineage$lin <- as.factor(haem@meta.data$lineage[haem@meta.data$lineage != "MPP" & colnames(haem) %ni% c("CMP.1","CMP.2")])
    
    data.withlineage[c("CLP.1","CLP.2"),]$lin <- as.factor("Lymphoid")
    
    data.withlineage[c("BEMP.2" , "CD9Hi.2", "MEP.1" , "MEP.2","CD9Hi.1", "BEMP.1" ),]$lin <- as.factor("Erythroid") 
    data.withlineage[c("GMP.1",          "GMP.2"  ,                           
                       "FcgRCD150.1",    "FcgRCD150.2" ,   "PreGMFlt3Neg.1" ,"PreGMFlt3Neg.2" ,"PreGMFlt3Pos.1" ,"PreGMFlt3Pos.2",
                       "GMP_IRF8lo.1" ,  "GMP_IRF8lo.2" ,  "GMP_IRF8int.1" , "GMP_IRF8int.2" , "GMP_IRF8hi.1" ,  "GMP_IRF8hi.2"   ),]$lin <- as.factor("Myeloid") 
    return(data.withlineage)
    
    
  }
  
  # prep the matrix that we will use for classification
  
  else{
    mat <- as.matrix(haem@assays$RNA@scale.data)
    dat <- mat[interesting.genes,haem@meta.data$lineage != "RPP" &
                 haem@meta.data$lineage != "MPP" ]
    data.withlineage <- as.data.frame(t(dat))
    
    #add lineage as an extra variable
    data.withlineage$lin <- as.factor(haem@meta.data$lineage[haem@meta.data$lineage 
                                                             != "RPP"& haem@meta.data$lineage != "MPP"])
    
    return(data.withlineage)
  }
}


# custom volcano plot of limma differentially expressed genes
LimmaVolcano <- function(res, main="", fct=2, pt=0.05){
  res$sign <- 0
  res$sign[which(res$adj.P.Val < pt & abs(res$logFC) > fct)] <- 1
  
  res$sign[which(res$adj.P.Val < pt & abs(res$logFC) > fct & res$genes %in% metabolic.genes)] <- 2
  p <- ggplot(data=res, aes(x=logFC, y=-log10(adj.P.Val), colour=as.factor(sign))) + geom_point( size=2) +
    
    scale_color_manual(name="", values=c("4" = "orange","3" = "blue" ,"2" = "red","1"=
                                           rgb(80/255,80/255, 80/255,0.2), "0"=rgb(220/255,220/255, 220/255,0.2))) +  
    ggtitle(paste0("Volcano Plot - Limma ", main)) +
    theme(legend.position = "none") + xlim(-12,12) + 
    xlab("log2 fold change") + ylab("-log10 adj pvalue") + 
    geom_vline(xintercept=c(-fct, fct), linetype=2) + 
    geom_hline(yintercept=-log10(pt), linetype=2)
  p
}



# convert the haemopedia data into a seurat object, update the metadata and scale 
# the countmatrix so that it is centered around zero
# this is useful for visualisation, and also for the machine learning pipeline. 
# Inputs: A count matrix for the haemopedia dataset
# Outputs: A seurat object
convertToSeurat <- function(countmatrix){
  #lets conver the data into a seurat object so we can use some of their visualisation functionality
  haem <- CreateSeuratObject(counts = countmatrix, min.cells = 0, min.features =0, 
                             project = "haemopedia")
  
  #add lineage as metadata
  haem@meta.data$lineage <- lin
  
  #tell seurat to use all genes when performing dimensionality reduction, and clustering
  haem@assays$RNA@var.features <- rownames(haem)
  
  #we scale the data prior to performing classification
  haem <- ScaleData(object = haem)
  
  return(haem)
  
}



# Preprocessing of the haemopedia datasets
# filter out lowly expressed genes, and also cells we do not want to analyse
# Inputs: A haemopedia dataset and a vector with the samples that you do not want to include in the analysis
# Outputs: return a filtered DGElist object
filterData <- function(haemopedia,cells.to.remove){
  # filter out lowly expressed genes, be careful not to change this filter or it can mess with the normalisation
  # as many erythroid cells are lowly expressed. We set this following a code review with nicolas servant
  logcounts <- edgeR::cpm(haemopedia)
  isexpr <- names(which(apply(logcounts, 1, function(x){length(which(x>=1))})>=3))
  haemopedia.genes_filtered  <- haemopedia$counts[isexpr,]
  
  #filter out cells that we do not want to include in our analysis
  haemopedia.genes_cells_filtered <- haemopedia.genes_filtered[, colnames(haemopedia.genes_filtered) %ni% cells.to.remove]
  
  #create a new DGE object from the filtered data
  haemopedia.filtered <- DGEList(haemopedia.genes_cells_filtered, genes=rownames(haemopedia.genes_cells_filtered))
  
  return(haemopedia.filtered)
}


# For the haemopedia dataset add the lineage annotation for each sample
# Inputs:  The processed Haemopedia dataset
# Outputs: returns a vector with each cell type, and the lineage it belongs to. Some of the restricted potential progenitors
#          are assigned a lineage as this helps with the supervised learning approach. 
performLineageAnnotation <- function(haemopedia){
  
  lymphoid <- c("B Cell Lineage", "T Cell Lineage" , "Innate Lymphocyte Lineage", "NK Cell Lineage")
  myeloid <- c("Dendritic Cell Lineage" , "Neutrophil Lineage", "Mast Cell Lineage",  "Macrophage Lineage",             
               "Basophil Lineage",  "Eosinophil Lineage"  )
  
  h.meta <-read.table("datasets/Haemopedia/preprocessing/Haemopedia-Mouse-RNASeq_samples.txt",sep="\t",header=TRUE)
  
  lin <- c()
  for(i in 1:nrow(haemopedia$samples)){
    lineage_column = 4
    metadata <- h.meta
    if(rownames(haemopedia$samples)[i] == "MEP.1" |rownames(haemopedia$samples)[i] == "MEP.2" ){
      lin <- c(lin,"RPP")
    }
    else{
      lin <- c(lin,paste(metadata[metadata$sampleId == rownames(haemopedia$samples)[i],lineage_column]))
    }
    
  }
  names(lin) <- rownames(haemopedia$sample)
  
  lin[lin %in% lymphoid] <- "Lymphoid"
  lin[lin %in% myeloid] <- "Myeloid"
  lin[lin ==  "Erythrocyte Lineage"] <- "Erythroid"
  lin[lin ==  "Megakaryocyte Lineage"] <- "Meg"
  lin[lin ==  "Multi Potential Progenitor"] <- "MPP"
  lin[lin ==  "Restricted Potential Progenitor"] <- "RPP"
  return(lin)
}





# This method creates a geneset library from a GMT file
# Inputs: the path to the GMT file
# Outputs: A geneset library
makeGeneSetFromGMT <- function(filename){
  gs <- GSA.read.gmt(filename)
  
  ## number of gene sets
  n <- length(gs$geneset.names)
  
  ## create environment
  env <- new.env(parent=globalenv())
  invisible(lapply(1:n,function(i) {
    genes <- as.character(unlist(gs$genesets[i]))
    name <- as.character(gs$geneset.names[i])
    assign(name, genes, envir = env)
  }))
  
  return(env)  
}



# Use the org.Mm.eg.db package for GO annotations for genes from a given dataset
# inputs: A dataset where you would like to extract some gene names
# outputs: A GO geneset library with all genes that are found in your dataset
makeGeneSetFromMM_eg_db <- function(dataset){
  # Translate gene names to ids
  ids <- unlist(lapply(mget(rownames(dataset), org.Mm.egALIAS2EG, ifnotfound = NA), function(x) x[1]))
  # Reverse map
  rids <- names(ids)
  names(rids) <- ids
  # Convert ids per GO category to gene names
  go.env <- eapply(org.Mm.egGO2ALLEGS, function(x) as.character(na.omit(rids[x])))
  go.env <- clean.gos(go.env) # Remove GOs with too few or too many genes
  go.env <- list2env(go.env)  # Convert to an environment
  
  return(go.env)
}



# This method allows you to convert an object from scanpy (in H5 format) into a seurat object
# this is useful as we took a lot of our datasets and preprocessing from the Theis lab github
# The method is rewritten from Seurat to read in a H5 file. The original method from Seurat was
# not working for me so i rewrote part of it. 

ReadH5AD.H5File <- function(
  filename,
  assay = 'RNA',
  layers = 'data',
  verbose = TRUE,
  ...
) {
  # Pull assay data
  # If X is an H5D, assume scaled
  # Otherwise, if file$exists(name = 'raw'), assume X is normalized
  # Otherwise, assume file[['X']] is raw counts
  file <- hdf5r::h5file(filename = filename, mode = 'r')
  
  if (verbose) {
    message("Pulling expression matrices and metadata")
  }
  #if (is(object = file[['X']], class2 = 'H5Group')) {
  # x <- as.sparse(x = file[['X']])
  # else {
  x <- file[['X']][, ]
  #}
  # x will be an S3 matrix if X was scaled, otherwise will be a dgCMatrix
  scaled <- is.matrix(x = x)
  if (verbose) {
    message("Data is ", ifelse(test = scaled, yes = 'scaled', no = 'unscaled'))
  }
  print("code is reached 1")
  
  # Pull cell- and feature-level metadata
  obs <- file[['obs']][]
  x.var <- file[['var']][]
  rownames(x = x) <- rownames(x = x.var) <- x.var$index
  colnames(x = x) <- rownames(x = obs) <- obs$index
  # Pull raw expression matrix and feature-level metadata
  if (file$exists(name = 'raw.X')) {
    print("code is reached point 2")
    raw <- as.sparse(x = file[['raw.X']][,])
    raw.var <- file[['raw.var']][]
    slot(object = raw, name = 'Dim') <- c(nrow(x = raw.var), nrow(x = obs))
    rownames(x = raw) <- rownames(x = raw.var) <- raw.var$index
    colnames(x = raw) <- obs$index
    raw.var <- raw.var[, -which(x = colnames(x = raw.var) == 'index'), drop = FALSE]
    x.slot <- ifelse(test = scaled, yes = 'scale.data', no = 'data')
  } else {
    # If X is scaled, we required normalized data present in raw
    if (scaled) {
      stop("Seurat requires normalized data present in the raw slot when X is scaled")
    } else {
      x.slot <- 'raw'
    }
  }
  
  print("code is reached")
  
  obs <- obs[, -which(x = colnames(x = obs) == 'index'), drop = FALSE]
  x.var <- x.var[, -which(x = colnames(x = x.var) == 'index'), drop = FALSE]
  # Merge raw.var and x.var
  # Only happens when we have a raw.X and raw.var in the h5ad file
  if (x.slot != 'raw') {
    if (verbose) {
      message("Merging feature-level metadata dataframes")
    }
    x.var <- x.var[, -which(x = colnames(x = x.var) %in% colnames(x = raw.var))]
    meta.features <- merge(x = raw.var, y = x.var, by = 0, all = TRUE)
    rownames(x = meta.features) <- meta.features$Row.names
    meta.features <- meta.features[, -which(x = colnames(x = meta.features) == 'Row.names'), drop = FALSE]
    rm(raw.var)
  } else {
    meta.features <- x.var
  }
  # Fix meta feature colnames
  colnames(x = meta.features) <- gsub(
    pattern = 'dispersions_norm',
    replacement = 'mvp.dispersion.scaled',
    x = colnames(x = meta.features)
  )
  colnames(x = meta.features) <- gsub(
    pattern = 'dispersions',
    replacement = 'mvp.dispersion',
    x = colnames(x = meta.features)
  )
  colnames(x = meta.features) <- gsub(
    pattern = 'means',
    replacement = 'mvp.mean',
    x = colnames(x = meta.features)
  )
  colnames(x = meta.features) <- gsub(
    pattern = '_',
    replacement = '.',
    x = colnames(x = meta.features)
  )
  if ('highly.variable' %in% colnames(x = meta.features)) {
    meta.features$highly.variable[is.na(x = meta.features$highly.variable)] <- FALSE
  }
  rm(x.var)
  #CheckGC()
  # Fix metadata colnames
  colnames(x = obs) <- gsub(
    pattern = '_',
    replacement = '.',
    x = colnames(x = obs)
  )
  colnames(x = obs) <- gsub(
    pattern = 'n.genes',
    replacement = paste0('nFeatures_', assay),
    x = colnames(x = obs)
  )
  colnames(x = obs) <- gsub(
    pattern = 'n.counts',
    replacement = paste0('nCount_', assay),
    x = colnames(x = obs)
  )
  # Assemble assay object
  if (verbose) {
    message("Creating assay object")
    message(
      "Storing X as ",
      x.slot,
      ifelse(
        test = x.slot != 'counts',
        yes = paste(" and raw as", ifelse(test = scaled, yes = 'data', no = 'counts')),
        no = ''
      )
    )
  }
  if (scaled) {
    assays <- list(CreateAssayObject(data = raw))
    assays[[1]] <- SetAssayData(
      object = assays[[1]],
      slot = 'scale.data',
      new.data = x
    )
    rm(raw)
  } else if (x.slot == 'data') {
    assays <- list(CreateAssayObject(counts = raw))
    assays[[1]] <- SetAssayData(
      object = assays[[1]],
      slot = 'data',
      new.data = x
    )
    rm(raw)
  } else {
    assays <- list(CreateAssayObject(counts = x))
  }
  names(x = assays) <- assay
  # Add meta feature information
  if (ncol(x = meta.features) > 0) {
    assays[[assay]][[names(x = meta.features)]] <- meta.features
  }
  # Add highly variable feature information
  if ('highly.variable' %in% colnames(x = assays[[assay]][[]])) {
    if (verbose) {
      message("Setting highly variable features")
    }
    hvf.info <- HVFInfo(object = assays[[assay]], selection.method = 'mvp')
    hvf.info <- hvf.info[order(hvf.info$dispersion, decreasing = TRUE), , drop = FALSE]
    means.use <- (hvf.info$mean > 0.1) & (hvf.info$mean < 8)
    dispersions.use <- (hvf.info$dispersion.scaled > 1) & (hvf.info$dispersion.scaled < Inf)
    top.features <- rownames(x = hvf.info)[which(x = means.use & dispersions.use)]
    VariableFeatures(object = assays[[assay]]) <- top.features
  } else if (verbose) {
    message("No variable feature expression found in h5ad file")
  }
  Key(object = assays[[assay]]) <- paste0(tolower(x = assay), '_')
  rm(x)
  #CheckGC()
  # Get dimensional reduction information
  # If data isn't scaled, don't bother
  if (scaled && file$exists(name = 'obsm')) {
    if (verbose) {
      message("Pulling dimensional reduction information")
      message("Pulling cell embeddings")
    }
    # Pull cell embeddings
    if (inherits(x = file[['obsm']], what = 'H5Group')) {
      embed.reduc <- names(x = file[['obsm']])
      embeddings <- sapply(
        X = embed.reduc,
        FUN = function(x) {
          return(t(x = file[['obsm']][[x]][, ]))
        },
        simplify = FALSE,
        USE.NAMES = TRUE
      )
    } else {
      embed.reduc <- file[['obsm']]$get_type()$get_cpd_labels()
      embed.n <- sapply(
        X = file[['obsm']]$get_type()$describe()$cpd_types,
        FUN = '[[',
        'array_dims'
      )
      names(x = embed.n) <- embed.reduc
      ncells <- file[['obsm']]$dims
      embeddings <- lapply(
        X = embed.reduc,
        FUN = function(r) {
          return(t(x = vapply(
            X = 1:ncells,
            FUN = function(i) {
              return(file[['obsm']][i][[r]])
            },
            FUN.VALUE = numeric(length = embed.n[[r]])
          )))
        }
      )
      names(x = embeddings) <- embed.reduc
    }
    
    
    
    
    # Set cell names for embeddings matrices
    for (i in 1:length(x = embeddings)) {
      rownames(x = embeddings[[i]]) <- colnames(x = assays[[assay]])
    }
    # Pull feature loadings
    if (file$exists(name = 'varm')) {
      if (verbose) {
        message("Pulling feature loadings")
      }
      if (inherits(x = file[['varm']], what = 'H5Group')) {
        load.reduc <- names(x = file[['varm']])
        loadings <- sapply(
          X = load.reduc,
          FUN = function(x) {
            return(t(x = file[['varm']][[x]][, ]))
          },
          simplify = FALSE,
          USE.NAMES = TRUE
        )
      } else {
        load.reduc <- file[['varm']]$get_type()$get_cpd_labels()
        load.n <- sapply(
          X = file[['varm']]$get_type()$describe()$cpd_types,
          FUN = '[[',
          'array_dims'
        )
        names(x = load.n) <- load.reduc
        nfeatures <- file[['varm']]$dims
        loadings <- lapply(
          X = load.reduc,
          FUN = function(r) {
            return(t(x = vapply(
              X = 1:nfeatures,
              FUN = function(i) {
                return(file[['varm']][i][[r]])
              },
              FUN.VALUE = numeric(length = load.n[[load.reduc]])
            )))
          }
        )
      }
      match.ind <- lapply(
        X = gsub(pattern = 's$', replacement = '', x = tolower(x = load.reduc)),
        FUN = grep,
        x = embed.reduc
      )
      no.match <- which(x = sapply(X = match.ind, FUN = length) != 1)
      if (length(x = no.match) >= 1) {
        warning(
          "Unable to determine where the following feature loadings belong: ",
          paste(load.reduc[no.match], collapse = ', '),
          call. = FALSE,
          immediate. = TRUE
        )
        loadings <- loadings[-no.match]
        load.reduc <- load.reduc[-no.match]
        match.ind <- match.ind[-no.match]
      }
      names(x = loadings) <- embed.reduc[unlist(x = match.ind)]
      for (i in 1:length(x = loadings)) {
        rownames(x = loadings[[i]]) <- rownames(x = GetAssayData(
          object = assays[[assay]],
          slot = 'scale.data'
        ))
      }
    } else {
      if (verbose) {
        message("No feature loadings found")
      }
      loadings <- list()
    }
    # Create DimReduc objects
    dim.reducs <- vector(mode = 'list', length = length(x = embed.reduc))
    for (i in 1:length(x = embed.reduc)) {
      r <- embed.reduc[i]
      key <- tolower(x = gsub(pattern = 'X_', replacement = '', x = r))
      key <- switch(
        EXPR = key,
        'pca' = 'PC',
        'tsne' = 'tSNE',
        toupper(x = key)
      )
      key <- paste0(key, '_')
      stdev <- if (r == 'X_pca' && file$exists(name = 'uns') && file$exists(name = 'uns/pca/variance')) {
        sqrt(x = file[['uns/pca/variance']][])
      } else {
        numeric(length = 0L)
      }
      dim.reducs[[i]] <- CreateDimReducObject(
        embeddings = embeddings[[r]],
        loadings = loadings[[r]] %||% new(Class = 'matrix'),
        assay = assay,
        stdev = stdev,
        key = key
      )
    }
    # Properly name dimensional reductions
    names(x = dim.reducs) <- gsub(
      pattern = 'X_',
      replacement = '',
      x = embed.reduc
    )
    # Clean up
    rm(embeddings, loadings)
    #CheckGC()
  } else {
    if (verbose) {
      message("No dimensional reduction information found")
    }
    dim.reducs <- list()
  }
  # Create the Seurat object
  if (verbose) {
    message("Assembling Seurat object")
  }
  # Create a project name, will be used as identity classes
  project <- gsub(
    pattern = '\\.h5ad',
    replacement = '',
    x = basename(path = file$filename)
  )
  object <- new(
    Class = 'Seurat',
    assays = assays,
    meta.data = obs,
    version = packageVersion(pkg = 'Seurat'),
    project.name = project
  )
  # Set default assay and identity information
  DefaultAssay(object = object) <- assay
  Idents(object = object) <- project
  # Add dimensional reduction infrom
  if (scaled && length(x = dim.reducs) >= 1) {
    for (r in names(x = dim.reducs)) {
      object[[r]] <- dim.reducs[[r]]
    }
  }
  # Get graph information
  if (scaled && file$exists(name = 'uns') && file$exists(name = 'uns/neighbors')) {
    if (verbose) {
      message("Finding nearest neighbor graph")
    }
    graph <- as.sparse(x = file[['uns/neighbors/distances']])
    colnames(x = graph) <- rownames(x = graph) <- colnames(x = object)
    method <- ifelse(
      test = file[['uns/neighbors/params']]$exists(name = 'method'),
      yes = file[['uns/neighbors/params/method']][],
      no = 'adata'
    )
    object[[paste(assay, method, sep = '_')]] <- as.Graph(x = graph)
  } else if (verbose) {
    message("No nearest-neighbor graph")
  }
  # Add layers
  if (isFALSE(x = layers)) {
    if (verbose) {
      message("Not pulling layers")
    }
  } else if (file$exists(name = 'layers')) {
    file.layers <- names(x = file[['layers']])
    layers <- rep_len(
      x = tolower(x = layers),
      length.out = length(x = file.layers)
    )
    if (!all(layers %in% c('counts', 'data'))) {
      stop("'layers' must be either 'counts' or 'data'", call. = FALSE)
    }
    names(x = layers) <- file.layers
    for (layer in file.layers) {
      layer.dest <- layers[[layer]]
      if (verbose) {
        message(
          "Reading ",
          layer,
          " into new assay, putting data into ",
          layer.dest
        )
      }
      layer.data <- if (inherits(x = file[['layers']][[layer]], what = 'H5Group')) {
        as.sparse(x = file[['layers']][[layer]])
      } else {
        file[['layers']][[layer]][, ]
      }
      dimnames(x = layer.data) <- dimnames(x = object)
      layer.assay <- switch(
        EXPR = layer.dest,
        'counts' = CreateAssayObject(
          counts = layer.data,
          min.cells = -1,
          min.features = -1
        ),
        'data' = CreateAssayObject(data = layer.data),
        stop("Unknown layer destination: ", layer.data, call. = FALSE)
      )
      object[[layer]] <- layer.assay
    }
  } else if (verbose) {
    message("No additional layers found")
  }
  return(object)
}


`%||%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(lhs)
  } else {
    return(rhs)
  }
}

