################# HELPER METHODS #################
#
# About: Helper methods for the lentiviral barcoding 
#        of WT vs G6PDTg MPPs
#
# Author: Jason Cosgrove (jason.cosgrove@curie.fr)
# Date: 05/04/2022
################################################


#perform column-wise normalisation of the barcode count matrix
columnNormalisation <- function(data,barcodes){
  
  data.normed <- apply(data,2, function(x) (x/sum(x))) # here we normalise for each sample
  rownames(data.normed) <- barcodes
  return(data.normed)
  
}

#remove clones that are not in any of the samples
removeUninformativeClones <- function(mat){
  m.filtered <- mat[rowSums(mat)> 0,]
  return(m.filtered)
}


# get the barcode expression matrix for a single mouse
# mat is the barcode expression matrix
# mouse.id is the mouse you want to analyse
getMouseMatrix <- function(mat, mouse.id){
  
  mouse.columns <- colnames(mat)[grepl(mouse.id,colnames(mat))]
  mouse <- data.frame(mat[,c(mouse.columns)])
  m <- removeUninformativeClones(mouse)
  colnames(m) <- c("B","M","E")
  return(m)
}



cell_counts <- read.csv('/Users/jasoncosgrove/Dropbox (Team_Perie)/Jason/Experiments/Wet_lab/JCW24_G6PD_barcoding/LentiviralBarcoding_Analysis/cell_counts.csv')

#normalise the counts by cell numbers in each sample
cellcountNormalisation <- function(mouse, mouse.mat){
  
  mat.cellnorm <- mouse.mat
  

  mat.cellnorm[,1] <- mouse.mat[,1] * cell_counts[cell_counts$Mouse == mouse  & cell_counts$Gate == 'B',6]
  mat.cellnorm[,2] <- mouse.mat[,2] * cell_counts[cell_counts$Mouse == mouse  & cell_counts$Gate == 'M',6]
  mat.cellnorm[,3] <- mouse.mat[,3] * cell_counts[cell_counts$Mouse == mouse  & cell_counts$Gate == 'E',6]
  
  colnames(mat.cellnorm) <- c("B","M","E")
  return(mat.cellnorm)
  
}


#merge CD62lhi and CD62llo into a single representative sample: LSKs
mergeMPPs <- function(mat){
  mat$lsk <- rowSums(mat[,1:2])
  return(mat)
}



#combine the expression matrices for each single mouse into a single metamouse sample
# so we have increased statistical power to look at rare barcodes. 
makeMetaMouse <- function(m1,m2,m3,m4){
  
  
  #rename columns so they are consistent across datasets
  colnames(m1) <- c("B","M","E")
  colnames(m2) <- c("B","M","E")
  colnames(m3) <- c("B","M","E")
  colnames(m4) <- c("B","M","E")
  
  #now make an integrated sample with all mice
  supermouse <- rbind(m1, m2,m3,m4)
  
  #clones must be in either the pos or the neg
  supermouse <- supermouse[rowSums(supermouse) > 0,]
  
  return(supermouse)
}



plotCloneSizes <- function(wt.m1,wt.m2,wt.m3,wt.m4,
                           tg.m1,tg.m2,tg.m3,tg.m4,lineage){
  
  
  wt.b <- c(median(wt.m1[,lineage]),
            median(wt.m2[,lineage]),
            median(wt.m3[,lineage]),
            median(wt.m4[,lineage]))
  
  tg.b <- c(median(tg.m1[,lineage]),
            median(tg.m2[,lineage]),
            median(tg.m3[,lineage]),
            median(tg.m4[,lineage]))
  
  p <- boxplot(wt.b, tg.b,
               main = paste("lineage: ",lineage,' \n p-value: ',
                            wilcox.test(wt.b,tg.b)$p.value),
               col = c("red","blue"), names= c('WT','Tg'))
  
}



performRowNormalisation <- function(matrix){
  
  matrix <- matrix[rowSums(matrix) > 0,]
  
  for(i in 1:nrow(matrix)){
    if(rowSums(matrix[i,c("B","M","E")]) > 0){
      matrix[i,c("B","M","E")] <-   
        matrix[i,c("B","M","E")]/sum(matrix[i,c("B","M","E")])
    }
  }
  
  return(matrix)
}



#classifier
multiOutcomeClassifier <- function(mat, threshold = 0.1){
  
  bias <- matrix("",nrow = nrow(mat), ncol = 4)
  colnames(bias) <- c("B","M","E","class")
  
  for(i in 1:nrow(mat)){
    
    if(mat[i,1] > threshold){
      bias[i,1] <- "B"
    }
    if(mat[i,2] > threshold){
      bias[i,2] <- "M"
    }
    if(mat[i,3] > threshold){
      bias[i,3] <- "E"
    }
    
    bias[i,4] <- paste(bias[i,1], bias[i,2] , bias[i,3], sep = "")
    
  }
  
  bias <- data.frame(bias)
  bias$class <- factor(bias$class, 
                          levels = c("B","BE",'BM',"BME","M","ME","E"))
  return(list(df = bias, table = table(bias$class)))
  
}



plotClassifierLineage <- function(lineage,threshold){
  
  wt <-c(multiOutcomeClassifier(m1.wt.cellnorm.norm, threshold = threshold)$table[[lineage]],
         multiOutcomeClassifier(m2.wt.cellnorm.norm, threshold = threshold)$table[[lineage]],
         multiOutcomeClassifier(m3.wt.cellnorm.norm, threshold = threshold)$table[[lineage]],
         multiOutcomeClassifier(m4.wt.cellnorm.norm, threshold = threshold)$table[[lineage]])
  
  tg <-c(multiOutcomeClassifier(m1.tg.cellnorm.norm, threshold = threshold)$table[[lineage]],
         multiOutcomeClassifier(m2.tg.cellnorm.norm, threshold = threshold)$table[[lineage]],
         multiOutcomeClassifier(m3.tg.cellnorm.norm, threshold = threshold)$table[[lineage]],
         multiOutcomeClassifier(m4.tg.cellnorm.norm, threshold = threshold)$table[[lineage]])
  
  
  boxplot(wt,tg, 
          main = paste('classifier bias: ',lineage,"\n p.val: ",
                       wilcox.test(wt,tg)$p.value),
          col = c("red","blue"),
          names = c("WT","Tg"),outline = T)
  
  
}





generateCumulativeProbability <- function(barcodematrix,sample){
  
  barcodematrix$data.normalised <- barcodematrix[,sample] 
  
  #rank each clone in terms of normalised abundance
  data.ranked <- barcodematrix[order(barcodematrix$data.normalised),]
  
  output <- data.frame(cbind((data.ranked[,sample]),
                             cumsum(data.ranked$data.normalised)))
  
  rownames(output) <- rownames(data.ranked)
  colnames(output) <- c("cell_counts","cumsum")
  return(output)
  
}



generateCumulativeProbabilityNorm <- function(barcodematrix,sample){
  
  #renormalise the clones to 1 for the myeloid lineage
  barcodematrix$data.normalised <- barcodematrix[,sample] / sum(barcodematrix[,sample])
  
  
  #rank each clone in terms of normalised abundance
  data.ranked <- barcodematrix[order(barcodematrix$data.normalised),]
  
  output <- data.frame(cbind((data.ranked[,sample]),
                             cumsum(data.ranked$data.normalised)))
  colnames(output) <- c("cell_counts","cumsum")
  rownames(output) <- rownames(data.ranked)
  return(output)
  
}
