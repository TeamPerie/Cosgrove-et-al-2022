########################################
# METABOLOMICS ANALYSIS
# Author: Jason Cosgrove (jason.cosgrove@curie.fr)
# Date: 10/04/2023 
###########################################


### Prepare the workspace
rm(list=ls())
library(pheatmap)
library("heatmaply")
library(dendextend)
library(RColorBrewer)
library(umap)
library(reshape2)

### helper methods

# Define a function called label_color_function that takes a label as input and returns a color code
label_color_function <- function(label) {
  # Check if the label contains the string "CD62L_HIGH"
  if (grepl("CD62L_HIGH", label)) {
    return("#56B4E9") # Set the color code to blue
  } 
  else if (grepl("CD62L_high", label)) {
    return("#56B4E9") # Set the color code to blue
  } 
  else if (grepl("CD62L_NEG", label)) {
    return("#999999") # Set the color code to gray
  } 
  else if (grepl("CMP", label)) {
    return("#E69F00") # Set the color code to orange
  } 
  else if (grepl("GMP", label)) {
    return("#009E73") # Set the color code to green
  } 
  else if (grepl("HSC", label)) {
    return("#CC79A7") # Set the color code to pink
  } 
  # If none of the above conditions are met, set the color code to black
  else {
    return("black")
  }
  
}



#Define a function called plotDendrogram that takes a dataframe, 
#distance metric, and clustering method as inputs
plotDendrogram <- function(df, dist = "euclidean", clust = "complete") {
  
  # Center and scale the rows of the dataframe to match the preprocessing of pheatmap
  df_scaled <- t(scale(t(df), center = TRUE, scale = TRUE))
  
  # Compute a dendrogram of the columns of the scaled dataframe using the specified distance metric and clustering method
  col_dendrogram <- as.dendrogram(hclust(dist(t(df_scaled), method = dist), method = clust))
  
  # Apply a label color function to the dendrogram labels
  label_colors <- sapply(labels(col_dendrogram), label_color_function)
  
  # Set the label colors of the dendrogram to the computed label colors
  labels_colors(col_dendrogram) <- paste(label_colors)
  
  # Set the margin parameters of the plot
  par(mar = c(12, 8, 8, 8))
  
  # Plot the dendrogram of the columns with a smaller font size
  plot(col_dendrogram, cex = 1.5, edgePar = list( lwd = 2.5))
  
}


# Define a function called PCA_ElbowPlot that takes a PCA result object as input
PCA_ElbowPlot <- function(pca_result) {
  
  # Compute the variance explained by each principal component
  variance_explained <- pca_result$sdev^2
  
  # Compute the total variance
  total_variance <- sum(variance_explained)
  
  # Compute the proportion of variance explained by each principal component
  prop_variance_explained <- variance_explained / total_variance
  
  # Compute the cumulative proportion of variance explained by each principal component
  cumulative_prop_variance_explained <- cumsum(prop_variance_explained)
  
  # Plot the cumulative proportion of variance explained by each principal component
  plot(cumulative_prop_variance_explained, xlab = "Number of Principal Components", 
       ylab = "Cumulative Proportion of Variance Explained", type = "b", main = "Elbow Plot")
  
}

# Define a function called plotPCA that takes a PCA result object as input
plotPCA <- function(pca) {
  
  # Convert the PCA result object to a data frame with the first two principal components
  pca_df <- as.data.frame(pca$x[, 1:2])
  
  # Add a grouping variable to the data frame based on sample type
  pca_df$group <- c(rep("HSC", 5), rep("CD62L_HIGH_MPP", 6), rep("CD62L_NEG_MPP", 6),
                    rep("MEP", 6), rep("GMP", 6), rep("CMP", 6))
  
  
  my_palette <- c("#56B4E9", "#999999", "#E69F00", "#009E73", "#CC79A7", "#000000", "#D55E00")
  
  # Create a scatter plot of the first two principal components, with points colored by group
  ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
    geom_point(size = 6) +
    scale_color_manual(values = my_palette) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),  # Remove light grey grid
      panel.background = element_rect(fill = "white", color = "black", size = 4),  # Set background color and border size
      #axis.line = element_line(color = "black", size = 2),  # Set axis line color and size
      axis.title = element_text(face = "bold", size = 12),
      axis.text = element_text(face = "bold", size = 14),
      legend.title = element_blank(),
      legend.text = element_text(face = "bold", size = 14)
    ) +
    labs(
      x = "Principal Component 1",
      y = "Principal Component 2",
      title = "PCA Plot Showing Similarity of LC-MS/MS profiles"
    )
  
}


plotUMAP <- function(pca,dims = 15){
  pca_scores <- pca$x[, 1:dims]
  
  umap_result <- umap(pca_scores, n_neighbors = 9, n_components = 2, metric = "euclidean")
  
  umap_df <- as.data.frame(umap_result$layout)
  
  umap_df$group <- c(rep("HSC",5), rep("CD62L_HIGH_MPP",6),rep("CD62L_NEG_MPP",6),
                     rep("MEP",6),rep("GMP",6),rep("CMP",6))
  
  my_palette <- c("#56B4E9", "#999999", "#E69F00", "#009E73", "#CC79A7", "#000000", "#D55E00")
  
  ggplot(umap_df, aes(x = V1, y = V2, color = group)) +
    geom_point(size = 6) +
    scale_color_manual(values = my_palette) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),  # remove light grey grid
      panel.background = element_rect(fill = "white", color = "black", size = 4),  # set background color and border size
      #axis.line = element_line(color = "black", size = 2),  # set axis line color and size
      axis.title = element_text(face = "bold", size = 12),
      axis.text = element_text(face = "bold", size = 14),
      legend.title = element_blank(),
      legend.text = element_text(face = "bold", size = 14)
    ) +
    labs(
      x = "UMAP 1",
      y = "UMAP 2",
      title = "UMAP Projection of LC-MS/MS Profiles",
    )
  
}


# Define a function called differentialExpression that takes a dataframe and two cell type names as inputs
differentialExpression <- function(df, celltype1 = "CD62L_high", celltype2 = "CD62L_NEG") {
  
  #HSCs dont have a mouse 3 so remove it from other celltypes or 
  # we cant do paired statistical tests
  if(celltype1 == "HSC" | celltype2 == "HSC"){
    df <- df[,!grepl("M3", colnames(df))]
  }
  
  
  # Get the column indices for the two specified cell types
  c1 <- which(grepl(celltype1, colnames(df)))
  c2 <- which(grepl(celltype2, colnames(df)))
  
  # Combine the column indices into a single vector
  c <- c(c1, c2)
  
  # Compute differential expression between the two cell types for each row of the dataframe
  out <- matrix(0, ncol = 1, nrow = nrow(df))
  colnames(out) <- "p_value"
  rownames(out) <- rownames(df)
  for (i in 1:nrow(df)) {
    if (shapiro.test(as.numeric(df[i, c]))$p.value < 0.05) {
      # Use a Wilcoxon signed-rank test if the data is not normally distributed
      out[i,] <- wilcox.test(as.numeric(df[i, c1]), as.numeric(df[i, c2]), paired = TRUE)$p.value
    } else {
      # Use a paired t-test if the data is normally distributed
      out[i,] <- t.test(as.numeric(df[i, c1]), as.numeric(df[i, c2]), paired = TRUE)$p.value
    }
  }
  
  # Return a matrix of p-values for differential expression between the two cell types
  return(as.data.frame(out))
}


setwd("/Users/jasoncosgrove/Dropbox (Team_Perie)/Jason/Experiments/Wet_lab/JCW45_mass_spec/ANALYSIS3")

#define a color palette for plotting
my_palette <- c( "#56B4E9", "#999999", "#E69F00","#009E73", "#CC79A7", "#000000" ,  "#D55E00" )


### read the data
df <- read.csv("DELTA_DEBRIS_ZEROED.csv",
               row.names =1)

df.median <- read.csv("DELTA_DEBRIS_ZEROED_MEDIAN.csv",
               row.names = 1)



### Preprocessing of the data

#remove any columns which are zero
#remove any rows which are zero
df <- df[, colSums(df) != 0]
df <- df[rowSums(df) != 0,]
df.median <- df.median[, colSums(df.median) != 0]
df.median <- df.median[rowSums(df.median) != 0,]


# remove samples with low detection or if they are a sorting control

df.median <- df.median[rownames(df.median) != "ATA",]
df.median <- df.median[rownames(df.median) != "CLF",]
df <- df[rownames(df) != "ATA",]
df <- df[rownames(df) != "CLF",]

#remove samples that are completely zero in one experiment

e1 <- df[,colnames(df)[grepl("M1|M2|M3",colnames(df))]]
e2 <- df[,colnames(df)[grepl("M4|M5|M6",colnames(df))]]

to.delete <- rownames(df)[ rowSums(e1) == 0 | rowSums(e2) == 0]

'%ni%' <- function(x,y)!('%in%'(x,y))

df.median <- df.median[rownames(df.median) %ni% to.delete,]
df <- df[rownames(df) %ni% to.delete,]

### Plot metabolite recovery

#Plot the total values for each sample to see if there are global differences 
#in the total amount of recovered metabolites
barplot(colSums(df),las = 2)
barplot(colSums(df.median),las = 2)



### hierarchical clustering of the data

plotDendrogram(df)
plotDendrogram(df.median)


### pheatmap visualisation
pheatmap(df.median,scale = "row",color = viridis(n = 50,  option = "magma"))


scale_limits <- seq(-2, 2, by = 0.1)
pheatmap(df,scale = "row",color = viridis(n = 50,  option = "magma")
         , breaks = scale_limits,clustering_method = "complete",fontsize =8)




# Create a heatmap using heatmaply with scaled rows, the viridis color palette, and clustering method
heatmaply(df, scale = "row", colors = viridis(n = 256, option = "magma"), 
          margins = c(50, 50, 50, 50),
          clustering_method = "complete",cexRow = 0.7, cexCol = 0.7)

# Create a heatmap using heatmaply with scaled rows, the viridis color palette, and custom row scaling factor
heatmaply(df.median, scale = "row", colors = viridis(n = 256, option = "magma"), 
          margins = c(50, 50, 50, 50),cexRow = 0.8,
          clustering_method = "complete")






### Dimensionality reduction plots
pca_result <- prcomp(t(df), scale. = TRUE)

PCA_ElbowPlot(pca_result)

plotPCA(pca_result)

plotUMAP(pca_result,dims =15)








######stats testing

#test each metabolite to see whether they are normally distributed or not

out <- matrix(0, ncol = 1, nrow = nrow(df))
for(i in 1:nrow(df)){
  out[i,] <- shapiro.test(as.numeric(df[i,6:17]))$p.value
}
table(out < 0.05)

#as we can see most metabolites are normally distributed when you look across both all samples, 
#but we are comparing just high and low MPPs and in this case many of them are not normally distributed


mpps <- differentialExpression(df,celltype1 = "CD62L_high", celltype2 = "CD62L_NEG")
mpps$p.adjust <- p.adjust(mpps$p_value,method = "BH")

hsc.mpp.neg <- differentialExpression(df,celltype1 = "HSC", celltype2 = "CD62L_NEG")
hsc.mpp.neg$p.adjust <- p.adjust(hsc.mpp.neg$p_value,method = "BH")

hsc.mpp.pos <- differentialExpression(df,celltype1 = "HSC", celltype2 = "CD62L_high")
hsc.mpp.pos$p.adjust <- p.adjust(hsc.mpp.pos$p_value,method = "BH")

#which metabolites are significantly different?
hit.list <- unique(c(rownames(mpps[mpps$p.adjust < 0.05,]),
                     rownames(hsc.mpp.neg[hsc.mpp.neg$p.adjust  < 0.05,]),
                     rownames(hsc.mpp.pos[hsc.mpp.pos$p.adjust  < 0.05,])))











