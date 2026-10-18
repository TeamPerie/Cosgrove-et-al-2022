rm(list=ls())

#library(reticulate)
library(SeuratDisk)
library(Seurat)
library(SeuratData)



setwd("/Users/jasoncosgrove/Dropbox (Team_Perie)/Jason/Experiments/Dry_Lab/JCB16_LoCelso_Plasmodium_Infection/wetransfer-2a31f3/")



#convert python format to seurat
#ad <- import("anndata", convert = FALSE)
#pbmc_ad <- ad$read_h5ad("malaria_control_data_genespector.h5ad")
#pbmc_ad$write("malaria_control_data_genespector.gzip.h5ad", compression="gzip")

#pbmc_ad <- ad$read_h5ad("malaria_infected_data_genespector.h5ad")
#pbmc_ad$write("malaria_infected_data_genespector.gzip.h5ad", compression="gzip")




Convert("malaria_control_data_genespector.gzip.h5ad", dest = "h5seurat", overwrite = TRUE)
data.control <- LoadH5Seurat("malaria_control_data_genespector.gzip.h5seurat")

Convert("malaria_infected_data_genespector.gzip.h5ad", dest = "h5seurat", overwrite = TRUE)
data.infected <- LoadH5Seurat("malaria_infected_data_genespector.gzip.h5seurat")

load("/Users/jasoncosgrove/Dropbox (Team Perié)/Jason/Experiments/JCB4_metabolism_bioinformatics_pipeline/genesets/geneset_Robjects/metabolic_signatures.Rda")






data.control@meta.data$condition <- rep("ctrl",ncol(data.control))
data.infected@meta.data$condition <- rep("infected",ncol(data.infected))

control.meta.data <- read.csv("control_metadata.csv")
infected.meta.data <- read.csv("infected_metadata.csv")

data.control@meta.data$clusters <- control.meta.data$clusters
data.control@meta.data$clusters <- control.meta.data$clusters

data.infected@meta.data$clusters <- infected.meta.data$clusters
data.infected@meta.data$sample <- infected.meta.data$sample



##### lets analyse the samples separately first

data.control <- SetIdent(data.control, value = data.control@meta.data$clusters)
data.infected <- SetIdent(data.infected, value = data.infected@meta.data$clusters)


sig <- read.csv(file.choose())
sig <- sig$x
data.control<- AddModuleScore(data.control, features =  list(sig),name = "metabolic_sig")
data.infected <- AddModuleScore(data.infected, features =  list(sig),name = "metabolic_sig")


VlnPlot(data.infected, "metabolic_sig1",group.by = "clusters")
VlnPlot(data.control, "metabolic_sig1",group.by = "clusters")

write.csv(data.control@meta.data$metabolic_sig1, "controlsig.csv")
write.csv(data.infected@meta.data$metabolic_sig1, "infectedsig.csv")

#you need to make sure that something is stored in the metadata befre you integrate or it doesnt like it so much...


#data.control.hspc <- subset(x = data.control, idents = "3 Myeloid")
#data.infected.hspc <- subset(x = data.infected, idents = "3 Myeloid")

dataset.list <- list()
dataset.list[[1]] <- data.control 
dataset.list[[2]] <- data.infected


dataset.list  <- lapply(X = dataset.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 5000)
})

features <- SelectIntegrationFeatures(object.list = dataset.list)

anchors <- FindIntegrationAnchors(object.list = dataset.list, anchor.features = features)


data.combined <- IntegrateData(anchorset = anchors)


DefaultAssay(data.combined) <- "integrated"

save(data.combined,file ="integratedData.Rda")

# Run the standard workflow for visualization and clustering
data.combined <- ScaleData(data.combined, verbose = FALSE)
data.combined <- RunPCA(data.combined, npcs = 30, verbose = FALSE)
ElbowPlot(data.combined,ndims = 12)
DimPlot(data.combined,group.by = "clusters",dims = c(1,2), cols = c("grey","grey","grey","grey",
                                                                   "grey","grey","grey"))

data.combined <- AddModuleScore(data.combined, features =  list(sig),name = "metab",assay = "RNA")
data.combined <- AddModuleScore(data.combined, features =  list(metabolic.signatures$sequential_myeloid),
                                name = "metab_old",assay = "RNA")

FeaturePlot(data.combined, features = "metab1",max.cutoff = "q90", min.cutoff = "q3" )
FeaturePlot(data.combined, features = "metab_old1",max.cutoff = "q90", min.cutoff = "q3" )

VlnPlot(data.combined, features = "metab1",group.by = "condition" ,pt.size = 0)
VlnPlot(data.combined, features = "metab_old1",group.by = "condition" ,pt.size = 0)

x <- VlnPlot(data.combined,c("metab1"),idents = "4 Primitive HSPC", group.by = "condition")

write.csv(x$data,"/users/jasoncosgrove/Desktop/metafateHSPCs.csv")



createDensityPlot <- function(sobj){
  
  tmp.all<-as.data.frame(Embeddings(object = sobj, reduction = "pca"))
  p <- ggplot(tmp.all, aes(x = PC_1, y = PC_2)) + geom_point(colour="#00000000") + stat_density_2d(aes(fill = stat(level)), geom = "polygon", bins=50) + scale_fill_gradientn(colors = c("#4169E100","royalblue", "darkolivegreen3","goldenrod1","red")) + theme_classic() + 
    theme(legend.position="bottom") + NoAxes()
  
  return(p)
}

Idents(data.combined) <- "condition"
ctrl <- subset(data.combined, idents = c("ctrl"))
infected <- subset(data.combined, idents = c("infected"))

createDensityPlot(ctrl)
createDensityPlot(infected)

VlnPlot(data.combined,c("Gsr","G6pdx","Pkm"),idents = "4 Primitive HSPC", group.by = "condition")


write.csv(VlnPlot(data.combined,"erythroid_sig1",idents = "4 Primitive HSPC", group.by = "condition")$data,"/Users/jasoncosgrove/Desktop/e_sig.csv")
write.csv(VlnPlot(data.combined,"myeloid_sig1",idents = "4 Primitive HSPC" ,group.by = "condition")$data,"/Users/jasoncosgrove/Desktop/m_sig.csv")


DefaultAssay(data.combined) <- "RNA"
markers <- FindConservedMarkers(data.combined, ident.1 = "3 Myeloid", grouping.var = "condition", verbose = FALSE)

FeaturePlot(data.combined, c("Mpo","Ltb","Klf1","Sell"),min.cutoff = "q3", max.cutoff = "q97",cols = c("grey90","grey20"))

FeaturePlot(data.combined, c("myeloid_sig1","erythroid_sig1","both_sig1"),min.cutoff = "q3", max.cutoff = "q97")
data.combined <- RunUMAP(data.combined, reduction = "pca", dims = 1:12)
DimPlot(data.combined,group.by = "condition",dims = c(1,2))
FeaturePlot(data.combined, "Sell",dims = c(1,2))



Idents(data.combined.subset ) <- "clusters"
data.combined.subset <- subset(x = data.combined, idents = "4 Primitive HSPC")
Idents(data.combined.subset ) <- "condition"
avg.dat<- as.data.frame(log1p(AverageExpression(data.combined.subset , verbose = FALSE)$RNA))
avg.dat$gene <- rownames(avg.dat)

avg.dat.subset <- avg.dat[avg.dat$ctrl > 0.1,]
avg.dat.subset$ratio <- avg.dat.subset$infected / avg.dat.subset$ctrl
avg.dat.subset <- avg.dat.subset[order(avg.dat.subset$ratio),]


genes.to.label = avg.dat.subset$gene[avg.dat.subset$ratio > 1.5]
p1 <- ggplot(avg.dat, aes(ctrl, infected)) + geom_point() + ggtitle("HSPCs")
p1 <- LabelPoints(plot = p1, points = genes.to.label,repel = TRUE)



Idents(data.combined) <- "clusters"
Idents(data.combined) <- factor(Idents(data.combined), 
                                  levels = c( "0 MEP","1 Basophil","3 Myeloid","4 Primitive HSPC", 
                                              "5 Megakaryocyte","2 Erythroid","0 Empty"))

DotPlot(data.combined, features = "Sell", cols = c("red", "red"), dot.scale = 8, split.by = "condition") +
  RotatedAxis()



VlnPlot(data.combined,"Sell",group.by = "condition",assay = "RNA",idents = "4 Primitive HSPC")

VlnPlot(data.combined,"Sell",group.by = "condition",assay = "RNA")


Idents(ctrl) <- "clusters"
Idents(infected) <- "clusters"
ctrl.hspc <- subset(ctrl, idents = c("4 Primitive HSPC"))
infected.hspc <- subset(infected, idents = c("4 Primitive HSPC"))





gene <- "G6pdx"
barplot(table(data.combined@meta.data$condition,
              data.combined@assays$RNA@data[gene,] > 0), legend = T)

df <- table(data.combined@meta.data$condition,
            data.combined@assays$RNA@data[gene,] > 0)

barplot(t(t(df) /colSums(df)), legend = T)



hspc <- subset(data.combined, idents = c("4 Primitive HSPC"))


gene <- "Gsr"
barplot(table(hspc@meta.data$condition,
              hspc@assays$RNA@data[gene,] > 0), legend = T)

df <- table(hspc@meta.data$condition,
            hspc@assays$RNA@data[gene,] > 0)

barplot(t(t(df) /colSums(df)), legend = T)

table(ctrl@assays$RNA@data[gene,] > 0)
table(infected@assays$RNA@data[gene,] > 0)

table(ctrl.hspc@assays$RNA@data[gene,] > 0)
table(infected.hspc@assays$RNA@data[gene,] > 0)

FeaturePlot(data.combined,"Sell")
VlnPlot(data.combined, "metabolic_sig1",group.by = "condition")




plots <- VlnPlot(data.combined, features = c("both_sig1"), split.by = "condition", group.by = "clusters",
                 pt.size = 0, combine = FALSE,assay = "RNA")
wrap_plots(plots = plots, ncol = 1)


RidgePlot(data.combined,features = "Sell",group.by = "condition",ident = "4 Primitive HSPC",assay = "RNA")

sig <- metabolic.signatures$sequential_myeloid
data.combined <- AddModuleScore(data.combined, features =  list(sig),name = "myeloid_sig",assay = "RNA")

sig <- metabolic.signatures$sequential_erythroid
data.combined <- AddModuleScore(data.combined, features =  list(sig),name = "erythroid_sig",assay = "RNA")

sig <- c(metabolic.signatures$sequential_erythroid,metabolic.signatures$sequential_myeloid)
data.combined <- AddModuleScore(data.combined, features =  list(sig),name = "both_sig",assay = "RNA")

FeaturePlot(data.infected, "metabolic_sig1",reduction = "draw_graph_",
            min.cutoff = "q3", max.cutoff = "q97")

data.control <- AddModuleScore(data.control, features =  list(sig),name = "metabolic_sig")

FeaturePlot(data.control, "metabolic_sig1",reduction = "draw_graph_",
            min.cutoff = "q3", max.cutoff = "q97")
