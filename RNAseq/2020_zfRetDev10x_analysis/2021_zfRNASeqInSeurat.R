Sys.setenv(RETICULATE_PYTHON = "/Users/angueyraaristjm/anaconda3/bin/python")

reticulate::use_python("/Users/angueyraaristjm/anaconda3/bin/python")
# reticulate::use_python("/usr/local/bin/python3")
reticulate::use_virtualenv()
reticulate::py_install(packages = 'umap-learn')


library(Seurat)
library(dplyr)
library(patchwork)
library(reticulate)
library(umap)
library(ggplot2)
library(gridExtra)
library(grid)

# Clear all plots -------------------------------------------------------------------
try(dev.off(dev.list()["RStudioGD"]),silent=TRUE)
try(dev.off(),silent=TRUE)

# Clear environment -------------------------------------------------------------------
rm(list=ls())
# Plot themes -------------------------------------------------------------------
eelTheme = function (base_size = 42, base_family = "") {
   theme_classic() %+replace% 
      theme(
         axis.line = element_line(colour = 'black', size = 1),
         axis.text = element_text(size=18),
         text = element_text(size=18)
      )
}
# Setup -------------------------------------------------------------------
setwd("/Users/angueyraaristjm/Documents/LiMolec/zfRNAseq/20190827/")
directory = "/Users/angueyraaristjm/Documents/LiMolec/zfRNAseq/20190827/"
exportDir = paste(directory,"SeuratAnalysis",sep="")
getwd()

gCounts<-read.table(file=paste0("/Users/angueyraaristjm/Documents/LiMolec/zfRNAseq/20190827/20190827_DESeq2/","gCount.csv"),sep=",", header = TRUE)
rownames(gCounts) = gCounts[,1]
gCounts = gCounts[,-1]
head(gCounts)

# manually picked photoreceptors
photoMP <- CreateSeuratObject(counts = gCounts, min.cells = 3, min.genes = 100, project = "zfRNAseq")
photoMP <- NormalizeData(photoMP, normalization.method = "LogNormalize", scale.factor = 10000)
photoMP <- ScaleData(photoMP, features = rownames(photoMP))
photoMP


# -------------------------------------------------------------------
# -------------------------------------------------------------------
photoMP$subtype = c("L","L","L","L","L","L",
                     "M","M","M","M","M","M",
                     "R","R","R","R","R","R",
                     "S","S","S","S","S","S","M", #correction for mislabeled sample during collection
                     "UV","UV","UV","UV","UV")
Idents(object = photoMP) = photoMP$subtype
Idents(photoMP) <- factor(x = Idents(photoMP), levels = c("R","UV","S","M","L"))

# -------------------------------------------------------------------
# Identification of highly variable features (feature selection)
photoMP <- FindVariableFeatures(object = photoMP, selection.method = "vst", nfeatures = 500)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(photoMP), 100)
top10

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(photoMP)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#run PCA
upDimLimit=20
photoMP <- RunPCA(photoMP, features = VariableFeatures(object = photoMP), npcs=upDimLimit)
# # Examine and visualize PCA results a few different ways
# print(photo[["pca"]], dims = 1:5, nfeatures = 5)
# VizDimLoadings(photo, dims = 1:2, reduction = "pca")
# DimPlot(photo, reduction = "pca")
# DimHeatmap(photo, dims = 1, cells = 500, balanced = TRUE)
# DimHeatmap(photo, dims = 1:15, cells = 500, balanced = TRUE)


# Determine number of clusters using Macosko, 2015 (random permutation of 1% of data and rerun PCA iteritavely)
# Or just use elbow plot (var explained)
ElbowPlot(photoMP, ndims=upDimLimit) + eelTheme()
# ggsave("larval_Elbow.png", path=exportDir, width = 140, height = 105, units = "mm")

# Cells don't need reclustering (one of the main reasons why we did manual picking)

DimPlot(photoMP, reduction = "pca", label = TRUE, pt.size = 1, label.size = 6) + eelTheme()
# ggsave("larval_PCAinitial.png", path=exportDir, width = 140, height = 105, units = "mm")

# UMAP using the same PCA dimensions for prettier visualization
photoMP <- RunUMAP(photoMP, dims = 1:upDimLimit)
DimPlot(photoMP, reduction = "umap", label = TRUE, pt.size = 6, label.size = 6) + eelTheme()
# ggsave("larval_UMAPinitial.png", path=exportDir, width = 140, height = 105, units = "mm")

# or use tSNE (clustering can't separate some M cones from L cones)
photoMP <- RunTSNE(photoMP, dims = 1:upDimLimit, perplexity = 8)
DimPlot(photoMP, reduction = "tsne", label = TRUE, pt.size = 6, label.size = 6) + eelTheme()
# ggsave("larval_TSNEinitial.png", path=exportDir, width = 140, height = 105, units = "mm")

# FeaturePlot(photo, reduction = 'tsne', features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1','efna1b'))
FeaturePlot(photoMP, reduction = 'tsne', features = c("rho","opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2'))
FeaturePlot(photoMP, reduction = 'tsne', features = c("rho","nrl","pde6a","nr2e3","crx",'gnat1',"saga","sagb","gucy2f","grk1a")) # rod markers
FeaturePlot(photoMP, reduction = 'tsne', features = c("rho","eno2","rom1a","lrrn1","unc119.2","pdca","cplx4c")) # rod genes from my own dataset
FeaturePlot(photoMP, reduction = 'tsne', features = c("gnat2",'arr3a','arr3b',"pde6c","pde6h","guca1d","grk7a","crx","neurod1","nr2f6b")) #cone markers
FeaturePlot(photoMP, reduction = 'tsne', features = c("nr2e3","nr2f6b","crx",'syt5a','syt5b','gnat2','arr3a','arr3b',"pde6g","pde6h", "neurod1")) #cone markers
FeaturePlot(photoMP, reduction = 'tsne', features = c("slc1a8b","rgs9a","slc25a24","ppa1a","sema7a","kera","nexn","dusp5","crhbp","prph2a")) #cone markers from my own dataset
FeaturePlot(photoMP, reduction = 'tsne', features = c("otx5","tbx2a","tbx2b","rxrga","rxrgb", "thrb", "sema7a","cnga3a","cnga3b"))
VlnPlot(photoMP, features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))
FeaturePlot(photoMP, reduction = 'tsne', features = c("nrl","mafa","mafb")) # rod markers

FeaturePlot(photoMP, reduction = 'tsne', features = c("crx","nr2e3","nr2f6a","nr2f6b",'syt5a','syt5b','gnat2',"pde6g","pde6h", "neurod1")) #dev markers


DotPlot(photoMP, features = c('tbx2a','tbx2b','efna1b','foxq2',"ntf3","nr2e3","nr2f6b","thrb"))

saveRDS(photoMP, file = "./zfAPhotoreceptors.rds")

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Running DESeq2 (subtype unique DEGs) and saving as csv

tempC = c("UV","S","M","L")
nC = length(tempC);
temp = subset(photoMP, idents = tempC)

# L cones
DESeq2_L = FindMarkers(temp, c("L"), c("UV","S","M"), test.use="DESeq2");
AvgExp_L = AverageExpression(temp, features=rownames(DESeq2_L))
AvgExp_L = AvgExp_L$RNA
DESeq2_L = merge.data.frame(DESeq2_L,AvgExp_L, by="row.names", sort=FALSE)
DESeq2_L$baseMean = rowMeans(DESeq2_L[,tempC])
rownames(DESeq2_L) = DESeq2_L$Row.names
DESeq2_L = DESeq2_L[,c("baseMean","avg_logFC","p_val","p_val_adj","pct.1","pct.2","UV","S","M","L")]
PctExp_L = DotPlot(temp, features=rownames(DESeq2_L))
PctNames = rownames(PctExp_L$data)[1:(length(rownames(PctExp_L$data))/nC)]
PctColNames = paste("pct",levels(PctExp_L$data$id),sep="")
PctExp_L = PctExp_L$data$pct.exp
PctExp_L = matrix(PctExp_L,nrow = length(PctExp_L)/nC, ncol=nC);
rownames(PctExp_L) = PctNames
colnames(PctExp_L) = PctColNames
DESeq2_L = merge.data.frame(DESeq2_L,PctExp_L, by="row.names", sort=FALSE)
rownames(DESeq2_L) = DESeq2_L$Row.names
DESeq2_L = DESeq2_L[,c("baseMean","avg_logFC","p_val","p_val_adj","pct.1","pct.2","UV","S","M","L","pctUV","pctS","pctM","pctL")]
head(DESeq2_L[order(-DESeq2_L$avg_logFC),],40)
write.csv(DESeq2_L,"./DESeq2_L.csv")

# M cones
DESeq2_M = FindMarkers(temp, c("M"), c("UV","S","L"), test.use="DESeq2");
AvgExp_M = AverageExpression(temp, features=rownames(DESeq2_M))
AvgExp_M = AvgExp_M$RNA
DESeq2_M = merge.data.frame(DESeq2_M,AvgExp_M, by="row.names", sort=FALSE)
DESeq2_M$baseMean = rowMeans(DESeq2_M[,tempC])
rownames(DESeq2_M) = DESeq2_M$Row.names
DESeq2_M = DESeq2_M[,c("baseMean","avg_logFC","p_val","p_val_adj","pct.1","pct.2","UV","S","M","L")]
PctExp_M = DotPlot(temp, features=rownames(DESeq2_M))
PctNames = rownames(PctExp_M$data)[1:(length(rownames(PctExp_M$data))/nC)]
PctColNames = paste("pct",levels(PctExp_M$data$id),sep="")
PctExp_M = PctExp_M$data$pct.exp
PctExp_M = matrix(PctExp_M,nrow = length(PctExp_M)/nC, ncol=nC);
rownames(PctExp_M) = PctNames
colnames(PctExp_M) = PctColNames
DESeq2_M = merge.data.frame(DESeq2_M,PctExp_M, by="row.names", sort=FALSE)
rownames(DESeq2_M) = DESeq2_M$Row.names
DESeq2_M = DESeq2_M[,c("baseMean","avg_logFC","p_val","p_val_adj","pct.1","pct.2","UV","S","M","L","pctUV","pctS","pctM","pctL")]
head(DESeq2_M[order(-DESeq2_M$avg_logFC),],40)
write.csv(DESeq2_M,"./DESeq2_M.csv")

# S cones
DESeq2_S = FindMarkers(temp, c("S"), c("UV","M","L"), test.use="DESeq2");
AvgExp_S = AverageExpression(temp, features=rownames(DESeq2_S))
AvgExp_S = AvgExp_S$RNA
DESeq2_S = merge.data.frame(DESeq2_S,AvgExp_S, by="row.names", sort=FALSE)
DESeq2_S$baseMean = rowMeans(DESeq2_S[,tempC])
rownames(DESeq2_S) = DESeq2_S$Row.names
DESeq2_S = DESeq2_S[,c("baseMean","avg_logFC","p_val","p_val_adj","pct.1","pct.2","UV","S","M","L")]
PctExp_S = DotPlot(temp, features=rownames(DESeq2_S))
PctNames = rownames(PctExp_S$data)[1:(length(rownames(PctExp_S$data))/nC)]
PctColNames = paste("pct",levels(PctExp_S$data$id),sep="")
PctExp_S = PctExp_S$data$pct.exp
PctExp_S = matrix(PctExp_S,nrow = length(PctExp_S)/nC, ncol=nC);
rownames(PctExp_S) = PctNames
colnames(PctExp_S) = PctColNames
DESeq2_S = merge.data.frame(DESeq2_S,PctExp_S, by="row.names", sort=FALSE)
rownames(DESeq2_S) = DESeq2_S$Row.names
DESeq2_S = DESeq2_S[,c("baseMean","avg_logFC","p_val","p_val_adj","pct.1","pct.2","UV","S","M","L","pctUV","pctS","pctM","pctL")]
head(DESeq2_S[order(-DESeq2_S$avg_logFC),],40)
write.csv(DESeq2_S,"./DESeq2_S.csv")


# UV cones
DESeq2_U = FindMarkers(temp, c("UV"), c("S","M","L"), test.use="DESeq2");
AvgExp_U = AverageExpression(temp, features=rownames(DESeq2_U))
AvgExp_U = AvgExp_U$RNA
DESeq2_U = merge.data.frame(DESeq2_U,AvgExp_U, by="row.names", sort=FALSE)
DESeq2_U$baseMean = rowMeans(DESeq2_U[,tempC])
rownames(DESeq2_U) = DESeq2_U$Row.names
DESeq2_U = DESeq2_U[,c("baseMean","avg_logFC","p_val","p_val_adj","pct.1","pct.2","UV","S","M","L")]
PctExp_U = DotPlot(temp, features=rownames(DESeq2_U))
PctNames = rownames(PctExp_U$data)[1:(length(rownames(PctExp_U$data))/nC)]
PctColNames = paste("pct",levels(PctExp_U$data$id),sep="")
PctExp_U = PctExp_U$data$pct.exp
PctExp_U = matrix(PctExp_U,nrow = length(PctExp_U)/nC, ncol=nC);
rownames(PctExp_U) = PctNames
colnames(PctExp_U) = PctColNames
DESeq2_U = merge.data.frame(DESeq2_U,PctExp_U, by="row.names", sort=FALSE)
rownames(DESeq2_U) = DESeq2_U$Row.names
DESeq2_U = DESeq2_U[,c("baseMean","avg_logFC","p_val","p_val_adj","pct.1","pct.2","UV","S","M","L","pctUV","pctS","pctM","pctL")]
head(DESeq2_U[order(-DESeq2_U$avg_logFC),],40)
write.csv(DESeq2_U,"./DESeq2_U.csv")
