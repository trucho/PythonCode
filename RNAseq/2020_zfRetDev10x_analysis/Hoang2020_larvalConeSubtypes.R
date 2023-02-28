# Attempt at separating cone subtypes in larval samples from Hoang 2020.
# Does not work well
# Requires some subsetting of original data which is part of Hoang2020_reanalysis.R

# Load Libraries ------------------------------------------------------------------------
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
# Setup -------------------------------------------------------------------
setwd(".code/R_SeuratAnalysis/")
directory = ".code/R_SeuratAnalysis/"
exportDir = paste(directory,"Hoang2020_10x_photoreceptors",sep="")
getwd()
# Plot themes -------------------------------------------------------------------
plotTheme = function (base_size = 42, base_family = "") {
   theme_classic() %+replace% 
      theme(
         axis.line = element_line(colour = 'black', size = 1),
         axis.text = element_text(size=18),
         text = element_text(size=18)
      )
}

# -------------------------------------------------------------------
# Load the dataset
photo_larval = readRDS("./cones_LarvalReClustered.rds");


# Subtypes in late stage cones: difficult, maybe can separate M-leaning cones vs L-leaning cones

Cll = subset(photo_larval, idents = "Cll")
# Identification of highly variable features (feature selection)
Cll <- FindVariableFeatures(Cll, selection.method = "vst", nfeatures = 50)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Cll), 200)
top10

#run PCA
Cll <- RunPCA(Cll, features = VariableFeatures(object = Cll))
# Use elbow plot (var explained)
ElbowPlot(Cll, ndims=50) + plotTheme()

# Cluster cells
# nFeatures = 100, ndim = 30, res = 1 -> no separation
# nFeatures = 250, ndim = 35, res = 1 -> no separation, except for rods in cluster 5; 0 might be M-cones?
# nFeatures = 1000, ndim = 20, res = 1 -> not good either
# nFeatures = 20, ndim = 15, res = 1 -> forces separation on opsins giving a mixed cluster and a pure L
# nFeatures = 50, ndim = 30, res = 1 -> no separation, except for rods in cluster 5
upDimLimit=30;
lsl <- FindNeighbors(lsl, dims = 1:upDimLimit)
lsl <- FindClusters(lsl, resolution = 1)

DimPlot(lsl, reduction = "pca", label=TRUE)

# UMAP using the same PCA dimensions for prettier visualization
lsl <- RunUMAP(lsl, dims = 1:upDimLimit)
DimPlot(lsl, reduction = "umap", label=TRUE, pt.size = 2, label.size = 8)

# or use tSNE (clustering can't separate some M cones from L cones)
lsl <- RunTSNE(lsl, dims = 1:upDimLimit)
DimPlot(lsl, reduction = "tsne", label=TRUE,  pt.size = 1, label.size = 6) + plotTheme()
# ------------------------------------------------------------------------------------------
# plot layouts for direct export
#OPSINS
p = FeaturePlot(lsl, reduction = 'tsne', features = c("rho","nr2e3","gnat1","gnat2","arr3b","arr3a","opn1sw1", "opn1sw2",'opn1mw1','opn1lw2'),
                pt.size=1, order=TRUE, combine=FALSE)
lay = rbind(c(1,2,3,NA),c(4,5,6,NA),c(7,8,9,10))
grid.arrange(grobs = p, layout_matrix = lay)
# ps = arrangeGrob(grobs = p, layout_matrix = lay)
# ggsave(ps, file="larval_lslMarkers.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")
# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
# Renaming to separate rods
new.cluster.ids <- c("lslC","lslC","lslC","lslC","lslR","lslC","lslC")
names(new.cluster.ids) <- levels(lsl)
lsl <- RenameIdents(lsl, new.cluster.ids)
DimPlot(lsl, reduction = "tsne", label = TRUE, pt.size = 1, label.size = 6) + plotTheme()

# ------------------------------------------------------------------------------------------
# Exploring msl briefly
# seems like L and M cones do have their own clusters, with a subdivision for cells that are farther along development with more expression of phototransduction proteins
# this is what PC2 is picking up on
msl =  subset(photo, idents = "mslPR")

# Identification of highly variable features (feature selection)
msl <- FindVariableFeatures(msl, selection.method = "vst", nfeatures = 200)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(msl), 200)
top10


#run PCA
msl <- RunPCA(msl, features = VariableFeatures(object = msl))
# Use elbow plot (var explained)
ElbowPlot(msl, ndims=50)

# Cluster cells
upDimLimit=30;
msl <- FindNeighbors(msl, dims = 1:upDimLimit)
msl <- FindClusters(msl, resolution = 1)

DimPlot(msl, reduction = "pca", label=TRUE)

# UMAP using the same PCA dimensions for prettier visualization
msl <- RunUMAP(msl, dims = 1:upDimLimit)
DimPlot(msl, reduction = "umap", label=TRUE, pt.size = 2, label.size = 8)

# or use tSNE (clustering can't separate some M cones from L cones)
msl <- RunTSNE(msl, dims = 1:upDimLimit)
DimPlot(msl, reduction = "tsne", label=TRUE, pt.size = 1, label.size = 6) + plotTheme()


FeaturePlot(msl, reduction = 'tsne', features = c("rho", "opn1sw1", "opn1sw2",'opn1mw1','opn1lw2','si:busm1-57f23.1'))
FeaturePlot(msl, reduction = 'tsne', features = c("gnat1","gnat2","opn1sw2","opn1lw2","arr3a","arr3b"))
# seems like L and M cones do have their own clusters, with a subdivision for cells that are farther along development with more expression of phototransduction proteins
FeaturePlot(msl, reduction = 'tsne', features = c("rho","nr2e3","gnat1","gnat2","arr3a","arr3b","opn1sw1", "opn1sw2",'opn1mw1','opn1lw2','thrb'))

# ------------------------------------------------------------------------------------------
# plot layouts for direct export
#Subtype markers
p = FeaturePlot(msl, reduction = 'tsne', features = c("rho","nr2e3","gnat1","gnat2","arr3b","arr3a","thrb","opn1sw1", "opn1sw2",'opn1mw1','opn1lw2'),
                pt.size=1, order=TRUE, combine=FALSE)
lay = rbind(c(1,2,3,NA),c(4,5,6,7),c(8,9,10,11))
grid.arrange(grobs = p, layout_matrix = lay)
# ps = arrangeGrob(grobs = p, layout_matrix = lay)
# ggsave(ps, file="larval_mslMarkers.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")
# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
# Exploring esl briefly
esl =  subset(photo_larval, idents = c("eslPR","eslS"))
# Identification of highly variable features (feature selection)
esl <- FindVariableFeatures(esl, selection.method = "vst", nfeatures = 500)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(esl), 200)
top10
#run PCA
esl <- RunPCA(esl, features = VariableFeatures(object = esl))
# Use elbow plot (var explained)
ElbowPlot(esl, ndims=50)

# Cluster cells
upDimLimit=20;
esl <- FindNeighbors(esl, dims = 1:upDimLimit)
esl <- FindClusters(esl, resolution = 1)

DimPlot(esl, reduction = "pca", label=TRUE)

# UMAP using the same PCA dimensions for prettier visualization
esl <- RunUMAP(esl, dims = 1:upDimLimit)
DimPlot(esl, reduction = "umap", label=TRUE, pt.size = 2, label.size = 8)
# or use tSNE (clustering can't separate some M cones from L cones)
esl <- RunTSNE(esl, dims = 1:upDimLimit)
DimPlot(esl, reduction = "tsne", label=TRUE, pt.size = 1, label.size = 6) + plotTheme()


FeaturePlot(esl, reduction = 'tsne', features = c("rho", "opn1sw1", "opn1sw2",'opn1mw1','opn1lw2','si:busm1-57f23.1'))
FeaturePlot(esl, reduction = 'tsne', features = c('tbx2a','tbx2b','efna1b','foxq2',"ntf3","nr2e3","thrb"))
FeaturePlot(esl, reduction = 'tsne', features = c('tbx2a','tbx2b','foxq2',"nr2e3","thrb","skor1a",'lrrfip1a','sall1a'))
DotPlot(esl, features = c('tbx2a','tbx2b','foxq2',"nr2e3","thrb","skor1a",'lrrfip1a','sall1a'))
# seems like L and M cones do have their own clusters still and there is a S cluster which is foxq2+
# no gnat yet, nr2e3+ everywhere, six7 is present almost in all clusters, 

# Ploting by age reveals that L cone cluster corresponds to 5dpf sample, only a few cones come from 2dpf (too early) and 2.5 - 3.5 dpf are intermingled
esl2 = esl
esl2 = SetIdent(esl2, value = "Sample")
DimPlot(esl2, reduction = "tsne", label = TRUE, pt.size = 2, label.size = 8)
# ------------------------------------------------------------------------------------------
# plot layouts for direct export
#Subtype markers
p = FeaturePlot(esl, reduction = 'tsne', features = c("rho","nr2e3","gnat1","gnat2","arr3b","arr3a","thrb","opn1sw1", "opn1sw2",'opn1mw1','opn1lw2'),
                pt.size=1, order=TRUE, combine=FALSE)
lay = rbind(c(1,2,3,NA),c(4,5,6,7),c(8,9,10,11))
grid.arrange(grobs = p, layout_matrix = lay)
# ps = arrangeGrob(grobs = p, layout_matrix = lay)
# ggsave(ps, file="larval_eslMarkers.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")
# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
# CONCLUSION: not enough information to pull out photoreceptor subtypes except early S cones. 
# USEFUL INFORMATION: some developmental timelines could be derived from clusters

# remerging to be able to keep rod and cone clusters in lsl
photo_larval = merge(x = lsl, y = c(msl,esl))

# 'pr' : "#E6B800"
# 'lslpr' : "#cca819"
# 'mslpr' : "#dcc360"
# 'eslpr' : "#dacd9a"
# 'esls' : '#8f9bcc',
# 'r' : '#747474',
# 'lr' : '#474747'

# ordering identities for plotting and assigning colors
# Idents(photoLarval) <- factor(x = Idents(photoLarval), levels = c("eslPR","eslL","eslM","eslS","mslPR","mslL","mslM","lslPR","R"))
# lcolors = c("#dacd9a","#d69f9e","#80bc89","#8f9bcc","#dcc360","#c67271","#57cb69","#cca819","#474747")

saveRDS(photo_larval, file = "./photo_LarvalSubtypes.rds")
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# -------------------------------------------------------------------
