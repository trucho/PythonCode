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

# Clear all plots -------------------------------------------------------------------
try(dev.off(dev.list()["RStudioGD"]),silent=TRUE)
try(dev.off(),silent=TRUE)

# Clear environment -------------------------------------------------------------------
rm(list=ls())

# Setup -------------------------------------------------------------------
setwd("/Users/angueyraaristjm/Documents/LiLab/RNAseq/zfRet_HoangBlackshaw2020/")
directory = "/Users/angueyraaristjm/Documents/LiLab/RNAseq/zfRet_HoangBlackshaw2020/"
getwd()
# -------------------------------------------------------------------
# Load the 10x dataset (28845 cells), after updating to Seurat_v3 (had to use biowulf)
pbmc = readRDS("~/Documents/LiLab/RNAseq/zfRet_HoangBlackshaw2020/zfDev_pbmc_v3.rds");
# pca, tsne and umap already done, so will keep their clusters and just separate photoreceptors
# ------------------------------------------------------------------------------------------

# Retinal progenitor: 0
# Photoreceptor precursor: 3
# Cones: adults = 4,  larval = 12
# Rods: 2
# HCs: 6
# BCs: 5, 15
# Immature MG: 16
# MG: 11, 9 , 1
# RGCs = 13 + 10
# ngng ACs = 7
# gly ACs = 14
# GABA ACs = 8

# Assigning cell type identity to clusters (it's basically just renaming)
new.cluster.ids <- c("RPC","MG1","R","PRPC","Ca","BC1","HC","AC","ACgaba","MG2","RGC1","MG3","Cl","RGC2","ACgly","BC2","MGi")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
# ------------------------------------------------------------------------------------------
# RODS
# There are 2 rod clusters, but it's not age since >95% on cells are from adult samples
rods <- subset(pbmc, idents = c("R"))
# saveRDS(rods, file = "./rods.rds") #preliminary saving to analyze later

DimPlot(rods, reduction = "tsne", label=TRUE)
zfAFlag = grepl("AdR",names(Idents(rods)));
rods = AddMetaData(rods,metadata=zfAFlag, col.name = "zfA");
rods = SetIdent(rods, value = "zfA")
FeaturePlot(rods, reduction = 'tsne', features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))
# smaller cluster has higher nr2e3 levels
FeaturePlot(rods, reduction = 'tsne', features = c("rho", "nr2e3",'nrl','saga','sagb'))

# ------------------------------------------------------------------------------------------
# CONES
# Subclustering cones (probably will then remerge with rods and photoreceptor progenitors to track gene expression)
# DimPlot(pbmc, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()
# almost impossible separation by subtypes combining all "cones"
photo <- subset(pbmc, idents = c("Ca","Cl"))
DimPlot(photo, reduction = "tsne", label=TRUE)
# clear separation between adult and larval cones
zfAFlag = grepl("AdR",names(Idents(photo)));
photo = AddMetaData(photo,metadata=zfAFlag, col.name = "zfA");
photo = SetIdent(photo, value = "zfA")
#Ca cluster still has a couple (n=2) cells from E120(indices 21 and 1577), so it's better to subset by flag
# photo = subset(pbmc, idents = "Ca") 
# ------------------------------------------------------------------------------------------

# Using only adult cone cluster
## Problem with this is that the Variable Features seem to be highly expressed genes in all subtypes, such that clustering is driven by this rather thatn cone identity

photo_larval = subset(photo, idents = FALSE)
# saveRDS(photo_larval, file = "./cones_Larval.rds") #preliminary saving to analyze later

photo = subset(photo, idents = TRUE)

# Identification of highly variable features (feature selection)
photo <- FindVariableFeatures(photo, selection.method = "vst", nfeatures = 150)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(photo), 100)
top10

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(photo)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#run PCA
photo <- RunPCA(photo, features = VariableFeatures(object = photo))
# # Examine and visualize PCA results a few different ways
# print(photo[["pca"]], dims = 1:5, nfeatures = 5)
# VizDimLoadings(photo, dims = 1:2, reduction = "pca")
# DimPlot(photo, reduction = "pca")
# DimHeatmap(photo, dims = 1, cells = 500, balanced = TRUE)
# DimHeatmap(photo, dims = 1:15, cells = 500, balanced = TRUE)


# Determine number of clusters using Macosko, 2015 (random permutation of 1% of data and rerun PCA iteritavely)
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
photo <- JackStraw(photo, num.replicate = 100)
photo <- ScoreJackStraw(photo, dims = 1:20)

JackStrawPlot(photo, dims = 1:20)
# Or just use elbow plot (var explained)
ElbowPlot(photo, ndims=50)

# Cluster cells
# nFeatures = 1000 or 2000 muddles everything, especially it mixes L and M cones
# nFeatures = 200 looks better
# nFeatures = 100 seems like not enough info
# nFeatures = 150 apparently good compromise
# dims = 8, resolution = 0.6 -> not bad but still some mixing
# dims = 23, resolution = 0.6 -> almost there; there are still some M-cones mixed with L-cones. UV and S have to be separated further

# Settling for nFeatures = 150, ndim = 23, res = 1
upDimLimit=23;
photo <- FindNeighbors(photo, dims = 1:upDimLimit)
photo <- FindClusters(photo, resolution = 1)

DimPlot(photo, reduction = "pca", label = TRUE, pt.size = 1)

# UMAP using the same PCA dimensions for prettier visualization
photo <- RunUMAP(photo, dims = 1:upDimLimit)
DimPlot(photo, reduction = "umap", label = TRUE, pt.size = 1)
# saveRDS(photo, file = "../photoreceptors_umap.rds")

# or use tSNE (clustering can't separate some M cones from L cones)
photo <- RunTSNE(photo, dims = 1:upDimLimit)
DimPlot(photo, reduction = "tsne", label = TRUE, pt.size = 1)
# saveRDS(photo, file = "../photoreceptors_tSNE.rds")

FeaturePlot(photo, reduction = 'tsne', features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))
FeaturePlot(photo, reduction = 'tsne', features = c("rho", "nr2e3","nr2f6","crx",'syt5a','syt5b'))
VlnPlot(photo, features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))


# Assigning cone subtype identity to clusters
new.cluster.ids <- c("L","L","L","L","L","L","L","L","US","M3","M4","M1","Ribo","L")
names(new.cluster.ids) <- levels(photo)
photo <- RenameIdents(photo, new.cluster.ids)
DimPlot(photo, reduction = "tsne", label = TRUE, pt.size = 2) + NoLegend()
DimPlot(photo, reduction = "umap", label = TRUE, pt.size = 1) + NoLegend()

FeaturePlot(photo, reduction = "tsne", features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))
VlnPlot(photo, features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))
# ------------------------------------------------------------------------------------------
# Attempting to separate UV and S cones

us = subset(photo, idents = "US")
# Identification of highly variable features (feature selection)
us <- FindVariableFeatures(us, selection.method = "vst", nfeatures = 100)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(us), 200)
top10


#run PCA
us <- RunPCA(us, features = VariableFeatures(object = us))
# Use elbow plot (var explained)
ElbowPlot(us, ndims=50)

# Cluster cells
# nFeatures = 100, ndim = 8, res = 1 -> 3 clusters; ok separation but UV cluster still has a lost of s cones
upDimLimit=12;
us <- FindNeighbors(us, dims = 1:upDimLimit)
us <- FindClusters(us, resolution = 1)

DimPlot(us, reduction = "pca", label=TRUE)

# UMAP using the same PCA dimensions for prettier visualization
us <- RunUMAP(us, dims = 1:upDimLimit)
DimPlot(us, reduction = "umap", label=TRUE)
# saveRDS(photo, file = "../photoreceptors_umap.rds")

# or use tSNE (clustering can't separate some M cones from L cones)
us <- RunTSNE(us, dims = 1:upDimLimit)
DimPlot(us, reduction = "tsne", label=TRUE)
# saveRDS(photo, file = "../photoreceptors_tSNE.rds")

FeaturePlot(us, reduction = 'tsne', features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))
VlnPlot(us, features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))
FeaturePlot(us, reduction = 'tsne', features = c('tbx2a','tbx2b','efna1b','foxq2'))

new.cluster.ids <- c("S","S","S","UV")
names(new.cluster.ids) <- levels(us)
us <- RenameIdents(us, new.cluster.ids)
DimPlot(us, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()

# pull out other clusters for DESeq2 analysis
lm = subset(photo, idents = c("L","M1","M3","M4","CPC"));
# remerge into single cone dataset (1894 cones)
photo = merge(x = lm, y = us)
photo

photo = ScaleData(photo)
photo <- FindVariableFeatures(photo, selection.method = "vst", nfeatures = 40)
photo <- RunPCA(photo, features = VariableFeatures(object = photo))
ElbowPlot(photo, ndims=50)
upDimLimit=16;
DimPlot(photo, reduction = "pca", label=TRUE)
# UMAP using the same PCA dimensions for prettier visualization
photo <- RunUMAP(photo, dims = 1:upDimLimit)
DimPlot(photo, reduction = "umap", label=TRUE, pt.size = 1) + NoLegend()
# or use tSNE (clustering can't separate some M cones from L cones)
photo <- RunTSNE(photo, dims = 1:upDimLimit)
DimPlot(photo, reduction = "tsne", label=TRUE, pt.size = 1) + NoLegend()

FeaturePlot(photo, reduction = 'umap', features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))
VlnPlot(photo, features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))
DotPlot(photo, features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2'))
saveRDS(photo, file = "./cones_AdultAll.rds")
# ------------------------------------------------------------------------------------------
# Merging M1 and M3 and kicking out Ribo cluster and M4 for DESeq2analysis. Will consider these 2 clusters separately
new.cluster.ids <- c("Ribo","L","M1","M3","M4","S","UV")
names(new.cluster.ids) <- levels(photoAll)
photo <- RenameIdents(photo, new.cluster.ids)
DimPlot(photo, reduction = "tsne", label=TRUE, pt.size = 1) + NoLegend()
photo = subset(photo, idents = c("UV","S","M","L"))


# photo = ScaleData(photo)
photo <- FindVariableFeatures(photo, selection.method = "vst", nfeatures = 80)
photo <- RunPCA(photo, features = VariableFeatures(object = photo))
ElbowPlot(photo, ndims=50)
upDimLimit=12;
DimPlot(photo, reduction = "pca", label=TRUE)
# UMAP using the same PCA dimensions for prettier visualization
photo <- RunUMAP(photo, dims = 1:upDimLimit)
DimPlot(photo, reduction = "umap", label=TRUE, pt.size = 1) + NoLegend()
# or use tSNE (clustering can't separate some M cones from L cones)
photo <- RunTSNE(photo, dims = 1:upDimLimit)
DimPlot(photo, reduction = "tsne", label=TRUE, pt.size = 1) + NoLegend()

FeaturePlot(photo, reduction = 'umap', features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))
VlnPlot(photo, features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))
DotPlot(photo, features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2'))
saveRDS(photo, file = "./cones_Adult.rds")


# ------------------------------------------------------------------------------------------
try(dev.off(dev.list()["RStudioGD"]),silent=TRUE)
try(dev.off(),silent=TRUE)
rm(list=ls())
# ------------------------------------------------------------------------------------------
# Can be restarted here
photo = readRDS("~/Documents/LiLab/RNAseq/zfRet_HoangBlackshaw2020/cones_Adult.rds");
# ------------------------------------------------------------------------------------------
# Running DESeq2 (all relevant comparisons) and saving as csv
# L-cones
DESeq2_L = FindMarkers(photo, c("L"), c("M","S","UV"), test.use="DESeq2");
AvgExp_L = AverageExpression(photo, features=rownames(DESeq2_L))
AvgExp_L = AvgExp_L$RNA
DESeq2_L = merge.data.frame(DESeq2_L,AvgExp_L, by="row.names", sort=FALSE)
DESeq2_L$baseMean = rowMeans(DESeq2_L[,c("L","M","S","UV")])
rownames(DESeq2_L) = DESeq2_L$Row.names
DESeq2_L = DESeq2_L[,c("baseMean","avg_logFC","p_val","p_val_adj","pct.1","pct.2","L","M","S","UV")]
PctExp_L = DotPlot(photo, features=rownames(DESeq2_L))
PctNames = rownames(PctExp_L$data)[1:(length(rownames(PctExp_L$data))/4)]
PctExp_L = PctExp_L$data$pct.exp
PctExp_L = matrix(PctExp_L,nrow = length(PctExp_L)/4, ncol=4);
rownames(PctExp_L) = PctNames
colnames(PctExp_L) = c("pctL","pctM","pctS","pctUV")
DESeq2_L = merge.data.frame(DESeq2_L,PctExp_L, by="row.names", sort=FALSE)
rownames(DESeq2_L) = DESeq2_L$Row.names
write.csv(DESeq2_L,"./DESeq2_L.csv")

# M-cones
DESeq2_M = FindMarkers(photo, c("M"), c("L","S","UV"), test.use="DESeq2");
AvgExp_M = AverageExpression(photo, features=rownames(DESeq2_M))
AvgExp_M = AvgExp_M$RNA
DESeq2_M = merge.data.frame(DESeq2_M,AvgExp_M, by="row.names", sort=FALSE)
DESeq2_M$baseMean = rowMeans(DESeq2_M[,c("L","M","S","UV")])
rownames(DESeq2_M) = DESeq2_M$Row.names
DESeq2_M = DESeq2_M[,c("baseMean","avg_logFC","p_val","p_val_adj","pct.1","pct.2","L","M","S","UV")]
PctExp_M = DotPlot(photo, features=rownames(DESeq2_M))
PctNames = rownames(PctExp_M$data)[1:(length(rownames(PctExp_M$data))/4)]
PctExp_M = PctExp_M$data$pct.exp
PctExp_M = matrix(PctExp_M,nrow = length(PctExp_M)/4, ncol=4);
rownames(PctExp_M) = PctNames
colnames(PctExp_M) = c("pctL","pctM","pctS","pctUV")
DESeq2_M = merge.data.frame(DESeq2_M,PctExp_M, by="row.names", sort=FALSE)
rownames(DESeq2_M) = DESeq2_M$Row.names
write.csv(DESeq2_M,"./DESeq2_M.csv")

# S-cones
DESeq2_S = FindMarkers(photo, c("S"), c("L","M","UV"), test.use="DESeq2");
AvgExp_S = AverageExpression(photo, features=rownames(DESeq2_S))
AvgExp_S = AvgExp_S$RNA
DESeq2_S = merge.data.frame(DESeq2_S,AvgExp_S, by="row.names", sort=FALSE)
DESeq2_S$baseMean = rowMeans(DESeq2_S[,c("L","M","S","UV")])
rownames(DESeq2_S) = DESeq2_S$Row.names
DESeq2_S = DESeq2_S[,c("baseMean","avg_logFC","p_val","p_val_adj","pct.1","pct.2","L","M","S","UV")]
PctExp_S = DotPlot(photo, features=rownames(DESeq2_S))
PctNames = rownames(PctExp_S$data)[1:(length(rownames(PctExp_S$data))/4)]
PctExp_S = PctExp_S$data$pct.exp
PctExp_S = matrix(PctExp_S,nrow = length(PctExp_S)/4, ncol=4);
rownames(PctExp_S) = PctNames
colnames(PctExp_S) = c("pctL","pctM","pctS","pctUV")
DESeq2_S = merge.data.frame(DESeq2_S,PctExp_S, by="row.names", sort=FALSE)
rownames(DESeq2_S) = DESeq2_S$Row.names
write.csv(DESeq2_S,"./DESeq2_S.csv")


# UV-cones
DESeq2_UV = FindMarkers(photo, c("UV"), c("L","M","S"), test.use="DESeq2");
AvgExp_UV = AverageExpression(photo, features=rownames(DESeq2_UV))
AvgExp_UV = AvgExp_UV$RNA
DESeq2_UV = merge.data.frame(DESeq2_UV,AvgExp_UV, by="row.names", sort=FALSE)
DESeq2_UV$baseMean = rowMeans(DESeq2_UV[,c("L","M","S","UV")])
rownames(DESeq2_UV) = DESeq2_UV$Row.names
DESeq2_UV = DESeq2_UV[,c("baseMean","avg_logFC","p_val","p_val_adj","pct.1","pct.2","L","M","S","UV")]
PctExp_UV = DotPlot(photo, features=rownames(DESeq2_UV))
PctNames = rownames(PctExp_UV$data)[1:(length(rownames(PctExp_UV$data))/4)]
PctExp_UV = PctExp_UV$data$pct.exp
PctExp_UV = matrix(PctExp_UV,nrow = length(PctExp_UV)/4, ncol=4);
rownames(PctExp_UV) = PctNames
colnames(PctExp_UV) = c("pctL","pctM","pctS","pctUV")
DESeq2_UV = merge.data.frame(DESeq2_UV,PctExp_UV, by="row.names", sort=FALSE)
rownames(DESeq2_UV) = DESeq2_UV$Row.names
write.csv(DESeq2_UV,"./DESeq2_UV.csv")


#UV&S vs L&M
DESeq2_USvsLM = FindMarkers(photo, c("M","L"), c("UV","S"), test.use="DESeq2");
AvgExp_USvsLM = AverageExpression(photo, features=rownames(DESeq2_USvsLM))
AvgExp_USvsLM = AvgExp_USvsLM$RNA
DESeq2_USvsLM = merge.data.frame(DESeq2_USvsLM,AvgExp_USvsLM, by="row.names", sort=FALSE)
DESeq2_USvsLM$baseMean = rowMeans(DESeq2_USvsLM[,c("L","M","S","UV")])
rownames(DESeq2_USvsLM) = DESeq2_USvsLM$Row.names
DESeq2_USvsLM = DESeq2_USvsLM[,c("baseMean","avg_logFC","p_val","p_val_adj","pct.1","pct.2","L","M","S","UV")]
PctExp_USvsLM = DotPlot(photo, features=rownames(DESeq2_USvsLM))
PctNames = rownames(PctExp_USvsLM$data)[1:(length(rownames(PctExp_USvsLM$data))/4)]
PctExp_USvsLM = PctExp_USvsLM$data$pct.exp
PctExp_USvsLM = matrix(PctExp_USvsLM,nrow = length(PctExp_USvsLM)/4, ncol=4);
rownames(PctExp_USvsLM) = PctNames
colnames(PctExp_USvsLM) = c("pctL","pctM","pctS","pctUV")
DESeq2_USvsLM = merge.data.frame(DESeq2_USvsLM,PctExp_USvsLM, by="row.names", sort=FALSE)
rownames(DESeq2_USvsLM) = DESeq2_USvsLM$Row.names
write.csv(DESeq2_USvsLM,"./DESeq2_USvsLM.csv")

#L vs M
DESeq2_LvsM = FindMarkers(photo, c("L"), c("M"), test.use="DESeq2");
AvgExp_LvsM = AverageExpression(photo, features=rownames(DESeq2_LvsM))
AvgExp_LvsM = AvgExp_LvsM$RNA
DESeq2_LvsM = merge.data.frame(DESeq2_LvsM,AvgExp_LvsM, by="row.names", sort=FALSE)
DESeq2_LvsM$baseMean = rowMeans(DESeq2_LvsM[,c("L","M","S","UV")])
rownames(DESeq2_LvsM) = DESeq2_LvsM$Row.names
DESeq2_LvsM = DESeq2_LvsM[,c("baseMean","avg_logFC","p_val","p_val_adj","pct.1","pct.2","L","M","S","UV")]
PctExp_LvsM = DotPlot(photo, features=rownames(DESeq2_LvsM))
PctNames = rownames(PctExp_LvsM$data)[1:(length(rownames(PctExp_LvsM$data))/4)]
PctExp_LvsM = PctExp_LvsM$data$pct.exp
PctExp_LvsM = matrix(PctExp_LvsM,nrow = length(PctExp_LvsM)/4, ncol=4);
rownames(PctExp_LvsM) = PctNames
colnames(PctExp_LvsM) = c("pctL","pctM","pctS","pctUV")
DESeq2_LvsM = merge.data.frame(DESeq2_LvsM,PctExp_LvsM, by="row.names", sort=FALSE)
rownames(DESeq2_LvsM) = DESeq2_LvsM$Row.names
write.csv(DESeq2_LvsM,"./DESeq2_LvsM.csv")

#UV vs S
DESeq2_UvsS = FindMarkers(photo, c("S"), c("UV"), test.use="DESeq2");
AvgExp_UvsS = AverageExpression(photo, features=rownames(DESeq2_UvsS))
AvgExp_UvsS = AvgExp_UvsS$RNA
DESeq2_UvsS = merge.data.frame(DESeq2_UvsS,AvgExp_UvsS, by="row.names", sort=FALSE)
DESeq2_UvsS$baseMean = rowMeans(DESeq2_UvsS[,c("L","M","S","UV")])
rownames(DESeq2_UvsS) = DESeq2_UvsS$Row.names
DESeq2_UvsS = DESeq2_UvsS[,c("baseMean","avg_logFC","p_val","p_val_adj","pct.1","pct.2","L","M","S","UV")]
PctExp_UvsS = DotPlot(photo, features=rownames(DESeq2_UvsS))
PctNames = rownames(PctExp_UvsS$data)[1:(length(rownames(PctExp_UvsS$data))/4)]
PctExp_UvsS = PctExp_UvsS$data$pct.exp
PctExp_UvsS = matrix(PctExp_UvsS,nrow = length(PctExp_UvsS)/4, ncol=4);
rownames(PctExp_UvsS) = PctNames
colnames(PctExp_UvsS) = c("pctL","pctM","pctS","pctUV")
DESeq2_UvsS = merge.data.frame(DESeq2_UvsS,PctExp_UvsS, by="row.names", sort=FALSE)
rownames(DESeq2_UvsS) = DESeq2_UvsS$Row.names
write.csv(DESeq2_UvsS,"./DESeq2_UvsS.csv")

# ------------------------------------------------------------------------------------------
# Can be restarted here
photoAll = readRDS("~/Documents/LiLab/RNAseq/zfRet_HoangBlackshaw2020/cones_AdultAll.rds");
# ------------------------------------------------------------------------------------------
# Ribo vs rest
DESeq2_Ribo = FindMarkers(photoAll, c("Ribo"), c("L","M1","M3","M4","S","UV"), test.use="DESeq2");
AvgExp_Ribo = AverageExpression(photoAll, features=rownames(DESeq2_Ribo))
AvgExp_Ribo = AvgExp_Ribo$RNA
DESeq2_Ribo = merge.data.frame(DESeq2_Ribo,AvgExp_Ribo, by="row.names", sort=FALSE)
DESeq2_Ribo$baseMean = rowMeans(DESeq2_Ribo[,c("Ribo","L","M1","M3","M4","S","UV")])
rownames(DESeq2_Ribo) = DESeq2_Ribo$Row.names
DESeq2_Ribo = DESeq2_Ribo[,c("baseMean","avg_logFC","p_val","p_val_adj","pct.1","pct.2","Ribo","L","M1","M3","M4","S","UV")]
PctExp_Ribo = DotPlot(photoAll, features=rownames(DESeq2_Ribo))
PctNames = rownames(PctExp_Ribo$data)[1:(length(rownames(PctExp_Ribo$data))/7)]
PctExp_Ribo = PctExp_Ribo$data$pct.exp
PctExp_Ribo = matrix(PctExp_Ribo,nrow = length(PctExp_Ribo)/7, ncol=7);
rownames(PctExp_Ribo) = PctNames
colnames(PctExp_Ribo) = c("pctRibo","pctL","pctM1","pctM3","pctM4","pctS","pctUV")
DESeq2_Ribo = merge.data.frame(DESeq2_Ribo,PctExp_Ribo, by="row.names", sort=FALSE)
rownames(DESeq2_Ribo) = DESeq2_Ribo$Row.names
write.csv(DESeq2_Ribo,"./DESeq2_Ribo.csv")


# opn1mw4 cluster vs others
DESeq2_mw4 = FindMarkers(photoAll, c("M4"), c("L","M1","M3","S","UV"), test.use="DESeq2");
AvgExp_mw4 = AverageExpression(photoAll, features=rownames(DESeq2_mw4))
AvgExp_mw4 = AvgExp_mw4$RNA
AvgExp_mw4 = AvgExp_mw4[c("L","M1","M3","M4","S","UV")]
DESeq2_mw4 = merge.data.frame(DESeq2_mw4,AvgExp_mw4, by="row.names", sort=FALSE)
DESeq2_mw4$baseMean = rowMeans(DESeq2_mw4[,c("L","M1","M3","M4","S","UV")])
rownames(DESeq2_mw4) = DESeq2_mw4$Row.names
DESeq2_mw4 = DESeq2_mw4[,c("baseMean","avg_logFC","p_val","p_val_adj","pct.1","pct.2","L","M1","M3","M4","S","UV")]
PctExp_mw4 = DotPlot(photoAll, features=rownames(DESeq2_mw4))
PctNames = rownames(PctExp_mw4$data)[1:(length(rownames(PctExp_mw4$data))/7)]
PctExp_mw4 = PctExp_mw4$data$pct.exp
PctExp_mw4 = matrix(PctExp_mw4,nrow = length(PctExp_mw4)/7, ncol=7);
rownames(PctExp_mw4) = PctNames
colnames(PctExp_mw4) = c("pctRibo","pctL","pctM1","pctM3","pctM4","pctS","pctUV")
PctExp_mw4 = PctExp_mw4[,c("pctL","pctM1","pctM3","pctM4","pctS","pctUV")]
DESeq2_mw4 = merge.data.frame(DESeq2_mw4,PctExp_mw4, by="row.names", sort=FALSE)
rownames(DESeq2_mw4) = DESeq2_mw4$Row.names
write.csv(DESeq2_mw4,"./DESeq2_mw4.csv")

# M1 vs M3
DESeq2_M1vsM3 = FindMarkers(photoAll, c("M1"), c("M3"), test.use="DESeq2");
AvgExp_M1vsM3 = AverageExpression(photoAll, features=rownames(DESeq2_M1vsM3))
AvgExp_M1vsM3 = AvgExp_M1vsM3$RNA
AvgExp_M1vsM3 = AvgExp_M1vsM3[c("L","M1","M3","M4","S","UV")]
DESeq2_M1vsM3 = merge.data.frame(DESeq2_M1vsM3,AvgExp_M1vsM3, by="row.names", sort=FALSE)
DESeq2_M1vsM3$baseMean = rowMeans(DESeq2_M1vsM3[,c("L","M1","M3","M4","S","UV")])
rownames(DESeq2_M1vsM3) = DESeq2_M1vsM3$Row.names
DESeq2_M1vsM3 = DESeq2_M1vsM3[,c("baseMean","avg_logFC","p_val","p_val_adj","pct.1","pct.2","L","M1","M3","M4","S","UV")]
PctExp_M1vsM3 = DotPlot(photoAll, features=rownames(DESeq2_M1vsM3))
PctNames = rownames(PctExp_M1vsM3$data)[1:(length(rownames(PctExp_M1vsM3$data))/7)]
PctExp_M1vsM3 = PctExp_M1vsM3$data$pct.exp
PctExp_M1vsM3 = matrix(PctExp_M1vsM3,nrow = length(PctExp_M1vsM3)/7, ncol=7);
rownames(PctExp_M1vsM3) = PctNames
colnames(PctExp_M1vsM3) = c("pctRibo","pctL","pctM1","pctM3","pctM4","pctS","pctUV")
PctExp_M1vsM3 = PctExp_M1vsM3[,c("pctL","pctM1","pctM3","pctM4","pctS","pctUV")]
DESeq2_M1vsM3 = merge.data.frame(DESeq2_M1vsM3,PctExp_M1vsM3, by="row.names", sort=FALSE)
rownames(DESeq2_M1vsM3) = DESeq2_M1vsM3$Row.names
write.csv(DESeq2_M1vsM3,"./DESeq2_M1vsM3.csv")
# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
# cell counter
temp = subset(photoAll, idents = c("M1", "M3", "M4"))
length(names(Idents(temp)))
