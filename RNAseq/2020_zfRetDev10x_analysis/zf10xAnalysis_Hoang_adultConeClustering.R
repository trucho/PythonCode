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

# Setup -------------------------------------------------------------------
setwd("/Users/angueyraaristjm/Documents/LiMolec/otherRNAseq/zfRet_HoangBlackshaw2020/")
directory = "/Users/angueyraaristjm/Documents/LiMolec/otherRNAseq/zfRet_HoangBlackshaw2020/"
exportDir = paste(directory,"eelAnalysis",sep="")
getwd()
# Plot themes -------------------------------------------------------------------
eelTheme = function (base_size = 42, base_family = "") {
   theme_classic() %+replace% 
      theme(
         axis.line = element_line(colour = 'black', size = 1),
         axis.text = element_text(size=18),
         text = element_text(size=18)
      )
}
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Load the 10x dataset (28845 cells), after updating to Seurat_v3 (had to use biowulf)
pbmc = readRDS("~/Documents/LiMolec/otherRNAseq/zfRet_HoangBlackshaw2020/zfDev_pbmc_v3.rds");
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
new.cluster.ids <- c("RetinalProgenitor","Glia","Rod","PhotoProgenitor","Cone(adult)","BC","HC","AC","AC","Glia","RGC(adult)","Glia","Cone(larval)","RGC(larval)","AC","BC","Glia(immature)")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

DimPlot(pbmc, reduction = "tsne", label =  TRUE, pt.size = 1, label.size = 6) + eelTheme()
ggsave("adult_allClusters03.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 1, label.size = 6, repel = TRUE) + eelTheme()
ggsave("adult_allClusters.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

# store naming information to object
pbmc$cellTypes <- Idents(object = pbmc)
# switch to Sample to breakdown by age
pbmc = SetIdent(pbmc, value = "Sample")
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 1, label.size = 6) + eelTheme()
ggsave("adult_allClusters02.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

# DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 1, label.size = 6, repel = TRUE) + NoLegend() + eelTheme()

DimPlot(pbmc, reduction = "tsne", label =  TRUE, pt.size = 1, label.size = 6, repel = TRUE) + eelTheme()
ggsave("adult_allClusters04.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

# And switch back
pbmc = SetIdent(pbmc, value = "cellTypes")
# can't trust tbx2 data from this
FeaturePlot(pbmc, reduction = 'tsne', features = c("crx","otx5", "tbx2a",'tbx2b'), pt.size = 2)

FeaturePlot(pbmc, reduction = 'umap', features = c("crx","otx5", "tbx2a",'tbx2b'), pt.size = 2)
VlnPlot(pbmc, features = c("crx","otx5", "tbx2a",'tbx2b'))
# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
# RODS
# There are 2 rod clusters, and small one corresponds to larval samples + larval looking rods.
rods <- subset(pbmc, idents = c("Rod"))
# saveRDS(rods, file = "./rods.rds") #preliminary saving to analyze later

DimPlot(rods, reduction = "tsne", label =  TRUE, pt.size = 1, label.size = 6) + eelTheme()
ggsave("adult_rods01.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")
rods2 = SetIdent(rods, value = "Sample")
new.cluster.ids <- c("2.5dpf","5dpf","Ad1","2.5dpf","Ad5","4dpf","Ad4","Ad3","Ad2" )
names(new.cluster.ids) <- levels(rods2)
rods2 <- RenameIdents(rods2, new.cluster.ids)
Idents(rods) <- factor(x = Idents(rods2), levels = c("2.5dpf","4dpf","5dpf","Ad1","Ad2","Ad3","Ad4","Ad5"))
DimPlot(rods2, reduction = "tsne", label = TRUE, pt.size = 1, label.size = 6, repel = TRUE, order = c("2.5dpf","4dpf","5dpf","AdR1","AdR2","AdR3","AdR4","AdR5")) + eelTheme()
ggsave("adult_rods02.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")


zfAFlag = grepl("AdR",names(Idents(rods)));
rods = AddMetaData(rods,metadata=zfAFlag, col.name = "zfA");
rods = SetIdent(rods, value = "zfA")
FeaturePlot(rods, reduction = 'tsne', features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))
# smaller cluster has higher nr2e3 levels because it belongs to larval samples mostly
FeaturePlot(rods, reduction = 'tsne', features = c("rho", "nr2e3",'nrl','saga','sagb'))
# ------------------------------------------------------------------------------------------
# plot layouts for direct export
#Subtype markers
p = FeaturePlot(rods, reduction = 'tsne', features = c("rho","nr2e3","gnat1","gnat2","arr3b","arr3a","thrb","opn1sw1", "opn1sw2",'opn1mw1','opn1lw1'),
                pt.size=1, order=FALSE, combine=FALSE)
lay = rbind(c(1,2,3,NA),c(4,5,6,7),c(8,9,10,11))
grid.arrange(grobs = p, layout_matrix = lay)
ps = arrangeGrob(grobs = p, layout_matrix = lay)
ggsave(ps, file="adult_rodsMarkers.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

#Subtype markers in whole dataset
p = FeaturePlot(pbmc, reduction = 'tsne', features = c("rho","nr2e3","gnat1","gnat2","arr3b","arr3a","thrb","opn1sw1", "opn1sw2",'opn1mw1','opn1lw1'),
                pt.size=1, order=FALSE, combine=FALSE)
lay = rbind(c(1,2,3,NA),c(4,5,6,7),c(8,9,10,11))
grid.arrange(grobs = p, layout_matrix = lay)
ps = arrangeGrob(grobs = p, layout_matrix = lay)
ggsave(ps, file="adult_allMarkers.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

# HERE
ps = DotPlot(pbmc, features = c("rho","nr2e3","gnat1","gnat2","arr3b","arr3a","thrb","opn1sw1", "opn1sw2",'opn1mw1','opn1lw1')) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ggsave(ps, file="adult_allMarkersDot.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")


p = FeaturePlot(pbmc, reduction = 'tsne', features = c("crx","sox2","syt5a","syt5b","cnga1","cnga3a","elovl4b"),
                pt.size=1, order=FALSE, combine=FALSE)
lay = rbind(c(1,2,NA,NA),c(3,4,NA,NA),c(5,6,7,NA))
grid.arrange(grobs = p, layout_matrix = lay)
ps = arrangeGrob(grobs = p, layout_matrix = lay)
ggsave(ps, file="adult_allMarkers02.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

ps = DotPlot(pbmc, features =  c("crx","sox2","syt5a","syt5b","cnga1","cnga3a",'elovl4b')) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ggsave(ps, file="adult_allMarkersDot02.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

FeaturePlot(pbmc, reduction = 'tsne', features = c("crx","tbx2a","sox2","syt5a","syt5b","gnat1","gnat2",'elovl4b'),
            pt.size=1, order=FALSE, combine=TRUE)

FeaturePlot(pbmc, reduction = 'tsne', features = c("myo7aa","pcdh15a","cacna1fa","cacna1da","cacna1db"),
            pt.size=0.1, order=FALSE, combine=TRUE)

FeaturePlot(pbmc, reduction = 'umap', features = c("foxq2"), pt.size=1, order=TRUE, combine=FALSE)
FeaturePlot(pbmc, reduction = 'umap', features = c("foxq2","prdm1b","pbx1a","rxrga","nr2f6b"), pt.size=1, order=TRUE, combine=TRUE)

DotPlot(pbmc, features = c("myo7aa","pcdh15a","cacna1fa", "cacna1da","cacna1db","il34","p2ry1","nr2e3","crx","otx5"))

DotPlot(pbmc, features = c("cpne1","cpne2","cpne3", "cpne4a","cpne4b","cpne5a","cpne5b","cpne7","cpne8","cpne9"))

# HC genes in Yamagata, 2021 chicken RNAseq
DotPlot(pbmc, features = c("onecut1","onecut2","onecut3a","lhx1a","lhx1b","isl1","ipcef1","oxt","ntrk1", "egfra","ltk"))

FeaturePlot(pbmc, reduction = 'tsne', features = c("onecut1","onecut2","onecut3a","lhx1a","lhx1b","isl1","ipcef1","oxt","ntrk1", "egfra","ltk","cpne8"),
            pt.size=0.2, order=TRUE, combine=TRUE)

FeaturePlot(pbmc, reduction = 'tsne', features = c("zgc:162612"), pt.size=0.2, order=TRUE, combine=TRUE)

# SAC genes in Yamagata, 2021 chicken RNAseq + tfap2a = INL and tfap2b = INL + GCL
# Exploring briefly, ACs that are tbx2+ in this dataset are clusters: 11,14,15,16,17(nmb),2,21,22,25 (OFF SAC),31(nts, penk),39,42 (penk),43,47,54,6,7 (ON SAC)
# In mouse amacrine dataset it's clearly in clusters 35 (PENK), 43 (GABA+), 51 (GHRH)
# less clearly in 16, 17 (SACs), 23, 24(nGnG-1), 26 (VIP), 3 (AII),30 (nGnG-3), 56 (VG1), 9 (Gly+)
FeaturePlot(pbmc, reduction = 'tsne', features = c("tbx2a","tbx2b","chata","sox10","tfap2a","tfap2b","slc5a7a","slc18a3a","slc18a3b"), pt.size=0.2, order=TRUE, combine=TRUE)

FeaturePlot(pbmc, reduction = 'tsne', features = c("tbx2a","tbx2b","nts","nmbb","slc18a3a","penka","penkb","ghrh","maff"), pt.size=0.2, order=TRUE, combine=TRUE)

# FeaturePlot(pbmc, reduction = 'tsne', features = c("opn1mw1","opn1mw2","opn1mw3","opn1mw4","opn1lw1","opn1lw2"),
#             pt.size=1, order=FALSE, combine=TRUE)
# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
# CONES
# Subclustering cones (probably will then remerge with rods and photoreceptor progenitors to track gene expression)
# DimPlot(pbmc, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()
# almost impossible separation by subtypes combining all "cones"
photo <- subset(pbmc, idents = c("Cone(adult)","Cone(larval)"))
ps = DimPlot(photo, reduction = "tsne", label=TRUE)+ eelTheme()
ggsave(ps, file="adult_coneClusters01.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")
# clear separation between adult and larval cones
photo2 = SetIdent(photo, value = "Sample")
ps = DimPlot(photo2, reduction = "tsne", label=TRUE)+ eelTheme()
ggsave(ps, file="adult_coneClusters02.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

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


# or use tSNE (clustering can't separate some M cones from L cones)
photo <- RunTSNE(photo, dims = 1:upDimLimit)
ps = DimPlot(photo, reduction = "tsne", label = TRUE, pt.size = 1)+ eelTheme()
ggsave(ps, file="adult_coneClusters03.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

photo2 = SetIdent(photo, value = "Sample")
ps = DimPlot(photo2, reduction = "tsne", label=TRUE)+ eelTheme()
ggsave(ps, file="adult_coneClusters04.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")


FeaturePlot(photo, reduction = 'tsne', features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))
FeaturePlot(photo, reduction = 'tsne', features = c("rho", "nr2e3","nr2f6","crx",'syt5a','syt5b'))
VlnPlot(photo, features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))

# ------------------------------------------------------------------------------------------
# plot layouts for direct export

#OPSINS
p = FeaturePlot(photo, reduction = 'tsne', features = c("rho","opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2'),
                pt.size=1, order=TRUE, combine=FALSE)
lay = rbind(c(1,NA,NA,NA),c(2,3,NA,NA),c(4,5,6,7),c(8,9,NA,NA))
grid.arrange(grobs = p, layout_matrix = lay)
ps = arrangeGrob(grobs = p, layout_matrix = lay)
ggsave(ps, file="adult_Opsins.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

#ROD STUFF (known)
p = FeaturePlot(photo, reduction = 'tsne', features = c("rho",'gnat1',"gnat2","saga",'pde6g',"guca1b","grk1a",'crx','nr2e3','nrl'),pt.size=1, order=FALSE, combine=FALSE)
lay = rbind(c(1,2,3,NA),c(4,5,6,7),c(8,9,10,NA))
grid.arrange(grobs = p, layout_matrix = lay)
ps = arrangeGrob(grobs = p, layout_matrix = lay)
ggsave(ps, file="adult_Rods.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")


ps = DotPlot(photo, features = c("rho","nr2e3","gnat1","gnat2","arr3a","arr3b","thrb","opn1sw1", "opn1sw2",'opn1mw4','opn1mw3','opn1mw2','opn1mw1','opn1lw1','opn1lw2')) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1)) 
ggsave(ps, file="adult_coneMarkersDot.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

# ------------------------------------------------------------------------------------------


# Assigning cone subtype identity to clusters
new.cluster.ids <- c("L","L","L","L","L","L","L","L","US","M3","M4","M1","Ribo","L")
names(new.cluster.ids) <- levels(photo)
photo <- RenameIdents(photo, new.cluster.ids)
Idents(photo) <- factor(x = Idents(photo), levels = c("Ribo","L","M1","M3","M4","US"))
ps = DimPlot(photo, reduction = "tsne", label = TRUE, pt.size = 2) + NoLegend()+ eelTheme()
ggsave(ps, file="adult_coneClusters05.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")
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

# or use tSNE (clustering can't separate some M cones from L cones)
us <- RunTSNE(us, dims = 1:upDimLimit)
ps = DimPlot(us, reduction = "tsne", label=TRUE, pt.size=4) + NoLegend()+ eelTheme()
ggsave(ps, file="adult_coneUS01.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")
# ------------------------------------------------------------------------------------------
# plot layouts for direct export

#OPSINS
p = FeaturePlot(us, reduction = 'tsne', features = c("rho","opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2'),
                pt.size=1, order=TRUE, combine=FALSE)
lay = rbind(c(1,NA,NA,NA),c(2,3,NA,NA),c(4,5,6,7),c(8,9,NA,NA))
grid.arrange(grobs = p, layout_matrix = lay)
ps = arrangeGrob(grobs = p, layout_matrix = lay)
ggsave(ps, file="adult_coneUS02.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")


ps = DotPlot(us, features = c("rho","nr2e3","gnat1","gnat2","arr3b","arr3a","thrb","opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2')) + eelTheme()+ theme(axis.text.x = element_text(angle = 45, hjust=1)) 
ggsave(ps, file="adult_coneUS03.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

# ------------------------------------------------------------------------------------------
FeaturePlot(us, reduction = 'tsne', features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))
VlnPlot(us, features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))
FeaturePlot(us, reduction = 'tsne', features = c('tbx2a','tbx2b','efna1b','foxq2'))

new.cluster.ids <- c("S","S","S","UV")
names(new.cluster.ids) <- levels(us)
us <- RenameIdents(us, new.cluster.ids)
ps = DimPlot(us, reduction = "tsne", label = TRUE, pt.size = 4) + NoLegend()+ eelTheme()
ggsave(ps, file="adult_coneUS04.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

# ------------------------------------------------------------------------------------------
# Probably better at this stage to remerge everything with rods also; will do that in next section after DESeq2 since subtype specific comparisons are still interesting
# pull out other clusters for DESeq2 analysis
lm = subset(photo, idents = c("L","M1","M3","M4","Ribo")); #decided to exclude Ribo cluster since I already explored it
# lm = subset(photo, idents = c("L","M1","M3","M4"));
# remerge into single cone dataset (1894 cones)
photo = merge(x = lm, y = us)
photo

photo = ScaleData(photo)
photo <- FindVariableFeatures(photo, selection.method = "vst", nfeatures = 40)
photo <- RunPCA(photo, features = VariableFeatures(object = photo))
ElbowPlot(photo, ndims=50)
upDimLimit=10; #16
DimPlot(photo, reduction = "pca", label=TRUE)
# UMAP using the same PCA dimensions for prettier visualization
photo <- RunUMAP(photo, dims = 1:upDimLimit)
DimPlot(photo, reduction = "umap", label=TRUE, pt.size = 1) + NoLegend()
# or use tSNE (clustering can't separate some M cones from L cones)
photo <- RunTSNE(photo, dims = 1:upDimLimit)
ps = DimPlot(photo, reduction = "tsne", label=TRUE, pt.size = 2) + NoLegend() +eelTheme()
ggsave(ps, file="adult_coneClusters06.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

# ------------------------------------------------------------------------------------------
# plot layouts for direct export

#OPSINS
p = FeaturePlot(photo, reduction = 'tsne', features = c("rho","opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2'),
                pt.size=1, order=TRUE, combine=FALSE)
lay = rbind(c(1,NA,NA,NA),c(2,3,NA,NA),c(4,5,6,7),c(8,9,NA,NA))
grid.arrange(grobs = p, layout_matrix = lay)
ps = arrangeGrob(grobs = p, layout_matrix = lay)
ggsave(ps, file="adult_coneClusters06.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")


ps = DotPlot(photo, features = c("rho","nr2e3","gnat1","gnat2","arr3b","arr3a","thrb","opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2')) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ggsave(ps, file="adult_coneMarkersDot02.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

# ------------------------------------------------------------------------------------------

FeaturePlot(photo, reduction = 'tsne', features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))
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
DimPlot(photo, reduction = "tsne", label=TRUE, pt.size = 2) + NoLegend()

FeaturePlot(photo, reduction = 'tsne', features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))
VlnPlot(photo, features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))
DotPlot(photo, features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2'))
saveRDS(photo, file = "./cones_Adult.rds")


# ------------------------------------------------------------------------------------------
try(dev.off(dev.list()["RStudioGD"]),silent=TRUE)
try(dev.off(),silent=TRUE)
rm(list=ls())
# ------------------------------------------------------------------------------------------
# Can be restarted here
photo = readRDS("~/Documents/LiMolec/otherRNAseq/zfRet_HoangBlackshaw2020/cones_Adult.rds");
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
photoAll = readRDS("~/Documents/LiMolec/otherRNAseq/zfRet_HoangBlackshaw2020/cones_AdultAll.rds");
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

# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
# remerge cones with rods still split by sutype
# remerge into single dataset (1894 cones + 3093 rods = 4987 PRs)
r = subset(pbmc, idents = c("Rod"));
photoAll = merge(x = lm, y = c(us,r))
photoAll

photoAll = ScaleData(photoAll)
photoAll <- FindVariableFeatures(photoAll, selection.method = "vst", nfeatures = 200)

photoAll <- RunPCA(photoAll, features = VariableFeatures(object = photoAll))
ElbowPlot(photoAll, ndims=50)
# if ever wanted to recluster again
# photo <- FindNeighbors(photoAll, dims = 1:upDimLimit)
# photo <- FindClusters(photoAll, resolution = 1)

# ordering identities for plotting and assigning colors
Idents(photoAll) <- factor(x = Idents(photoAll), levels = c("Rod","UV","S","M1","M3","M4","L"))
prcolors = c("#a3a3a3","#B540B7","#4669F2","#04CD22","#04CD22","#04CD22","#CC2C2A")

#small change for rods, using just "R" as label for consistency
new.cluster.ids <- c("R","UV","S","M1","M3","M4","L")
names(new.cluster.ids) <- levels(photoAll)
photoAll <- RenameIdents(photoAll, new.cluster.ids)
DimPlot(photoAll, reduction = "tsne", label = TRUE, pt.size = 1, label.size = 6) + eelTheme()

# 'r' : '#a3a3a3',
# 'lr' : '#474747'
# 'u' : '#B540B7',

# 's' : '#4669F2',
# 'esls' : '#8f9bcc',

# 'm' : '#04CD22',
# 'lslm' : '#18e236',
# 'mslm' : '#57cb69',
# 'eslm' : '#80bc89',


# 'l' : '#CC2C2A',
# 'lsll' : '#ce4340',
# 'msll' : '#c67271',
# 'esl' : '#d69f9e',

# 'pr' : "#E6B800"
# 'lslpr' : "#cca819"
# 'mslpr' : "#dcc360"
# 'eslpr' : "#dacd9a"
saveRDS(photoAll, file = "./photoreceptors_AdultAll.rds")
# ------------------------------------------------------------------------------------------
# Can be restarted here
photoAll = readRDS("~/Documents/LiMolec/otherRNAseq/zfRet_HoangBlackshaw2020/photoreceptors_AdultAll.rds");
prcolors = c("#a3a3a3","#B540B7","#4669F2","#04CD22","#04CD22","#04CD22","#CC2C2A")
# ------------------------------------------------------------------------------------------
upDimLimit=30;
DimPlot(photoAll, reduction = "pca", label=TRUE, cols=prcolors)
# UMAP using the same PCA dimensions for prettier visualization
photoAll <- RunUMAP(photoAll, dims = 1:upDimLimit)
DimPlot(photoAll, reduction = "umap", label=TRUE, pt.size = 1) + NoLegend()
# or use tSNE (clustering can't separate some M cones from L cones)
photoAll <- RunTSNE(photoAll, dims = 1:upDimLimit)
ps = DimPlot(photoAll, reduction = "tsne", label=FALSE, pt.size = 2, cols=prcolors) +eelTheme()
ggsave(ps, file="adult_coneClusters07.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

FeaturePlot(photoAll, reduction = 'tsne', features = c("rho","opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2'),
            pt.size=1, order=FALSE, combine=TRUE)
# ------------------------------------------------------------------------------------------
# plot layouts for direct export

#OPSINS
p = FeaturePlot(photoAll, reduction = 'tsne', features = c("rho","opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2'),
                pt.size=1, order=FALSE, combine=FALSE)
lay = rbind(c(1,NA,NA,NA),c(2,3,NA,NA),c(4,5,6,7),c(8,9,NA,NA))
grid.arrange(grobs = p, layout_matrix = lay)
ps = arrangeGrob(grobs = p, layout_matrix = lay)
ggsave(ps, file="adult_coneClusters08.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

# OPSINS + KNOWN MARKERS
ps = DotPlot(photoAll, features = c("rho","nr2e3","gnat1","gnat2","arr3b","arr3a","thrb","opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2')) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ggsave(ps, file="adult_coneMarkersDot03.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

#EXPLORING GENES ID'D IN OUR OWN RNAseq
ps = DotPlot(photoAll, features = c("sema3aa","sema3ab","sema3b","sema3bl","sema3c","sema3d","sema3e","sema3fa","sema3fb","sema3ga","sema3gb","sema3h","sema4aa","sema4ab","sema4ba","sema4bb","sema4c","sema4d","sema4e","sema4f","sema4ga","sema4gb","sema5a","sema5ba","sema5bb","sema6a","sema6ba","sema6bb","sema6d","sema6dl","sema6e","sema7a")) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ggsave(ps, file="adult_coneMarkersDot04_sema.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

ps = DotPlot(photoAll, features = c("vamp1","vamp2","vamp3","vamp4","vamp5","vamp8")) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ps
ggsave(ps, file="adult_coneMarkersDot04_vamp.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

ps = DotPlot(photoAll, features = c("syt1a","syt1b","syt2a","syt3","syt4","syt5a","syt5b","syt6a","syt6b","syt7a","syt7b")) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ps
ggsave(ps, file="adult_coneMarkersDot04_syt.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

ps = DotPlot(photoAll, features = c("efna1a","efna1b","efna2a","efna2b","efna3a","efna3b","efna5a","efna5b")) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ps
ggsave(ps, file="adult_coneMarkersDot04_efna.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")


ps = DotPlot(photoAll, features = c("tbx1","tbx15","tbx16","tbx18","tbx20","tbx21","tbx22","tbx2a","tbx2b","tbx3a","tbx3b","tbx4","tbx5a","tbx5b","tbx6")) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ps
ggsave(ps, file="adult_coneMarkersDot04_tbx.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")


ps = DotPlot(photoAll, features = c("nr2c1","nr2c2","nr2e1","nr2e3","nr2f1a","nr2f1b","nr2f2","nr2f5","nr2f6a","nr2f6b")) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ps
ggsave(ps, file="adult_coneMarkersDot04_nr2.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

ps = DotPlot(photoAll, features = c("ntn1a","ntn1b","ntn2","ntn4","ntn5","ntng1a","ntng2a","ntng2b")) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ps
ggsave(ps, file="adult_coneMarkersDot04_ntn.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

# Rod genes
ps = DotPlot(photoAll, features = c("rho","gnat1","dhrs13l1","lingo1b","cabp4","saga","sagb","gnb1b","rgs9bp","eno2","guca1b","grk1a","rom1b","cplx4c","lrrn1","ncs1a","sgce","kcnv2a","pdca")) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ps
# ggsave(ps, file="adult_coneMarkersDot04_nr2.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")
# Cone Genes
ps = DotPlot(photoAll, features = c("slc1a8b","rgs9a","slc25a24","drd4b","si:busm1-57f23.1","kera","nexn","clic1","gnat2","tgif1","crhbp","kcnv2b","anks1b","clul1")) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ps

# UV+S Genes
ps = DotPlot(photoAll, features = c("grk7b","skor1a","pcdh11","srgap3","ablim1a","nav2a","jam2a","chl1a","s100z","foxq2","sh3bp5b","kcnk1a","ntng2b","nxph1","pik3r3b","myl4","pacrg","fah")) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ps


# M+L Genes
ps = DotPlot(photoAll, features = c("slc32a1","apln","lactbl1b","nrtn","pcdh10a","myo7aa")) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ps

# Mup/Ldown Genes
ps = DotPlot(photoAll, features = c("dok6","spock3","lrrfip1a","sema3fb","itga1","auts2a","esama","esamb","plxnb1a","cgnb","phf19","lrrc20")) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ps

# Lup/Mdown Genes
ps = DotPlot(photoAll, features = c("slc32a1","s100v2","smad5","snap25a","ttyh2l","arhgap11a","ggctb","fbxo32","pcdh10a","abracl")) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ps

# UVup/Sdown Genes
ps = DotPlot(photoAll, features = c("myl4","ttyh2l","lhx1a","rx3","itgb1bp1")) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ps

# Sup/UVdown Genes
ps = DotPlot(photoAll, features = c("fkbp5","prss23","chkb","mpzl2b","nr1d1","camk2d1","frmpd2","foxo1a","prom1a")) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ps
   
   
# TFs in DEGs
ps = DotPlot(photoAll, features = c("tbx2a","tbx2b","otx5","tgif1","crhbp","six7","six6b","ndrg1b","egr1","nr1d1","sall1a","skor1a","foxq2","lhx1a","ntf3","rxrga","thrb","fgf1b","eya2","sox6","sox4b","pbx1a","hmgb3a","fbxo21","tfe3a")) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ps

# ------------------------------------------------------------------------------------------

FeaturePlot(photoAll, reduction = 'tsne', features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))
VlnPlot(photoAll, features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))
DotPlot(photoAll, features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2'))
saveRDS(photoAll, file = "./photoreceptors_AdultAll.rds")

FeaturePlot(photoAll, reduction = 'tsne', features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','cyp26c1','cyp26a1'), order=TRUE)

# cell counter
table(Idents(photoAll))

# Rod   UV    S   M1   M3   M4    L 
# 3097   20   82   56   73   70 1589 

# ------------------------------------------------------------------------------------------
# getting % and Mean expression for all genes into csv.
# DoHeatmap and DotPlot use the @scale.data slot for average expression display, which z-scored expression values (for example, as those used in PCA).
# Cells with a value > 0 represent cells with expression above the population mean (a value of 1 would represent cells with expression 1SD away from the population mean). Hope that helps!
# This is very detailed explanation: https://github.com/satijalab/seurat/issues/2798

# Average expression in non-log scale:
avgExp = AverageExpression(photoAll, slot="counts")
avgExp = avgExp$RNA
# avgExp = head(avgExp,10)
nI = length(levels(Idents(photoAll)))
colnames(avgExp) = paste("avg", levels(Idents(photoAll)), sep = "")

# Calculate mean expression across all subtypes
avgExp$avgMean = rowMeans(avgExp)
avgExp[,c(ncol(avgExp),1:(ncol(avgExp)-1))]


# Percent expression can be obtained through DotPlot, but not Average expression, as this is log z-cored normalized data # avgExpZ = dPlot$data$avg.exp
# Percent expression in larvae:
dPlot = DotPlot(photoAll, features=rownames(avgExp))
pctExp = dPlot$data$pct.exp
pctExp = as.data.frame(matrix(pctExp,nrow = length(pctExp)/nI, ncol=nI));

rownames(pctExp) = rownames(avgExp)
colnames(pctExp) = paste("pct", levels(Idents(photoAll)), sep = "")

# Calculate mean percent expression across all ages
pctExp$pctMean = rowMeans(pctExp)
pctExp[,c(ncol(pctExp),1:(ncol(pctExp)-1))]


photoreceptors_Hoang = cbind(pctExp,avgExp)



# these are duplicated genes with lower vs upper case
dups=toupper(c("arid5b","asph", "aste1", "crip2", "ctbp1", "dab2", "eif1b", "flnb", "frmd7", "galnt10", "grxcr1", "gse1", "hist1h4l", "hspb11", "kif1c", "lamp1", "maf1", "pamr1", "pcdh20", "pde6h", "phlpp2", "psmb10", "ptp4a3", "reep6", "rfesd", "rgs9bp", "rps17", "shank2", "slc16a7", "slc25a10", "slc6a13", "slc9a1", "srbd1", "tatdn3", "tenm3", "tmem178b", "tmem241", "tom1l2", "tp53inp2", "trappc9", "tsc22d3", "ube2o", "zc3h12a", "znf423"))
for(i in 1:length(dups)) {
   rownames(photoreceptors_Hoang)[rownames(photoreceptors_Hoang) == dups[i]] = paste(dups[i],"_ii",sep="")
}
# now everything can be lower cased
rownames(photoreceptors_Hoang) = tolower(rownames(photoreceptors_Hoang))




photoreceptors_Hoang = photoreceptors_Hoang[c("pctMean","pctRod","pctUV","pctS","pctM1","pctM3","pctM4","pctL","avgMean","avgRod","avgUV","avgS","avgM1","avgM3","avgM4","avgL")]



scangenes=c('syt5a','syt5b','nr2e3','tbx2a','tbx2b','arr3a','arr3b','opn1lw1','opn1lw2');
photoreceptors_Hoang[scangenes,]

conesAdult_Hoang[scangenes,]

# still need to manually add symbol to first column
write.csv(photoreceptors_Hoang,"./photoreceptors_Hoang.csv",quote=FALSE)
