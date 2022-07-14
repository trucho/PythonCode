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
setwd("/Users/angueyraaristjm/Documents/eelMolec/zfRNAseq/2020_Hoang_zfRet10x/")
directory = "/Users/angueyraaristjm/Documents/eelMolec/zfRNAseq/2020_Hoang_zfRet10x/"
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


# 'r' : '#7d7d7d',
# 'lr' : '#a3a3a3'
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
# 'prp' : "#dfdac8"
# -------------------------------------------------------------------
# Load the previously saved larval cone dataset (flagged by sample name != AdR)
photo = readRDS(file = "./cones_Larval.rds")


# Identification of highly variable features (feature selection)
photo <- FindVariableFeatures(photo, selection.method = "vst", nfeatures = 200)
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
ElbowPlot(photo, ndims=50) + eelTheme()
ggsave("larval_Elbow.png", path=exportDir, width = 140, height = 105, units = "mm")

# Cluster cells
# nFeatures = 500
#        ndim = 24, res = 1 -> c = 10; 
               # one big cluster expresses syt5a and gnat2 with a few rods inside and it's basically al cells from 5dpf
               # another cluster is nr2e3+
# nFeatures = 100
#        ndim = 33, res = 1 -> c = 9;
               # same 2 big clusters but some are misasigned
# nFeatures = 200
#        ndim = 26, res = 0.6 -> c = 9;
               # seems like clustering is pretty much the same with all of these.
#        ndim = 12, res = 0.4 -> c = 8; not good; a lot of empty space in tSNE and umap
#        ndim = 12, res = 1 -> c = 12; not good; too many clusters
#        ndim = 45, res = 1 -> c = 12; good elongated plot that shows dev age in diagonal
               # I think that 3 groups can be parsed apart: late-stage (syt5a+, gnat+,nr2f6b+), early-stage (nr2e3+) and mid-stage(syt5a-, syt5b+, arr3+,nr2e3-)

# Settling for nFeatures = 200; ndim = 45, res = 1
upDimLimit=45;
photo <- FindNeighbors(photo, dims = 1:upDimLimit)
photo <- FindClusters(photo, resolution = 1)

DimPlot(photo, reduction = "pca", label = TRUE, pt.size = 1, label.size = 6) + eelTheme()
ggsave("larval_PCAinitial.png", path=exportDir, width = 140, height = 105, units = "mm")

# UMAP using the same PCA dimensions for prettier visualization
photo <- RunUMAP(photo, dims = 1:upDimLimit)
DimPlot(photo, reduction = "umap", label = TRUE, pt.size = 1, label.size = 6) + eelTheme()
ggsave("larval_UMAPinitial.png", path=exportDir, width = 140, height = 105, units = "mm")

# or use tSNE (clustering can't separate some M cones from L cones)
photo <- RunTSNE(photo, dims = 1:upDimLimit)
DimPlot(photo, reduction = "tsne", label = TRUE, pt.size = 1, label.size = 6) + eelTheme()
ggsave("larval_TSNEinitial.png", path=exportDir, width = 140, height = 105, units = "mm")

# Plot by age
photo2 = SetIdent(photo, value = "Sample")
DimPlot(photo2, reduction = "tsne", label = TRUE, pt.size = 1, label.size = 6) + eelTheme()
ggsave("larval_TSNEdpf.png", path=exportDir, width = 140, height = 105, units = "mm")

FeaturePlot(photo, reduction = 'tsne',features = c("foxq2")) + eelTheme()
ggsave("larval_TSNEfoxq2.png", path=exportDir, width = 140, height = 105, units = "mm")

FeaturePlot(photo, reduction = 'tsne',features = c("opn1sw2")) + eelTheme()
ggsave("larval_TSNEsws2.png", path=exportDir, width = 140, height = 105, units = "mm")

# FeaturePlot(photo, reduction = 'tsne', features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1','efna1b'))
FeaturePlot(photo, reduction = 'tsne', features = c("rho","opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2'))
FeaturePlot(photo, reduction = 'tsne', features = c("rho","nrl","pde6a","nr2e3","crx",'gnat1',"saga","sagb","gucy2f","grk1a")) # rod markers
FeaturePlot(photo, reduction = 'tsne', features = c("rho","eno2","rom1a","lrrn1","unc119.2","pdca","cplx4c")) # rod genes from my own dataset
FeaturePlot(photo, reduction = 'tsne', features = c("gnat2",'arr3a','arr3b',"pde6c","pde6h","guca1d","grk7a","crx","neurod1","nr2f6b")) #cone markers
FeaturePlot(photo, reduction = 'tsne', features = c("nr2e3","nr2f6b","crx",'syt5a','syt5b','gnat2','arr3a','arr3b',"pde6g","pde6h", "neurod1")) #cone markers
FeaturePlot(photo, reduction = 'tsne', features = c("slc1a8b","rgs9a","slc25a24","ppa1a","sema7a","kera","nexn","dusp5","crhbp","prph2a")) #cone markers from my own dataset
FeaturePlot(photo, reduction = 'tsne', features = c("otx5","tbx2a","tbx2b","rxrga","rxrgb", "thrb", "sema7a","cnga3a","cnga3b"))
VlnPlot(photo, features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))
FeaturePlot(photo, reduction = 'tsne', features = c("nrl","mafa","mafb")) # rod markers

FeaturePlot(photo, reduction = 'tsne', features = c("crx","nr2e3","nr2f6a","nr2f6b",'syt5a','syt5b','gnat2',"pde6g","pde6h", "neurod1")) #dev markers
# ------------------------------------------------------------------------------------------
# plot layouts for direct export
#OPSINS
p = FeaturePlot(photo, reduction = 'tsne', features = c("rho","opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2'),
                pt.size=1, order=TRUE, combine=FALSE)
lay = rbind(c(1,NA,NA,NA),c(2,3,NA,NA),c(4,5,6,7),c(8,9,NA,NA))
grid.arrange(grobs = p, layout_matrix = lay)
ps = arrangeGrob(grobs = p, layout_matrix = lay)
ggsave(ps, file="larval_Opsins.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

#ROD STUFF (known)
p = FeaturePlot(photo, reduction = 'tsne', features = c("rho",'gnat1',"saga",'pde6g',"guca1b","grk1a",'crx','nr2e3','nrl'),pt.size=1, order=TRUE, combine=FALSE)
lay = rbind(c(1,2,NA,NA),c(3,4,5,6),c(7,8,9,NA))
grid.arrange(grobs = p, layout_matrix = lay)
ps = arrangeGrob(grobs = p, layout_matrix = lay)
ggsave(ps, file="larval_Rods.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

#JUAN'S DEV MARKERS STUFF
p = FeaturePlot(photo, reduction = 'tsne', features = c("syt5a","syt5b","gnat2","pde6g","pde6h","crx","neurod1","nr2e3","nr2f6b","nr2f1b"),pt.size=1, order=TRUE, combine=FALSE)
lay = rbind(c(1,2,NA,NA,NA),c(3,4,5,NA,NA),c(6,7,8,9,10))
grid.arrange(grobs = p, layout_matrix = lay)
ps = arrangeGrob(grobs = p, layout_matrix = lay)
ggsave(ps, file="larval_DevClusters.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")


# ------------------------------------------------------------------------------------------

# FeaturePlot(photo, reduction = 'tsne', features = c("cyp26c1","cyp26a1"))




# Assigning cone subtype identity to clusters
# new.cluster.ids <- c("lslPR","mslPR","lslPR","mslPR","eslPR","eslPR","eslPR","sw2eslPR","mslPR")
new.cluster.ids <- c("lslPR","lslPR","eslPR","mslPR","mslPR","eslPR","mslPR","eslPR","eslS","mslPR")
names(new.cluster.ids) <- levels(photo)
photo <- RenameIdents(photo, new.cluster.ids)

Idents(photo) <- factor(x = Idents(photo), levels = c("eslPR","mslPR","lslPR","eslS"))
lcolors = c("#dacd9a","#dcc360","#cca819","#8f9bcc")
levels(Idents(photo))


DimPlot(photo, reduction = "tsne", label = TRUE, pt.size = 1, label.size = 6, cols=lcolors) +eelTheme()
ggsave("larval_TSNEdev.png", path=exportDir, width = 140, height = 105, units = "mm")
DimPlot(photo, reduction = "umap", label = TRUE, pt.size = 1) + NoLegend()

FeaturePlot(photo, reduction = "tsne", features = c("rho","opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))
VlnPlot(photo, features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))
VlnPlot(photo, features = c("nr2e3","nr2f6b","crx",'syt5a','syt5b','gnat2','arr3a','arr3b',"pde6g","pde6h", "neurod1"))

saveRDS(photo, file = "./cones_LarvalReClustered.rds")
# ------------------------------------------------------------------------------------------
# Attempting to separate subtypes in late-stage larval photoreceptors: FAILED except rods vs. cones

# there is a small group of cell that are cyp26c1/a1 (enzymes involved in RetAcid degradation and retinal patterning)

lsl = subset(photo, idents = "lslPR")
# Identification of highly variable features (feature selection)
lsl <- FindVariableFeatures(lsl, selection.method = "vst", nfeatures = 50)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(lsl), 200)
top10


#run PCA
lsl <- RunPCA(lsl, features = VariableFeatures(object = lsl))
# Use elbow plot (var explained)
ElbowPlot(lsl, ndims=50) + eelTheme()
ggsave("larval_lslElbow.png", path=exportDir, width = 140, height = 105, units = "mm")

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
DimPlot(lsl, reduction = "tsne", label=TRUE,  pt.size = 1, label.size = 6) + eelTheme()
ggsave("larval_lslTSNEinitial.png", path=exportDir, width = 140, height = 105, units = "mm")

FeaturePlot(lsl, reduction = 'tsne', features = c("rho", "opn1sw1", "opn1sw2",'opn1mw1','opn1lw2','si:busm1-57f23.1'))
FeaturePlot(lsl, reduction = 'tsne', features = c("gnat1","gnat2","syt5a","rho","sagb"))
VlnPlot(lsl, features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))
FeaturePlot(lsl, reduction = 'tsne', features = c('tbx2a','tbx2b','efna1b','foxq2',"ntf3","nr2e3","nr2f6b","thrb"))
DotPlot(lsl, features = c('tbx2a','tbx2b','efna1b','foxq2',"ntf3","nr2e3","nr2f6b","thrb"))

FeaturePlot(lsl, reduction = 'tsne', features = c("rho","gnat1","gnat2","opn1sw1", "opn1sw2",'opn1mw1','opn1lw2'))
# ------------------------------------------------------------------------------------------
# plot layouts for direct export
#OPSINS
p = FeaturePlot(lsl, reduction = 'tsne', features = c("rho","nr2e3","gnat1","gnat2","arr3b","arr3a","opn1sw1", "opn1sw2",'opn1mw1','opn1lw2'),
                pt.size=1, order=TRUE, combine=FALSE)
lay = rbind(c(1,2,3,NA),c(4,5,6,NA),c(7,8,9,10))
grid.arrange(grobs = p, layout_matrix = lay)
ps = arrangeGrob(grobs = p, layout_matrix = lay)
ggsave(ps, file="larval_lslMarkers.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")
# ------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------
# Removing rods and trying clustering again
new.cluster.ids <- c("C","C","C","C","R","C")
names(new.cluster.ids) <- levels(lsl)
lsl <- RenameIdents(lsl, new.cluster.ids)
DimPlot(lsl, reduction = "tsne", label = TRUE, pt.size = 1, label.size = 6) + eelTheme()
ggsave("larval_lslTSNEfinal.png", path=exportDir, width = 140, height = 105, units = "mm")

lslRods = subset(lsl, idents = "R")

lsl = subset(lsl, idents = "C")
# ------------------------------------------------------------------------------------------
lsl <- FindVariableFeatures(lsl, selection.method = "vst", nfeatures = 40)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(lsl), 200)
top10


#run PCA
lsl <- RunPCA(lsl, features = VariableFeatures(object = lsl))
# Use elbow plot (var explained)
ElbowPlot(lsl, ndims=50)

# Cluster cells
# nFeatures = 100, ndim = 30, res = 1 -> c=6, 
# nFeatures = 50, ndim = 30, res = 1 -> no separation, except for rods in cluster 5
# nFeatures = 250, ndim = 35, res = 1 -> no separation, except for rods in cluster 5; 0 might be M-cones?
# nFeatures = 1000, ndim = 20, res = 1 -> not good either
# nFeatures = 20, ndim = 15, res = 1 -> forces separation on opsins giving a mixed cluster and a pure L
upDimLimit=30;
lsl <- FindNeighbors(lsl, dims = 1:upDimLimit)
lsl <- FindClusters(lsl, resolution = 1)

DimPlot(lsl, reduction = "pca", label=TRUE)

# UMAP using the same PCA dimensions for prettier visualization
lsl <- RunUMAP(lsl, dims = 1:upDimLimit)
DimPlot(lsl, reduction = "umap", label=TRUE, pt.size = 2, label.size = 8)

# or use tSNE (clustering can't separate some M cones from L cones)
lsl <- RunTSNE(lsl, dims = 1:upDimLimit)
DimPlot(lsl, reduction = "tsne", label=TRUE, pt.size = 2, label.size = 8)


FeaturePlot(lsl, reduction = 'tsne', features = c("rho", "opn1sw1", "opn1sw2",'opn1mw1','opn1lw2','si:busm1-57f23.1'))
FeaturePlot(lsl, reduction = 'tsne', features = c("gnat1","gnat2","syt5a",'syt5b',"sox2","atoh7"))
VlnPlot(lsl, features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1lw1','opn1lw2','si:busm1-57f23.1'))
FeaturePlot(lsl, reduction = 'tsne', features = c('tbx2a','tbx2b','efna1b','foxq2',"ntf3","nr2e3","thrb"))
DotPlot(lsl, features = c('tbx2a','tbx2b','efna1b','foxq2',"ntf3","nr2e3","thrb"))


# VARIABILITY IN THIS DATASET DOES NOT DEPEND readily on cone subtype except for sws1 expression but same cells have reads for many opsins
print(lsl[["pca"]], dims = 1:5, nfeatures = 5)

print(lsl[["pca"]], dims = 1:2, nfeatures = 20)

FeaturePlot(lsl, reduction = 'tsne', features = c("opn1sw1","gngt2a","gngt2b", "gnao1b","aanat2", "calb2a","mt2","nme2b.1","cyp26c1","cyp26a1","gad1b"))
FeaturePlot(photo, reduction = 'tsne', features = c("opn1sw1","gngt2a","gngt2b", "gnao1b","aanat2", "calb2a","mt2","nme2b.1","cyp26c1","cyp26a1","gad1b"))
FeaturePlot(lsl, reduction = 'tsne', features = c("marcksl1b","RPL41", "her6","hmgb2b", "junba","stmn1a","fosab","hmgn2","her15.1"))



# correlation between feature
FeaturePlot(object = lsl, reduction="tsne" ,features = c('opn1sw1', 'opn1lw2'), blend = TRUE)
FeaturePlot(object = lsl, reduction="tsne" ,features = c('opn1sw1', 'opn1sw2'), blend = TRUE)
FeaturePlot(object = lsl, reduction="tsne" ,features = c('opn1mw1', 'opn1lw2'), blend = TRUE)
FeatureScatter(lsl,"opn1sw1","opn1sw2")
FeatureScatter(lsl,"arr3a","opn1lw2")

# expression by known markers
lsl2 = lsl
uvsMarker <- list(c("opn1sw1", "opn1sw2", "arr3b"))
lsl2 <- AddModuleScore(object = lsl2, features = uvsMarker, name = "UVS")
FeaturePlot(object = lsl2, reduction="tsne", features = "UVS1")
lMarker <- list(c('opn1lw2','si:busm1-57f23.1',"thrb","arr3b"))
lsl2 <- AddModuleScore(object = lsl2, features = lMarker, name = "L")
FeaturePlot(object = lsl2, reduction="tsne", features = "L1")

# renaming clusters as lsl for remerge later
new.cluster.ids <- c("lslPR","lslPR","lslPR","lslPR","lslPR")
names(new.cluster.ids) <- levels(lsl)
lsl <- RenameIdents(lsl, new.cluster.ids)
DimPlot(lsl, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()
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
# nFeatures = 100, ndim = 30, res = 1 -> no separation
# nFeatures = 50, ndim = 30, res = 1 -> no separation, except for rods in cluster 5
# nFeatures = 250, ndim = 35, res = 1 -> no separation, except for rods in cluster 5; 0 might be M-cones?
# nFeatures = 1000, ndim = 20, res = 1 -> not good either
# nFeatures = 20, ndim = 15, res = 1 -> forces separation on opsins giving a mixed cluster and a pure L
upDimLimit=30;
msl <- FindNeighbors(msl, dims = 1:upDimLimit)
msl <- FindClusters(msl, resolution = 1)

DimPlot(msl, reduction = "pca", label=TRUE)

# UMAP using the same PCA dimensions for prettier visualization
msl <- RunUMAP(msl, dims = 1:upDimLimit)
DimPlot(msl, reduction = "umap", label=TRUE, pt.size = 2, label.size = 8)

# or use tSNE (clustering can't separate some M cones from L cones)
msl <- RunTSNE(msl, dims = 1:upDimLimit)
DimPlot(msl, reduction = "tsne", label=TRUE, pt.size = 1, label.size = 6) + eelTheme()
ggsave("larval_mslTSNEinitial.png", path=exportDir, width = 140, height = 105, units = "mm")

FeaturePlot(msl, reduction = 'tsne', features = c("rho", "opn1sw1", "opn1sw2",'opn1mw1','opn1lw2','si:busm1-57f23.1'))
FeaturePlot(msl, reduction = 'tsne', features = c("gnat1","gnat2","syt5a","syt5b","arr3a","arr3b"))
VlnPlot(msl, features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))
FeaturePlot(msl, reduction = 'tsne', features = c('tbx2a','tbx2b','efna1b','foxq2',"ntf3","nr2e3","thrb"))
DotPlot(msl, features = c('tbx2a','tbx2b','efna1b','foxq2',"ntf3","nr2e3","thrb"))
# seems like L and M cones do have their own clusters, with a subdivision for cells that are farther along development with more expression of phototransduction proteins

FeaturePlot(msl, reduction = 'tsne', features = c("rho","nr2e3","gnat1","gnat2","arr3a","arr3b","opn1sw1", "opn1sw2",'opn1mw1','opn1lw2','thrb'))

# ------------------------------------------------------------------------------------------
# plot layouts for direct export
#Subtype markers
p = FeaturePlot(msl, reduction = 'tsne', features = c("rho","nr2e3","gnat1","gnat2","arr3b","arr3a","thrb","opn1sw1", "opn1sw2",'opn1mw1','opn1lw2'),
                pt.size=1, order=TRUE, combine=FALSE)
lay = rbind(c(1,2,3,NA),c(4,5,6,7),c(8,9,10,11))
grid.arrange(grobs = p, layout_matrix = lay)
ps = arrangeGrob(grobs = p, layout_matrix = lay)
ggsave(ps, file="larval_mslMarkers.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")
# ------------------------------------------------------------------------------------------


# this is what PC2 is picking up on
FeaturePlot(msl, reduction = 'tsne', features = c("opn1mw1","opn1lw2","gnat2","arr3a","rgs9a","aanat2"))
FeaturePlot(object = msl, reduction="tsne" ,features = c('gnat2', 'rgs9a'), blend = TRUE)


print(msl[["pca"]], dims = 1:2, nfeatures = 20)
#PC1 genes
FeaturePlot(msl, reduction = 'tsne', features = c("mdka","rho","gngt1","gstp1","stmn1b","gnat1","sepp1a","sncgb","rgs9a"))
#PC2 genes
FeaturePlot(msl, reduction = 'tsne', features = c("rcvrn3", "gnb3b", "si:ch211-81a5.8", "gnat2", "gnb5b", "cplx4a", "nme2a", "pde6h", "ppa1a"))
FeaturePlot(msl, reduction = 'tsne', features = c("pdcb", "prph2b", "prph2a", "si:dkey-17e16.15", "si:busm1-57f23.1", "opn1lw2", "si:ch211-285j22.3", "aqp9b", "pde6c", "mdka"))

# Ploting by age reveals that there is just one 2.5 dpf cell and maybe 6 cells from 5 dpf; rest are intermingled 3 and 4 dpf cones (without clear batch effects either).
msl2 = msl
msl2 = SetIdent(msl2, value = "Sample")
DimPlot(msl2, reduction = "tsne", label = TRUE, pt.size = 2, label.size = 8)

# Trying to id what small cluster 3 is -> DEGs are for low counts of phototransduction genes (arr3a, gngt2a/b,rs1a,gnat2 + hmgn6,...) and high nr2e3
FeaturePlot(msl, reduction = 'tsne', features = c("arr3a","gngt2a","gngt2b","rs1a","hmgn6","gnat2","si:ch211-81a5.8","nr2e3","rgs9a"))
# DESeq2_temp = FindMarkers(msl, c("3"), c("0","1","2"), test.use="DESeq2");
# AvgExp_temp = AverageExpression(msl, features=rownames(DESeq2_temp))
# AvgExp_temp = AvgExp_temp$RNA
# DESeq2_temp = merge.data.frame(DESeq2_temp,AvgExp_temp, by="row.names", sort=FALSE)
# DESeq2_temp$baseMean = rowMeans(DESeq2_temp[,c("0","1","2","3")])
# rownames(DESeq2_temp) = DESeq2_temp$Row.names
# DESeq2_temp = DESeq2_temp[,c("baseMean","avg_logFC","p_val","p_val_adj","pct.1","pct.2","0","1","2","3")]
# PctExp_temp = DotPlot(msl, features=rownames(DESeq2_temp))
# PctNames = rownames(PctExp_temp$data)[1:(length(rownames(PctExp_temp$data))/4)]
# PctExp_temp = PctExp_temp$data$pct.exp
# PctExp_temp = matrix(PctExp_temp,nrow = length(PctExp_temp)/4, ncol=4);
# rownames(PctExp_temp) = PctNames
# colnames(PctExp_temp) = c("pct0","pct1","pct2","pct3")
# DESeq2_temp = merge.data.frame(DESeq2_temp,PctExp_temp, by="row.names", sort=FALSE)
# rownames(DESeq2_temp) = DESeq2_temp$Row.names
# write.csv(DESeq2_temp,"./DESeq2_temp.csv")





new.cluster.ids <- c("mslPR","mslM","mslM","mslM","mslL","mslL")
names(new.cluster.ids) <- levels(msl)
msl <- RenameIdents(msl, new.cluster.ids)
DimPlot(msl, reduction = "tsne", label = TRUE, pt.size = 1, label.size = 6) + eelTheme()
ggsave("larval_mslTSNEmid.png", path=exportDir, width = 140, height = 105, units = "mm")

FeaturePlot(msl, reduction = 'tsne', features = c("rho", "opn1sw1", "opn1sw2",'opn1mw1','opn1lw2','si:busm1-57f23.1'))
VlnPlot(msl, features = c("rho", "opn1sw1", "opn1sw2",'opn1mw1','opn1lw2','si:busm1-57f23.1'))
# ------------------------------------------------------------------------------------------
# Exploring esl briefly
esl =  subset(photo, idents = c("eslPR","eslS"))

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
# nFeatures = 100, ndim = 30, res = 1 -> no separation
# nFeatures = 50, ndim = 30, res = 1 -> no separation, except for rods in cluster 5
# nFeatures = 250, ndim = 35, res = 1 -> no separation, except for rods in cluster 5; 0 might be M-cones?
# nFeatures = 1000, ndim = 20, res = 1 -> not good either
# nFeatures = 20, ndim = 15, res = 1 -> forces separation on opsins giving a mixed cluster and a pure L
upDimLimit=20;
esl <- FindNeighbors(esl, dims = 1:upDimLimit)
esl <- FindClusters(esl, resolution = 1)

DimPlot(esl, reduction = "pca", label=TRUE)

# UMAP using the same PCA dimensions for prettier visualization
esl <- RunUMAP(esl, dims = 1:upDimLimit)
DimPlot(esl, reduction = "umap", label=TRUE, pt.size = 2, label.size = 8)


# or use tSNE (clustering can't separate some M cones from L cones)
esl <- RunTSNE(esl, dims = 1:upDimLimit)
DimPlot(esl, reduction = "tsne", label=TRUE, pt.size = 1, label.size = 6) + eelTheme()
ggsave("larval_eslTSNEinitial.png", path=exportDir, width = 140, height = 105, units = "mm")

FeaturePlot(esl, reduction = 'tsne', features = c("rho", "opn1sw1", "opn1sw2",'opn1mw1','opn1lw2','si:busm1-57f23.1'))
FeaturePlot(esl, reduction = 'tsne', features = c("gnat1","gnat2","syt5a","syt5b","arr3a","arr3b"))
VlnPlot(esl, features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))
FeaturePlot(esl, reduction = 'tsne', features = c('tbx2a','tbx2b','efna1b','foxq2',"ntf3","nr2e3","thrb"))
DotPlot(esl, features = c('tbx2a','tbx2b','efna1b','foxq2',"ntf3","nr2e3","thrb"))
# seems like L and M cones do have their own clusters still and there is a S cluster which is foxq2+
# no gnat yet, nr2e3+ everywhere, six7 is present almost in all clusters, 
# one cluster is dscamb+ (not found in msl or lsl)
# another cluster is zgc:101100+
# some cells are also aanat2+

FeaturePlot(esl, reduction = 'tsne', features = c("opn1lw2","si:busm1-57f23.1","thrb","arr3a"))

FeaturePlot(esl, reduction = 'tsne', features = c("opn1mw1","opn1lw2","gnat2","arr3a","rgs9a"))
FeaturePlot(object = esl, reduction="tsne" ,features = c('opn1sw2', 'foxq2'), blend = TRUE)
FeaturePlot(object = esl, reduction="tsne" ,features = c('opn1mw1', 'arr3a'), blend = TRUE)
FeaturePlot(object = esl, reduction="tsne" ,features = c('opn1lw2', 'arr3a'), blend = TRUE)


print(esl[["pca"]], dims = 1:3, nfeatures = 20)
#PC1 genes
FeaturePlot(esl, reduction = 'tsne', features = c("nr2e3", "si:ch73-28h20.1", "six7", "neurod1", "foxq2", "sall1a", "opn6a", "dscamb", "tulp1a", "zgc:101100"))
FeaturePlot(esl, reduction = 'tsne', features = c("zbtb18", "myclb", "rcvrn2", "fabp7a", "rs1a", "thrb", "nusap1", "pcdh8", "rpl39", "cuedc1a"))
# FeaturePlot(msl, reduction = 'tsne', features = c("zbtb18", "myclb", "rcvrn2", "fabp7a", "rs1a", "thrb", "nusap1", "pcdh8", "rpl39", "cuedc1a"))
# FeaturePlot(esl, reduction = 'tsne', features = c("mdka", "rho", "gngt1", "gnat1", "sepp1a", "pde6h", "sncgb", "mt2", "nrgna", "id1"))
# FeaturePlot(esl, reduction = 'tsne', features = c("si:ch211-113d22.2", "ndrg4", "pcp4a", "efhd1", "tox", "bcl2l13", "cnga1", "mdkb", "pax6b", "gad1b"))

#PC2 genes
FeaturePlot(esl, reduction = 'tsne', features = c("foxq2", "zgc:101100", "rs1a", "opn1sw2", "si:ch211-81a5.8", "arr3b", "hsp70l", "gngt2b", "tdh2", "ppp1r14c"))
FeaturePlot(esl, reduction = 'tsne', features = c("opn1mw1", "tulp1a", "prom1b", "apln", "opn6a", "arl3l2", "hspb1", "opn1sw1", "rpl39", "prph2a","aanat2"))

#PC3 genes
FeaturePlot(esl, reduction = 'tsne', features = c("opn6a", "tulp1a", "opn1mw1", "dscamb", "arr3a", "rcvrn2", "guk1b", "si:ch211-81a5.8", "prom1b", "tdh2"))
FeaturePlot(esl, reduction = 'tsne', features = c("rs1a", "cuedc1a", "gngt2b", "rbp4l", "apln", "atp6v0cb", "si:ch211-285j22.3", "si:ch211-207l14.1", "tulp1b", "arl3l2"))

# FeaturePlot(msl, reduction = 'tsne', features = c("si:ch211-81a5.8","opn1mw1"))
# Ploting by age reveals that L cone cluster corresponds to 5dpf sample, only a few cones come from 2dpf (too early) and 2.5 - 3.5 dpf are intermingled
esl2 = esl
esl2 = SetIdent(esl2, value = "Sample")
DimPlot(esl2, reduction = "tsne", label = TRUE, pt.size = 2, label.size = 8)

FeaturePlot(esl, reduction = 'tsne', features = c("opn1sw2","foxq2",'opn1mw1',"six7","samd7",'opn1lw2','thrb',"rxrga","aanat2",'dscamb',"zgc:101100"))

FeaturePlot(esl, reduction = 'tsne', features = c('tbx2a','tbx2b','gdf6a','cnbpa','egr1', 'ndrg1a','ndrg1b','pbx1a','nfia','rx3','lin9'))
FeaturePlot(esl, reduction = 'tsne', features = c('prdm1b','nbeaa','flrt3'))
            
# ------------------------------------------------------------------------------------------
# plot layouts for direct export
#Subtype markers
p = FeaturePlot(esl, reduction = 'tsne', features = c("rho","nr2e3","gnat1","gnat2","arr3b","arr3a","thrb","opn1sw1", "opn1sw2",'opn1mw1','opn1lw2'),
                pt.size=1, order=TRUE, combine=FALSE)
lay = rbind(c(1,2,3,NA),c(4,5,6,7),c(8,9,10,11))
grid.arrange(grobs = p, layout_matrix = lay)
ps = arrangeGrob(grobs = p, layout_matrix = lay)
ggsave(ps, file="larval_eslMarkers.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")


#Clustering markers
p = FeaturePlot(esl, reduction = 'tsne', features = c("opn1sw2","foxq2",'opn1mw1',"six7","samd7",'opn1lw2','thrb',"rxrga","aanat2",'dscamb',"zgc:101100"),
                pt.size=1, order=TRUE, combine=FALSE)
lay = rbind(c(1,2,NA,NA),c(3,4,5,NA),c(6,7,8,9),c(10,11,NA,NA))
grid.arrange(grobs = p, layout_matrix = lay)
ps = arrangeGrob(grobs = p, layout_matrix = lay)
ggsave(ps, file="larval_eslMarkers.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")
# ------------------------------------------------------------------------------------------

# I think conclusion is that I should do DESeq2 comparison between foxq2+ eslCluster; can also separate M-cones and L-cones

new.cluster.ids <- c("eslL","eslM","eslS","eslPR","eslPR","eslL","eslL")
names(new.cluster.ids) <- levels(esl)
esl <- RenameIdents(esl, new.cluster.ids)
DimPlot(esl, reduction = "tsne", label = TRUE, pt.size = 1, label.size = 6) + eelTheme()
ggsave("larval_eslTSNEsubtypes.png", path=exportDir, width = 140, height = 105, units = "mm")

VlnPlot(esl, features = c("opn1sw1", "opn1sw2","foxq2", "opn1mw1","six7","arr3a", "opn1lw2","thrb","rxrga"))

# Explorign other DEGs in S-cone/foxq2 cluster
#DEGs that have adj_pval < 0.8
VlnPlot(esl, features = c("foxq2","opn1sw2","gngt2a","rcvrn2","thrb","si:ch73-1a9.3","samd7"))
# Upregulated
VlnPlot(esl, features = c("foxq2", "opn1sw2","pcdh11","si:ch211-106m9.1","gngt2a","skor1a","scinla","hspbp1","atp6ap1lb","glula","adipor2","krt91"))
VlnPlot(esl, features = c("spock3","atf3","pcdh17","CU693494.2","arl13a","tmx3","RDH13 (1 of many).2","pcdh19","ubb","hopx","si:rp71-56i13.6","jun","acp5a","ntf3","hunk"))
# Downregulated
VlnPlot(esl, features = c("rcvrn2","hmgn2","thrb","si:ch73-281n10.2","hmgb2a","si:ch73-1a9.3","tubb4b","h2afvb","hmga1a","hmgb2b","h1f0","cirbpa","fabp7a"))
VlnPlot(esl, features = c("opn1lw2","ptmab","histh1l","h3f3b.1","arr3a","tubb2b","aanat2","zbtb18","cfl1","hnrnpa0b","RPL41","nr2e3","serbp1a","nusap1","hnrnpaba","USMG5"))

# DESeq2_temp = FindMarkers(esl, c("eslS"), c("eslPR","eslL","eslM"), test.use="DESeq2");
# AvgExp_temp = AverageExpression(esl, features=rownames(DESeq2_temp))
# AvgExp_temp = AvgExp_temp$RNA
# DESeq2_temp = merge.data.frame(DESeq2_temp,AvgExp_temp, by="row.names", sort=FALSE)
# DESeq2_temp$baseMean = rowMeans(DESeq2_temp[,c("eslPR","eslL","eslM","eslS")])
# rownames(DESeq2_temp) = DESeq2_temp$Row.names
# DESeq2_temp = DESeq2_temp[,c("baseMean","avg_logFC","p_val","p_val_adj","pct.1","pct.2","eslPR","eslL","eslM","eslS")]
# PctExp_temp = DotPlot(msl, features=rownames(DESeq2_temp))
# PctNames = rownames(PctExp_temp$data)[1:(length(rownames(PctExp_temp$data))/4)]
# PctExp_temp = PctExp_temp$data$pct.exp
# PctExp_temp = matrix(PctExp_temp,nrow = length(PctExp_temp)/4, ncol=4);
# rownames(PctExp_temp) = PctNames
# colnames(PctExp_temp) = c("pctPR","pctL","pctM","pctS")
# DESeq2_temp = merge.data.frame(DESeq2_temp,PctExp_temp, by="row.names", sort=FALSE)
# rownames(DESeq2_temp) = DESeq2_temp$Row.names
# write.csv(DESeq2_temp,"./DESeq2_temp.csv")

# ------------------------------------------------------------------------------------------
# correlation plot for 2 genes
FeatureScatter(lsl,feature1 = "opn1sw2",feature2 = "foxq2")
FeatureScatter(lsl,feature1 = "rho",feature2 = "opn1lw2")
FeatureScatter(lsl,feature1 = "opn1lw2",feature2 = "syt5a")
# Look into rx1 and neuroD
# ------------------------------------------------------------------------------------------
# Will remerge with 2 strategies:
# have fine grained version with eslS - eslM - eslL - eslPR / mslM - mslL - mslPR / lslCones - lslRods
# then rename as just esl/msl/lslCones/lslRods as general one and probably kick out rods for most DESeq2 comparisons again

photoLarval = merge(x = lsl, y = c(msl,esl,lslRods), merge.data = TRUE,)
photoLarval = ScaleData(photoLarval)
# ordering identities for plotting and assigning colors
Idents(photoLarval) <- factor(x = Idents(photoLarval), levels = c("eslPR","eslL","eslM","eslS","mslPR","mslL","mslM","lslPR","R"))
lcolors = c("#dacd9a","#d69f9e","#80bc89","#8f9bcc","#dcc360","#c67271","#57cb69","#cca819","#474747")

# 'r' : '#747474',
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


# Identification of highly variable features (feature selection)
photoLarval <- FindVariableFeatures(photoLarval, selection.method = "vst", nfeatures = 200)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(photoLarval), 100)
top10


#run PCA
photoLarval <- RunPCA(photoLarval, features = VariableFeatures(object = photoLarval))

# photoLarval <- JackStraw(photoLarval, num.replicate = 100)
# photo <- ScoreJackStraw(photoLarval, dims = 1:20)
# JackStrawPlot(photoLarval, dims = 1:20)

# Or just use elbow plot (var explained)
ElbowPlot(photoLarval, ndims=50)

# Trying to find a pleasing tsne plot
# nFeatures = 250; ndim = 50 -> lays dev age as tsne1, so not bad
# nFeatures = 250; ndim = 10 -> curves developmental age and clusters start smoothing out
# nFeatures = 250; ndim = 25 -> almost back to original
# nFeatures = 250; ndim = 32|45 -> like it less than 25 because it inverts time. Could just flip it

# nFeatures = 100; ndim = 10|32|50 -> bad for tsne; rods get lost
# nFeatures = 400; ndim = 10|50 -> not great
# nFeatures = 400; ndim = 25|30 -> better but eslS is spread out and hard to id as clustered or rods get lost
# nFeatures = 400; ndim = 20 -> good but not better than original values

# In original clustering: nFeatures = 200; ndim = 45, res = 1; sticking to these, maybe doing -tsne_2
upDimLimit=45;
photoLarval <- FindNeighbors(photoLarval, dims = 1:upDimLimit)

# Omit this step to keep original identities before merging
# photoLarval <- FindClusters(photoLarval, resolution = 2)

DimPlot(photoLarval, reduction = "pca", label = TRUE, pt.size = 2, cols=lcolors)

# UMAP using the same PCA dimensions for prettier visualization
photoLarval <- RunUMAP(photoLarval, dims = 1:upDimLimit)
DimPlot(photoLarval, reduction = "umap", label = TRUE, pt.size = 2, cols=lcolors)

# or use tSNE (clustering can't separate some M cones from L cones)
photoLarval <- RunTSNE(photoLarval, dims = 1:upDimLimit)
DimPlot(photoLarval, reduction = "tsne", label = TRUE, pt.size = 1, label.size = 6, cols=lcolors) + eelTheme()
ggsave("larval_lslTSNEsubtypes.png", path=exportDir, width = 140, height = 105, units = "mm")

# stashing semi-manual clustering in metadata to be able to swap in and out
photoLarval$subClusters <- Idents(object = photoLarval)

# Plot by age (and getting rid of replicate information for clarity)
photoLarval2 = SetIdent(photoLarval, value = "Sample")
new.cluster.ids <- c("5dpf","4dpf","3dpf","4dpf","4dpf","2.5dpf","2dpf","2.5dpf")
names(new.cluster.ids) <- levels(photoLarval2)
photoLarval2 <- RenameIdents(photoLarval2, new.cluster.ids)
Idents(photoLarval2) <- factor(x = Idents(photoLarval2), levels = c("2dpf","2.5dpf","3dpf","4dpf","5dpf"))
DimPlot(photoLarval2, reduction = "tsne", label = TRUE, pt.size = 1, label.size = 6) + eelTheme()
ggsave("larval_lslTSNEsubtypesAge.png", path=exportDir, width = 140, height = 105, units = "mm")

# give this information to merged object
photoLarval$devAge <- Idents(object = photoLarval2)

# Remove crappy subtyping and just stick with early/mid/late
levels(Idents(photoLarval))
new.cluster.ids <- c("eslPR","eslPR","eslPR","eslPR","mslPR","mslPR","mslPR","lslPR","lslR")
names(new.cluster.ids) <- levels(photoLarval)
photoLarval <- RenameIdents(photoLarval, new.cluster.ids)
DimPlot(photoLarval, reduction = "tsne", label = TRUE, pt.size = 1, label.size = 6) + eelTheme()
ggsave("larval_lslTSNEnosubtypes.png", path=exportDir, width = 140, height = 105, units = "mm")
# give this information to merged object
photoLarval$devStage <- Idents(object = photoLarval)

# swap to subclusters again
Idents(object = photoLarval) = photoLarval$subClusters
DimPlot(photoLarval, reduction = "tsne", label = TRUE, pt.size = 2, label.size = 8, cols = lcolors) + NoLegend()




# exploring some alternative plots
# this one is hard to read because 0 peak is always so so big
# RidgePlot(photoLarval, features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))
# this does not add much here and it's too crowded
# VlnPlot(photoLarval, features = c("nr2e3","nr2f6b","crx",'syt5a','syt5b','gnat2','arr3a','arr3b',"pde6g","pde6h", "neurod1"), split.by = "devAge") #cone markers
# if clustering wasn't already basically age, this could be helpful
# FeaturePlot(photoLarval, reduction = 'tsne', features = c("syt5a"), split.by = "devStage")
# Cool but labels need to be white!
# DimPlot(photoLarval2, reduction = "tsne", label = TRUE, pt.size = 2, label.size = 8) + NoLegend() + DarkTheme()


                                                          
# Most discriminating features from previous analysis
# opsins
FeaturePlot(photoLarval, reduction = 'tsne', features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))
# l/m/e markers
FeaturePlot(photoLarval, reduction = 'tsne', features = c("nr2e3","nr2f6b","crx",'syt5a','syt5b','gnat1','gnat2','arr3a','arr3b',"pde6g","pde6h", "neurod1")) #cone markers
# subclusters: rods/S+foxq2
FeaturePlot(photoLarval, reduction = 'tsne', features = c("rho","saga","nr2e3","foxq2","opn1sw2","pcdh11",'opn1lw2',"thrb","rxrga","rxrgb","opn1mw1"))

#other genes of interest
FeaturePlot(photoLarval, reduction = 'tsne', features = c("prdm1a","prdm1b","tbx2a","tbx2b","nr2e3","nr2f1b","nr2f6b",'foxq2',"lmo4a","apln","aplnra","aplnrb","skor1a"))
FeaturePlot(photoLarval, reduction = 'tsne', features = c("fabp7a","ntf3","slc1a8b","ndrg1a","ndrg1b","nme3","samd7","ankrd33ab","cplx4a","cplx4c"))


VlnPlot(photoLarval, features = c("crx","nr2e3","nr2f6b",'syt5a','syt5b',"gnat1",'gnat2','arr3a','arr3b',"pde6g","pde6h", "neurod1"), cols = lcolors, ncol = 3) #cone markers
VlnPlot(photoLarval, features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'), cols = lcolors) #subtype markers
VlnPlot(photoLarval, features = c("saga", "nr2e3",'foxq2','opn1sw2','opn1lw2','thrb','rxrga'), cols = lcolors) #Subtype markers

VlnPlot(photoLarval, features = c("tbx2a", "tbx2b",'sema7a','sema3fa','sema6d',"elovl4a","elovl4b",'efna1b','syt5a','syt5b','pcdh10a','pcdh11'), cols = lcolors) #zfA_DEGs


# cell counter
table(Idents(photoLarval))
# 2021
# eslPR  eslL  eslM  eslS mslPR  mslL  mslM lslPR     lslR 
# 167   109    63    73   141    96    91   372    41 

# eslPR mslPR lslPR     lslR 
# 412   328   372    41 

# Jul 2022
# eslPR  eslL  eslM  eslS mslPR  mslL  mslM lslPR     R 
# 85   178    87    62   110    73   219   297    42 

# eslPR mslPR lslPR  lslR 
# 412   402   297    42 


saveRDS(photoLarval, file = "./cones_LarvalReAnalyzed.rds")
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# -------------------------------------------------------------------
rm(list=ls())
# -------------------------------------------------------------------
# CAN BE RESTARTED HERE 
# Load the previously saved larval cone dataset after trying to identify photoreceptor subtypes and not succeeding well
photoLarval = readRDS(file = "./cones_LarvalReAnalyzed.rds")
# Load the previously saved adult cone dataset after trying to more or less clustering by photoreceptor subtypes
photo = readRDS("./cones_AdultAll.rds")

# make sure unreliable subtype info is removed
Idents(object = photoLarval) = photoLarval$devStage

# getting % and Mean expression for all genes into csv.
# DoHeatmap and DotPlot use the @scale.data slot for average expression display, which z-scored expression values (for example, as those used in PCA).
# Cells with a value > 0 represent cells with expression above the population mean (a value of 1 would represent cells with expression 1SD away from the population mean). Hope that helps!
# This is very detailed explanation: https://github.com/satijalab/seurat/issues/2798

# Average expression in non-log scale for larvae:
avgExp = AverageExpression(photoLarval, slot="counts")
avgExp = avgExp$RNA
# avgExp = head(avgExp,10)
nI = length(levels(Idents(photoLarval)))
colnames(avgExp) = paste("avg", levels(Idents(photoLarval)), sep = "")


# Average expression in non-log scale for adults:
avgExpAdult = AverageExpression(photo, slot="counts", features=rownames(avgExp))
avgExpAdult = avgExpAdult$RNA
# avgExpAdult = head(avgExpAdult,10)
nIAd = length(levels(Idents(photo)))
colnames(avgExpAdult) = paste("avg", levels(Idents(photo)), sep = "")
avgExpAdult$avgMean = rowMeans(avgExpAdult)
avgExpAdult[,c(ncol(avgExpAdult),1:(ncol(avgExpAdult)-1))]

# Add Adult column into larval data frame
avgExp$avgadPR = avgExpAdult$avgMean
# Calculate mean expression across all ages
avgExp$avgMean = rowMeans(avgExp)
avgExp[,c(ncol(avgExp),1:(ncol(avgExp)-1))]


# Percent expression can be obtained through DotPlot, but not Average expression, as this is log z-cored normalized data # avgExpZ = dPlot$data$avg.exp
# Percent expression in larvae:
dPlot = DotPlot(photoLarval, features=rownames(avgExp))
pctExp = dPlot$data$pct.exp
pctExp = as.data.frame(matrix(pctExp,nrow = length(pctExp)/nI, ncol=nI));

rownames(pctExp) = rownames(avgExp)
colnames(pctExp) = paste("pct", levels(Idents(photoLarval)), sep = "")


# Percent expression in adults:
dPlotAd = DotPlot(photo, features=rownames(avgExp))
pctExpAdult = dPlotAd$data$pct.exp
pctExpAdult = as.data.frame(matrix(pctExpAdult,nrow = length(pctExpAdult)/nIAd, ncol=nIAd));

rownames(pctExpAdult) = rownames(avgExp)
colnames(pctExpAdult) = paste("pct", levels(Idents(photo)), sep = "")
pctExpAdult$pctMean = rowMeans(pctExpAdult)
pctExpAdult[,c(ncol(pctExpAdult),1:(ncol(pctExpAdult)-1))]

# Add Adult column into larval data frame
pctExp$pctadPR = pctExpAdult$pctMean
# Calculate mean percent expression across all ages
pctExp$pctMean = rowMeans(pctExp)
pctExp[,c(ncol(pctExp),1:(ncol(pctExp)-1))]


# cones_Hoang <- cones_Hoang[ order((cones_Hoang$pctMean)), ]

conesLarval_Hoang = cbind(pctExp,avgExp)
conesAdult_Hoang = cbind(pctExpAdult,avgExpAdult)



# these are duplicated genes with lower vs upper case
dups=toupper(c("arid5b","asph", "aste1", "crip2", "ctbp1", "dab2", "eif1b", "flnb", "frmd7", "galnt10", "grxcr1", "gse1", "hist1h4l", "hspb11", "kif1c", "lamp1", "maf1", "pamr1", "pcdh20", "pde6h", "phlpp2", "psmb10", "ptp4a3", "reep6", "rfesd", "rgs9bp", "rps17", "shank2", "slc16a7", "slc25a10", "slc6a13", "slc9a1", "srbd1", "tatdn3", "tenm3", "tmem178b", "tmem241", "tom1l2", "tp53inp2", "trappc9", "tsc22d3", "ube2o", "zc3h12a", "znf423"))
for(i in 1:length(dups)) {
   rownames(conesLarval_Hoang)[rownames(conesLarval_Hoang) == dups[i]] = paste(dups[i],"_ii",sep="")
   rownames(conesAdult_Hoang)[rownames(conesAdult_Hoang) == dups[i]] = paste(dups[i],"_ii",sep="")
}
# now everything can be lower cased
rownames(conesLarval_Hoang) = tolower(rownames(conesLarval_Hoang))
rownames(conesAdult_Hoang) = tolower(rownames(conesAdult_Hoang))



conesLarval_Hoang = conesLarval_Hoang[c("pctMean","pcteslPR","pctmslPR","pctlslPR","pctadPR","pctlslR","avgMean","avgeslPR","avgmslPR","avglslPR","avgadPR","avglslR")]
conesAdult_Hoang = conesAdult_Hoang[c("pctMean","pctUV","pctS","pctM","pctL","avgMean","avgUV","avgS","avgM","avgL")]


scangenes=c('syt5a','syt5b','nr2e3','tbx2a','tbx2b','arr3a','arr3b','opn1lw1','opn1lw2');
conesLarval_Hoang[scangenes,]

conesAdult_Hoang[scangenes,]


# this is pretty good, but adult data set is missing rods
# then larval rods are contaminated by cones; and adult cones are contaminated by rods
# also could integrate photoreceptor progenitor cluster for completion
# still need to manually add symbol to first column
write.csv(conesLarval_Hoang,"./conesLarval_Hoang.csv",quote=FALSE)
write.csv(conesAdult_Hoang,"./conesAdult_Hoang.csv",quote=FALSE)
# -------------------------------------------------------------------
# Adding photoreceptor progenitor cluster (and exploring briefly)

# make sure unreliable subtype info in larval samples is removed
Idents(object = photoLarval) = photoLarval$devStage
# or include subtype info?
# Idents(object = photoLarval) = photoLarval$subClusters
# check labels
levels(Idents(photoLarval))
levels(Idents(photoAll))
levels(Idents(pbmc))

#removing subM clusters, they have been DeSEQ'd anyway
new.cluster.ids <- c("R","UV","S","M","M","M","L")
names(new.cluster.ids) <- levels(photoAll)
photoAll <- RenameIdents(photoAll, new.cluster.ids)
# merge clusters
prp = subset(pbmc, idents = c("PhotoProgenitor"));
new.cluster.ids <- c("PRP")
names(new.cluster.ids) <- levels(prp)
prp <- RenameIdents(prp, new.cluster.ids)

photoDev = merge(x = prp, y = c(photoLarval,photoAll))
levels(Idents(photoDev))

# ordering identities for plotting and assigning colors
Idents(photoDev) <- factor(x = Idents(photoDev), levels = c("PRP","eslPR","mslPR","lslPR","lslR","R","UV","S","M","L"))
prcolors = c("#dfdac8","#dacd9a","#dcc360","#cca819","#a3a3a3","#7d7d7d","#B540B7","#4669F2","#04CD22","#CC2C2A")
levels(Idents(photoDev))


photoDev = ScaleData(photoDev)
photoDev <- FindVariableFeatures(photoDev, selection.method = "vst", nfeatures = 200)

photoDev <- RunPCA(photoDev, features = VariableFeatures(object = photoDev))
ElbowPlot(photoDev, ndims=50)
# if ever wanted to recluster again
# photo <- FindNeighbors(photoDev, dims = 1:upDimLimit)
# photo <- FindClusters(photoDev, resolution = 1)



upDimLimit = 30;
photoDev <- RunUMAP(photoDev, dims = 1:upDimLimit)
DimPlot(photoDev, reduction = "umap", label=TRUE, pt.size = 1, cols=prcolors) + NoLegend()
# or use tSNE (clustering can't separate some M cones from L cones)
photoDev <- RunTSNE(photoDev, dims = 1:upDimLimit)
ps = DimPlot(photoDev, reduction = "tsne", label=TRUE, repel=TRUE,pt.size = 2, cols=prcolors) +eelTheme() +NoLegend()
ps
ggsave(ps, file="dev_Clusters01.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

FeaturePlot(photoDev, reduction = 'tsne', features = c("rho","opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2'),
            pt.size=1, order=FALSE, combine=TRUE)

# "Known" PRP genes
FeaturePlot(photoDev, reduction = 'tsne', features = c("notch1a","ccnd1","cdk1","fgf19","vsx2","hopx","lhx2a","lhx2b","rlbp1a","rlbp1b"),
            pt.size=1, order=FALSE, combine=TRUE)

# Me trying to map out things
FeaturePlot(photoDev, reduction = 'tsne', features = c("crx","nr2e3","neurod1","otx5","sox2","prdm1a","prdm1b","syt5a","syt5b","foxq2","skor1a"),
            pt.size=1, order=FALSE, combine=TRUE)

FeaturePlot(photoDev, reduction = 'tsne', features = c("opn1sw2","foxq2","skor1a","nrtn"),
            pt.size=1, order=TRUE, combine=TRUE)

# Also checking "known rod markers"
FeaturePlot(photoDev, reduction = 'tsne', features = c("prkca","sebox","guca1b","nrl","rho","rom1a","rom1b","nr2e3","saga"),
            pt.size=1, order=FALSE, combine=TRUE)


DotPlot(photoDev,features = c("ntf3","spock3","hs3st3l","sema3fb","si:busm1-57f23.1","eya2","c1galt1a","itga1","slc32a1","dok6"))



saveRDS(photoDev, file = "./photoreceptors_Development.rds")

# -------------------------------------------------------------------
# CAN BE RESTARTED HERE  DEVELOPMENTLINK
photoDev = readRDS(file = "./photoreceptors_Development.rds")

# getting % and Mean expression for all genes into csv.
# DoHeatmap and DotPlot use the @scale.data slot for average expression display, which z-scored expression values (for example, as those used in PCA).
# Cells with a value > 0 represent cells with expression above the population mean (a value of 1 would represent cells with expression 1SD away from the population mean). Hope that helps!
# This is very detailed explanation: https://github.com/satijalab/seurat/issues/2798

# Average expression in non-log scale for larvae:
avgExp = AverageExpression(photoDev, slot="counts")
avgExp = avgExp$RNA
# avgExp = head(avgExp,10)
nI = length(levels(Idents(photoDev)))
colnames(avgExp) = paste("avg", levels(Idents(photoDev)), sep = "")
# Calculate mean expression across all ages
avgExp$avgMean = rowMeans(avgExp)

# Average expression in non-log scale for adults to complete dev profile without subtype and also have by subtype:
avgExpAdult = AverageExpression(photoAll, slot="counts", features=rownames(avgExp))
avgExpAdult = avgExpAdult$RNA
# avgExpAdult = head(avgExpAdult,10)
nIAd = length(levels(Idents(photoAll)))

# Add average adult column into data frame
avgExp$avgadPR = rowMeans(avgExpAdult)
head(avgExp,10)

# Percent expression can be obtained through DotPlot, but not Average expression, as this is log z-cored normalized data # avgExpZ = dPlot$data$avg.exp
# Percent expression in larvae:
dPlot = DotPlot(photoDev, features=rownames(avgExp))
pctExp = dPlot$data$pct.exp
pctExp = as.data.frame(matrix(pctExp,nrow = length(pctExp)/nI, ncol=nI));

rownames(pctExp) = rownames(avgExp)
colnames(pctExp) = paste("pct", levels(Idents(photoDev)), sep = "")

# Calculate mean percent expression across all ages
pctExp$pctMean = rowMeans(pctExp)
head(pctExp,10)

# Percent expression in adults:
dPlotAd = DotPlot(photoAll, features=rownames(avgExpAdult))
pctExpAdult = dPlotAd$data$pct.exp
pctExpAdult = as.data.frame(matrix(pctExpAdult,nrow = length(pctExpAdult)/nIAd, ncol=nIAd));

rownames(pctExpAdult) = rownames(avgExpAdult)
colnames(pctExpAdult) = paste("pct", levels(Idents(photoAll)), sep = "")

# Add Adult column into larval data frame
pctExp$pctadPR =  rowMeans(pctExpAdult)
head(pctExp,10)

# cones_Hoang <- cones_Hoang[ order((cones_Hoang$pctMean)), ]

conesDev_Hoang = cbind(pctExp,avgExp)
head(conesDev_Hoang,10)


# these are duplicated genes with lower vs upper case
dups=toupper(c("arid5b","asph", "aste1", "crip2", "ctbp1", "dab2", "eif1b", "flnb", "frmd7", "galnt10", "grxcr1", "gse1", "hist1h4l", "hspb11", "kif1c", "lamp1", "maf1", "pamr1", "pcdh20", "pde6h", "phlpp2", "psmb10", "ptp4a3", "reep6", "rfesd", "rgs9bp", "rps17", "shank2", "slc16a7", "slc25a10", "slc6a13", "slc9a1", "srbd1", "tatdn3", "tenm3", "tmem178b", "tmem241", "tom1l2", "tp53inp2", "trappc9", "tsc22d3", "ube2o", "zc3h12a", "znf423"))
for(i in 1:length(dups)) {
   rownames(conesDev_Hoang)[rownames(conesDev_Hoang) == dups[i]] = paste(dups[i],"_ii",sep="")
}
# now everything can be lower cased
rownames(conesDev_Hoang) = tolower(rownames(conesDev_Hoang))
# also replacing (1 of many) with .1
rownames(conesDev_Hoang) = gsub("(1 of many)",".1",rownames(conesDev_Hoang))


conesDev_Hoang = conesDev_Hoang[c("avgMean","avgPRP","avgeslPR","avgmslPR","avglslPR","avgadPR","avglslR","avgR","avgUV","avgS","avgM","avgL",
                                  "pctMean","pctPRP","pcteslPR","pctmslPR","pctlslPR","pctadPR","pctlslR","pctR","pctUV","pctS","pctM","pctL")]
# sort by avgExpression
conesDev_Hoang=conesDev_Hoang[order(-conesDev_Hoang$avgMean),]
head(conesDev_Hoang,100)

scangenes=c('syt5a','syt5b','nr2e3','tbx2a','tbx2b','arr3a','arr3b','opn1lw1','opn1lw2');
conesDev_Hoang[scangenes,]

# this is as far as I can take this
# take aways: 
# larval rods and cones are all contaminated by L cones
# and adult cones are contaminated by rods

write.csv(conesDev_Hoang,"./conesDev_Hoang.csv",quote=FALSE)
            
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Some plots for lab meeting
#EXPLORING GENES ID'D IN OUR OWN RNAseq
ps = DotPlot(photoLarval, features = c("sema3aa","sema3ab","sema3b","sema3bl","sema3c","sema3d","sema3e","sema3fa","sema3fb","sema3ga","sema3gb","sema3h","sema4aa","sema4ab","sema4ba","sema4bb","sema4c","sema4d","sema4e","sema4f","sema4ga","sema4gb","sema5a","sema5ba","sema5bb","sema6a","sema6ba","sema6bb","sema6d","sema6dl","sema6e","sema7a")) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ps
ggsave(ps, file="larval_coneMarkersDot04_sema.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

ps = DotPlot(photoLarval, features = c("vamp1","vamp2","vamp3","vamp4","vamp5","vamp8")) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ps
ggsave(ps, file="larval_coneMarkersDot04_vamp.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

ps = DotPlot(photoLarval, features = c("syt1a","syt1b","syt2a","syt3","syt4","syt5a","syt5b","syt6a","syt6b","syt7a","syt7b")) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ps
ggsave(ps, file="larval_coneMarkersDot04_syt.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

ps = DotPlot(photoLarval, features = c("efna1a","efna1b","efna2a","efna2b","efna3a","efna3b","efna5a","efna5b")) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ps
ggsave(ps, file="larval_coneMarkersDot04_efna.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")


ps = DotPlot(photoLarval, features = c("tbx1","tbx15","tbx16","tbx18","tbx20","tbx21","tbx22","tbx2a","tbx2b","tbx3a","tbx3b","tbx4","tbx5a","tbx5b","tbx6")) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ps
ggsave(ps, file="larval_coneMarkersDot04_tbx.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")


ps = DotPlot(photoLarval, features = c("nr2c1","nr2c2","nr2e1","nr2e3","nr2f1a","nr2f1b","nr2f2","nr2f5","nr2f6a","nr2f6b")) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ps
ggsave(ps, file="larval_coneMarkersDot04_nr2.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

ps = DotPlot(photoLarval, features = c("ntn1a","ntn1b","ntn2","ntn4","ntn5","ntng1a","ntng2a","ntng2b")) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ps
ggsave(ps, file="larval_coneMarkersDot04_ntn.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

# Rod genes
ps = DotPlot(photoLarval, features = c("rho","gnat1","dhrs13l1","lingo1b","cabp4","saga","sagb","gnb1b","rgs9bp","eno2","guca1b","grk1a","rom1b","cplx4c","lrrn1","ncs1a","sgce","kcnv2a","pdca")) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ps
ggsave(ps, file="larval_coneMarkersDot04_RodGenes.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")
# Cone Genes
ps = DotPlot(photoLarval, features = c("slc1a8b","rgs9a","slc25a24","drd4b","si:busm1-57f23.1","kera","nexn","clic1","gnat2","tgif1","crhbp","kcnv2b","anks1b","clul1")) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ps
ggsave(ps, file="larval_coneMarkersDot04_ConeGenes.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

# UV+S Genes
ps = DotPlot(photoLarval, features = c("grk7b","skor1a","pcdh11","srgap3","ablim1a","nav2a","jam2a","chl1a","s100z","foxq2","sh3bp5b","kcnk1a","ntng2b","nxph1","pik3r3b","myl4","pacrg","fah")) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ps
ggsave(ps, file="larval_coneMarkersDot04_USGenes.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

# M+L Genes
ps = DotPlot(photoLarval, features = c("slc32a1","apln","lactbl1b","nrtn","pcdh10a","myo7aa")) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ps
ggsave(ps, file="larval_coneMarkersDot04_MLGenes.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")
# Mup/Ldown Genes
ps = DotPlot(photoLarval, features = c("dok6","spock3","lrrfip1a","sema3fb","itga1","auts2a","esama","esamb","plxnb1a","cgnb","phf19","lrrc20")) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ps
ggsave(ps, file="larval_coneMarkersDot04_MnotLGenes.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")
# Lup/Mdown Genes
ps = DotPlot(photoLarval, features = c("slc32a1","s100v2","smad5","snap25a","ttyh2l","arhgap11a","ggctb","fbxo32","pcdh10a","abracl")) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ps
ggsave(ps, file="larval_coneMarkersDot04_LnotMGenes.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")
# UVup/Sdown Genes
ps = DotPlot(photoLarval, features = c("myl4","ttyh2l","lhx1a","rx3","itgb1bp1")) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ps
ggsave(ps, file="larval_coneMarkersDot04_UnotSGenes.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")
# Sup/UVdown Genes
ps = DotPlot(photoLarval, features = c("fkbp5","prss23","chkb","mpzl2b","nr1d1","camk2d1","frmpd2","foxo1a","prom1a")) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ps
ggsave(ps, file="larval_coneMarkersDot04_SnotUGenes.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

# TFs in DEGs
ps = DotPlot(photoLarval, features = c("tbx2a","tbx2b","otx5","tgif1","crhbp","six7","six6b","ndrg1b","egr1","nr1d1","sall1a","skor1a","foxq2","lhx1a","ntf3","rxrga","thrb","fgf1b","eya2","sox6","sox4b","pbx1a","hmgb3a","fbxo21","tfe3a")) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ps
ggsave(ps, file="larval_coneMarkersDot04_TFs.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

# Exploring all celltypes
ps = DotPlot(pbmc, features = c("tbx2a","tbx2b","otx5","tgif1","crhbp","six7","six6b","ndrg1b","egr1","nr1d1","sall1a","skor1a","foxq2","lhx1a","ntf3","rxrga","thrb","fgf1b","eya2","sox6","sox4b","pbx1a","hmgb3a","fbxo21","tfe3a")) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ps

ps = DotPlot(pbmc, features = c("nr2e3","crx","syt5b")) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ps

ps = DotPlot(pbmc, features = c("vamp1","vamp2")) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ps

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Running DESeq2 (some comparisons) and saving as csv
# rods vs cones
tempC = c("R","UV","S","M","L")
nC = length(tempC);
temp = subset(photoDev, idents = tempC)
DESeq2_R = FindMarkers(temp, c("R"), c("UV","S","M","L"), test.use="DESeq2");
AvgExp_R = AverageExpression(temp, features=rownames(DESeq2_R))
AvgExp_R = AvgExp_R$RNA
DESeq2_R = merge.data.frame(DESeq2_R,AvgExp_R, by="row.names", sort=FALSE)
DESeq2_R$baseMean = rowMeans(DESeq2_R[,tempC])
rownames(DESeq2_R) = DESeq2_R$Row.names
DESeq2_R = DESeq2_R[,c("baseMean","avg_logFC","p_val","p_val_adj","pct.1","pct.2","R","UV","S","M","L")]
PctExp_R = DotPlot(temp, features=rownames(DESeq2_R))
PctNames = rownames(PctExp_R$data)[1:(length(rownames(PctExp_R$data))/nC)]
PctColNames = paste("pct",levels(PctExp_R$data$id),sep="")
PctExp_R = PctExp_R$data$pct.exp
PctExp_R = matrix(PctExp_R,nrow = length(PctExp_R)/nC, ncol=nC);
rownames(PctExp_R) = PctNames
colnames(PctExp_R) = PctColNames
DESeq2_R = merge.data.frame(DESeq2_R,PctExp_R, by="row.names", sort=FALSE)
rownames(DESeq2_R) = DESeq2_R$Row.names
DESeq2_R = DESeq2_R[,c("baseMean","avg_logFC","p_val","p_val_adj","pct.1","pct.2","R","UV","S","M","L","pctR","pctUV","pctS","pctM","pctL")]
write.csv(DESeq2_R,"./DESeq2_RvsC.csv")

# eslsw2 vs eslPR
Idents(object = photoLarval) = photoLarval$subClusters
tempC = c("eslS","eslPR","eslM","eslL")
nC = length(tempC);
temp = subset(photoLarval, idents = tempC)
DESeq2_eS = FindMarkers(temp, c("eslS"), c("eslPR","eslM","eslL"), test.use="DESeq2");
AvgExp_eS = AverageExpression(temp, features=rownames(DESeq2_eS))
AvgExp_eS = AvgExp_eS$RNA
DESeq2_eS = merge.data.frame(DESeq2_eS,AvgExp_eS, by="row.names", sort=FALSE)
DESeq2_eS$baseMean = rowMeans(DESeq2_eS[,tempC])
rownames(DESeq2_eS) = DESeq2_eS$Row.names
DESeq2_eS = DESeq2_eS[,c("baseMean","avg_logFC","p_val","p_val_adj","pct.1","pct.2","eslS","eslM","eslL","eslPR")]
PctExp_eS = DotPlot(temp, features=rownames(DESeq2_eS))
PctNames = rownames(PctExp_eS$data)[1:(length(rownames(PctExp_eS$data))/nC)]
PctColNames = paste("pct",levels(PctExp_eS$data$id),sep="")
PctExp_eS = PctExp_eS$data$pct.exp
PctExp_eS = matrix(PctExp_eS,nrow = length(PctExp_eS)/nC, ncol=nC);
rownames(PctExp_eS) = PctNames
colnames(PctExp_eS) = PctColNames
DESeq2_eS = merge.data.frame(DESeq2_eS,PctExp_eS, by="row.names", sort=FALSE)
rownames(DESeq2_eS) = DESeq2_eS$Row.names
DESeq2_eS = DESeq2_eS[,c("baseMean","avg_logFC","p_val","p_val_adj","pct.1","pct.2","eslS","eslM","eslL","eslPR","pcteslS","pcteslM","pcteslL","pcteslPR")]
write.csv(DESeq2_eS,"./DESeq2_eSvsePR.csv")
