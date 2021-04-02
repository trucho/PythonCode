# Enter commands in R (or R studio, if installed)
# install.packages('Seurat')
# install.packages("reticulate")
# install.packages("umap")
# install.packages("openssl")
# install.packages('BiocManager')
# BiocManager::install('limma')

Sys.setenv(RETICULATE_PYTHON = "/Users/angueyraaristjm/anaconda3/bin/python")

# reticulate::use_python("/Users/angueyraaristjm/anaconda3/bin/python")
reticulate::use_python("/usr/local/bin/python3")
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
setwd("/Users/angueyraaristjm/Documents/LiMolec/otherRNAseq/zfRetDev10x_Xu2020/ret72hpf/")
directory <- "/Users/angueyraaristjm/Documents/LiMolec/otherRNAseq/zfRetDev10x_Xu2020/ret72hpf/"
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
# Load the 10x dataset (2826 cells)
pbmcA.data <- Read10X(data.dir = "/Users/angueyraaristjm/Documents/LiMolec/otherRNAseq/zfRetDev10x_Xu2020/ret72hpf/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmcA <- CreateSeuratObject(counts = pbmcA.data, project = "zf10X_72hpfA", min.cells = 3, min.features = 200)
pbmcA

# Load the replciate 10x dataset (405 cells)
pbmcB.data <- Read10X(data.dir = "/Users/angueyraaristjm/Documents/LiMolec/otherRNAseq/zfRetDev10x_Xu2020/ret72hpfB/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmcB <- CreateSeuratObject(counts = pbmcB.data, project = "zf10X_72hpfB", min.cells = 3, min.features = 200)
pbmcB

# merge the 2 datasets (6882)
pbmc = merge(pbmcA, y = pbmcB, add.cell.ids = c("72hA","72hB"), project = "zf10X_72hpf")
pbmc

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
# pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2

# remove unwanted cells
# pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200)
# normalize data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000) # this is to get natural log-transformed data using log10
# pbmc <- NormalizeData(pbmc, normalization.method = "RC", scale.factor = 1e6) #this is counts per million

# Identification of highly variable features (feature selection)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", dispersion.function = LogVMR, nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
top10
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#scale data for PCA
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

#run PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))


# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca")

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)


# Determine number of clusters using Macosko, 2015 (random permutation of 1% of data and rerun PCA iteritavely)
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

JackStrawPlot(pbmc, dims = 1:20)

# Or just use elbow plot (var explained)
ElbowPlot(pbmc)

# Cluster cells
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.6)

# Look at cluster IDs of the first 5 cells
# head(Idents(pbmc), 5)

# or use tSNE or UMAP using the same PCA dimensions for prettier visualization
# pbmc <- RunUMAP(pbmc, dims = 1:10)
# DimPlot(pbmc, reduction = "umap")
# saveRDS(pbmc, file = "../output/umap.rds")
pbmc <- RunTSNE(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "tsne", label=TRUE)
saveRDS(pbmc, file = "./tSNE.rds")

# 
# # Finding differentially expressed features (cluster biomarkers)
# # find all markers of cluster 1
# cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
# head(cluster1.markers, n = 5)
# 
# # find all markers distinguishing cluster 5 from clusters 0 and 3
# cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
# head(cluster5.markers, n = 5)


# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)



# check counts in each cluster
VlnPlot(pbmc, features = c("gnat2", "rho",'elovl4b','gnat1','si:busm1-57f23.1')) #log norm values
VlnPlot(pbmc, features = c("gnat2", "rho",'elovl4b','gnat1','si:busm1-57f23.1'), slot = "counts", log = TRUE) #raw counts
RidgePlot(pbmc, features = c("gnat2", "rho",'elovl4b','gnat1','si:busm1-57f23.1'))
DotPlot(pbmc, features = c("gnat2", "rho",'elovl4b','gnat1','si:busm1-57f23.1'))

DotPlot(pbmc, features = c("sema7a","efna1b","ntng2b","syt1a","syt1b","syt5a","syt5b","tbx2a","tbx2b","eml1"))
# check counts in tSNE space
# FeaturePlot(pbmc, features = c("gnat2",'gngt2b', "gnat1",'rho','si:busm1-57f23.1','opn1lw2', 'crx'))

# heatmap of top genes
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()


# Using published DEGs in FigS1. to redefine clusters
# RGCs = 11 & 4 (onlyA dataset 5 & 7)
FeaturePlot(pbmc, features = c("islr2","rbpms2a",'pou4f1'))

# ACs = 1 & 6 & 10 (onlyA dataset 0 & 1 & 3)
FeaturePlot(pbmc, features = c("slc32a1","pax10",'tfap2a','gad2','slc6a1b'))

# HCs = 13 & 15 & 16 (why so separated? HC subtypes?) (onlyA dataset 9)
FeaturePlot(pbmc, features = c("rem1","CR361564.1",'prkar2ab','rprmb','opn4.1','atp1b1a'))

# BCs = 7? & 9 (onlyA dataset 4)
FeaturePlot(pbmc, features = c("samsn1a","vsx1",'gnb3a','pcp4l1','neurod4','gnao1b'))

# MG = 14 (onlyA dataset 12)
FeaturePlot(pbmc, features = c("cavin2a","rhbg",'slc1a2b','atp1a1b'))

# CMZ = 0 & 3 & 2 (2 must be photoreceptor progenitors) (onlyA dataset 11, 8, 2 (11 must be photoreceptor progenitors))
FeaturePlot(pbmc, features = c("hmgb2b","stmn1a",'hmgn2','rrm2.1'))

# Photoreceptors = 5 & 8 & 12 (onlyA dataset 6 & 10)
FeaturePlot(pbmc, features = c("tulp1a","elovl4b",'ckmt2a','gngt2a'))
FeaturePlot(pbmc, features = c("opn1sw1", "opn1sw2",'opn1mw1','gnat2','opn1lw2','si:busm1-57f23.1','rho','saga'))

# 0   "CMZ"
# 1   "ACs"
# 2   "CMZ" (or "Photo")probably Photoreceptor progenitors)
# 3   "CMZ"
# 4   "RGCs"
# 5   "Photo"
# 6   "ACs"
# 7   "BCs"
# 8   "Photo"
# 9   "BCs"
# 10  "ACs"
# 11  "RGCs"
# 12  "Photo"
# 13  "HCs"
# 14  "MG"
# 15  "HCs"
# 16  "HCs"

# Assigning cell type identity to clusters if markers (it's basically just renaming)
# only A dataset 
# new.cluster.ids <- c("ACs","ACs","CMZ","ACs","BCs","RGCs","Photo","RGCs","CMZ","HCs","Photo","CMZ","MG")
# merged datasets
new.cluster.ids <- c("CMZ", "ACs", "CMZ", "CMZ", "RGCs", "Photo", "ACs", "BCs", "Photo", "BCs", "ACs", "RGCs", "Photo", "HCs", "MG", "HCs", "HCs")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "tsne", label = TRUE, pt.size = 2) + NoLegend()
saveRDS(pbmc, file = "./allCells_reClsutered.rds")
# ------------------------------------------------------------------------------------------

FeaturePlot(pbmc, reduction = 'tsne', features = c("crx","otx5", "tbx2a",'tbx2b'), pt.size = 2)
# FeaturePlot(pbmc, reduction = 'umap', features = c("crx","otx5", "tbx2a",'tbx2b'), pt.size = 2)
VlnPlot(pbmc, features = c("crx","otx5", "tbx2a",'tbx2b'))

# ------------------------------------------------------------------------------------------

# Clear all plots -------------------------------------------------------------------
try(dev.off(),silent=TRUE)
# Clear environment -------------------------------------------------------------------
rm(list=ls())

# ------------------------------------------------------------------------------------------
pbmc <- readRDS(file = "./allCells_reClsutered.rds")
# Subclustering photoreceptors
photo <- subset(pbmc, idents = "Photo")

# Identification of highly variable features (feature selection)
photo <- FindVariableFeatures(photo, selection.method = "vst", nfeatures = 500)
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
# photo <- JackStraw(photo, num.replicate = 100)
# photo <- ScoreJackStraw(photo, dims = 1:20)
# 
# JackStrawPlot(photo, dims = 1:20)
# Or just use elbow plot (var explained)
ElbowPlot(photo, ndims=100, reduction = "pca")

# Cluster cells
upDimLimit=15
photo <- FindNeighbors(photo, dims = 1:upDimLimit)
photo <- FindClusters(photo, resolution = 1)

DimPlot(photo, reduction = "pca", label = TRUE, pt.size = 1, label.size = 6) + eelTheme()

# UMAP using the same PCA dimensions for prettier visualization
photo <- RunUMAP(photo, dims = 1:upDimLimit)
DimPlot(photo, reduction = "umap", label = TRUE, pt.size = 1, label.size = 6) + eelTheme()

# or use tSNE (clustering can't separate some M cones from L cones)
photo <- RunTSNE(photo, dims = 1:upDimLimit)
DimPlot(photo, reduction = "tsne", label = TRUE, pt.size = 1, label.size = 6) + eelTheme()

new.cluster.ids <- c("mslPR","mslL","mslL","mslL","mslPR","mslS","mslUV","mslR","mslM")
names(new.cluster.ids) <- levels(photo)
photo <- RenameIdents(photo, new.cluster.ids)
Idents(photo) <- factor(x = Idents(photo), levels = c("mslR","mslUV","mslS","mslM","mslL","mslPR"))

DimPlot(photo, reduction = "tsne", label = TRUE, pt.size = 1, label.size = 6) + eelTheme()


DotPlot(photo, features = c('rho',"opn1sw1", "opn1sw2",'opn1mw1','opn1lw2','si:busm1-57f23.1'))
DotPlot(photo, features = c('rho',"opn1sw1", "opn1sw2",'opn1mw1','opn1lw2','si:busm1-57f23.1',"nr2e3","nr2f1b","nr2f6b",'syt5a','syt5b','gnat2','arr3a','arr3b',"foxq2")) + theme(axis.text.x = element_text(angle = 45, hjust=1))

FeaturePlot(photo, reduction = 'tsne', features = c("rho","opn1sw1", "opn1sw2",'opn1mw1','opn1lw2'))
FeaturePlot(photo, reduction = 'tsne', features = c("rho","opn1sw1", "opn1sw2","foxq2",'opn1mw1',"pcdh10a",'opn1lw2',"thrb","skor1a")) #markers to parse clusters apart
FeaturePlot(photo, reduction = 'tsne', features = c("rho","nrl","pde6a","nr2e3","crx",'gnat1',"saga","sagb","gucy2f","grk1a")) # rod markers
FeaturePlot(photo, reduction = 'tsne', features = c("rho","eno2","rom1a","lrrn1","unc119.2","pdca","cplx4c")) # rod genes from my own dataset
FeaturePlot(photo, reduction = 'tsne', features = c("gnat2",'arr3a','arr3b',"pde6c","pde6h","guca1d","grk7a","crx","neurod1","nr2f6b")) #cone markers
FeaturePlot(photo, reduction = 'tsne', features = c("nr2e3","nr2f6b","crx",'syt5a','syt5b','gnat2','arr3a','arr3b',"pde6g","pde6h", "neurod1")) #cone markers
FeaturePlot(photo, reduction = 'tsne', features = c("slc1a8b","rgs9a","slc25a24","ppa1a","sema7a","kera","nexn","dusp5","crhbp","prph2a")) #cone markers from my own dataset
FeaturePlot(photo, reduction = 'tsne', features = c("otx5","tbx2a","tbx2b","rxrga","rxrgb", "thrb", "sema7a","cnga3a","cnga3b"))

FeaturePlot(photo, reduction = 'tsne', features = c("nrl","mafa","mafb")) # rod markers

FeaturePlot(photo, reduction = 'tsne', features = c("tbx2a","tbx2b"))
FeaturePlot(photo, reduction = 'tsne', features = c("sema7a","efna1b","ntng2a","syt1a","syt5a","syt5b","tbx2a","tbx2b","eml1"))

FeaturePlot(photo, reduction = 'tsne', features = c("nr2e3","nr2f1b","nr2f6b",'syt5a','syt5b','gnat2','arr3a','arr3b',"foxq2")) #markers to parse clusters apart


DotPlot(photo, features = c("sema7a","efna1b","ntng2b","syt1a","syt1b","syt5a","syt5b","tbx2a","tbx2b","eml1"))
# plotting semaphorins
p = DotPlot(photo, features = c("sema3aa","sema3ab","sema3b","sema3bl","sema3c","sema3d","sema3e","sema3fa","sema3fb","sema3ga","sema3gb","sema3h","sema4aa","sema4ab","sema4ba","sema4bb","sema4c","sema4d","sema4e","sema4f","sema4ga","sema4gb","sema5a","sema5ba","sema5bb","sema6a","sema6ba","sema6bb","sema6d","sema6dl","sema6e","sema7a")) 
p + eelTheme() + theme(axis.text.x=element_text(angle=45, hjust=1))
# ephrins
p = DotPlot(photo, features = c("efna1a","efna1b","efna2a","efna2b","efna3a","efna3b","efna5a","efna5b","efnb1","efnb2a","efnb2b","efnb3a","efnb3b")) 
p + eelTheme() + theme(axis.text.x=element_text(angle=45, hjust=1))
# vamps
p = DotPlot(photo, features = c("vamp1","vamp2","vamp3","vamp4","vamp5","vamp8")) 
p + eelTheme() + theme(axis.text.x=element_text(angle=45, hjust=1))
# syts
p = DotPlot(photo, features = c("syt1a","syt1b","syt2a","syt3","syt4","syt5a","syt5b","syt6a","syt6b","syt7a","syt7b")) 
p + eelTheme() + theme(axis.text.x=element_text(angle=45, hjust=1))
# netrins
p = DotPlot(photo, features = c("ntn1a","ntn1b","ntn2","ntn4","ntn5","ntng1a","ntng2a","ntng2b")) 
p + eelTheme() + theme(axis.text.x=element_text(angle=45, hjust=1))
# nr2
p = DotPlot(photo, features = c("nr2c1","nr2c2","nr2e1","nr2e3","nr2f1a","nr2f1b","nr2f2","nr2f5","nr2f6a","nr2f6b")) 
p + eelTheme() + theme(axis.text.x=element_text(angle=45, hjust=1))
# tbx
p = DotPlot(photo, features = c("tbx1","tbx15","tbx16","tbx18","tbx20","tbx21","tbx22","tbx2a","tbx2b","tbx3a","tbx3b","tbx4","tbx5a","tbx5b","tbx6")) 
p + eelTheme() + theme(axis.text.x=element_text(angle=45, hjust=1))
# rod genes
p = DotPlot(photo, features = c("rho","gnat1","dhrs13l1","lingo1b","cabp4","saga","sagb","gnb1b","rgs9bp","eno2","guca1b","grk1a","rom1b","cplx4c","lrrn1","ncs1a","sgce","kcnv2a","pdca")) 
p + eelTheme() + theme(axis.text.x=element_text(angle=45, hjust=1))
# cone genes
p = DotPlot(photo, features = c("slc1a8b","rgs9a","slc25a24","drd4b","si:busm1-57f23.1","kera","nexn","clic1","gnat2","tgif1","crhbp","kcnv2b","anks1b","clul1")) 
p + eelTheme() + theme(axis.text.x=element_text(angle=45, hjust=1))
# UV+S enriched
p = DotPlot(photo, features = c("pcdh11.1","grk7b","skor1a","pcdh11","srgap3","ablim1a","nav2a","jam2a","chl1a","s100z","foxq2","sh3bp5b","kcnk1a","ntng2b","nxph1","pik3r3b","myl4","pacrg","fah")) 
p + eelTheme() + theme(axis.text.x=element_text(angle=45, hjust=1))
# M+L enriched
p = DotPlot(photo, features = c("slc32a1","apln","lactbl1b","nrtn","pcdh10a","myo7aa")) 
p + eelTheme() + theme(axis.text.x=element_text(angle=45, hjust=1))
# Mup/Ldown
p = DotPlot(photo, features = c("dok6","spock3","lrrfip1a","sema3fb","itga1","auts2a","esama","esamb","plxnb1a","cgnb","phf19","lrrc20")) 
p + eelTheme() + theme(axis.text.x=element_text(angle=45, hjust=1))
# Lup/Mdown
p = DotPlot(photo, features = c("slc32a1","s100v2","smad5","snap25a","ttyh2l","arhgap11a","ggctb","fbxo32","pcdh10a","abracl")) 
p + eelTheme() + theme(axis.text.x=element_text(angle=45, hjust=1))
# UVup/Sdown Genes
p = DotPlot(photo, features = c("myl4","ttyh2l","lhx1a","rx3","itgb1bp1")) 
p + eelTheme() + theme(axis.text.x=element_text(angle=45, hjust=1))
# Sup/UVdown Genes
p = DotPlot(photo, features = c("fkbp5","prss23","chkb","mpzl2b","nr1d1","camk2d1","frmpd2","foxo1a","prom1a")) 
p + eelTheme() + theme(axis.text.x=element_text(angle=45, hjust=1))
# TFs in DEGs
p = DotPlot(photo, features = c("tbx2a","tbx2b","otx5","tgif1","crhbp","six7","six6b","ndrg1b","egr1","nr1d1","sall1a","skor1a","foxq2","lhx1a","ntf3","rxrga","thrb","fgf1b","eya2","sox6","sox4b","pbx1a","hmgb3a","fbxo21","tfe3a","prdm1a")) 
p + eelTheme() + theme(axis.text.x=element_text(angle=45, hjust=1))


# heatmap of top genes
# # find markers for every cluster compared to all remaining cells, report only the positive ones
photo.markers <- FindAllMarkers(photo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
photo.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
top10 <- photo.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(photo, features = top10$gene) + NoLegend()


saveRDS(photo, file = "./photoreceptors_3dpf.rds")
## -----------------------------------------------------------------------------
# Exploring all celltypes
ps = DotPlot(pbmc, features = c("tbx2a","tbx2b","otx5","tgif1","crhbp","six7","six6b","ndrg1b","egr1","nr1d1","sall1a","skor1a","foxq2","lhx1a","ntf3","rxrgb","thrb","fgf1b","eya2","sox6","sox4b","pbx1a","hmgb3a","fbxo21","tfe3a")) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ps

ps = DotPlot(pbmc, features = c("nr2e3","crx","syt5b")) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ps

ps = DotPlot(pbmc, features = c("vamp1","vamp2")) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ps

ps = DotPlot(pbmc, features = c("tbx2a","tbx2b","foxq2","skor1a","nrtn","jam2a","pcdh10a","pcdh11","ntng2a","vamp2","syt5a","syt5b")) + eelTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ps
# # ------------------------------------------------------------------------------------------
# -------------------------------------------------------------------
# CAN BE RESTARTED HERE
photo = readRDS(file = "./photoreceptors_3dpf.rds")

# getting % and Mean expression for all genes into csv.
# DoHeatmap and DotPlot use the @scale.data slot for average expression display, which z-scored expression values (for example, as those used in PCA).
# Cells with a value > 0 represent cells with expression above the population mean (a value of 1 would represent cells with expression 1SD away from the population mean). Hope that helps!
# This is very detailed explanation: https://github.com/satijalab/seurat/issues/2798

# Average expression in non-log scale for larvae:
avgExp = AverageExpression(photo, slot="counts")
avgExp = avgExp$RNA
# avgExp = head(avgExp,10)
nI = length(levels(Idents(photo)))
colnames(avgExp) = paste("avg", levels(Idents(photo)), sep = "")
# Calculate mean expression across all ages
avgExp$avgMean = rowMeans(avgExp)

# Percent expression can be obtained through DotPlot, but not Average expression, as this is log z-cored normalized data # avgExpZ = dPlot$data$avg.exp
# Percent expression in larvae:
dPlot = DotPlot(photo, features=rownames(avgExp))
pctExp = dPlot$data$pct.exp
pctExp = as.data.frame(matrix(pctExp,nrow = length(pctExp)/nI, ncol=nI));

rownames(pctExp) = rownames(avgExp)
colnames(pctExp) = paste("pct", levels(Idents(photo)), sep = "")

# Calculate mean percent expression across all ages
pctExp$pctMean = rowMeans(pctExp)
head(pctExp,10)

conesXu = cbind(pctExp,avgExp)
head(conesXu,10)


# these are duplicated genes with lower vs upper case
dups=toupper(c('brd8','cct2','dtx4','eif1b','flrt2','frmd5','kcnb2','lamp1','lyrm4','mepce','mras','ncam2','pcdh8','pcp4l1','rps17','serp1','slc25a22','tatdn3','tenm3','tmem151a','tppp','trappc9','ubb','ublcp1','vps72','znf423'))
for(i in 1:length(dups)) {
   rownames(conesXu)[rownames(conesXu) == dups[i]] = paste(dups[i],"_ii",sep="")
}
# now everything can be lower cased
rownames(conesXu) = tolower(rownames(conesXu))
# also replacing (1 of many) with .1
rownames(conesXu) = gsub("(1 of many)",".1",rownames(conesXu))


conesXu = conesXu[c("avgMean","avgmslR","avgmslUV","avgmslS","avgmslM","avgmslL","avgmslPR",
                    "pctMean","pctmslR","pctmslUV","pctmslS","pctmslM","pctmslL","pctmslPR")]
# sort by avgExpression
conesXu=conesXu[order(-conesXu$avgMean),]
head(conesXu,100)

scangenes=c('syt5a','syt5b','nr2e3','tbx2a','tbx2b','arr3a','arr3b','opn1lw1','opn1lw2');
conesXu[scangenes,]

# this is as far as I can take this
# take aways: 
# ribosomal proteins are the most expressed genes
# separation worked quite well and there is consistency between datasets it seems

write.csv(conesXu,"./conesXu.csv",quote=FALSE)

# -------------------------------------------------------------------
# -------------------------------------------------------------------



# # ------------------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------------------
# # Subclustered photoreceptors (Vincent's dataset)
# # This has 473 cells, so my guess is that he's filtered out anything that doesn't have opsin counts
# photo2 <- readRDS(file = "./VPK_ZF-Cones4-Subset.rds")
# 
# DimPlot(photo2, reduction = "umap")
# 
# FeaturePlot(photo2, features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw3','opn1lw2','si:busm1-57f23.1','rho','saga'))


# check counts in each cluster
VlnPlot(photo, features = c("gnat2", "rho",'elovl4b','gnat1','si:busm1-57f23.1')) #log norm values

VlnPlot(photo, features = c("sema7a","efna1a","efna1b","ntng2b","syt1a","syt1b","syt5a","syt5b","tbx2a","tbx2b","eml1","ntng2a"))

DotPlot(photo, features = c("sema7a","efna1b","ntng2b","syt1a","syt1b","syt5a","syt5b","tbx2a","tbx2b","eml1"))

# find markers for every cluster compared to all remaining cells, report only the positive ones
photo.markers <- FindAllMarkers(photo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
photo.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
# heatmap of top genes
top20 <- photo.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
DoHeatmap(photo, features = top20$gene) + NoLegend()

FeaturePlot(photo, features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw3','opn1lw2','si:busm1-57f23.1','rho','saga'))

# Manually checking adult DEGs
#enriched in UV and S
# "skor1a" "grk7b" "pcdh11" "arr3b" "nav2a" "jam2a" "pip5k1bb"
VlnPlot(photo, features = c("skor1a", "grk7b", "pcdh11", "arr3b", "nav2a", "jam2a", "pip5k1bb"))
# "ablim1a" #UV and S + M
# "abcg4b" #UV, S and M
# "il34" #UV, S and M
VlnPlot(photo, features = c("ablim1a", "abcg4b", "il34"))
# enriched in M and L
VlnPlot(photo, features = c("wdr47b", "nrtn", "hs6st1a", "phf19", "selenom",'fox6'))
# enriched in L
VlnPlot(photo, features = c("ntng2a", "slc32a1", "lactbl1b", "arhgap11a", "thrb", "si:busm1-57f23.1", "s100v2", "rxrga", "eya2", "c1galt1a"))
# down in L
VlnPlot(photo, features = c("ntf3", "dok6", "hs3st3l", "lrrfip1a"))
# enriched in M
VlnPlot(photo, features = c("sema3fb", "cgnb"))
#neriched in UV (ttyh2l in UV and L)
VlnPlot(photo, features = c("tbx2a", "myl4", "ttyh2l", "galnt12"))
# enriched in S
VlnPlot(photo, features = c("foxq2", "prss23", "mpzl2b", "prss23", "tefb", "pacrg", "slc41a2b"))
# enriched in rods
VlnPlot(photo, features = c("alpl","cabp4","arhgap20","dhrs13l1","lingo1b","comtb","mapkapk3","nr2e3"))
# enriched in cones
VlnPlot(photo, features = c("stxbp4","si:ch211-22d5.2","slc1a8b","rgs9a","gnat2"))
# phototransduction

# transducin
VlnPlot(photo, features = c("gnat1","gnat2"))

VlnPlot(photo, features = c('gngt1','gngt2a','gngt2b'))

# plotName = 'gnb'
VlnPlot(photo, features = c('gnb1a','gnb1b','gnb3a','gnb3b','gnb4b','gnb5a','gnb5b'))

# plotName = 'pde'
VlnPlot(photo, features = c('pde6a','pde6b','pde6c','pde6d','pde6ga','pde6gb','pde6c','pde6ha'))
# genelist = []

# plotName = 'GC'
VlnPlot(photo, features = c('gucy1a1','gucy1b1','gucy1b2','gucy2c','gucy2d','gucy2f','gucy2g'))

# plotName = 'GC_retinaSpecific'
VlnPlot(photo, features = c('gucy2d','gucy2f'))

# plotName = 'GCAP'
VlnPlot(photo, features = c('guca1a','guca1b','guca1c','guca1d','guca1e','guca1g'))

# plotName = 'CNGalpha'
VlnPlot(photo, features = c('cnga1a','cnga1b','cnga3a','cnga3b','cnga4','cnga2a','cnga2b'))

# plotName = 'CNGalpha_retinaSpecific'
VlnPlot(photo, features = c('cnga1a','cnga1b','cnga3a','cnga3b'))

# plotName = 'CNGbeta'
VlnPlot(photo, features = c('cngb1a','cngb3.1','cngb3.2'))

# plotName = 'RGS'
VlnPlot(photo, features = c('rgs9a','rgs9b','rgs9bp','rgs3a','rgs6','rgs11','rgs12a','rgs16','rgs20'))

# plotName = 'GRK'
VlnPlot(photo, features = c('grk1a','grk1b','grk3','grk4','grk5','grk5l', 'grk6','grk7a','grk7b'))

# genelist = []
# plotName = 'rcvrn'
VlnPlot(photo, features = c('rcvrna','rcvrnb','rcvrn2','rcvrn3','ncalda','ncaldb'))

# plotName = 'arrestins'
VlnPlot(photo, features = c('saga','sagb','arr3a','arr3b'))

#special
VlnPlot(photo, features = c("pcdh10a", "slc1a8b","frmpd2","prdm1b","nme2a","lrrn1","nr2f6b"))
VlnPlot(photo, features = c("elovl8a", "elovl4b","syt5a","syt5b"))
VlnPlot(photo, features = c("ankrd33aa", "ankrd33ba","foxg1a","foxg1b","lrit3a"))

VlnPlot(photo, features = c("ntng2b", "ntng2a","nlgn1","nlgn2a","nlgn2b","nlgn3a","nlgn3b","nlgn4xa","nlgn4xb"))

# Expression in all retinal cells
VlnPlot(pbmc, features = c("ntng2b", "ntng2a","nlgn1","nlgn2a","nlgn2b","nlgn3a","nlgn3b","nlgn4xa","nlgn4xb"))

VlnPlot(pbmc, features = c("syt5a","syt5b"))

VlnPlot(pbmc, features = c("sema7a","efna1a","efna1b","ntng2b","syt1a","syt1b","syt5a","syt5b","tbx2a","tbx2b","eml1","ntng2a"))
VlnPlot(pbmc, features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw3','opn1lw2','si:busm1-57f23.1','rho','saga'))

# complete agreement with manual picking 
VlnPlot(photo, features = c("cadm1a", "cadm1b",'cadm2a','cadm2b','cadm3','cadm4'))
# but all classes expresses cadm3
VlnPlot(pbmc, features = c("cadm1a", "cadm1b",'cadm2a','cadm2b','cadm3','cadm4'))


VlnPlot(photo, features = c("prdm1a","prdm1b"))





