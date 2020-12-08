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
setwd("/Users/angueyraaristjm/Documents/LiLab/RNAseq/zfRetDev10x_Xu2020/ret72hpf/")
directory <- "/Users/angueyraaristjm/Documents/LiLab/RNAseq/zfRetDev10x_Xu2020/ret72hpf/"
getwd()

# -------------------------------------------------------------------
# Load the 10x dataset (2826 cells)
pbmcA.data <- Read10X(data.dir = "/Users/angueyraaristjm/Documents/LiLab/RNAseq/zfRetDev10x_Xu2020/ret72hpf/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmcA <- CreateSeuratObject(counts = pbmcA.data, project = "zf10X_72hpfA", min.cells = 3, min.features = 200)
pbmcA

# Load the replciate 10x dataset (405 cells)
pbmcB.data <- Read10X(data.dir = "/Users/angueyraaristjm/Documents/LiLab/RNAseq/zfRetDev10x_Xu2020/ret72hpfB/")
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
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000) # this is to get natural log-transformed data using log1p
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
DimPlot(pbmc, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()
saveRDS(pbmc, file = "./allCells_reClsutered.rds")
# ------------------------------------------------------------------------------------------

# Clear all plots -------------------------------------------------------------------
try(dev.off(),silent=TRUE)
# Clear environment -------------------------------------------------------------------
rm(list=ls())

# ------------------------------------------------------------------------------------------
# Setup -------------------------------------------------------------------
setwd("/Users/angueyraaristjm/Documents/LiLab/RNAseq/zfRetDev10x_Xu2020/ret72hpf/")
directory <- "/Users/angueyraaristjm/Documents/LiLab/RNAseq/zfRetDev10x_Xu2020/ret72hpf/"
getwd()

# ------------------------------------------------------------------------------------------
pbmc <- readRDS(file = "./allCells_reClsutered.rds")
# Subclustering photoreceptors
photo <- subset(pbmc, idents = "Photo")

# Identification of highly variable features (feature selection)
photo <- FindVariableFeatures(photo, selection.method = "vst", nfeatures = 200)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(photo), 20)
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
ElbowPlot(photo)

# Cluster cells
photo <- FindNeighbors(photo, dims = 1:3)
photo <- FindClusters(photo, resolution = 0.2)

# UMAP using the same PCA dimensions for prettier visualization
photo <- RunUMAP(photo, dims = 1:3)
DimPlot(photo, reduction = "umap", label=TRUE)
# saveRDS(photo, file = "../photoreceptors_umap.rds")

# or use tSNE (clustering can't separate some M cones from L cones)
photo <- RunTSNE(photo, dims = 1:7)
DimPlot(photo, reduction = "tsne", label=TRUE)
# saveRDS(photo, file = "../photoreceptors_tSNE.rds")

# find markers for every cluster compared to all remaining cells, report only the positive ones
photo.markers <- FindAllMarkers(photo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
photo.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)

# dataset A only:
# 0 = L-cones too or junk?
# 1 = UV-cones + S-cones
# 2 = L-cones
# 3 = M-cones
# 4 = L-cones
# 5 = rods

# combined datasets:
# 0 = not sure at all
# 1 = L-cones
# 2 = L-cones and M-cones
# 3 = UV-cones
# 4 = rods
# 5 = rods


FeaturePlot(photo, features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw3','opn1lw2','si:busm1-57f23.1','rho','saga'))
VlnPlot(photo, features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw3','opn1lw2','si:busm1-57f23.1','rho','saga')) #raw counts
# VlnPlot(photo, features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw3','opn1lw2','si:busm1-57f23.1','rho','saga'), slot = "counts", log = TRUE) #raw counts

# check counts in each cluster
VlnPlot(photo, features = c("gnat2", "rho",'elovl4b','gnat1','si:busm1-57f23.1')) #log norm values
VlnPlot(photo, features = c("gnat2", "rho",'elovl4b','gnat1','si:busm1-57f23.1'), slot = "counts", log = TRUE) #raw counts



RidgePlot(photo, features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1','rho'))
DotPlot(photo, features = c("gnat2", "rho",'elovl4b','gnat1','si:busm1-57f23.1'))

DotPlot(photo, features = c("sema7a","efna1b","ntng2b","syt1a","syt1b","syt5a","syt5b","tbx2a","tbx2b","eml1"))
# plotting semaphorins
p = DotPlot(photo, features = c("sema3aa","sema3ab","sema3b","sema3bl","sema3c","sema3d","sema3e","sema3fa","sema3fb","sema3ga","sema3gb","sema3h","sema4aa","sema4ab","sema4ba","sema4bb","sema4c","sema4d","sema4e","sema4f","sema4ga","sema4gb","sema5a","sema5ba","sema5bb","sema6a","sema6ba","sema6bb","sema6d","sema6dl","sema6e","sema7a")) 
p + theme(axis.text.x=element_text(angle=45, hjust=1))
# ephrins
p = DotPlot(photo, features = c("efna1a","efna1b","efna2a","efna2b","efna3a","efna3b","efna5a","efna5b","efnb1","efnb2a","efnb2b","efnb3a","efnb3b")) 
p + theme(axis.text.x=element_text(angle=45, hjust=1))

# check counts in tSNE space
# FeaturePlot(photo, features = c("gnat2",'gngt2b', "gnat1",'rho','si:busm1-57f23.1','opn1lw2', 'crx'))

# heatmap of top genes
top10 <- photo.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(photo, features = top10$gene) + NoLegend()

FeaturePlot(photo, features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw3','opn1lw2','si:busm1-57f23.1','rho','saga'))

# Clear all plots -------------------------------------------------------------------
try(dev.off(),silent=TRUE)
# Clear environment -------------------------------------------------------------------
rm(list=ls())

# ------------------------------------------------------------------------------------------
# Setup -------------------------------------------------------------------
setwd("/Users/angueyraaristjm/Documents/LiLab/RNAseq/zfRetDev10x_Xu2020/ret72hpf/")
directory <- "/Users/angueyraaristjm/Documents/LiLab/RNAseq/zfRetDev10x_Xu2020/ret72hpf/"
getwd()
# ------------------------------------------------------------------------------------------
# Subclustered photoreceptors (Vincent's dataset)
photo <- readRDS(file = "./VPK_ZF-Cones4-Subset.rds")

DimPlot(photo, reduction = "umap")

FeaturePlot(photo, features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw3','opn1lw2','si:busm1-57f23.1','rho','saga'))
VlnPlot(photo, features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw3','opn1lw2','si:busm1-57f23.1','rho','saga'))
VlnPlot(photo, features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw3','opn1lw2','si:busm1-57f23.1','rho','saga'), slot = "counts", log = TRUE) #raw counts

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




# Exporting photoreceptor count matrix to make comparisons with manual data
table(Idents(photo))
rods <- subset(photo, idents = "Rods")
GetAssayData(object = rods, slot = "counts")
