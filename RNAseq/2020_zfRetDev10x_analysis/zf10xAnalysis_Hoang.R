# Enter commands in R (or R studio, if installed)
# install.packages('Seurat')
# install.packages("reticulate")
# install.packages("umap")
# install.packages("openssl")
# install.packages('BiocManager')
# BiocManager::install('limma')
# BiocManager::install('DESeq2')

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
getwd()
# -------------------------------------------------------------------
# Load the 10x dataset (28845 cells), after updating to Seurat_v3 (had to use biowulf)
pbmc = readRDS("~/Documents/LiMolec/otherRNAseq/zfRet_HoangBlackshaw2020/zfDev_pbmc_v3.rds");
# Load the adult cone dataset (~2000 cells), after reclustering (separate R script)
photo = readRDS("~/Documents/LiMolec/otherRNAseq/zfRet_HoangBlackshaw2020/cones_Adult.rds");
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
# -------------------------------------------------------------------
# check counts in each cluster
VlnPlot(pbmc, features = c("gnat2", "rho",'elovl4b','gnat1','si:busm1-57f23.1')) #log norm values
VlnPlot(pbmc, features = c("gnat2", "rho",'elovl4b','gnat1','si:busm1-57f23.1'), slot = "counts", log = TRUE) #raw counts
RidgePlot(pbmc, features = c("gnat2", "rho",'elovl4b','gnat1','si:busm1-57f23.1'))
DotPlot(pbmc, features = c("gnat2", "rho",'elovl4b','gnat1','si:busm1-57f23.1'))

VlnPlot(pbmc, features = c("sema7a","efna1b","ntng2b","syt1a","syt1b","syt5a","syt5b","tbx2a","tbx2b","eml1")) #log norm values
DotPlot(pbmc, features = c("sema7a","efna1b","ntng2b","syt1a","syt1b","syt5a","syt5b","tbx2a","tbx2b","eml1"))

VlnPlot(pbmc, features = c("opn1sw1", "opn1sw2",'opn1mw1','gnat2','opn1lw2','si:busm1-57f23.1','rho','saga')) #log norm values
# check counts in tSNE space
# FeaturePlot(pbmc, features = c("gnat2",'gngt2b', "gnat1",'rho','si:busm1-57f23.1','opn1lw2', 'crx'))

# heatmap of top genes: THIS IS NOT HELPFUL IN HEATMAP.
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
# DoHeatmap(pbmc, features = top10$gene) + NoLegend() # THIS IS NOT HELPFUL IN HEATMAP.

#Using https://proteinpaint.stjude.org/?genome=danRer10&singlecell=files/danRer10/NEI.AGI.retina/singlecell/dev/view.json to map
# Retinal progenitor: 0
# Photoreceptor precursor: 3
FeaturePlot(pbmc, features = c("hmgb2b","stmn1a",'hmgn2','rrm2.1'))

# Cones: adults = 4,  larval = 12
# Rods: 2
FeaturePlot(pbmc, features = c("tulp1a","elovl4b",'ckmt2a','gngt2a'))
FeaturePlot(pbmc, features = c("opn1sw1", "opn1sw2",'opn1mw1','gnat2','opn1lw2','si:busm1-57f23.1','rho','saga'))



# HCs: 6
FeaturePlot(pbmc, features = c("rem1","CR361564.1",'prkar2ab','rprmb','opn4.1','atp1b1a'))

# BCs: 5, 15
FeaturePlot(pbmc, features = c("samsn1a","vsx1",'gnb3a','pcp4l1','neurod4','gnao1b'))

# Immature MG: 16
# MG: 11, 9 , 1
FeaturePlot(pbmc, features = c("cavin2a","rhbg",'slc1a2b','atp1a1b'))

# RGCs = 13 + 10
FeaturePlot(pbmc, features = c("islr2","rbpms2a",'pou4f1'))

# ngng ACs = 7
# gly ACs = 14
# GABA ACs = 8
FeaturePlot(pbmc, features = c("slc32a1","pax10",'tfap2a','gad2','slc6a1b'))
# Clear all plots -------------------------------------------------------------------
try(dev.off(),silent=TRUE)


# ------------------------------------------------------------------------------------------

#some genes that overlap in this dataset and ours
AverageExpression(photo, features=c("tbx2a", "tbx2b","sema7a","mpzl2b","foxq2","tefb","cnga3b","cnga3a","ccdc85a","ntf3","lrrfip1a","sema3fb","cgnb","s100v2","rxrga","eya2","mstnb","pbx1a","cxxc4","thrb","smad5","snap25a","snap25b","si:ch211-160o17.6","fbxo32","pcdh10a","myo7aa"))

# Load DESeq2 results
# L-cones
DESeq2_L = read.csv("./DESeq2_L.csv", header = TRUE, row.names=1, sep = ",", quote = "\"",dec = ".", fill = TRUE, comment.char = "")
# M-cones
DESeq2_M = read.csv("./DESeq2_M.csv", header = TRUE, row.names=1, sep = ",", quote = "\"",dec = ".", fill = TRUE, comment.char = "")
# S-cones
DESeq2_S = read.csv("./DESeq2_S.csv", header = TRUE, row.names=1, sep = ",", quote = "\"",dec = ".", fill = TRUE, comment.char = "")
# UV-cones
DESeq2_UV = read.csv("./DESeq2_UV.csv", header = TRUE, row.names=1, sep = ",", quote = "\"",dec = ".", fill = TRUE, comment.char = "")
#UV&S vs L&M
DESeq2_USvsLM = read.csv("./DESeq2_USvsLM.csv", header = TRUE, row.names=1, sep = ",", quote = "\"",dec = ".", fill = TRUE, comment.char = "")
#L vs M
DESeq2_LvsM = read.csv("./DESeq2_LvsM.csv", header = TRUE, row.names=1, sep = ",", quote = "\"",dec = ".", fill = TRUE, comment.char = "")
#UV vs S
DESeq2_UvsS = read.csv("./DESeq2_UvsS.csv", header = TRUE, row.names=1, sep = ",", quote = "\"",dec = ".", fill = TRUE, comment.char = "")

#run PCA
photo = ScaleData(photo)

# Identification of highly variable features (feature selection)
photo <- FindVariableFeatures(photo, selection.method = "vst", nfeatures = 150)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(photo), 100)
top10


photo <- RunPCA(photo, features = VariableFeatures(object = photo))
# # Examine and visualize PCA results a few different ways
# print(photo[["pca"]], dims = 1:5, nfeatures = 5)
# VizDimLoadings(photo, dims = 1:2, reduction = "pca")
# DimPlot(photo, reduction = "pca")
# DimHeatmap(photo, dims = 1, cells = 500, balanced = TRUE)
# DimHeatmap(photo, dims = 1:15, cells = 500, balanced = TRUE)


FeaturePlot(photo, reduction = 'tsne', features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))
VlnPlot(photo, features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))


DotPlot(photo, features = c("foxq2", "ntf3", "rx1","rxrga", "tbx2a", "tbx2b")) #log norm values


FeatureScatter(photo, feature1 = "opn1sw2", feature2 = "foxq2")
DoHeatmap(photo, top10)

# find markers for every cluster compared to all remaining cells, report only the positive ones
photo.markers <- FindAllMarkers(photo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
photo.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)


# check counts in each cluster
FeaturePlot(photo, reduction = 'tsne', features = c('tbx2a','opn1mw1','opn1mw2','opn1mw3','opn1mw4'))

FeaturePlot(photo, reduction = 'tsne', features = c('tbx2a','tbx2b','syt5a','syt5b','arr3a','arr3b'))
VlnPlot(photo, features = c("gnat2", "rho",'elovl4b','gnat1','si:busm1-57f23.1')) #log norm values
# VlnPlot(photo, features = c("gnat2", "rho",'elovl4b','gnat1','si:busm1-57f23.1'), slot = "counts", log = TRUE) #raw counts

RidgePlot(photo, features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1','rho'))
DotPlot(photo, features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1','rho'))
DotPlot(photo, features = c("gnat2", "rho",'elovl4b','gnat1','si:busm1-57f23.1'))

DotPlot(photo, features = c("syt5a","syt5b","tbx2a","tbx2b",'si:busm1-57f23.1'))

DotPlot(photo, features = c("sema7a","efna1b","ntng2b","syt1a","syt1b","syt5a","syt5b","tbx2a","tbx2b","eml1"))
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
setwd("/Users/angueyraaristjm/Documents/LiMolec/otherRNAseqzfRetDev10x_Xu2020/ret72hpf/")
directory <- "/Users/angueyraaristjm/Documents/LiMolec/otherRNAseqzfRetDev10x_Xu2020/ret72hpf/"
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
