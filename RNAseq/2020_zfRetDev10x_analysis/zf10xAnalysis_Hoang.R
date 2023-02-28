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
setwd("/Users/angueyraaristjm/Documents/eelMolec/zfRNAseq/2020_Hoang_zfRet10x/")
directory = "/Users/angueyraaristjm/Documents/eelMolec/zfRNAseq/2020_Hoang_zfRet10x/"
getwd()
# Plot themes -------------------------------------------------------------------
eelTheme = function (base_size = 42, base_family = "Avenir") {
   theme_classic() %+replace%
      theme(
         axis.line = element_line(colour = 'black', size = 1),
         axis.text = element_text(size=28, family = "Avenir"),
         text = element_text(size=42, family = "Avenir"),
      )
}
# -------------------------------------------------------------------
# Load the 10x dataset (28845 cells), after updating to Seurat_v3 (had to use biowulf)
pbmc = readRDS("~/Documents/eelMolec/zfRNAseq/2020_Hoang_zfRet10x/zfDev_pbmc_v3.rds");
# Load the adult cone dataset (~2000 cells), after reclustering (separate R script)
# photo = readRDS("~/Documents/LiMolec/otherRNAseq/zfRet_HoangBlackshaw2020/cones_Adult.rds");
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
new.cluster.ids <- c("RPC","MG1","Rods","PRPC","Cones_adult","BC1","HC","AC","ACgaba","MG2","RGC1","MG3","Cones_larval","RGC2","ACgly","BC2","MGi")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
Idents(pbmc) <- factor(x = Idents(pbmc), levels = c("RPC","PRPC","Cones_larval","Cones_adult","Rods","HC","BC1","BC2","AC","ACgaba","ACgly","RGC1","RGC2","MGi","MG1","MG2","MG3"))
# retcolors = c("#4d493c","#dfdac8","#cca819","#ffd024","#7d7d7d","#3db500","#871308","#871308","#e7851d","#e7851d","#8c561b","#2a81de","#2a81de","#d277ff","#ff13c8","#ff13c8","#ff13c8")

# UMAP plot
ps = DimPlot(pbmc, reduction = "umap", label=TRUE, repel = TRUE, pt.size = 1, label.size = 10) + #, cols=retcolors) +
   eelTheme() + NoLegend() + xlim(-12,12) + ylim(-10,10)
ps
ggsave(ps, file="00_UMAPClusters.png", path="./zfConeRNAseqFigure/UMAP/", width = 140*2, height = 105*2, units = "mm")
# tSNE plot
ps = DimPlot(pbmc, reduction = "tsne", label=TRUE, repel = TRUE, pt.size = 1, label.size = 10) +
   eelTheme() + NoLegend() + xlim(-50,50) + ylim(-60,60)
ps
ggsave(ps, file="00_TSNEClusters_labels.png", path="./zfConeRNAseqFigure/TSNE/", width = 140*2, height = 105*2, units = "mm")

# -----------------------------------------------------------------------------
max(pbmc$nFeature_RNA) # number of genes in each cell
max(pbmc$nCount_RNA) # number of UMIs in each cell
min(pbmc$nCount_RNA)
mean(pbmc$nCount_RNA)
median(pbmc$nCount_RNA)

metadata = pbmc@meta.data
# distribution of unique genes per cell
metadata %>%
   ggplot(aes(x=nFeature_RNA)) +
   geom_density(alpha = 0.2) +
   scale_x_log10() +
   theme_classic() +
   ylab("Cell density") +
   geom_vline(xintercept = 500)

# distribution of number of genes per cell
metadata %>%
   ggplot(aes(x=nCount_RNA)) +
   geom_density(alpha = 0.2) +
   scale_x_log10() +
   theme_classic() +
   ylab("Cell density") +
   geom_vline(xintercept = 500)

write.csv(x=pbmc$nFeature_RNA, file="Hoang2020_allCells_nUniqueGenes.csv", row.names = FALSE)
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# getting % and Mean expression for all genes into csv.
# DoHeatmap and DotPlot use the @scale.data slot for average expression display, which z-scored expression values (for example, as those used in PCA).
# Cells with a value > 0 represent cells with expression above the population mean (a value of 1 would represent cells with expression 1SD away from the population mean). Hope that helps!
# This is very detailed explanation: https://github.com/satijalab/seurat/issues/2798

# Average expression in non-log scale for larvae:
avgExp = AverageExpression(pbmc, slot="counts")
avgExp = avgExp$RNA
# avgExp = head(avgExp,10)
nI = length(levels(Idents(pbmc)))
colnames(avgExp) = paste("avg", levels(Idents(pbmc)), sep = "")
avgExp = as.data.frame(avgExp)
# Calculate mean expression across all ages
avgExp$avgMean = rowMeans(avgExp)
avgExp = avgExp[,c(ncol(avgExp),1:(ncol(avgExp)-1))]


# Percent expression can be obtained through DotPlot, but not Average expression, as this is log z-cored normalized data # avgExpZ = dPlot$data$avg.exp
# dPlot = DotPlot(pbmc, features=rownames(avgExp)) # not working; too many warnings
# dPlot <- FetchData(pbmc, rownames(avgExp)) # dot plot is not working. Trying this instead; does not work either
# ------------------------------------------------------------------------------------------------
# breaking down into pieces to try to find row that creates error
length(rownames(avgExp))
# let's break it in half: nope!; a quarter? works
dPlot = DotPlot(pbmc, features=rownames(avgExp[c(1:5000),]))
# second quarter: fails; keep splitting.
dPlot21 = DotPlot(pbmc, features=rownames(avgExp[c(5001:6000),])) #pass
dPlot22 = DotPlot(pbmc, features=rownames(avgExp[c(6001:7000),])) #pass
dPlot23 = DotPlot(pbmc, features=rownames(avgExp[c(7001:8000),])) #pass
dPlot24 = DotPlot(pbmc, features=rownames(avgExp[c(8001:9000),])) #pass
dPlot25 = DotPlot(pbmc, features=rownames(avgExp[c(9001:10000),])) #pass
# second half? works
dPlot3 = DotPlot(pbmc, features=rownames(avgExp[c(10001:15000),]))
dPlot4 = DotPlot(pbmc, features=rownames(avgExp[c(15001:20000),]))
dPlot5 = DotPlot(pbmc, features=rownames(avgExp[c(20001:24711),]))
# rebind everything
# ------------------------------------------------------------------------------------------------
# testing how to give rownames appropriately
nLevels = length(levels(Idents(pbmc)))

temp = dPlot
tempnames = rownames(temp$data[seq(1,dim(temp$data)[1]/nLevels,1),])
temppct = temp$data$pct.exp
tempout1 = as.data.frame(matrix(temppct, nrow = length(temppct)/nLevels, ncol = nLevels, dimnames = list(c(tempnames),c(levels(Idents(pbmc))))))
temp = dPlot21
tempnames = rownames(temp$data[seq(1,dim(temp$data)[1]/nLevels,1),])
temppct = temp$data$pct.exp
tempout2 = as.data.frame(matrix(temppct, nrow = length(temppct)/nLevels, ncol = nLevels, dimnames = list(c(tempnames),c(levels(Idents(pbmc))))))
temp = dPlot22
tempnames = rownames(temp$data[seq(1,dim(temp$data)[1]/nLevels,1),])
temppct = temp$data$pct.exp
tempout3 = as.data.frame(matrix(temppct, nrow = length(temppct)/nLevels, ncol = nLevels, dimnames = list(c(tempnames),c(levels(Idents(pbmc))))))
temp = dPlot23
tempnames = rownames(temp$data[seq(1,dim(temp$data)[1]/nLevels,1),])
temppct = temp$data$pct.exp
tempout4 = as.data.frame(matrix(temppct, nrow = length(temppct)/nLevels, ncol = nLevels, dimnames = list(c(tempnames),c(levels(Idents(pbmc))))))
temp = dPlot24
tempnames = rownames(temp$data[seq(1,dim(temp$data)[1]/nLevels,1),])
temppct = temp$data$pct.exp
tempout5 = as.data.frame(matrix(temppct, nrow = length(temppct)/nLevels, ncol = nLevels, dimnames = list(c(tempnames),c(levels(Idents(pbmc))))))
temp = dPlot25
tempnames = rownames(temp$data[seq(1,dim(temp$data)[1]/nLevels,1),])
temppct = temp$data$pct.exp
tempout6 = as.data.frame(matrix(temppct, nrow = length(temppct)/nLevels, ncol = nLevels, dimnames = list(c(tempnames),c(levels(Idents(pbmc))))))
temp = dPlot3
tempnames = rownames(temp$data[seq(1,dim(temp$data)[1]/nLevels,1),])
temppct = temp$data$pct.exp
tempout7 = as.data.frame(matrix(temppct, nrow = length(temppct)/nLevels, ncol = nLevels, dimnames = list(c(tempnames),c(levels(Idents(pbmc))))))
temp = dPlot4
tempnames = rownames(temp$data[seq(1,dim(temp$data)[1]/nLevels,1),])
temppct = temp$data$pct.exp
tempout8 = as.data.frame(matrix(temppct, nrow = length(temppct)/nLevels, ncol = nLevels, dimnames = list(c(tempnames),c(levels(Idents(pbmc))))))
temp = dPlot5
tempnames = rownames(temp$data[seq(1,dim(temp$data)[1]/nLevels,1),])
temppct = temp$data$pct.exp
tempout9 = as.data.frame(matrix(temppct, nrow = length(temppct)/nLevels, ncol = nLevels, dimnames = list(c(tempnames),c(levels(Idents(pbmc))))))

# ------------------------------------------------------------------------------------------------
tempout = rbind(tempout1,tempout2)
tempout = rbind(tempout,tempout3)
tempout = rbind(tempout,tempout4)
tempout = rbind(tempout,tempout5)
tempout = rbind(tempout,tempout6)
tempout = rbind(tempout,tempout7)
tempout = rbind(tempout,tempout8)
tempout = rbind(tempout,tempout9)
# fixing names replaced by DotPlot
rownames(tempout) = gsub(pattern = 'rna_', replacement = '',x = rownames(tempout))
pctExp = tempout
colnames(pctExp) = paste("pct", colnames(pctExp), sep = "")
pctExp$pctMean = rowMeans(pctExp)
pctExp = pctExp[,c(ncol(pctExp),1:(ncol(pctExp)-1))]
# merging data frames
pbmcData = merge(avgExp,pctExp, by="row.names", all=TRUE)
rownames(pbmcData) = pbmcData[,1]
pbmcData = pbmcData[-c(1)]
# ------------------------------------------------------------------------------------------------
# quick check
scangenes=c('rho','nr2e3','opn1lw1','opn1lw2','thrb','skor1a','tbx2a','tbx2b','arr3a','arr3b');
pbmcData[scangenes,]


write.csv(pbmcData,"./retCells_Hoang.csv",quote=FALSE)
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
# -------------------------------------------------------------------
# check counts in each cluster
VlnPlot(pbmc, features = c("vamp1","vamp2","vamp3","vamp4","vamp5","vamp8")) #log norm values
VlnPlot(pbmc, features = c("gnat2", "rho",'elovl4b','gnat1','si:busm1-57f23.1')) #log norm values
VlnPlot(pbmc, features = c("gnat2", "rho",'elovl4b','gnat1','si:busm1-57f23.1'), slot = "counts", log = TRUE) #raw counts
RidgePlot(pbmc, features = c("gnat2", "rho",'elovl4b','gnat1','si:busm1-57f23.1'))
DotPlot(pbmc, features = c("gnat2", "rho",'elovl4b','gnat1','si:busm1-57f23.1'))

VlnPlot(pbmc, features = c("sema7a","efna1b","ntng2b","syt1a","syt1b","syt5a","syt5b","tbx2a","tbx2b","eml1")) #log norm values
DotPlot(pbmc, features = c("sema7a","efna1b","ntng2b","syt1a","syt1b","syt5a","syt5b","tbx2a","tbx2b","eml1"))
DotPlot(pbmc, features = c("vamp1","vamp2","vamp3","vamp4","vamp5","vamp8"))
FeaturePlot(pbmc, reduction = 'tsne', features = c("vamp1","vamp2","vamp3","vamp4","vamp5","vamp8"))
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
FeaturePlot(pbmc,reduction = 'tsne', features = c("hmgb2b","stmn1a",'hmgn2','rrm2.1'))

# Cones: adults = 4,  larval = 12
# Rods: 2
FeaturePlot(pbmc,reduction = 'tsne', features = c("tulp1a","elovl4b",'ckmt2a','gngt2a'))
FeaturePlot(pbmc,reduction = 'tsne', features = c("opn1sw1", "opn1sw2",'opn1mw1','gnat2','opn1lw2','si:busm1-57f23.1','rho','saga'))
FeaturePlot(pbmc,reduction = 'tsne', features = c("crx","otx5","nr2e3",'rx1','meis1b'),pt.size=1, order=TRUE, combine=TRUE, ncol = 3)
# --------------------------------------------------------------------------------
# Photoreceptor precursor
genelist = c("hmgb2b","stmn1a",'hmgn2','rrm2.1')
for(i in 1:length(genelist)) {
   genename = genelist[i]
   ps = FeaturePlot(pbmc, reduction = 'tsne', features = c(genename),pt.size=1, order=TRUE, combine=TRUE) + eelTheme() + NoLegend() + xlim(-50,50) + ylim(-60,60)
   ggsave(ps, file=paste("00_TSNE",genename,".png",sep=''), path="./zfConeRNAseqFigure/TSNE/", width = 140*2, height = 105*2, units = "mm")
}
FeaturePlot(pbmc, reduction = 'tsne', features = genelist,pt.size=1, order=TRUE, combine=TRUE, ncol = 3)
# --------------------------------------------------------------------------------
# opsins and others
genelist = c("tulp1a","elovl4b",'ckmt2a','gngt2a',"rho","opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','gnat2','gnat1','saga','arr3a','arr3b')
for(i in 1:length(genelist)) {
   genename = genelist[i]
   ps = FeaturePlot(pbmc, reduction = 'tsne', features = c(genename),pt.size=1, order=TRUE, combine=TRUE) + eelTheme() + NoLegend() + xlim(-50,50) + ylim(-60,60)
   ggsave(ps, file=paste("00_TSNE",genename,".png",sep=''), path="./zfConeRNAseqFigure/TSNE/", width = 140*2, height = 105*2, units = "mm")
}
FeaturePlot(pbmc, reduction = 'tsne', features = genelist,pt.size=1, order=TRUE, combine=TRUE, ncol = 3)
# --------------------------------------------------------------------------------
# opsins and others
genelist = c("tulp1a","elovl4b",'ckmt2a','gngt2a',"rho","opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','gnat2','gnat1','saga','arr3a','arr3b')
for(i in 1:length(genelist)) {
   genename = genelist[i]
   ps = FeaturePlot(pbmc, reduction = 'tsne', features = c(genename),pt.size=1, order=TRUE, combine=TRUE) + eelTheme() + NoLegend() + xlim(-50,50) + ylim(-60,60)
   ggsave(ps, file=paste("00_TSNE",genename,".png",sep=''), path="./zfConeRNAseqFigure/TSNE/", width = 140*2, height = 105*2, units = "mm")
}
FeaturePlot(pbmc, reduction = 'tsne', features = genelist,pt.size=1, order=TRUE, combine=TRUE, ncol = 3)
# --------------------------------------------------------------------------------



# HCs: 6
FeaturePlot(pbmc,reduction = 'tsne', features = c("rem1","CR361564.1",'prkar2ab','rprmb','opn4.1','atp1b1a'))

# BCs: 5, 15
FeaturePlot(pbmc,reduction = 'tsne', features = c("samsn1a","vsx1",'gnb3a','pcp4l1','neurod4','gnao1b','rho'))

# Immature MG: 16
# MG: 11, 9 , 1
FeaturePlot(pbmc,reduction = 'tsne', features = c("cavin2a","rhbg",'slc1a2b','atp1a1b'))

# RGCs = 13 + 10
FeaturePlot(pbmc,reduction = 'tsne', features = c("islr2","rbpms2a",'pou4f1'))

# ngng ACs = 7 -> these are actually larval ACs
# gly ACs = 14
# GABA ACs = 8
FeaturePlot(pbmc,reduction = 'tsne', features = c("slc32a1","pax10",'tfap2a','gad2','slc6a1b','mafa'))
FeaturePlot(pbmc,reduction = 'tsne', features = c("slc32a1","pax10",'tfap2a','gad2','slc6a1b', "chata", "megf10","megf11","tbx2a","tbx2b",'rho'))

# chicken SAC markers:
FeaturePlot(pbmc,reduction = 'tsne', features = c("rho","chata","slc5a7a","slc18a3a","slc18a3b","fezf1","tenm3","zfhx3","foxq2"))


# SAC genes in Yamagata, 2021 chicken RNAseq + tfap2a = INL and tfap2b = INL + GCL
# Exploring briefly, ACs that are tbx2+ in this dataset are clusters: 11,14,15,16,17(nmb),2,21,22,25 (OFF SAC),31(nts, penk),39,42 (penk),43,47,54,6,7 (ON SAC)
# In mouse amacrine dataset it's clearly in clusters 35 (PENK), 43 (GABA+), 51 (GHRH)
# less clearly in 16, 17 (SACs), 23, 24(nGnG-1), 26 (VIP), 3 (AII),30 (nGnG-3), 56 (VG1), 9 (Gly+)
FeaturePlot(pbmc, reduction = 'tsne', features = c("tbx2a","tbx2b","chata","sox10","tfap2a","tfap2b","slc5a7a","slc18a3a","slc18a3b","fezf1"), pt.size=0.2, order=TRUE, combine=TRUE)

FeaturePlot(pbmc, reduction = 'tsne', features = c("tbx2a","tbx2b","nts","nmbb","slc18a3a","penka","penkb","ghrh","maff"), pt.size=0.2, order=TRUE, combine=TRUE)

FeaturePlot(pbmc, reduction = 'tsne', features = c("tbx2a","tbx2b","slc5a7a","slc18a3a","fezf1","megf10"), pt.size=0.2, order=TRUE, combine=TRUE)

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


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Subcluster HCs
# This is a great resource: https://www.frontiersin.org/articles/10.3389/fnana.2016.00077/full
# Thesis is even better: http://uu.diva-portal.org/smash/get/diva2:168823/FULLTEXT01.pdf
# H1s should be lhx1+, while H2 and H3 should be isl1+. H2s? (Candelabrum) should also be TrkA+(ntrk1)
# HC generation relies on ptf1a and foxn4 and also on sall3 (Blackshaw Spalt paper). Foxn4 seems to be important for mid-RPCs so it affects RGCs, PRs and HCs


#  c("rem1","CR361564.1",'prkar2ab','rprmb','opn4.1','atp1b1a')) -> used to id HC cluster


# hcs <- subset(pbmc, idents = c("HC"))
# saveRDS(hcs, file = "./hcs.rds") #preliminary saving to analyze later

hcs = readRDS("~/Documents/eelMolec/zfRNAseq/2020_Hoang_zfRet10x/hcs.rds");


DimPlot(hcs, reduction = "tsne", label =  TRUE, pt.size = 1, label.size = 6) + eelTheme()
hcs2 = SetIdent(hcs, value = "Sample")
DimPlot(hcs2, reduction = "tsne", label =  FALSE, pt.size = 1, label.size = 6) + eelTheme()

FeaturePlot(hcs, reduction = 'tsne', features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))
# smaller cluster has higher nr2e3 levels because it belongs to larval samples mostly
FeaturePlot(hcs, reduction = 'tsne', features = c("rho", "nr2e3",'nrl','saga','sagb'))
FeaturePlot(hcs, reduction = 'tsne', features = c("sema7a", "efna1b",'efna1a','ntng2b','ntng2a'))

#run PCA
hcs <- FindVariableFeatures(hcs, selection.method = "vst", nfeatures = 1000)
hcs <- RunPCA(hcs, features = VariableFeatures(object = hcs))
# # Examine and visualize PCA results a few different ways
# print(photo[["pca"]], dims = 1:5, nfeatures = 5)
# VizDimLoadings(photo, dims = 1:2, reduction = "pca")
# DimPlot(photo, reduction = "pca")
# DimHeatmap(photo, dims = 1, cells = 500, balanced = TRUE)
# DimHeatmap(photo, dims = 1:15, cells = 500, balanced = TRUE)

# Determine number of clusters using elbow plot (var explained)
ElbowPlot(hcs, ndims=50)

# Cluster cells
# nFeatures = 1448 gives too many clusters
# nFeatures = 100 seems like not enough info; just separates by sample
# nFeatures = 200 seems like not enough info either

# nFeatures = 150 apparently good compromise
# dims = 8, resolution = 0.6 -> not bad but still some mixing
# dims = 23, resolution = 0.6 -> almost there; there are still some M-cones mixed with L-cones. UV and S have to be separated further

# Settling for nFeatures = 1000, ndim = 12, res = .3
upDimLimit=24;
hcs <- FindNeighbors(hcs, dims = 1:upDimLimit)
hcs <- FindClusters(hcs, resolution = .5)

DimPlot(hcs, reduction = "pca", label = TRUE, pt.size = 1)

# UMAP using the same PCA dimensions for prettier visualization
hcs <- RunUMAP(hcs, dims = 1:upDimLimit)
DimPlot(hcs, reduction = "umap", label = TRUE, pt.size = 1)

hcs2 = SetIdent(hcs, value = "Sample")
DimPlot(hcs2, reduction = "umap", label = TRUE, pt.size = 1)

# or use tSNE (clustering can't separate some M cones from L cones)
hcs <- RunTSNE(hcs, dims = 1:upDimLimit)
DimPlot(hcs, reduction = "tsne", label = TRUE, pt.size = 1)+ eelTheme()
# ps = DimPlot(hcs, reduction = "tsne", label = TRUE, pt.size = 1)+ eelTheme()
# ggsave(ps, file="adult_coneClusters03.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

hcs2 = SetIdent(hcs, value = "Sample")
DimPlot(hcs2, reduction = "tsne", label = TRUE, pt.size = 1)



FeaturePlot(hcs, reduction = 'tsne', features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))
FeaturePlot(hcs, reduction = 'tsne', features = c("rho", "nr2e3","nr2f6","crx",'syt5a','syt5b'))

FeaturePlot(hcs, reduction = 'tsne', features = c("rho", "onecut1","onecut2","lhx1a","efna1a","isl1","ntrk3a"))
FeaturePlot(hcs,reduction = 'tsne', features = c("rem1","CR361564.1",'prkar2ab','rprmb','opn4.1','atp1b1a'))

# trying to find DEGs for clusters
# In Yamagata/Sanes 2021:
#     lhx1+/isl1- have 2 clusters, one that is oxt+ during early embryo then ipcef1+ later
#     lhx1-/isl1+ have 3 clusters diff by ntrk1, egfr and ltk
# seems like opn4.1 and opn9 are good; isl1/efna1a splits things in half
# In mouse, HCs are lhx1+
# sall3b
# good dividers for clusters
FeaturePlot(hcs,reduction = 'tsne', features = c("lhx1a",'vsx1', 'isl1','ntrk3a',"opn4.1",'opn9','egfra','ltk'))
FeaturePlot(hcs,reduction = 'tsne', features = c("opn4.1","lhx1a",'vsx1', 'isl1','ntrk3a','slc6a1l','gad1b','opn9','ret'))
FeaturePlot(hcs,reduction = 'tsne', features = c("hmx3a","plekhd1",'tmtc1','nrxn1a','barhl2','sal3b','irs1','cacng5a'))
FeaturePlot(hcs,reduction = 'tsne', features = c("efna1a","arhgef18b",'rpap1','sgk494a','barhl2','plekhd1','shisa7b','grip1'))
FeaturePlot(hcs,reduction = 'tsne', features = c("sept7a","sema3fb",'cdh6','cd9a','ndnf','msna','foxp4'))

# some controls
FeaturePlot(hcs,reduction = 'tsne', features = c('rho','syt5b','fabp7a',"neurod1",'gnat1','gnat2','syt5a','prdm1a'))


FeaturePlot(hcs,reduction = 'tsne', features = c("opn4.1","lhx1a",'nrgna','neurod4','vsx1', 'isl1','ntrk3a'))
FeaturePlot(hcs,reduction = 'tsne', features = c('opn9','nrxn1a','ntrk3a',"sall3b","PDE6H",'isl1','efna1a','vsx1'))
FeaturePlot(hcs,reduction = 'tsne', features = c('ntrk1','ntrk2a',"ntrk2b","ntrk3a","ntrk3b"))

VlnPlot(hcs, features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))


FeaturePlot(hcs,reduction = 'tsne', features = c('rassf4','calm2b',"cd9a","znf385b","ndnf",'msna','foxp4','ppdpfb','prdm13'))

FeaturePlot(hcs,reduction = 'tsne', features = c('onecut1', 'onecutl', 'onecut2', 'onecut3a'))
FeaturePlot(hcs,reduction = 'tsne', features = c('epha2a', 'epha3', 'epha7', 'epha4b'))
FeaturePlot(hcs,reduction = 'tsne', features = c('cx52.7', 'cx52.6', 'cx52.9','gad2','lrit1a','tfap2a'))


# Clusters:
# there is clear separation between cells derived from AdR1 (less rod contamination) and AdR2 (more rod contamination)
# 0: H?s -> isl1-, lhx1a-, tmtc1+, slc6a1l+, gad1b+
# 1: H2/3s -> isl1+
# 2: H2/3s -> isl1+, opn9+,
# 3: ipHCs -> isl1+, opn4.1+, sall3b, nrxn1a+, arhgef18b+, rpap1+, sept7a+, sema3fb+, cdh6+, foxp4+ (subsplit comes from AdR1 and AdR2)
# 4: larval HCs
# 5: H2/3s -> isl1+, hmx3a+, irs1+
# 6: AdR1 cells, opn9+
# 7: H1s -> lhx1a+, barhl2+, slc6a1l+, gad1b+


FeaturePlot(hcs,reduction = 'tsne', features = c("lhx1a",'vsx1', 'isl1','ntrk3a',"opn4.1",'opn9','gad1b','hmx3a'))

FeaturePlot(hcs,reduction = 'tsne', features = c("lhx1a","rom1a"))
# Putting all things together and embodying a joiner and not a splitter, it seems that loosely
# H1s = lxh1a+
# H2s = isl1+, ntrk3a+ (mostly), opn9+ (mostly)
# H3s = isl1+ opn4.1+
# larvalHcs
# H?s = lhx1a-, isl1-

new.cluster.ids <- c("HCx","HC2","HC2","HC3","lHC","HC2","HC2","HC1")
names(new.cluster.ids) <- levels(hcs)
hcs <- RenameIdents(hcs, new.cluster.ids)
Idents(hcs) <- factor(x = Idents(hcs), levels = c("lHC","HC1","HC2","HC3","HCx"))

DimPlot(hcs, reduction = "tsne", label = TRUE, pt.size = 1)+ eelTheme()


FeaturePlot(hcs,reduction = 'tsne', features = c("rem1","CR361564.1",'prkar2ab','rprmb','opn4.1','atp1b1a','efna1a'))
VlnPlot(hcs, features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))
VlnPlot(hcs, features = c('rho','syt5b','fabp7a',"neurod1",'gnat1','gnat2','syt5a','prdm1a'))
VlnPlot(hcs, features = c("lhx1a",'vsx1', 'isl1','ntrk3a',"opn4.1",'opn9','gad1b','hmx3a','efna1a'))

# Trying out Mike Country's review receptors for Ca (https://pubmed.ncbi.nlm.nih.gov/6385010/)
#NMDA receptors: nothing great. grin1a in all HCs
FeaturePlot(hcs,reduction = 'tsne', features = c('grin1a','grin1b','grin2aa','grin2ab','grin2ba','grin2bb','grin2ca','grin2da','grin2db','grin3ba','grinaa','grinab'))
# Ryanodine receptors: nothing
FeaturePlot(hcs,reduction = 'tsne', features = c('ryr1a','ryr1b','ryr2a','ryr2b','ryr3'))
# voltage gated calcium channels
FeaturePlot(hcs,reduction = 'tsne', features = c('cacna1aa','cacna1ab','cacna1ba','cacna1bb','cacna1c','cacna1da','cacna1db','cacna1ea','cacna1fa','cacna1fb'))
FeaturePlot(hcs,reduction = 'tsne', features = c('itpr1a','itpr1b','itpr2','itpr3','itprip','itprid1','itprid2'))
FeaturePlot(hcs,reduction = 'tsne', features = c('panx1a','panx1b','panx2','panx3'))


# PC_ 3
# Positive:  slc6a1l, sall3b, ret, gad1b, opn4.1, barhl2, tmtc1, rpap1, arhgef18b, lhx1a
# shisa7b, sgk494a, valopa, vipr2, grip1, stk17a, rassf4, calm2b, irs1, si:dkeyp-41f9.3
# hmp19, sept7a, cd9a, nrxn1a, znf385b, ndnf, msna, foxp4, ppdpfb, prdm13
# Negative:  opn9, ntrk3a, hmx3a, isl1, efna1a, plekhd1, cryba4, nme2b.2, saga, vamp1
# CU861453.2, lmo1, elmo1, zgc:162144, lrrtm1, guca1a, si:ch211-113d22.2, mt2, BX004774.2, PDE6H
# guca1b, pdca, gnat1, kcnv2a, sagb, tnfaip8l3, gngt1, ptn, rho, tuba1a
# PC_ 4
# Positive:  si:ch211-56a11.2, lhx1a, ret, slc6a1l, opn9, ntrk3a, barhl2, gad1b, RAPGEF2 (1 of many), slc2a3a
# nme2b.2, hmx3a, lmo1, si:dkey-1d7.3, cd82a, palm1b, ip6k2b, slc16a9b, tnfaip8l3, syt5b
# cspg5b, h3f3b.1.2, histh1l, cry1bb, ppdpfb, sept7a, vsx1, si:ch73-1a9.3, rai14, lin7a
# Negative:  rpap1, arhgef18b, shisa7b, sgk494a, cacng5a, grip1, valopa, vipr2, si:dkeyp-41f9.3, irs1
# sall3b, PDE6H, saga, gnat1, si:ch211-113d22.2, guca1a, guca1b, zgc:162144, gngt1, rho
# rassf4, sagb, pdca, foxp4, pde6a, rbp4l, BX004774.2, msna, CNGA1 (1 of many), cnga1
# PC_ 5
# Positive:  cabp2a, rs1a, saga, guca1a, nrn1lb, cabp5a, si:ch73-256g18.2, zgc:162144, PDE6H, rgs16
# efna1b, guca1b, sagb, syt5a, BX004774.2, gnao1b, eml1, rpap1, ptmab, pdca
# pde6a, si:ch1073-83n3.2, opn4.1, sypb, lrit1a, gnat1, si:ch211-113d22.2, vsx1, pcp4a, atp1b2a
# Negative:  plekhd1, ptn, nme2b.2, CU861453.2, slc6a1l, tmtc1, barhl2, ret, onecut1, lhx1a
# gad1b, zgc:158463, nme2b.1, efna1a, sgol1, hmx3a, cdca8, kif4, knstrn, top2a
# ube2c, lbr, rps12, h2afx, rps29, anp32b, fam60al.1, kiaa0101, rpl13, rpl36a
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# Subcluster larval HCs
# Does not work well, not enough cells in any grouping to do this properly
lhcs <- subset(hcs, idents = c("lHC"))
# saveRDS(hcs, file = "./hcs.rds") #preliminary saving to analyze later

DimPlot(lhcs, reduction = "tsne", label =  TRUE, pt.size = 1, label.size = 6) + eelTheme()
lhcs2 = SetIdent(lhcs, value = "Sample")
DimPlot(lhcs2, reduction = "tsne", label =  FALSE, pt.size = 1, label.size = 6) + eelTheme()

FeaturePlot(lhcs, reduction = 'tsne', features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))
# smaller cluster has higher nr2e3 levels because it belongs to larval samples mostly
FeaturePlot(lhcs, reduction = 'tsne', features = c("rho", "nr2e3",'nrl','saga','sagb'))
FeaturePlot(lhcs, reduction = 'tsne', features = c("sema7a", "efna1b",'efna1a','ntng2b','ntng2a'))

#run PCA
lhcs <- FindVariableFeatures(lhcs, selection.method = "vst", nfeatures = 1000)
lhcs <- RunPCA(lhcs, features = VariableFeatures(object = lhcs))
# # Examine and visualize PCA results a few different ways
# print(photo[["pca"]], dims = 1:5, nfeatures = 5)
# VizDimLoadings(photo, dims = 1:2, reduction = "pca")
# DimPlot(photo, reduction = "pca")
# DimHeatmap(photo, dims = 1, cells = 500, balanced = TRUE)
# DimHeatmap(photo, dims = 1:15, cells = 500, balanced = TRUE)

# Determine number of clusters using elbow plot (var explained)
ElbowPlot(lhcs, ndims=50)

# Cluster cells
# nFeatures = 1448 gives too many clusters
# nFeatures = 100 seems like not enough info; just separates by sample
# nFeatures = 200 seems like not enough info either

# nFeatures = 150 apparently good compromise
# dims = 8, resolution = 0.6 -> not bad but still some mixing
# dims = 23, resolution = 0.6 -> almost there; there are still some M-cones mixed with L-cones. UV and S have to be separated further

# Settling for nFeatures = 1000, ndim = 12, res = .3
upDimLimit=12;
lhcs <- FindNeighbors(lhcs, dims = 1:upDimLimit)
lhcs <- FindClusters(lhcs, resolution = 2)

DimPlot(lhcs, reduction = "pca", label = TRUE, pt.size = 1)

# UMAP using the same PCA dimensions for prettier visualization
lhcs <- RunUMAP(lhcs, dims = 1:upDimLimit)
DimPlot(lhcs, reduction = "umap", label = TRUE, pt.size = 1)

lhcs2 = SetIdent(lhcs, value = "Sample")
DimPlot(lhcs2, reduction = "umap", label = TRUE, pt.size = 1)

# or use tSNE (clustering can't separate some M cones from L cones)
lhcs <- RunTSNE(lhcs, dims = 1:upDimLimit)
DimPlot(lhcs, reduction = "tsne", label = TRUE, pt.size = 1)+ eelTheme()
# ps = DimPlot(hcs, reduction = "tsne", label = TRUE, pt.size = 1)+ eelTheme()
# ggsave(ps, file="adult_coneClusters03.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

lhcs2 = SetIdent(lhcs, value = "Sample")
DimPlot(lhcs2, reduction = "tsne", label = TRUE, pt.size = 1)



FeaturePlot(lhcs, reduction = 'tsne', features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))
FeaturePlot(lhcs,reduction = 'tsne', features = c("lhx1a",'vsx1', 'isl1','ntrk3a',"opn4.1",'opn9','gad1b','hmx3a'))

# ----------------------------------------------------------------------------------------------------------------
# getting % and Mean expression for all genes into csv.
# DoHeatmap and DotPlot use the @scale.data slot for average expression display, which z-scored expression values (for example, as those used in PCA).
# Cells with a value > 0 represent cells with expression above the population mean (a value of 1 would represent cells with expression 1SD away from the population mean). Hope that helps!
# This is very detailed explanation: https://github.com/satijalab/seurat/issues/2798

# Average expression in non-log scale:
avgExp = AverageExpression(hcs, slot="counts")
avgExp = avgExp$RNA
# avgExp = head(avgExp,10)
nI = length(levels(Idents(hcs)))
colnames(avgExp) = paste("avg", levels(Idents(hcs)), sep = "")

# Calculate mean expression across all subtypes
meanAvg = as.data.frame(rowMeans(avgExp))
colnames(meanAvg) = "avgExp"
head(meanAvg)
avgExp = cbind(meanAvg,avgExp)
# avgExp$avgMean = rowMeans(avgExp)
# avgExp[,c(ncol(avgExp),1:(ncol(avgExp)-1))]


# Percent expression can be obtained through DotPlot, but not Average expression, as this is log z-cored normalized data # avgExpZ = dPlot$data$avg.exp
# Percent expression in larvae:
dPlot = DotPlot(hcs, features=rownames(avgExp))
pctExp = dPlot$data$pct.exp
pctExp = as.data.frame(matrix(pctExp,nrow = length(pctExp)/nI, ncol=nI));

rownames(pctExp) = rownames(avgExp)
colnames(pctExp) = paste("pct", levels(Idents(hcs)), sep = "")

# Calculate mean percent expression across all ages
meanPct = as.data.frame(rowMeans(pctExp))
colnames(meanPct) = "pctExp"
head(meanPct)
pctExp = cbind(meanPct,pctExp)
# pctExp$pctMean = rowMeans(pctExp)
# pctExp[,c(ncol(pctExp),1:(ncol(pctExp)-1))]


hcs_Hoang = cbind(avgExp,pctExp)



# these are duplicated genes with lower vs upper case
dups=toupper(c("arid5b","asph", "aste1", "crip2", "ctbp1", "dab2", "eif1b", "flnb", "frmd7", "galnt10", "grxcr1", "gse1", "hist1h4l", "hspb11", "kif1c", "lamp1", "maf1", "pamr1", "pcdh20", "pde6h", "phlpp2", "psmb10", "ptp4a3", "reep6", "rfesd", "rgs9bp", "rps17", "shank2", "slc16a7", "slc25a10", "slc6a13", "slc9a1", "srbd1", "tatdn3", "tenm3", "tmem178b", "tmem241", "tom1l2", "tp53inp2", "trappc9", "tsc22d3", "ube2o", "zc3h12a", "znf423"))
for(i in 1:length(dups)) {
   rownames(hcs_Hoang)[rownames(hcs_Hoang) == dups[i]] = paste(dups[i],"_ii",sep="")
}
# now everything can be lower cased
rownames(hcs_Hoang) = tolower(rownames(hcs_Hoang))




# hcs_Hoang = hcs_Hoang[c("pctMean","pctRod","pctUV","pctS","pctM1","pctM3","pctM4","pctL","avgMean","avgRod","avgUV","avgS","avgM1","avgM3","avgM4","avgL")]



scangenes=c('syt5a','syt5b','nr2e3','lhx1a','isl1','gad1b','vsx1','ntrk3a','opn9');
hcs_Hoang[scangenes,]



# still need to manually add symbol to first column
write.csv(hcs_Hoang,"./hcs_Hoang.csv",quote=FALSE)
