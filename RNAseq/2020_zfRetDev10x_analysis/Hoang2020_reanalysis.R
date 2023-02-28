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
# -------------------------------------------------------------------
# Reanalysis of Hoang 2020 zebrafish developmental dataset
# Seurat_v2 version available here: http://bioinfo.wilmer.jhu.edu/jiewang/scRNAseq/Seurat_objects/
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Load the 10x dataset (developemntal data, no treatment, 28845 cells), after updating to Seurat_v3
pbmc = readRDS("./zfDev_pbmc_v3.rds");
# pca, tsne and umap already done; will use the original clustering
# ------------------------------------------------------------------------------------------
# Here are the correspondences:
   # Retinal progenitor cells (RPC): 0
   # Photoreceptor precursors (PRP): 3
   # Cones_adults: 4
   # Cones_larval: 12
   # Rods: 2
   # Horizontal Cells (HC): 6
   # Bipolar Cells (BC): 5
   # Bipolar Cells larval (BC): 15
   # Immature Muller Glia (MGi): 16
   # Muller Glia (MG): 11, 9 , 1
   # Retinal Ganglion Cells (RGC): 13, 10
   # Amacrine Cells Glycinergic (ACgly): 18
   # Amacrine Cells GABAergic (ACgaba): 8
   # Amacrine Cells ngng (ACngng): 7

# Assigning cell type identity to clusters (it's basically just renaming)
new.cluster.ids <- c("RPC","MG","R","PRP","C","BC","HC","AClarval","ACgaba","MG","RGC","MG","Clarval","RGClarval","ACgly","BClarval","MGi")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

# forcing order to be able to define colors
Idents(pbmc) <- factor(x = Idents(pbmc), levels = c("RPC","PRP","Clarval","C","R","HC","BClarval","BC","AClarval","ACgaba","ACgly","RGClarval","RGC","MGi","MG"))
pbmc_colors = c('#DADADA','#dfdac8','#dcc360','#ffd429','#7d7d7d','#FC7715','#ccf2ff','#1B98C3','#3DF591','#3DF5C3','#56F53D','#F53D59','#BB0622','#EA9D81','#A2644E','#7E4835','#613728')

DimPlot(pbmc, reduction = "tsne", label =  FALSE, pt.size = 0.2, label.size = 6, repel = TRUE, cols = pbmc_colors) + plotTheme() +NoLegend()
# ggsave("adult_allClusters_tSNE.png", path=exportDir, width = 140, height = 105, units = "mm")
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 1, label.size = 6, repel = TRUE, cols = pbmc_colors) + plotTheme()
# ggsave("adult_allClusters_UMAP.png", path=exportDir, width = 140, height = 105, units = "mm")

# store naming information to object
pbmc$cellTypes <- Idents(object = pbmc)
# switch to Sample to breakdown by age
pbmc = SetIdent(pbmc, value = "Sample")
DimPlot(pbmc, reduction = "tsne", label =  TRUE, pt.size = 1, label.size = 6, repel = TRUE) + plotTheme()
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 1, label.size = 6, repel = TRUE) + plotTheme()



# And switch back to cluster names
pbmc = SetIdent(pbmc, value = "cellTypes")
# plot some genes of interest
FeaturePlot(pbmc, reduction = 'tsne', features = c("crx","otx5", "onecut1","vsx2","atoh7","ptf1a"), pt.size = 0.5)


# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
# RODS
# CONCLUSION: There are 2 rod clusters; larval cells are in small cluster with adult cells (maybe immature rods in adult?).
rods <- subset(pbmc, idents = c("R"))
# saveRDS(rods, file = "./rods.rds") #preliminary saving to analyze later

DimPlot(rods, reduction = "tsne", label =  TRUE, pt.size = 1, label.size = 6) + plotTheme()
# replot using age
rods2 = SetIdent(rods, value = "Sample")
new.cluster.ids <- c("2.5dpf","5dpf","Ad1","2.5dpf","Ad5","4dpf","Ad4","Ad3","Ad2" )
names(new.cluster.ids) <- levels(rods2)
rods2 <- RenameIdents(rods2, new.cluster.ids)
# reorganization
Idents(rods) <- factor(x = Idents(rods2), levels = c("2.5dpf","4dpf","5dpf","Ad1","Ad2","Ad3","Ad4","Ad5"))
DimPlot(rods2, reduction = "tsne", label = TRUE, pt.size = 1, label.size = 6, repel = TRUE, order = c("2.5dpf","4dpf","5dpf","AdR1","AdR2","AdR3","AdR4","AdR5")) + plotTheme()

# flag cells that belong to adult samples
zfAFlag = grepl("AdR",names(Idents(rods)));
rods = AddMetaData(rods,metadata=zfAFlag, col.name = "zfA");
rods = SetIdent(rods, value = "zfA")
FeaturePlot(rods, reduction = 'tsne', features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))
# smaller cluster has higher nr2e3 levels because it belongs to larval samples mostly
FeaturePlot(rods, reduction = 'tsne', features = c("rho", "nr2e3",'nrl','saga','sagb'))
# ------------------------------------------------------------------------------------------
# plot layouts for direct export
# Photoreceptor subtype markers
p = FeaturePlot(rods, reduction = 'tsne', features = c("rho","nr2e3","gnat1","gnat2","arr3b","arr3a","thrb","opn1sw1", "opn1sw2",'opn1mw1','opn1lw1'),
                pt.size=1, order=FALSE, combine=FALSE)
lay = rbind(c(1,2,3,NA),c(4,5,6,7),c(8,9,10,11))
grid.arrange(grobs = p, layout_matrix = lay)
ps = arrangeGrob(grobs = p, layout_matrix = lay)
# ggsave(ps, file="adult_rodsMarkers.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

# Photoreceptor subtype markers in whole dataset (rod genes are high in adult cells)
p = FeaturePlot(pbmc, reduction = 'tsne', features = c("rho","nr2e3","gnat1","gnat2","arr3b","arr3a","thrb","opn1sw1", "opn1sw2",'opn1mw1','opn1lw1'),
                pt.size=1, order=FALSE, combine=FALSE)
lay = rbind(c(1,2,3,NA),c(4,5,6,7),c(8,9,10,11))
grid.arrange(grobs = p, layout_matrix = lay)
ps = arrangeGrob(grobs = p, layout_matrix = lay)
# ggsave(ps, file="adult_allMarkers.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

# Same using dot plot
ps = DotPlot(pbmc, features = c("rho","nr2e3","gnat1","gnat2","arr3b","arr3a","thrb","opn1sw1", "opn1sw2",'opn1mw1','opn1lw1')) + plotTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ps

# Transcription Factors (TF) of interest
FeaturePlot(pbmc, reduction = 'tsne', features = c("foxq2"), pt.size=0.5, order=FALSE, combine=FALSE)
FeaturePlot(pbmc, reduction = 'tsne', features = c("nr2e3","foxq2","skor1a","lrrfip1a","xbp1","tbx2a","tbx2b"), pt.size=0.5, order=FALSE, combine=TRUE)

# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
# CONES
# Subclustering cones (will subsequently remerge with rods, photoreceptor progenitors and retinal progenitors to track gene expression)

# CONCLUSION: No clear separation by subtypes combining adult and larval cones; need to do separately for adults and larval; data mainly varies according to age and not subtypes
photo <- subset(pbmc, idents = c("C","Clarval"))
ps = DimPlot(photo, reduction = "tsne", label=TRUE)+ plotTheme()
ps

# clear separation between adult and larval cones
photo2 = SetIdent(photo, value = "Sample")
ps = DimPlot(photo2, reduction = "tsne", label=TRUE)+ plotTheme()
ps

# flag cells form Adult samples
zfAFlag = grepl("AdR",names(Idents(photo)));
photo = AddMetaData(photo,metadata=zfAFlag, col.name = "zfA");
photo = SetIdent(photo, value = "zfA")
# ------------------------------------------------------------------------------------------
# Adult cone cluster
## CONCLUSION: Variable Features contain many highly expressed genes in ALL subtypes, such that clustering is driven by variablity in these rather that subtype identity

# separate larval cones
photo_larval = subset(photo, idents = FALSE)
# saveRDS(photo_larval, file = "./photo_Larval.rds") #preliminary saving to analyze later
# separate adult cones
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

# photo <- JackStraw(photo, num.replicate = 100)
# photo <- ScoreJackStraw(photo, dims = 1:20)
# JackStrawPlot(photo, dims = 1:20)

# Or just use elbow plot (var explained)
ElbowPlot(photo, ndims=50)

# Cluster cones
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
DimPlot(photo, reduction = "tsne", label = TRUE, pt.size = 1)+ plotTheme()

# clustering is influenced by sample
photo2 = SetIdent(photo, value = "Sample")
DimPlot(photo2, reduction = "tsne", label=TRUE)+ plotTheme()

FeaturePlot(photo, reduction = 'tsne', features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))
FeaturePlot(photo, reduction = 'tsne', features = c("rho", "nr2e3","nr2f6","crx"))
VlnPlot(photo, features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))

# ------------------------------------------------------------------------------------------
# plot layouts for direct export

#OPSINS
p = FeaturePlot(photo, reduction = 'tsne', features = c("rho","opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2'),
                pt.size=1, order=TRUE, combine=FALSE)
lay = rbind(c(1,NA,NA,NA),c(2,3,NA,NA),c(4,5,6,7),c(8,9,NA,NA))
grid.arrange(grobs = p, layout_matrix = lay)
# ps = arrangeGrob(grobs = p, layout_matrix = lay)
# ggsave(ps, file="adult_Opsins.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

#ROD STUFF (known)
p = FeaturePlot(photo, reduction = 'tsne', features = c("rho",'gnat1',"gnat2","saga",'pde6g',"guca1b","grk1a",'crx','nr2e3','nrl'),pt.size=1, order=FALSE, combine=FALSE)
lay = rbind(c(1,2,3,NA),c(4,5,6,7),c(8,9,10,NA))
grid.arrange(grobs = p, layout_matrix = lay)


DotPlot(photo, features = c("rho","nr2e3","gnat1","gnat2","arr3a","arr3b","thrb","opn1sw1", "opn1sw2",'opn1mw4','opn1mw3','opn1mw2','opn1mw1','opn1lw1','opn1lw2')) + plotTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1)) 
# ------------------------------------------------------------------------------------------
# Assigning cone subtype identity to clusters
# Cluster 12 (most cells from AdR1) has markers for immature photoreceptors but also high counts of ribosomal genes; will discard it from further analysis
# Can't separate UV and S cones at this stage
new.cluster.ids <- c("L","L","L","L","L","L","L","L","US","M3","M4","M1","Ribo","L")
names(new.cluster.ids) <- levels(photo)
photo <- RenameIdents(photo, new.cluster.ids)
Idents(photo) <- factor(x = Idents(photo), levels = c("Ribo","L","M1","M3","M4","US"))
photo_colors = c('#AAAAAA','#CC2C2A','#04CD22','#04CD22','#CDCD04','#4669F2')
DimPlot(photo, reduction = "tsne", label = TRUE, pt.size = 2, cols = photo_colors) + NoLegend()+ plotTheme()
DimPlot(photo, reduction = "umap", label = TRUE, pt.size = 1, cols = photo_colors) + NoLegend()

VlnPlot(photo, features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))
# ------------------------------------------------------------------------------------------
# Attempting to separate UV and S cones
# CONCLUSION: some separation but the number of cells is quite low
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
DimPlot(us, reduction = "tsne", label=TRUE, pt.size=4) + NoLegend()+ plotTheme()
# ------------------------------------------------------------------------------------------
# plot layouts for direct export

#OPSINS
p = FeaturePlot(us, reduction = 'tsne', features = c("rho","opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2'),
                pt.size=1, order=TRUE, combine=FALSE)
lay = rbind(c(1,NA,NA,NA),c(2,3,NA,NA),c(4,5,6,7),c(8,9,NA,NA))
grid.arrange(grobs = p, layout_matrix = lay)
# ps = arrangeGrob(grobs = p, layout_matrix = lay)
# ggsave(ps, file="adult_coneUS02.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

DotPlot(us, features = c("rho","nr2e3","gnat1","gnat2","arr3b","arr3a","thrb","opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2')) + plotTheme()+ theme(axis.text.x = element_text(angle = 45, hjust=1)) 

# ------------------------------------------------------------------------------------------
FeaturePlot(us, reduction = 'tsne', features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))
VlnPlot(us, features = c("opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2','si:busm1-57f23.1'))
FeaturePlot(us, reduction = 'tsne', features = c('tbx2a','tbx2b','foxq2'))

new.cluster.ids <- c("S","UV","S","S")
names(new.cluster.ids) <- levels(us)
us <- RenameIdents(us, new.cluster.ids)

Idents(us) <- factor(x = Idents(us), levels = c("UV","S"))
us_colors = c('#B540B7','#4669F2')
DimPlot(us, reduction = "tsne", label = TRUE, pt.size = 4, cols = us_colors) + NoLegend()+ plotTheme()

# ------------------------------------------------------------------------------------------
# Remerging adult cones to perform DESeq2 analysis between subtypes

lm = subset(photo, idents = c("L","M1","M3","M4")); #excluded Ribo cluster
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
DimPlot(photo, reduction = "tsne", label=TRUE, pt.size = 2) + NoLegend() +plotTheme()
# ------------------------------------------------------------------------------------------
# plot layouts for direct export

#OPSINS
p = FeaturePlot(photo, reduction = 'tsne', features = c("rho","opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2'),
                pt.size=1, order=TRUE, combine=FALSE)
lay = rbind(c(1,NA,NA,NA),c(2,3,NA,NA),c(4,5,6,7),c(8,9,NA,NA))
grid.arrange(grobs = p, layout_matrix = lay)
# ps = arrangeGrob(grobs = p, layout_matrix = lay)
# ggsave(ps, file="adult_coneClusters06.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

DotPlot(photo, features = c("rho","nr2e3","gnat1","gnat2","arr3b","arr3a","thrb","opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2')) + plotTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
# ------------------------------------------------------------------------------------------
# save result to skip analysis up to this point
# saveRDS(photo, file = "./photo_AdultAll.rds")
# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
# Analysis can be restarted here
# try(dev.off(dev.list()["RStudioGD"]),silent=TRUE)
# try(dev.off(),silent=TRUE)
# rm(list=ls())
# photo = readRDS("~/Documents/LiMolec/otherRNAseq/zfRet_HoangBlackshaw2020/photo_Adult.rds");
# ------------------------------------------------------------------------------------------
# DESeq2 example: L-cones
DESeq2_L = FindMarkers(photo, c("L"), c("M1","M3","M4","S","UV"), test.use="DESeq2");
AvgExp_L = AverageExpression(photo, features=rownames(DESeq2_L))
AvgExp_L = AvgExp_L$RNA
DESeq2_L = merge.data.frame(DESeq2_L,AvgExp_L, by="row.names", sort=FALSE)
DESeq2_L$baseMean = rowMeans(DESeq2_L[,c("L","M1","M3","M4","S","UV")])
rownames(DESeq2_L) = DESeq2_L$Row.names
DESeq2_L = DESeq2_L[,c("baseMean","avg_logFC","p_val","p_val_adj","pct.1","pct.2","L","M1","M3","M4","S","UV")]
PctExp_L = DotPlot(photo, features=rownames(DESeq2_L))
PctNames = rownames(PctExp_L$data)[1:(length(rownames(PctExp_L$data))/4)]
PctExp_L = PctExp_L$data$pct.exp
PctExp_L = matrix(PctExp_L,nrow = length(PctExp_L)/4, ncol=4);
rownames(PctExp_L) = PctNames
colnames(PctExp_L) = c("pctL","pctM1","pctM3","pctM4","pctS","pctUV")
DESeq2_L = merge.data.frame(DESeq2_L,PctExp_L, by="row.names", sort=FALSE)
rownames(DESeq2_L) = DESeq2_L$Row.names
write.csv(DESeq2_L,"./DESeq2_L.csv")
# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
# cell counter
table(Idents(photo))
# L      M1   M3   M4    S   UV 
# 1587   53   80   68   80   22

# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
# remerge rods with cones still split by subtype
# remerge into single dataset (1890 cones + 3097 rods = 4987 PRs)
r = subset(pbmc, idents = c("R"));

photoAll = merge(x = lm, y = c(us,r))
photoAll = subset(photoAll, idents = c("R","UV","S","M1","M3","M4","L"))
photoAll

photoAll = ScaleData(photoAll)
photoAll <- FindVariableFeatures(photoAll, selection.method = "vst", nfeatures = 200)

photoAll <- RunPCA(photoAll, features = VariableFeatures(object = photoAll))
ElbowPlot(photoAll, ndims=50)
# if ever wanted to recluster again
# photo <- FindNeighbors(photoAll, dims = 1:upDimLimit)
# photo <- FindClusters(photoAll, resolution = 1)

# ordering identities for plotting and assigning colors
Idents(photoAll) <- factor(x = Idents(photoAll), levels = c("R","UV","S","M1","M3","M4","L"))
prcolors = c("#7D7D7D","#B540B7","#4669F2","#04CD22","#04CD22","#CDCD04","#CC2C2A")

# save result to skip analysis up to this point
# saveRDS(photoAll, file = "./photo_AdultAll.rds")
# ------------------------------------------------------------------------------------------
# Analysis can be restarted here
# try(dev.off(dev.list()["RStudioGD"]),silent=TRUE)
# try(dev.off(),silent=TRUE)
# rm(list=ls())
# photoAll = readRDS("./photo_AdultAll");
# prcolors = c("#7D7D7D","#B540B7","#4669F2","#04CD22","#04CD22","#CDCD04","#CC2C2A")
# ------------------------------------------------------------------------------------------

DimPlot(photoAll, reduction = "pca", label=TRUE, cols=prcolors)

photoAll <- RunUMAP(photoAll, dims = 1:upDimLimit)
DimPlot(photoAll, reduction = "umap", label=TRUE, pt.size = 1, cols=prcolors) + NoLegend()

photoAll <- RunTSNE(photoAll, dims = 1:upDimLimit)
DimPlot(photoAll, reduction = "tsne", label=FALSE, pt.size = 2, cols=prcolors) +plotTheme()

# final results 
   # rod contamination in adult samples pulls identified subtypes into rod cluster
   # uneven numbers per cluster with overrepresentation of rods and L cones.
FeaturePlot(photoAll, reduction = 'tsne', features = c("rho","opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2'),
            pt.size=0.2, order=FALSE, combine=TRUE)
# ------------------------------------------------------------------------------------------
# plot layouts for direct export

#OPSINS
p = FeaturePlot(photoAll, reduction = 'tsne', features = c("rho","opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2'),
                pt.size=1, order=FALSE, combine=FALSE)
lay = rbind(c(1,NA,NA,NA),c(2,3,NA,NA),c(4,5,6,7),c(8,9,NA,NA))
grid.arrange(grobs = p, layout_matrix = lay)

# OPSINS + KNOWN MARKERS
 DotPlot(photoAll, features = c("rho","nr2e3","gnat1","gnat2","arr3b","arr3a","thrb","opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2')) + plotTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
# TFs in DEGs
ps = DotPlot(photoAll, features = c("nr2e3", "foxq2", "skor1a","sall1a","lrrfip1a","xbp1","tbx2a","tbx2b")) + plotTheme() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ps
# ------------------------------------------------------------------------------------------
# cell counter
table(Idents(photoAll))

# R      UV    S   M1   M3   M4   L 
# 3097   22   80   53   80   68   1587 

# ------------------------------------------------------------------------------------------
# getting % and Mean expression for all genes into csv.
# DoHeatmap and DotPlot use the @scale.data slot for average expression display, which z-scored expression values (for example, as those used in PCA).
# Cells with a value > 0 represent cells with expression above the population mean (a value of 1 would represent cells with expression 1SD away from the population mean). Hope that helps!
# This is very detailed explanation: https://github.com/satijalab/seurat/issues/2798

# Average expression in non-log scale:
avgExp = AverageExpression(photoAll, slot="counts")
avgExp = avgExp$RNA

nI = length(levels(Idents(photoAll)))
colnames(avgExp) = paste("avg", levels(Idents(photoAll)), sep = "")

# Calculate mean expression across all subtypes
meanAvg = as.data.frame(rowMeans(avgExp))
colnames(meanAvg) = "avgExp"
head(meanAvg)
avgExp = cbind(meanAvg,avgExp)



# Percent expression can be obtained through DotPlot, but not Average expression, as this is log z-cored normalized data # avgExpZ = dPlot$data$avg.exp
# Percent expression:
dPlot = DotPlot(photoAll, features=rownames(avgExp))
pctExp = dPlot$data$pct.exp
pctExp = as.data.frame(matrix(pctExp,nrow = length(pctExp)/nI, ncol=nI));

rownames(pctExp) = rownames(avgExp)
colnames(pctExp) = paste("pct", levels(Idents(photoAll)), sep = "")

# Calculate mean percent expression across all ages
meanPct = as.data.frame(rowMeans(pctExp))
colnames(meanPct) = "pctExp"
head(meanPct)
pctExp = cbind(meanPct,pctExp)

photoreceptors_Hoang = cbind(pctExp,avgExp)

# these are duplicated genes with lower vs upper case
dups=toupper(c("arid5b","asph", "aste1", "crip2", "ctbp1", "dab2", "eif1b", "flnb", "frmd7", "galnt10", "grxcr1", "gse1", "hist1h4l", "hspb11", "kif1c", "lamp1", "maf1", "pamr1", "pcdh20", "pde6h", "phlpp2", "psmb10", "ptp4a3", "reep6", "rfesd", "rgs9bp", "rps17", "shank2", "slc16a7", "slc25a10", "slc6a13", "slc9a1", "srbd1", "tatdn3", "tenm3", "tmem178b", "tmem241", "tom1l2", "tp53inp2", "trappc9", "tsc22d3", "ube2o", "zc3h12a", "znf423"))
for(i in 1:length(dups)) {
   rownames(photoreceptors_Hoang)[rownames(photoreceptors_Hoang) == dups[i]] = paste(dups[i],"_ii",sep="")
}
# now everything can be lower cased
rownames(photoreceptors_Hoang) = tolower(rownames(photoreceptors_Hoang))




photoreceptors_Hoang = photoreceptors_Hoang[c("pctExp","pctR","pctUV","pctS","pctM1","pctM3","pctM4","pctL","avgExp","avgR","avgUV","avgS","avgM1","avgM3","avgM4","avgL")]



scangenes=c('nr2e3','tbx2a','tbx2b','arr3a','arr3b','opn1lw1','opn1lw2');
photoreceptors_Hoang[scangenes,]

# still need to manually add "symbol" label to first column
write.csv(photoreceptors_Hoang,"./photoreceptors_Hoang.csv",quote=FALSE)

# -----------------------------------------------------------------------------
# some raw numbers for paper
max(photoAll$nFeature_RNA) # number of genes in each cell
max(photoAll$nCount_RNA) # number of UMIs in each cell
min(photoAll$nCount_RNA)
mean(photoAll$nCount_RNA)
median(photoAll$nCount_RNA)

metadata = photoAll@meta.data
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

# save list of unique genes
# write.csv(x=photoAll$nFeature_RNA, file="Hoang2020_photoAdult_nUniqueGenes.csv", row.names = FALSE)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# LARVAL PHOTORECEPTOR ANALYSIS
# -------------------------------------------------------------------
# Load the previously saved larval cone dataset (flagged by sample name != AdR)
# photo_larval = readRDS(file = "./photo_Larval.rds")


# Identification of highly variable features (feature selection)
photo_larval <- FindVariableFeatures(photo_larval, selection.method = "vst", nfeatures = 200)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(photo_larval), 100)
top10

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(photo_larval)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#run PCA
photo_larval <- RunPCA(photo_larval, features = VariableFeatures(object = photo_larval))

# Elbow plot (var explained)
ElbowPlot(photo, ndims=50) + plotTheme()

# Cluster cells
# nFeatures = 500
#        ndim = 24, res = 1 -> c = 10; 
# one big cluster expresses nr2f6b and gnat2 with a few rods inside and it's basically all cells from 5dpf
# another cluster is nr2e3+
# nFeatures = 100
#        ndim = 33, res = 1 -> c = 9;
# same 2 big clusters but some are misassigned
# nFeatures = 200
#        ndim = 26, res = 0.6 -> c = 9;
# seems like clustering is pretty much the same with all of these.
#        ndim = 12, res = 0.4 -> c = 8; not good; a lot of empty space in tSNE and umap
#        ndim = 12, res = 1 -> c = 12; not good; too many clusters
#        ndim = 45, res = 1 -> c = 12; good elongated plot that shows dev age in diagonal
# CONCLUSIONS: 
   # 3 groups can be parsed apart: late-stage (gnat2+,nr2f6b+), early-stage (nr2e3+) and mid-stage(arr3+,nr2e3-)
   # Results are somewhat similar if just developmental age is used
   # Subtype identity is not the main driver of variability in this dataset

# Settling for nFeatures = 200; ndim = 45, res = 2
upDimLimit=45;
photo_larval <- FindNeighbors(photo_larval, dims = 1:upDimLimit)
photo_larval <- FindClusters(photo_larval, resolution = 2)

DimPlot(photo_larval, reduction = "pca", label = TRUE, pt.size = 1, label.size = 6) + plotTheme()

# UMAP using the same PCA dimensions for prettier visualization
photo_larval <- RunUMAP(photo_larval, dims = 1:upDimLimit)
DimPlot(photo_larval, reduction = "umap", label = TRUE, pt.size = 1, label.size = 6) + plotTheme()


# or use tSNE (clustering can't separate some M cones from L cones)
photo_larval <- RunTSNE(photo_larval, dims = 1:upDimLimit)
DimPlot(photo_larval, reduction = "tsne", label = TRUE, pt.size = 1, label.size = 6) + plotTheme()

# Plot by age
photo_larval2 = SetIdent(photo_larval, value = "Sample")
DimPlot(photo_larval2, reduction = "tsne", label = TRUE, pt.size = 1, label.size = 6) + plotTheme()
# ggsave("larval_TSNEdpf.png", path=exportDir, width = 140, height = 105, units = "mm")

# Plot by age (and getting rid of replicate information for simplicity)
new.cluster.ids <- c("5dpf","4dpf","4dpf","2dpf","2.5dpf","2.5dpf","3dpf","4dpf")
names(new.cluster.ids) <- levels(photo_larval2)
photo_larval2 <- RenameIdents(photo_larval2, new.cluster.ids)
Idents(photo_larval2) <- factor(x = Idents(photo_larval2), levels = c("2dpf","2.5dpf","3dpf","4dpf","5dpf"))
DimPlot(photo_larval2, reduction = "tsne", label = TRUE, pt.size = 1, label.size = 6) + plotTheme()

# separated cluster is foxq2+ from most developmental ages but nr2e3+, so clearly S cones
FeaturePlot(photo_larval, reduction = 'tsne', features = c('foxq2','opnslw2','nr2e3', 'nr2f6b'))

# opsins only detectable in later stages and not clear clustering by it
FeaturePlot(photo_larval, reduction = 'tsne', features = c("rho","opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2'))
# rod genes: larval rods are cluster 12
FeaturePlot(photo_larval, reduction = 'tsne', features = c("rho","nrl","pde6a","nr2e3","crx",'gnat1',"saga","sagb","gucy2f","grk1a")) # rod markers
# cone phototransdcution genes parallel opsins: mostly in later stages
FeaturePlot(photo_larval, reduction = 'tsne', features = c("gnat2",'arr3a','arr3b',"pde6c","pde6h","guca1d","grk7a","crx","neurod1","nr2f6b")) #cone markers
# clusters can be separated based on these genes (nr2e3 is in early cones!!!!)
FeaturePlot(photo_larval, reduction = 'tsne', features = c("nr2e3","nr2f6b","pde6g","pde6h", "neurod1")) #cone markers
# checking some TFs
FeaturePlot(photo_larval, reduction = 'tsne', features = c("otx5","tbx2a","tbx2b","rxrga","rxrgb", "thrb", "foxq2","sall1a"))

# ------------------------------------------------------------------------------------------
# plot layouts for direct export
#OPSINS
p = FeaturePlot(photo_larval, reduction = 'tsne', features = c("rho","opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2'),
                pt.size=1, order=TRUE, combine=FALSE)
lay = rbind(c(1,NA,NA,NA),c(2,3,NA,NA),c(4,5,6,7),c(8,9,NA,NA))
grid.arrange(grobs = p, layout_matrix = lay)
# ps = arrangeGrob(grobs = p, layout_matrix = lay)
# ggsave(ps, file="larval_Opsins.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

#ROD STUFF (known)
p = FeaturePlot(photo_larval, reduction = 'tsne', features = c("rho",'gnat1',"saga",'pde6g',"guca1b","grk1a",'crx','nr2e3','nrl'),pt.size=1, order=TRUE, combine=FALSE)
lay = rbind(c(1,2,NA,NA),c(3,4,5,6),c(7,8,9,NA))
grid.arrange(grobs = p, layout_matrix = lay)
# ps = arrangeGrob(grobs = p, layout_matrix = lay)
# ggsave(ps, file="larval_Rods.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

#BEST'S DEVELOPMENT MARKERS TO NAME CLUSTERS
# early = nr2e3+/nr2f1b-; mid = nr2e3+/nr2f1b+; late = nr2f6b+
p = FeaturePlot(photo_larval, reduction = 'tsne', features = c("gnat2","pde6g","pde6h","crx","neurod1","nr2e3","nr2f6b","nr2f1b"),pt.size=1, order=TRUE, combine=FALSE)
lay = rbind(c(1,2,3),c(4,5,NA),c(6,7,8))
grid.arrange(grobs = p, layout_matrix = lay)
# ps = arrangeGrob(grobs = p, layout_matrix = lay)
# ggsave(ps, file="larval_DevClusters.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")
# ------------------------------------------------------------------------------------------

# Assigning cone subtype identity to clusters
# foxq2+ cluster corresponds to larval S cones in early stage; won't be able to cleanly separate any other subtype, so getting rid of this information
new.cluster.ids <- c("Cle","Cll","Cll","Cll","Cle","Clm","Cle","Clm","Cle","Clm","Clm","Clm","Rll","Clm")
names(new.cluster.ids) <- levels(photo_larval)
photo_larval <- RenameIdents(photo_larval, new.cluster.ids)

Idents(photo_larval) <- factor(x = Idents(photo_larval), levels = c("Cle","Clm","Cll","Rll"))
Photol_colors = c("#DACD9A","#DCC360","#CCA819","#A3A3A3")
levels(Idents(photo_larval))


DimPlot(photo_larval, reduction = "tsne", label = FALSE, pt.size = 1, label.size = 6, cols=Photol_colors) +plotTheme() +NoLegend()
# ggsave("larvalPR_TSNEdev.png", path=exportDir, width = 140, height = 105, units = "mm")
# ggsave("Hoang2020_02_larvalPR_TSNEdev.png", path=exportDir, width = 140, height = 105, units = "mm")
DimPlot(photo_larval, reduction = "umap", label = TRUE, pt.size = 1, cols=Photol_colors) + NoLegend()


# compare reclustering vs. just using age of samples
VlnPlot(photo_larval, features = c("gnat2","pde6g","pde6h","crx","neurod1","nr2e3","nr2f6b","nr2f1b"), cols = Photol_colors)
VlnPlot(photo_larval2, features = c("gnat2","pde6g","pde6h","crx","neurod1","nr2e3","nr2f6b","nr2f1b"))


VlnPlot(photo_larval, features = c("nr2e3", "foxq2", "skor1a","sall1a","lrrfip1a","xbp1","tbx2a","tbx2b"), cols = Photol_colors)
DotPlot(photo_larval, features = c("nr2e3", "foxq2", "skor1a","sall1a","lrrfip1a","xbp1","tbx2a","tbx2b"))
DotPlot(photo_larval2, features = c("nr2e3", "foxq2", "skor1a","sall1a","lrrfip1a","xbp1","tbx2a","tbx2b"))

saveRDS(photo_larval, file = "./photo_LarvalReClustered.rds")
saveRDS(photo_larval2, file = "./photo_LarvalByDPF.rds")
# ------------------------------------------------------------------------------------------
# For attempt to separate subtypes in larval photoreceptors, check out Hoang2020_larvalConeSubtypes.R 
# -------------------------------------------------------------------
# CAN BE RESTARTED HERE 
# rm(list=ls())

# # reload all cells
# pbmc = readRDS("./zfDev_pbmc_v3.rds");
# new.cluster.ids <- c("RPC","MG","R","PRP","C","BC","HC","AClarval","ACgaba","MG","RGC","MG","Clarval","RGClarval","ACgly","BClarval","MGi")
# names(new.cluster.ids) <- levels(pbmc)
# pbmc <- RenameIdents(pbmc, new.cluster.ids)
# Idents(pbmc) <- factor(x = Idents(pbmc), levels = c("RPC","PRP","Clarval","C","R","HC","BClarval","BC","AClarval","ACgaba","ACgly","RGClarval","RGC","MGi","MG"))
# pbmc_colors = c('#DADADA','#dfdac8','#dcc360','#ffd429','#7d7d7d','#FC7715','#ccf2ff','#1B98C3','#3DF591','#3DF5C3','#56F53D','#F53D59','#BB0622','#EA9D81','#A2644E','#7E4835','#613728')
# # reload adult photoreceptors
# photoAll = readRDS(file = "./photo_AdultAll.rds")
# # reload larval photoreceptors
# photo_larval = readRDS(file = "./photo_LarvalReClustered.rds")
# -------------------------------------------------------------------
# check labels
levels(Idents(photo_larval))
levels(Idents(photoAll))
levels(Idents(pbmc))

# Subtype information only reliable is adults and already saved as photoreceptors_Hoang.csv
# For developmental timelines, will not use subtype information other than rods vs. cones
new.cluster.ids <- c("R","C","C","C","C","C","C")
names(new.cluster.ids) <- levels(photoAll)
photoAll <- RenameIdents(photoAll, new.cluster.ids)

# obtain relevant clusters and remerge into single object
rpc = subset(pbmc, idents = c("RPC")); #retinal progenitors
prp = subset(pbmc, idents = c("PRP")); #photoreceptor progenitors

photoDev = merge(x = rpc, y = c(prp,photo_larval,photoAll))
levels(Idents(photoDev))

# ordering identities for plotting and assigning colors
Idents(photoDev) <- factor(x = Idents(photoDev), levels = c("RPC","PRP","Cle","Clm","Cll","C","Rll","R"))
prcolors = c("#DADADA","#dfdac8","#dacd9a","#DCC360",'#CCA819','#FFD429','#A3A3A3','#7D7D7D')
levels(Idents(photoDev))
rm(rpc,prp)

photoDev = ScaleData(photoDev)
photoDev <- FindVariableFeatures(photoDev, selection.method = "vst", nfeatures = 200)

photoDev <- RunPCA(photoDev, features = VariableFeatures(object = photoDev))
ElbowPlot(photoDev, ndims=50)


upDimLimit = 10;
photoDev <- RunUMAP(photoDev, dims = 1:upDimLimit)
DimPlot(photoDev, reduction = "umap", label=TRUE, pt.size = 1, cols=prcolors) + NoLegend()

photoDev <- RunTSNE(photoDev, dims = 1:upDimLimit, check_duplicates = FALSE)
ps = DimPlot(photoDev, reduction = "tsne", label=FALSE, repel=TRUE, pt.size = 0.5, cols=prcolors) +plotTheme() + NoLegend()
ps
# ggsave(ps, file="photoDev_tSNE.png", path=exportDir, width = 140, height = 105, units = "mm")

FeaturePlot(photoDev, reduction = 'tsne', features = c("rho","opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2'),
            pt.size=1, order=TRUE, combine=TRUE)

# "Known" PRP genes
FeaturePlot(photoDev, reduction = 'tsne', features = c("notch1a","ccnd1","cdk1","fgf19","vsx2","hopx","lhx2a","lhx2b"),
            pt.size=1, order=TRUE, combine=TRUE)

# Trying to further map out things
FeaturePlot(photoDev, reduction = 'tsne', features = c("sox2","crx","otx5","neurod1","prdm1a","prdm1b","nr2e3","nr2f1b","nr2f6b"),
            pt.size=1, order=TRUE, combine=TRUE)

saveRDS(photoDev, file = "./photo_Development.rds")

# cell counter
table(Idents(photoDev))

# RPC  PRP   Cle  Clm  Cll    C  Rll   R 
# 3802 2732  407  346  364 1890  36    3097 
# -------------------------------------------------------------------
# CAN BE RESTARTED HERE  DEVELOPMENTLINK
# rm(list=ls())
# photoDev = readRDS(file = "./photo_Development.rds")

# getting % and Mean expression for all genes into csv.
# DoHeatmap and DotPlot use the @scale.data slot for average expression display, which z-scored expression values (for example, as those used in PCA).
# Cells with a value > 0 represent cells with expression above the population mean (a value of 1 would represent cells with expression 1SD away from the population mean). Hope that helps!
# This is very detailed explanation: https://github.com/satijalab/seurat/issues/2798
# Percent expression can be obtained through DotPlot, but not Average expression, as this is log z-cored normalized data # avgExpZ = dPlot$data$avg.exp


# Average expression in non-log scale:
avgExp = AverageExpression(photoDev, slot="counts")
avgExp = avgExp$RNA
nI = length(levels(Idents(photoDev)))
colnames(avgExp) = paste("avg", levels(Idents(photoDev)), sep = "")
# Calculate mean expression across all clusters
meanAvg = as.data.frame(rowMeans(avgExp))
colnames(meanAvg) = "avgExp"
head(meanAvg)
avgExp = cbind(meanAvg,avgExp)

# Percent expression:
dPlot = DotPlot(photoDev, features=rownames(avgExp))
pctExp = dPlot$data$pct.exp
pctExp = as.data.frame(matrix(pctExp,nrow = length(pctExp)/nI, ncol=nI));

rownames(pctExp) = rownames(avgExp)
colnames(pctExp) = paste("pct", levels(Idents(photoDev)), sep = "")

# Calculate mean expression across all clusters
meanPct = as.data.frame(rowMeans(pctExp))
colnames(meanPct) = "pctExp"
head(meanPct)
pctExp = cbind(meanPct,pctExp)

photoDev_Hoang = cbind(pctExp,avgExp)
head(photoDev_Hoang,10)


# these are duplicated genes with lower vs upper case
dups=toupper(c("arid5b","asph", "aste1", "crip2", "ctbp1", "dab2", "eif1b", "flnb", "frmd7", "galnt10", "grxcr1", "gse1", "hist1h4l", "hspb11", "kif1c", "lamp1", "maf1", "pamr1", "pcdh20", "pde6h", "phlpp2", "psmb10", "ptp4a3", "reep6", "rfesd", "rgs9bp", "rps17", "shank2", "slc16a7", "slc25a10", "slc6a13", "slc9a1", "srbd1", "tatdn3", "tenm3", "tmem178b", "tmem241", "tom1l2", "tp53inp2", "trappc9", "tsc22d3", "ube2o", "zc3h12a", "znf423"))
for(i in 1:length(dups)) {
   rownames(photoDev_Hoang)[rownames(photoDev_Hoang) == dups[i]] = paste(dups[i],"_ii",sep="")
}
# now everything can be lower cased
rownames(photoDev_Hoang) = tolower(rownames(photoDev_Hoang))
# also replacing (1 of many) with .1
rownames(photoDev_Hoang) = gsub("(1 of many)",".1",rownames(photoDev_Hoang))


photoDev_Hoang = photoDev_Hoang[c("avgExp","avgRPC","avgPRP","avgCle","avgClm","avgCll","avgC","avgRll","avgR",
                                  "pctExp","pctRPC","pctPRP","pctCle","pctClm","pctCll","pctC","pctRll","pctR")]

# Huge problem for dev timeline visualization is that Rll (larval rods) are very contaminated with cones genes and skew all plots to highlight these cells
scangenes=c('gnat1','gnat2','nr2e3','tbx2a','tbx2b','arr3a','arr3b','opn1lw1','opn1lw2');
photoDev_Hoang[scangenes,]

scangenes=c('nr2e3','skor1a');
photoDev_Hoang[scangenes,]

# removing larval rods
photoDev_Hoang = photoDev_Hoang[c("avgExp","avgRPC","avgPRP","avgCle","avgClm","avgCll","avgC","avgR",
                                  "pctExp","pctRPC","pctPRP","pctCle","pctClm","pctCll","pctC","pctR")]
# sort by avgExpression
photoDev_Hoang=photoDev_Hoang[order(-photoDev_Hoang$avgExp),]
head(photoDev_Hoang,100)

scangenes=c('gnat1','gnat2','nr2e3','tbx2a','tbx2b','arr3a','arr3b','opn1lw1','opn1lw2');
photoDev_Hoang[scangenes,]

scangenes=c('nr2e3','skor1a');
photoDev_Hoang[scangenes,]


write.csv(photoDev_Hoang,"./Hoang2020_10x_photoDev.csv",quote=FALSE)
# remember to manually add "symbol" to first column name!
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# -------------------------------------------------------------------