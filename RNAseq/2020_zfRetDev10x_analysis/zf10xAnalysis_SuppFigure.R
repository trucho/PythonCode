
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
directory = "/Users/angueyraaristjm/Documents/LiLab/Presentations/revealjs/resources/20210312_LLL/"
exportDir = paste(directory,"eelAnalysis",sep="")
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

# 

# Load data -------------------------------------------------------------------
# Final reclustering inclduing photoProgenitors, e/m/l photoLarval and photoAdult
photoDev = readRDS(file = "./photoreceptors_Development.rds")
levels(Idents(photoDev))
prcolors = c("#dfdac8","#dacd9a","#dcc360","#cca819","#a3a3a3","#7d7d7d","#B540B7","#4669F2","#04CD22","#CC2C2A")

prcolors2 = c("#7d7d7d","#CC2C2A","#cca819","#dcc360","#dacd9a","#dfdac8","#a3a3a3","#04CD22","#4669F2","#B540B7")
# UMAP plot
ps = DimPlot(photoDev, reduction = "umap", label=FALSE, repel = TRUE, pt.size = 1, cols=prcolors2, label.size = 10, 
             order = c("UV","S","M","lslR","PRP","eslPR","mslPR","lslPR","L","R")) + 
   eelTheme() + NoLegend() + xlim(-12,12) + ylim(-10,10)
ps
ggsave(ps, file="01_UMAPClusters.png", path="./zfConeRNAseqFigure/UMAP/", width = 140*2, height = 105*2, units = "mm")
# tSNE plot
ps = DimPlot(photoDev, reduction = "tsne", label=FALSE, repel = TRUE, pt.size = 1, cols=prcolors2, label.size = 10, 
             order = c("UV","S","M","lslR","PRP","eslPR","mslPR","lslPR","L","R")) + 
   eelTheme() + NoLegend() + xlim(-50,50) + ylim(-60,60)
ps
ggsave(ps, file="01_TSNEClusters.png", path="./zfConeRNAseqFigure/TSNE/", width = 140*2, height = 105*2, units = "mm")

# --------------------------------------------------------------------------------
# Plot by age (and getting rid of replicate information for clarity)
photoDev2 = SetIdent(photoDev, value = "Sample")
new.cluster.ids <- c("48hpf","60hpf","60hpf","72hpf","96hpf","96hpf","96hpf","120hpf","Adult","Adult","Adult","Adult","Adult")
names(new.cluster.ids) <- levels(photoDev2)
photoDev2 <- RenameIdents(photoDev2, new.cluster.ids)
ps = DimPlot(photoDev2, reduction = "tsne", label=TRUE, repel = TRUE, pt.size = 1,  label.size = 10, 
             order = c("48hpf")) + 
   eelTheme() + NoLegend() + xlim(-50,50) + ylim(-60,60)
ps
ggsave(ps, file="01_TSNEClusters_AgeLabels.png", path="./zfConeRNAseqFigure/TSNE/", width = 140*2, height = 105*2, units = "mm")
ps = DimPlot(photoDev2, reduction = "tsne", label=FALSE, repel = TRUE, pt.size = 1,  label.size = 10, 
             order = c("48hpf")) + 
   eelTheme() + NoLegend() + xlim(-50,50) + ylim(-60,60)
ps
ggsave(ps, file="01_TSNEClusters_Age.png", path="./zfConeRNAseqFigure/TSNE/", width = 140*2, height = 105*2, units = "mm")

# --------------------------------------------------------------------------------


# --------------------------------------------------------------------------------
# opsins
genelist = c("rho","opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2')
for(i in 1:length(genelist)) {
   genename = genelist[i]
   ps = FeaturePlot(photoDev, reduction = 'tsne', features = c(genename),pt.size=1, order=TRUE, combine=TRUE) + eelTheme() + NoLegend() + xlim(-50,50) + ylim(-60,60)
   ggsave(ps, file=paste("02_TSNE",genename,".png",sep=''), path="./zfConeRNAseqFigure/TSNE/", width = 140*2, height = 105*2, units = "mm")
}
FeaturePlot(photoDev, reduction = 'tsne', features = genelist,pt.size=1, order=TRUE, combine=TRUE, ncol = 3)
# --------------------------------------------------------------------------------
# PRP markers
genelist = c("notch1a","ccnd1","cdk1","neurod1","hes2.2","sox2")
for(i in 1:length(genelist)) {
   genename = genelist[i]
   ps = FeaturePlot(photoDev, reduction = 'tsne', features = c(genename),pt.size=1, order=TRUE, combine=TRUE) + eelTheme() + NoLegend() + xlim(-50,50) + ylim(-60,60)
   ggsave(ps, file=paste("02_TSNE",genename,".png",sep=''), path="./zfConeRNAseqFigure/TSNE/", width = 140*2, height = 105*2, units = "mm")
}
FeaturePlot(photoDev, reduction = 'tsne', features = genelist,pt.size=1, order=TRUE, combine=TRUE, ncol = 3)

# --------------------------------------------------------------------------------
# larval markers
genelist = c("crx","otx2","otx5", "prdm1a",'nr2e3','nr2f1a','nr2f1b','nr2f6a','nr2f6b','syt5a','syt5b')
for(i in 1:length(genelist)) {
   genename = genelist[i]
   ps = FeaturePlot(photoDev, reduction = 'tsne', features = c(genename),pt.size=1, order=TRUE, combine=TRUE) + eelTheme() + NoLegend() + xlim(-50,50) + ylim(-60,60)
   ggsave(ps, file=paste("02_TSNE",genename,".png",sep=''), path="./zfConeRNAseqFigure/TSNE/", width = 140*2, height = 105*2, units = "mm")
}
FeaturePlot(photoDev, reduction = 'tsne', features = genelist,pt.size=1, order=TRUE, combine=TRUE, ncol = 3)

# --------------------------------------------------------------------------------
# subtype TFs
genelist = c('nrl',"thrb","roraa", "rorab",'rxrga','rxrgb','tbx2a','tbx2b','foxq2','skor1a','xbp1','gdf6a')
for(i in 1:length(genelist)) {
   genename = genelist[i]
   ps = FeaturePlot(photoDev, reduction = 'tsne', features = c(genename),pt.size=1, order=TRUE, combine=TRUE) + eelTheme() + NoLegend() + xlim(-50,50) + ylim(-60,60)
   ggsave(ps, file=paste("02_TSNE",genename,".png",sep=''), path="./zfConeRNAseqFigure/TSNE/", width = 140*2, height = 105*2, units = "mm")
}
FeaturePlot(photoDev, reduction = 'tsne', features = genelist,pt.size=1, order=TRUE, combine=TRUE, ncol = 3)
# --------------------------------------------------------------------------------
# subtype other genes: TRANSDUCIN & ARRESTINS
genelist = c('gnat1',"gnat2",'gnb1a','gnb1b','gnb3a','gnb3b','gngt1','gngt2a','gngt2b',"saga", "sagb",'arr3a','arr3b')
for(i in 1:length(genelist)) {
   genename = genelist[i]
   ps = FeaturePlot(photoDev, reduction = 'tsne', features = c(genename),pt.size=1, order=TRUE, combine=TRUE) + eelTheme() + NoLegend() + xlim(-50,50) + ylim(-60,60)
   ggsave(ps, file=paste("02_TSNE",genename,".png",sep=''), path="./zfConeRNAseqFigure/TSNE/", width = 140*2, height = 105*2, units = "mm")
}
FeaturePlot(photoDev, reduction = 'tsne', features = genelist,pt.size=1, order=TRUE, combine=TRUE, ncol = 3)
# --------------------------------------------------------------------------------
# subtype other genes: cGMP-activated channel
genelist = c('cnga1','cnga1(.1)','cnga3a','cnga3b','cngb1','cngb3.1','cngb3.2')
for(i in 1:length(genelist)) {
   genename = genelist[i]
   ps = FeaturePlot(photoDev, reduction = 'tsne', features = c(genename),pt.size=1, order=TRUE, combine=TRUE) + eelTheme() + NoLegend() + xlim(-50,50) + ylim(-60,60)
   ggsave(ps, file=paste("02_TSNE",genename,".png",sep=''), path="./zfConeRNAseqFigure/TSNE/", width = 140*2, height = 105*2, units = "mm")
}
FeaturePlot(photoDev, reduction = 'tsne', features = genelist,pt.size=1, order=TRUE, combine=TRUE, ncol = 3)
# --------------------------------------------------------------------------------
# contamination controls
genelist = c('grm6a',"grm6b",'vsx2','prkca','syt2a','pou4f1','pou4f2','pou4f3','rbpms',"gfap", "rpe65a")
for(i in 1:length(genelist)) {
   genename = genelist[i]
   ps = FeaturePlot(photoDev, reduction = 'tsne', features = c(genename),pt.size=1, order=TRUE, combine=TRUE) + eelTheme() + NoLegend() + xlim(-50,50) + ylim(-60,60)
   ggsave(ps, file=paste("02_TSNE",genename,".png",sep=''), path="./zfConeRNAseqFigure/TSNE/", width = 140*2, height = 105*2, units = "mm")
}
FeaturePlot(photoDev, reduction = 'tsne', features = genelist,pt.size=1, order=TRUE, combine=TRUE, ncol = 3)
# --------------------------------------------------------------------------------
# CAMs 01
genelist = c('dscama','dscamb','dag1','lrit1a','lrit1b','lrit2','lrit3a','cadm1a','cadm1b','cadm2a','cadm2b','cadm3','cadm4')
for(i in 1:length(genelist)) {
   genename = genelist[i]
   ps = FeaturePlot(photoDev, reduction = 'tsne', features = c(genename),pt.size=1, order=TRUE, combine=TRUE) + eelTheme() + NoLegend() + xlim(-50,50) + ylim(-60,60)
   ggsave(ps, file=paste("02_TSNE",genename,".png",sep=''), path="./zfConeRNAseqFigure/TSNE/", width = 140*2, height = 105*2, units = "mm")
}
FeaturePlot(photoDev, reduction = 'tsne', features = genelist,pt.size=1, order=TRUE, combine=TRUE, ncol = 3)
# --------------------------------------------------------------------------------
# CAMs 02
genelist = c('sema3fa','sema3fb','sema4ab','sema6d','sema7a','efna1b')
for(i in 1:length(genelist)) {
   genename = genelist[i]
   ps = FeaturePlot(photoDev, reduction = 'tsne', features = c(genename),pt.size=1, order=TRUE, combine=TRUE) + eelTheme() + NoLegend() + xlim(-50,50) + ylim(-60,60)
   ggsave(ps, file=paste("02_TSNE",genename,".png",sep=''), path="./zfConeRNAseqFigure/TSNE/", width = 140*2, height = 105*2, units = "mm")
}
FeaturePlot(photoDev, reduction = 'tsne', features = genelist,pt.size=1, order=TRUE, combine=TRUE, ncol = 3)
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# Violin plots to make things simpler? or just heatmap?
photoDev3 = photoDev
new.cluster.ids <- c("Progenitor","Larval(early)","Larval(mid)","Larval(late)","Larval(late)","Adult","Adult","Adult","Adult","Adult")
names(new.cluster.ids) <- levels(photoDev3)
photoDev3 <- RenameIdents(photoDev3, new.cluster.ids)
prcolors3 = c("#dfdac8","#dacd9a","#dcc360","#cca819","#ffd429")

# new.cluster.ids <- c("Progenitor","Larval(early)","Larval(mid)","Larval(late)","Larval(late)","Adult(rods)","Adult(cones)","Adult(cones)","Adult(cones)","Adult(cones)")
# prcolors3 = c("#dfdac8","#dacd9a","#dcc360","#cca819",'#7d7d7d',"#ffd429")
# --------------------------------------------------------------------------------
# Dev Markers
genelist = c("otx2", "nr2f6a",'nr2e3','prdm1a','nr2f1b','nr2f6b')
for(i in 1:length(genelist)) {
   genename = genelist[i]
   ps = VlnPlot(photoDev3, features = c(genename),cols=prcolors3,pt.size=3) + eelTheme() + NoLegend()
   ggsave(ps, file=paste("03_DEV",genename,".png",sep=''), path="./zfConeRNAseqFigure/", width = 140*2.5, height = 105*2.5, units = "mm")
}
VlnPlot(photoDev3, features = genelist, cols=prcolors3,log=FALSE)
# DotPlot(photoDev3, features = c("otx2", "nr2f6a",'nr2e3','prdm1a','nr2f1b','nr2f6b'))
# --------------------------------------------------------------------------------

genelist = c("otx2", "nr2f6a",'nr2e3','prdm1a','nr2f1b','nr2f6b')
genelist = c("gnat1", "gnat2",'arr3a','arr3b')
VlnPlot(photoDev3, features = genelist, cols=prcolors3,log=FALSE)


VlnPlot(photoDev3, features = "neurod1", cols=prcolors3,log=FALSE) + geom_boxplot(width=0.1, fill="white")

# RidgePlot(photoDev3, features = genelist, cols=prcolors3)
# DoHeatmap(photoDev3, features = genelist)
# VlnPlot(photoDev, features = c("saga", "nr2e3",'foxq2','opn1sw2','opn1lw2','thrb','rxrga'), cols=prcolors) + eelTheme() + NoLegend() #Subtype markers
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# subtypes TFs
genelist = c('nr2e3','samd7','samd11','nrl','six6b','six7','thrb','ahr1b','neurod1','rxrga','prdm1a','prdm1b','thraa','tbx2a','tbx2b','skor1a','rx1','tefa','tefb','foxq2','xbp1')
DotPlot(photoDev, features = genelist)

VlnPlot(photoDev, features = genelist, cols=prcolors,log=FALSE)

genelist = c('samd7','samd11','thrb','rxrga','thraa','tbx2a','tbx2b','skor1a','rx1','tefa','tefb','foxq2','xbp1','lrrfip1a')

genelist = c('samd11','thrb','opn1lw2','tbx2a','tbx2b','skor1a','foxq2','lrrfip1a','opn1mw2')
genelist = c('opn1mw4','vax1','vax2')
FeaturePlot(photoDev, reduction = 'tsne', features = genelist,pt.size=1, order=TRUE, combine=TRUE, ncol = 3)






# test plot --------------------------------------------------------------------------------
genelist = c('efn1a','efna1b','ntng2b')
genelist = c('otx2',"neurod1","rho","nrl","guca1b","rom1a","rom1b","gnat1","nr2e3") #Hoang's ROD markers for zf all fail
FeaturePlot(photoDev, reduction = 'tsne', features = genelist,pt.size=1, order=TRUE, combine=TRUE)

genelist = c("gnat1","gnat2") #Hoang's ROD markers zf
FeaturePlot(photoDev, reduction = 'tsne', features = genelist,pt.size=1, order=FALSE, combine=TRUE)

# --------------------------------------------------------------------------------
# "Known" PRP genes
FeaturePlot(photoDev, reduction = 'tsne', features = c("notch1a","ccnd1","cdk1","fgf19","vsx2","hopx","lhx2a","lhx2b","rlbp1a","rlbp1b"),
            pt.size=1, order=FALSE, combine=TRUE)
# Me trying to map out things
FeaturePlot(photoDev, reduction = 'tsne', features = c("crx","nr2e3","neurod1","otx5","sox2","prdm1a","prdm1b","syt5a","syt5b","foxq2","skor1a"),
            pt.size=1, order=FALSE, combine=TRUE)

FeaturePlot(photoDev, reduction = 'umap', features = c("opn1sw2","foxq2","skor1a","thrb",'xbp1'),
            pt.size=1, order=TRUE, combine=TRUE)







# Also checking "known rod markers"
FeaturePlot(photoDev, reduction = 'tsne', features = c("prkca","sebox","guca1b","nrl","rho","rom1a","rom1b","nr2e3","saga"),
            pt.size=1, order=FALSE, combine=TRUE)


DotPlot(photoDev,features = c("ntf3","spock3","hs3st3l","sema3fb","si:busm1-57f23.1","eya2","c1galt1a","itga1","slc32a1","dok6"))



saveRDS(photoDev, file = "./photoreceptors_Development.rds")

# -------------------------------------------------------------------
# CAN BE RESTARTED HERE  DEVELOPMENTLINK
