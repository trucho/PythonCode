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
setwd("/Users/angueyraaristjm/Documents/LiMolec/otherRNAseq/zfRet_Ogawa2021/")
directory = "/Users/angueyraaristjm/Documents/LiMolec/otherRNAseq/zfRet_Ogawa2021/"
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
photo = readRDS("~/Documents/LiMolec/otherRNAseq/zfRet_Ogawa2021/GSM5351368_v432_dradult_Final.rds");

# ------------------------------------------------------------------------------------------
new.cluster.ids <- c("M","L","R","S","M4/L1","onBC","offBC","UV")
names(new.cluster.ids) <- levels(photo)
photo <- RenameIdents(photo, new.cluster.ids)
Idents(photo) <- factor(x = Idents(photo), levels = c("R","UV","S","M","L","M4/L1","onBC","offBC"))

# cell counter
table(Idents(photo))

# R     UV     S     M     L     M4/L1  onBC    offBC
# 332   105   301   506   490   195     136     121


# ---------------------------------------------------------------
DimPlot(photo, reduction = "pca", label =  TRUE, pt.size = 1, label.size = 6) + eelTheme()
ggsave("01_originalPCA.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

DimPlot(photo, reduction = "umap", label = TRUE, pt.size = 1, label.size = 6, repel = TRUE) + eelTheme()
ggsave("01_originalUMAP.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

photo <- RunTSNE(photo, dims = 1:33)
DimPlot(photo, reduction = "tsne", label = TRUE, pt.size = 1, label.size = 6, repel = TRUE) + eelTheme()
ggsave("01_originalTSNE.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")


FeaturePlot(photo, reduction = 'umap', features = c("crx","otx5", "tbx2a",'tbx2b','foxq2'), pt.size = 1, order=TRUE)

FeaturePlot(photo, reduction = 'umap', features = c("nr2e3","rorb", "samd7",'ar','skia','neurod1'), pt.size = 2)

FeaturePlot(photo, reduction = 'umap', features = c("nr2e3","nr1d1", "nr2f1a",'nr2f1b','nr2f6a','nr2f6b'), pt.size = 1, order=TRUE)

#L/M TFs
FeaturePlot(photo, reduction = 'umap', features = c("nfe2l1b","ahr1b", "churc1",'eloca','taf13','eaf1','stat1a','six7','neurod1'), pt.size = 0.5, order=FALSE)
#L TFs
FeaturePlot(photo, reduction = 'umap', features = c("samd7","thrb", "rxrga",'smad5','prdm1a','prdm1b','rx2','klf13'), pt.size = 1, order=TRUE)
#M TFs
FeaturePlot(photo, reduction = 'umap', features = c("cux2b", "epas1b",'thraa','lrrfip1a','pbx1a'), pt.size = 1, order=TRUE)
#S/Uv TFs
FeaturePlot(photo, reduction = 'umap', features = c("skor1a","tbx2b", "rx1",'crema','etv5b','tefb','foxq2','erfl3','tefa','tbx2a','xbp1'), pt.size = 0.5, order=FALSE)

FeaturePlot(photo, reduction = 'umap', features = c("crx","otx5", "tbx2a",'tbx2b'), pt.size = 2)
VlnPlot(photo, features = c("crx","otx5", "tbx2a",'tbx2b'))
# ------------------------------------------------------------------------------------------


VlnPlot(photo, features = c("vamp1","vamp2","vamp3","vamp4","vamp5","vamp8")) #log norm values
VlnPlot(photo, features = c("gnat2", "rho",'elovl4b','gnat1','si:busm1-57f23.1')) #log norm values
VlnPlot(photo, features = c("gnat2", "rho",'elovl4b','gnat1','si:busm1-57f23.1'), slot = "counts", log = TRUE) #raw counts
RidgePlot(photo, features = c("gnat2", "rho",'elovl4b','gnat1','si:busm1-57f23.1'))
DotPlot(photo, features = c("gnat2", "rho",'elovl4b','gnat1','si:busm1-57f23.1'))

VlnPlot(photo, features = c("sema7a","efna1b","ntng2b","syt1a","syt1b","syt5a","syt5b","tbx2a","tbx2b","eml1")) #log norm values
DotPlot(photo, features = c("sema7a","efna1b","ntng2b","syt1a","syt1b","syt5a","syt5b","tbx2a","tbx2b","eml1"))
DotPlot(photo, features = c("vamp1","vamp2","vamp3","vamp4","vamp5","vamp8")) 
FeaturePlot(photo, reduction = 'tsne', features = c("vamp1","vamp2","vamp3","vamp4","vamp5","vamp8")) 
VlnPlot(photo, features = c("opn1sw1", "opn1sw2",'opn1mw1','gnat2','opn1lw2','si:busm1-57f23.1','rho','saga')) #log norm values
# ------------------------------------------------------------------------------------------
# plot layouts for direct export
#Subtype markers
p = FeaturePlot(photo, reduction = 'tsne', features = c("rho","nr2e3","gnat1","gnat2","arr3b","arr3a","thrb","opn1sw1", "opn1sw2",'opn1mw1','opn1lw1'),
                pt.size=1, order=FALSE, combine=FALSE)
lay = rbind(c(1,2,3,NA),c(4,5,6,7),c(8,9,10,11))
grid.arrange(grobs = p, layout_matrix = lay)
ps = arrangeGrob(grobs = p, layout_matrix = lay)
ggsave(ps, file="02_PRMarkers.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

#OPSINS
p = FeaturePlot(photo, reduction = 'tsne', features = c("rho","opn1sw1", "opn1sw2",'opn1mw1','opn1mw2','opn1mw3','opn1mw4','opn1lw1','opn1lw2'),
                pt.size=1, order=TRUE, combine=FALSE)
lay = rbind(c(1,NA,NA,NA),c(2,3,NA,NA),c(4,5,6,7),c(8,9,NA,NA))
grid.arrange(grobs = p, layout_matrix = lay)
ps = arrangeGrob(grobs = p, layout_matrix = lay)
ggsave(ps, file="02_Opsins.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

#ROD STUFF (known)
p = FeaturePlot(photo, reduction = 'tsne', features = c("rho",'gnat1',"gnat2","saga",'pde6ga',"guca1b","grk1a",'crx','nr2e3','nrl'),pt.size=1, order=FALSE, combine=FALSE)
lay = rbind(c(1,2,3,NA),c(4,5,6,7),c(8,9,10,NA))
grid.arrange(grobs = p, layout_matrix = lay)
ps = arrangeGrob(grobs = p, layout_matrix = lay)
ggsave(ps, file="02_RodMarkers.png", path=exportDir, width = 140*2, height = 105*2, units = "mm")

# -----------------------------------------------------------------------------
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
avgExp = as.data.frame(avgExp)
# Calculate mean expression across all ages
avgExp$avgMean = rowMeans(avgExp)
avgExp = avgExp[,c(ncol(avgExp),1:(ncol(avgExp)-1))]


# Percent expression can be obtained through DotPlot, but not Average expression, as this is log z-cored normalized data # avgExpZ = dPlot$data$avg.exp
# dPlot = DotPlot(photo, features=rownames(avgExp)) # not working
# dPlot <- FetchData(photo, rownnames(avgExp)) # dot plot is not working. Trying this instead; does not work either
# ------------------------------------------------------------------------------------------------
# breaking down into pieces to try to find row that creates error
length(rownames(avgExp))
# let's break it in half: nope!; a quarter? works
dPlot = DotPlot(photo, features=rownames(avgExp[c(1:5000),]))
# second quarter: fails; keep splitting.
dPlot21 = DotPlot(photo, features=rownames(avgExp[c(5001:6000),])) #pass
dPlot22 = DotPlot(photo, features=rownames(avgExp[c(6001:7000),])) #pass
dPlot23 = DotPlot(photo, features=rownames(avgExp[c(7001:8000),])) #pass
dPlot24 = DotPlot(photo, features=rownames(avgExp[c(8001:9000),])) #pass
dPlot25 = DotPlot(photo, features=rownames(avgExp[c(9001:10000),])) #pass
# second half? works
dPlot3 = DotPlot(photo, features=rownames(avgExp[c(10001:22129),]))
# rebind everything
# ------------------------------------------------------------------------------------------------
# testing how to give rownames appropriately
temp = dPlot
tempnames = rownames(temp$data[seq(1,dim(temp$data)[1]/8,1),])
temppct = temp$data$pct.exp
tempout1 = as.data.frame(matrix(temppct, nrow = length(temppct)/8, ncol = 8, dimnames = list(c(tempnames),c(levels(Idents(photo))))))
temp = dPlot21
tempnames = rownames(temp$data[seq(1,dim(temp$data)[1]/8,1),])
temppct = temp$data$pct.exp
tempout2 = as.data.frame(matrix(temppct, nrow = length(temppct)/8, ncol = 8, dimnames = list(c(tempnames),c(levels(Idents(photo))))))
temp = dPlot22
tempnames = rownames(temp$data[seq(1,dim(temp$data)[1]/8,1),])
temppct = temp$data$pct.exp
tempout3 = as.data.frame(matrix(temppct, nrow = length(temppct)/8, ncol = 8, dimnames = list(c(tempnames),c(levels(Idents(photo))))))
temp = dPlot23
tempnames = rownames(temp$data[seq(1,dim(temp$data)[1]/8,1),])
temppct = temp$data$pct.exp
tempout4 = as.data.frame(matrix(temppct, nrow = length(temppct)/8, ncol = 8, dimnames = list(c(tempnames),c(levels(Idents(photo))))))
temp = dPlot24
tempnames = rownames(temp$data[seq(1,dim(temp$data)[1]/8,1),])
temppct = temp$data$pct.exp
tempout5 = as.data.frame(matrix(temppct, nrow = length(temppct)/8, ncol = 8, dimnames = list(c(tempnames),c(levels(Idents(photo))))))
temp = dPlot25
tempnames = rownames(temp$data[seq(1,dim(temp$data)[1]/8,1),])
temppct = temp$data$pct.exp
tempout6 = as.data.frame(matrix(temppct, nrow = length(temppct)/8, ncol = 8, dimnames = list(c(tempnames),c(levels(Idents(photo))))))
temp = dPlot3
tempnames = rownames(temp$data[seq(1,dim(temp$data)[1]/8,1),])
temppct = temp$data$pct.exp
tempout7 = as.data.frame(matrix(temppct, nrow = length(temppct)/8, ncol = 8, dimnames = list(c(tempnames),c(levels(Idents(photo))))))

temp = dPlot
tempnames = rownames(temp$data[seq(1,dim(temp$data)[1]/8,1),])
temppct = temp$data$pct.exp
tempout1 = as.data.frame(matrix(temppct, nrow = length(temppct)/8, ncol = 8, dimnames = list(c(tempnames),c(levels(Idents(photo))))))
# ------------------------------------------------------------------------------------------------
tempout = rbind(tempout1,tempout2)
tempout = rbind(tempout,tempout3)
tempout = rbind(tempout,tempout4)
tempout = rbind(tempout,tempout5)
tempout = rbind(tempout,tempout6)
tempout = rbind(tempout,tempout7)
# fixing names replaced by DotPlot
rownames(tempout) = gsub(pattern = 'rna_', replacement = '',x = rownames(tempout))
pctExp = tempout
colnames(pctExp) = paste("pct", colnames(pctExp), sep = "")
pctExp$pctMean = rowMeans(pctExp)
pctExp = pctExp[,c(ncol(pctExp),1:(ncol(pctExp)-1))]
# merging data frames
photoData = merge(avgExp,pctExp, by="row.names", all=TRUE)
rownames(photoData) = photoData[,1]
photoData = photoData[-c(1)]
# ------------------------------------------------------------------------------------------------
# quick check
scangenes=c('rho','nr2e3','opn1lw1','opn1lw2','thrb','skor1a','tbx2a','tbx2b','arr3a','arr3b');
photoData[scangenes,]


write.csv(photoData,"./photoreceptors_Ogawa.csv",quote=FALSE)
