# Package installation -------------------------------------------------------------------
if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install()

BiocManager::install(c("DESeq2"))
BiocManager::install(c("GenomicFeatures"))
BiocManager::install(c("AnnotationDbi"))
BiocManager::install(c("org.Dr.eg.db"))
BiocManager::install(c("apeglm"))
BiocManager::install(c("pheatmap"))
BiocManager::install("ReportingTools")
BiocManager::install("pcaExplorer")
BiocManager::install("purr")
BiocManager::install("dplyr")



# Library loading -------------------------------------------------------------------
library("DESeq2")
library("apeglm")
library("ggplot2")
library("org.Dr.eg.db")
library("dplyr")
library("ReportingTools")
library("pcaExplorer")
library("EnhancedVolcano")

# Use this if loading of Dr.db fails
# options(connectionObserver = NULL)

# L-cone specific -------------------------------------------------------------------
rm(list=ls())

setwd("/Users/angueyraaristjm/Documents/LiMolec/zfRNAseq/20190827/20190827_DESeq2/")
directory <- "/Users/angueyraaristjm/Documents/LiMolec/zfRNAseq/20190827/20190827_DESeq2/"
getwd()

# Removing merging with count data since this is not fpkm
# Flipping fold-chage to make enriched genes positive
# NEED TO SAVE VST TRANSFORMED DATA FOR HEATMAPS FROM RvC: vsd <- vst(dds)

# Diagnostic plots: 
# MAplot: 
# plotMA( res, ylim = c(-1, 1) )
# Ratio of small p values for groups of genes binned by mean normalized count: 
# hist( res$pvalue, breaks=20, col="grey" )


# remove subtrancript info from gCount.csv (sed 's/\.[0-9]//' gCount_ensdart.csv > gCount.csv)
# Read Data and run DESeq2 -------------------------------------------------------------------
countData <- as.matrix(read.csv(paste0(directory,"gCount.csv"), row.names="gene_id"))
#prepDE.py does not observe the order of the provided gtf_list
# countData <- countData[,c(13,14,15,16,17,18,25,26,27,28,29,30,19,20,21,22,23,24,25,7,8,9,10,11,12,1,2,3,4,5,6)]
colData <- read.csv(paste0(directory,"PHENO_DATA_L.csv"), sep="\t", row.names=1)

# use colData to reorganize countData
countDataSubset <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countDataSubset))
# save total number of mapped reads
# write.csv(colSums(countData), file = "nMappedReads.csv")

dds <- DESeqDataSetFromMatrix(countData = countDataSubset, colData = colData, design = ~subtype)
dds <- DESeq(dds)
res <- results(dds)
summary(res)
sprintf('n(DEGenes) = %g (p<0.05) ', sum(res$padj < 0.05, na.rm=TRUE))

# vst transformation ------------------------------------------------------------------------------------------------
vsd <- vst(dds)

# Log fold change shrinkage for visualization and ranking -------------------------------------------------------------------
resLFC <- lfcShrink(dds, coef="subtype_USM_vs_L", type="apeglm")
head(resLFC)


# Include Genename (descriptive) -------------------------------------------------
# Run this only if things have changed:
# resdata <- res # # # not merging with normalized counts
# resdata <- resdata[order(resdata$padj),]
# # For some reason this breaks the csv writing and symbol column is replaced by numbers
# # resdata$symbol <- tolower(resdata$symbol)
# genenames <- mapIds(org.Dr.eg.db, keys=rownames(resdata), column=c("GENENAME"), keytype="SYMBOL", multivals='first')
# # write.csv(genenames, file = "genenames.csv", col.names=c("symbol","genename"))
# write.csv(genenames, file = "genenames_L.csv", col.names=NA)


genenames <- read.csv("genenames_L.csv", sep=",")
colnames(genenames) <- c("symbol","genename")
genenames$genename <- gsub(",","",genenames$genename)
head(genenames)

# Save Results ------------------------------------------------------------
# results 
resdata <- res # # # not merging with normalized counts: resdata = merge(as.data.frame(resLFC), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
resdata <- resdata[order(resdata$padj),]
resdata$log2FoldChange = -resdata$log2FoldChange #inverting fold-chage
head(resdata)

# save raw results for python plotting
if (all(genenames$symbol == resdata$symbol)) {
  print("data frames DO match")
  resdata$genename = genenames$genename
  resdata <- resdata[c(ncol(resdata),1:ncol(resdata)-1)] #not sure why it's adding symbol again
  head(resdata)
  write.csv(resdata, file = "00_bySubtype/LvsUSM_rawFC.csv", row.names=TRUE, quote=FALSE)
  print("saved data frames to csv files")
} else {
  print("data frames do NOT match")
}


resdataLFC <- resLFC # # # not merging with normalized counts: resdata = merge(as.data.frame(resLFC), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
resdataLFC <- resdataLFC[order(resdataLFC$padj),]
# resdata$log2FoldChange = -resdata$log2FoldChange #inverting fold-chage
head(resdataLFC)

# save lfc results for python heatmaps
if (all(genenames$symbol == resdataLFC$symbol)) {
  print("data frames DO match")
  resdataLFC$genename = genenames$genename
  resdataLFC <- resdataLFC[c(ncol(resdataLFC),1:ncol(resdataLFC)-1)] #not sure why it's adding symbol again
  head(resdataLFC)
  write.csv(resdataLFC, file = "00_bySubtype/LvsUSM_shrinkFC.csv", row.names=TRUE, quote=FALSE)
  print("saved data frames to csv files")
} else {
  print("data frames do NOT match")
}


# resdata[rownames(resdata)=='opn1sw1',]
# -------------------------------------------------------------------------
# volcano plot
resDF = resdata
# resDF = resdataLFC
EnhancedVolcano(resDF,
                lab = rownames(resDF),
                x = 'log2FoldChange',y = 'pvalue',
                xlim = c(-30, 30),
                title = 'L versus UV/S/M',
                pCutoff = 0.05,
                FCcutoff = 0.6,col=c('grey75', 'grey50', 'grey25', 'red','blue'),
                colAlpha = 1)




# M-cone specific -------------------------------------------------------------------
rm(list=ls())

setwd("/Users/angueyraaristjm/Documents/LiMolec/zfRNAseq/20190827/20190827_DESeq2/")
directory <- "/Users/angueyraaristjm/Documents/LiMolec/zfRNAseq/20190827/20190827_DESeq2/"
getwd()

# remove subtrancript info from gCount.csv (sed 's/\.[0-9]//' gCount_ensdart.csv > gCount.csv)
# Read Data and run DESeq2 -------------------------------------------------------------------
countData <- as.matrix(read.csv(paste0(directory,"gCount.csv"), row.names="gene_id"))
#prepDE.py does not observe the order of the provided gtf_list
# countData <- countData[,c(13,14,15,16,17,18,25,26,27,28,29,30,19,20,21,22,23,24,25,7,8,9,10,11,12,1,2,3,4,5,6)]
colData <- read.csv(paste0(directory,"PHENO_DATA_M.csv"), sep="\t", row.names=1)

# use colData to reorganize countData
countDataSubset <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countDataSubset))
# save total number of mapped reads
# write.csv(colSums(countData), file = "nMappedReads.csv")

dds <- DESeqDataSetFromMatrix(countData = countDataSubset, colData = colData, design = ~subtype)
dds <- DESeq(dds)
res <- results(dds)
summary(res)
sprintf('n(DEGenes) = %g (p<0.05) ', sum(res$padj < 0.05, na.rm=TRUE))

# Log fold change shrinkage for visualization and ranking -------------------------------------------------------------------
resLFC <- lfcShrink(dds, coef="subtype_USL_vs_M", type="apeglm")
head(resLFC)

# Include Genename (descriptive) -------------------------------------------------
# Run this only if things have changed:
# resdata <- res # # # not merging with normalized counts
# resdata <- resdata[order(resdata$padj),]
# # For some reason this breaks the csv writing and symbol column is replaced by numbers
# # resdata$symbol <- tolower(resdata$symbol)
# genenames <- mapIds(org.Dr.eg.db, keys=rownames(resdata), column=c("GENENAME"), keytype="SYMBOL", multivals='first')
# # write.csv(genenames, file = "genenames.csv", col.names=c("symbol","genename"))
# write.csv(genenames, file = "genenames_M.csv", col.names=NA)


genenames <- read.csv("genenames_M.csv", sep=",")
colnames(genenames) <- c("symbol","genename")
genenames$genename <- gsub(",","",genenames$genename)
head(genenames)

# Save Results ------------------------------------------------------------
# results 
resdata <- res # # # not merging with normalized counts: resdata = merge(as.data.frame(resLFC), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
resdata <- resdata[order(resdata$padj),]
resdata$log2FoldChange = -resdata$log2FoldChange #inverting fold-chage
head(resdata)

# save raw results for python plotting
if (all(genenames$symbol == resdata$symbol)) {
  print("data frames DO match")
  resdata$genename = genenames$genename
  resdata <- resdata[c(ncol(resdata),1:ncol(resdata)-1)] #not sure why it's adding symbol again
  head(resdata)
  write.csv(resdata, file = "00_bySubtype/MvsUSL_rawFC.csv", row.names=TRUE, quote=FALSE)
  print("saved data frames to csv files")
} else {
  print("data frames do NOT match")
}

resdataLFC <- resLFC # # # not merging with normalized counts: resdata = merge(as.data.frame(resLFC), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
resdataLFC <- resdataLFC[order(resdataLFC$padj),]
# resdata$log2FoldChange = -resdata$log2FoldChange #inverting fold-chage
head(resdataLFC)

# save raw results for python plotting
if (all(genenames$symbol == resdataLFC$symbol)) {
  print("data frames DO match")
  resdataLFC$genename = genenames$genename
  resdataLFC <- resdataLFC[c(ncol(resdataLFC),1:ncol(resdataLFC)-1)] #not sure why it's adding symbol again
  head(resdataLFC)
  write.csv(resdataLFC, file = "00_bySubtype/MvsUSL_shrinkFC.csv", row.names=TRUE, quote=FALSE)
  print("saved data frames to csv files")
} else {
  print("data frames do NOT match")
}
# -------------------------------------------------------------------------
# volcano plot
resDF = resdata
# resDF = resdataLFC
EnhancedVolcano(resDF,
                lab = rownames(resDF),
                x = 'log2FoldChange',y = 'pvalue',
                xlim = c(-30, 30),
                title = 'M versus UV/S/L',
                pCutoff = 0.05,
                FCcutoff = 0.6,col=c('grey75', 'grey50', 'grey25', 'red','blue'),
                colAlpha = 1)


# S-cone specific -------------------------------------------------------------------
rm(list=ls())

setwd("/Users/angueyraaristjm/Documents/LiMolec/zfRNAseq/20190827/20190827_DESeq2/")
directory <- "/Users/angueyraaristjm/Documents/LiMolec/zfRNAseq/20190827/20190827_DESeq2/"
getwd()

# remove subtrancript info from gCount.csv (sed 's/\.[0-9]//' gCount_ensdart.csv > gCount.csv)
# Read Data and run DESeq2 -------------------------------------------------------------------
countData <- as.matrix(read.csv(paste0(directory,"gCount.csv"), row.names="gene_id"))
#prepDE.py does not observe the order of the provided gtf_list
# countData <- countData[,c(13,14,15,16,17,18,25,26,27,28,29,30,19,20,21,22,23,24,25,7,8,9,10,11,12,1,2,3,4,5,6)]
colData <- read.csv(paste0(directory,"PHENO_DATA_S.csv"), sep="\t", row.names=1)

# use colData to reorganize countData
countDataSubset <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countDataSubset))
# save total number of mapped reads
# write.csv(colSums(countData), file = "nMappedReads.csv")

dds <- DESeqDataSetFromMatrix(countData = countDataSubset, colData = colData, design = ~subtype)
dds <- DESeq(dds)
res <- results(dds)
summary(res)
sprintf('n(DEGenes) = %g (p<0.05) ', sum(res$padj < 0.05, na.rm=TRUE))

# Log fold change shrinkage for visualization and ranking -------------------------------------------------------------------
resLFC <- lfcShrink(dds, coef="subtype_UML_vs_S", type="apeglm")
head(resLFC)

# Include Genename (descriptive) -------------------------------------------------
# # Run this only if things have changed:
# resdata <- res # # # not merging with normalized counts
# resdata <- resdata[order(resdata$padj),]
# # For some reason this breaks the csv writing and symbol column is replaced by numbers
# # resdata$symbol <- tolower(resdata$symbol)
# genenames <- mapIds(org.Dr.eg.db, keys=rownames(resdata), column=c("GENENAME"), keytype="SYMBOL", multivals='first')
# # write.csv(genenames, file = "genenames.csv", col.names=c("symbol","genename"))
# write.csv(genenames, file = "genenames_S.csv", col.names=NA)


genenames <- read.csv("genenames_S.csv", sep=",")
colnames(genenames) <- c("symbol","genename")
genenames$genename <- gsub(",","",genenames$genename)
head(genenames)

# Save Results ------------------------------------------------------------
# results 
resdata <- res # # # not merging with normalized counts: resdata = merge(as.data.frame(resLFC), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
resdata <- resdata[order(resdata$padj),]
resdata$log2FoldChange = -resdata$log2FoldChange #inverting fold-chage
head(resdata)

# save raw results for python plotting
if (all(genenames$symbol == resdata$symbol)) {
  print("data frames DO match")
  resdata$genename = genenames$genename
  resdata <- resdata[c(ncol(resdata),1:ncol(resdata)-1)] #not sure why it's adding symbol again
  head(resdata)
  write.csv(resdata, file = "00_bySubtype/SvsUML_rawFC.csv", row.names=TRUE, quote=FALSE)
  print("saved data frames to csv files")
} else {
  print("data frames do NOT match")
}

resdataLFC <- resLFC # # # not merging with normalized counts: resdata = merge(as.data.frame(resLFC), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
resdataLFC <- resdataLFC[order(resdataLFC$padj),]
# resdata$log2FoldChange = -resdata$log2FoldChange #inverting fold-chage
head(resdataLFC)

# save raw results for python plotting
if (all(genenames$symbol == resdataLFC$symbol)) {
  print("data frames DO match")
  resdataLFC$genename = genenames$genename
  resdataLFC <- resdataLFC[c(ncol(resdataLFC),1:ncol(resdataLFC)-1)] #not sure why it's adding symbol again
  head(resdataLFC)
  write.csv(resdataLFC, file = "00_bySubtype/SvsUML_shrinkFC.csv", row.names=TRUE, quote=FALSE)
  print("saved data frames to csv files")
} else {
  print("data frames do NOT match")
}

# -------------------------------------------------------------------------
# volcano plot
resDF = resdata
# resDF = resdataLFC
EnhancedVolcano(resDF,
                lab = rownames(resDF),
                x = 'log2FoldChange',y = 'pvalue',
                xlim = c(-30, 30),
                title = 'S versus UV/M/L',
                pCutoff = 0.05,
                FCcutoff = 0.6,col=c('grey75', 'grey50', 'grey25', 'red','blue'),
                colAlpha = 1)




# UV-cone specific -------------------------------------------------------------------
rm(list=ls())

setwd("/Users/angueyraaristjm/Documents/LiMolec/zfRNAseq/20190827/20190827_DESeq2/")
directory <- "/Users/angueyraaristjm/Documents/LiMolec/zfRNAseq/20190827/20190827_DESeq2/"
getwd()

# remove subtrancript info from gCount.csv (sed 's/\.[0-9]//' gCount_ensdart.csv > gCount.csv)
# Read Data and run DESeq2 -------------------------------------------------------------------
countData <- as.matrix(read.csv(paste0(directory,"gCount.csv"), row.names="gene_id"))
#prepDE.py does not observe the order of the provided gtf_list
# countData <- countData[,c(13,14,15,16,17,18,25,26,27,28,29,30,19,20,21,22,23,24,25,7,8,9,10,11,12,1,2,3,4,5,6)]
colData <- read.csv(paste0(directory,"PHENO_DATA_U.csv"), sep="\t", row.names=1)

# use colData to reorganize countData
countDataSubset <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countDataSubset))
# save total number of mapped reads
# write.csv(colSums(countData), file = "nMappedReads.csv")

dds <- DESeqDataSetFromMatrix(countData = countDataSubset, colData = colData, design = ~subtype)
dds <- DESeq(dds, 'minReplicatesForReplace' = 5)
res <- results(dds)
summary(res)
sprintf('n(DEGenes) = %g (p<0.05) ', sum(res$padj < 0.05, na.rm=TRUE))

# Log fold change shrinkage for visualization and ranking -------------------------------------------------------------------
resLFC <- lfcShrink(dds, coef="subtype_U_vs_SML", type="apeglm")
head(resLFC)

# Include Genename (descriptive) -------------------------------------------------
# # Run this only if things have changed:
# resdata <- res # # # not merging with normalized counts
# resdata <- resdata[order(resdata$padj),]
# # For some reason this breaks the csv writing and symbol column is replaced by numbers
# # resdata$symbol <- tolower(resdata$symbol)
# genenames <- mapIds(org.Dr.eg.db, keys=rownames(resdata), column=c("GENENAME"), keytype="SYMBOL", multivals='first')
# # write.csv(genenames, file = "genenames.csv", col.names=c("symbol","genename"))
# write.csv(genenames, file = "genenames_U.csv", col.names=NA)


genenames <- read.csv("genenames_U.csv", sep=",")
colnames(genenames) <- c("symbol","genename")
genenames$genename <- gsub(",","",genenames$genename)
head(genenames)

# Save Results ------------------------------------------------------------
# results 
resdata <- res # # # not merging with normalized counts: resdata = merge(as.data.frame(resLFC), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
resdata <- resdata[order(resdata$padj),]
# resdata$log2FoldChange = -resdata$log2FoldChange #inverting fold-chage
head(resdata)

# save raw results for python plotting
if (all(genenames$symbol == resdata$symbol)) {
  print("data frames DO match")
  resdata$genename = genenames$genename
  resdata <- resdata[c(ncol(resdata),1:ncol(resdata)-1)] #not sure why it's adding symbol again
  head(resdata)
  write.csv(resdata, file = "00_bySubtype/UvsSML_rawFC.csv", row.names=TRUE, quote=FALSE)
  print("saved data frames to csv files")
} else {
  print("data frames do NOT match")
}

resdataLFC <- resLFC # # # not merging with normalized counts: resdata = merge(as.data.frame(resLFC), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
resdataLFC <- resdataLFC[order(resdataLFC$padj),]
# resdata$log2FoldChange = -resdata$log2FoldChange #inverting fold-chage
head(resdataLFC)

# save raw results for python plotting
if (all(genenames$symbol == resdataLFC$symbol)) {
  print("data frames DO match")
  resdataLFC$genename = genenames$genename
  resdataLFC <- resdataLFC[c(ncol(resdataLFC),1:ncol(resdataLFC)-1)] #not sure why it's adding symbol again
  head(resdataLFC)
  write.csv(resdataLFC, file = "00_bySubtype/UvsSML_shrinkFC.csv", row.names=TRUE, quote=FALSE)
  print("saved data frames to csv files")
} else {
  print("data frames do NOT match")
}

# -------------------------------------------------------------------------
# volcano plot
resDF = resdata
# resDF = resdataLFC
EnhancedVolcano(resDF,
                lab = rownames(resDF),
                x = 'log2FoldChange',y = 'pvalue',
                xlim = c(-30, 30),
                title = 'U versus S/M/L',
                pCutoff = 0.05,
                FCcutoff = 0.6,col=c('grey75', 'grey50', 'grey25', 'red','blue'),
                colAlpha = 1)


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Redoing PCA axis to have LFC and non-LFC versions
# ---------------------------------------------------------------------------------
# Rods vs Cones -------------------------------------------------------------------
rm(list=ls())

setwd("/Users/angueyraaristjm/Documents/LiMolec/zfRNAseq/20190827/20190827_DESeq2/")
directory <- "/Users/angueyraaristjm/Documents/LiMolec/zfRNAseq/20190827/20190827_DESeq2/"
getwd()

# remove subtrancript info from gCount.csv (sed 's/\.[0-9]//' gCount_ensdart.csv > gCount.csv)
# Read Data and run DESeq2 -------------------------------------------------------------------
countData <- as.matrix(read.csv(paste0(directory,"gCount.csv"), row.names="gene_id"))
#prepDE.py does not observe the order of the provided gtf_list
# countData <- countData[,c(13,14,15,16,17,18,25,26,27,28,29,30,19,20,21,22,23,24,25,7,8,9,10,11,12,1,2,3,4,5,6)]
colData <- read.csv(paste0(directory,"PHENO_DATA.csv"), sep="\t", row.names=1)
# use colData to reorganize countData
countDataSubset <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countDataSubset))
# save total number of mapped reads
# write.csv(colSums(countData), file = "nMappedReads.csv")

dds <- DESeqDataSetFromMatrix(countData = countDataSubset, colData = colData, design = ~type)
dds <- DESeq(dds, 'minReplicatesForReplace' = 5)
res <- results(dds)
summary(res)
sprintf('n(DEGenes) = %g (p<0.05) ', sum(res$padj < 0.05, na.rm=TRUE))

# vst transformation ------------------------------------------------------------------------------------------------
vsd <- vst(dds, blind=FALSE) # not blind because it's expected that differences in samples can cause large differences in gene counts
write.csv(assay(vsd), file = "00_bySubtype/RvsC_vstData.csv", row.names=TRUE, quote=FALSE)
# Log fold change shrinkage for visualization and ranking -------------------------------------------------------------------
resLFC <- lfcShrink(dds, coef="type_Rod_vs_Cone", type="apeglm")
head(resLFC)

# Include Genename (descriptive) -------------------------------------------------
# Run this only if things have changed:
resdata <- res # # # not merging with normalized counts
resdata <- resdata[order(resdata$padj),]
# For some reason this breaks the csv writing and symbol column is replaced by numbers
# resdata$symbol <- tolower(resdata$symbol)
genenames <- mapIds(org.Dr.eg.db, keys=rownames(resdata), column=c("GENENAME"), keytype="SYMBOL", multivals='first')
# write.csv(genenames, file = "genenames.csv", col.names=c("symbol","genename"))
write.csv(genenames, file = "genenames.csv", col.names=NA)


genenames <- read.csv("genenames.csv", sep=",")
colnames(genenames) <- c("symbol","genename")
genenames$genename <- gsub(",","",genenames$genename)
head(genenames)

# Save Results ------------------------------------------------------------
# results 
resdata <- res # # # not merging with normalized counts: resdata = merge(as.data.frame(resLFC), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
resdata <- resdata[order(resdata$padj),]
# resdata$log2FoldChange = -resdata$log2FoldChange #inverting fold-chage
head(resdata)

# save raw results for python plotting
if (all(genenames$symbol == resdata$symbol)) {
  print("data frames DO match")
  resdata$genename = genenames$genename
  resdata <- resdata[c(ncol(resdata),1:ncol(resdata)-1)] #not sure why it's adding symbol again
  head(resdata)
  write.csv(resdata, file = "00_bySubtype/RvsC_rawFC.csv", row.names=TRUE, quote=FALSE)
  print("saved data frames to csv files")
} else {
  print("data frames do NOT match")
}

resdataLFC <- resLFC # # # not merging with normalized counts: resdata = merge(as.data.frame(resLFC), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
resdataLFC <- resdataLFC[order(resdataLFC$padj),]
# resdata$log2FoldChange = -resdata$log2FoldChange #inverting fold-chage
head(resdataLFC)

# save raw results for python plotting
if (all(genenames$symbol == resdataLFC$symbol)) {
  print("data frames DO match")
  resdataLFC$genename = genenames$genename
  resdataLFC <- resdataLFC[c(ncol(resdataLFC),1:ncol(resdataLFC)-1)] #not sure why it's adding symbol again
  head(resdataLFC)
  write.csv(resdataLFC, file = "00_bySubtype/RvsC_shrinkFC.csv", row.names=TRUE, quote=FALSE)
  print("saved data frames to csv files")
} else {
  print("data frames do NOT match")
}

resdata[rownames(resdata)=='rho',]
# -------------------------------------------------------------------------
# volcano plot
resDF = resdata
resDF = resdataLFC
EnhancedVolcano(resDF,
                lab = rownames(resDF),
                x = 'log2FoldChange',y = 'pvalue',
                xlim = c(-30, 30),
                title = 'Rods versus Cones',
                pCutoff = 0.05,
                FCcutoff = 0.6,col=c('grey75', 'grey50', 'grey25', 'red','blue'),
                colAlpha = 1)




# ---------------------------------------------------------------------------------
#  US vs LM -------------------------------------------------------------------
rm(list=ls())

setwd("/Users/angueyraaristjm/Documents/LiMolec/zfRNAseq/20190827/20190827_DESeq2/")
directory <- "/Users/angueyraaristjm/Documents/LiMolec/zfRNAseq/20190827/20190827_DESeq2/"
getwd()

# remove subtrancript info from gCount.csv (sed 's/\.[0-9]//' gCount_ensdart.csv > gCount.csv)
# Read Data and run DESeq2 -------------------------------------------------------------------
countData <- as.matrix(read.csv(paste0(directory,"gCount.csv"), row.names="gene_id"))
#prepDE.py does not observe the order of the provided gtf_list
# countData <- countData[,c(13,14,15,16,17,18,25,26,27,28,29,30,19,20,21,22,23,24,25,7,8,9,10,11,12,1,2,3,4,5,6)]
colData <- read.csv(paste0(directory,"PHENO_DATA_Groups.csv"), sep="\t", row.names=1)
# use colData to reorganize countData
countDataSubset <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countDataSubset))
# save total number of mapped reads
# write.csv(colSums(countData), file = "nMappedReads.csv")

dds <- DESeqDataSetFromMatrix(countData = countDataSubset, colData = colData, design = ~subtype)
dds <- DESeq(dds, 'minReplicatesForReplace' = 5)
res <- results(dds)
summary(res)
sprintf('n(DEGenes) = %g (p<0.05) ', sum(res$padj < 0.05, na.rm=TRUE))

# Log fold change shrinkage for visualization and ranking -------------------------------------------------------------------
resLFC <- lfcShrink(dds, coef="subtype_US_vs_ML", type="apeglm")
head(resLFC)

# Include Genename (descriptive) -------------------------------------------------
# Run this only if things have changed:
resdata <- res # # # not merging with normalized counts
resdata <- resdata[order(resdata$padj),]
# For some reason this breaks the csv writing and symbol column is replaced by numbers
# resdata$symbol <- tolower(resdata$symbol)
genenames <- mapIds(org.Dr.eg.db, keys=rownames(resdata), column=c("GENENAME"), keytype="SYMBOL", multivals='first')
# write.csv(genenames, file = "genenames.csv", col.names=c("symbol","genename"))
write.csv(genenames, file = "genenames_Groups.csv", col.names=NA)


genenames <- read.csv("genenames_Groups.csv", sep=",")
colnames(genenames) <- c("symbol","genename")
genenames$genename <- gsub(",","",genenames$genename)
head(genenames)

# Save Results ------------------------------------------------------------
# results 
resdata <- res # # # not merging with normalized counts: resdata = merge(as.data.frame(resLFC), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
resdata <- resdata[order(resdata$padj),]
# resdata$log2FoldChange = -resdata$log2FoldChange #inverting fold-chage
head(resdata)

# save raw results for python plotting
if (all(genenames$symbol == resdata$symbol)) {
  print("data frames DO match")
  resdata$genename = genenames$genename
  resdata <- resdata[c(ncol(resdata),1:ncol(resdata)-1)] #not sure why it's adding symbol again
  head(resdata)
  write.csv(resdata, file = "00_bySubtype/USvsML_rawFC.csv", row.names=TRUE, quote=FALSE)
  print("saved data frames to csv files")
} else {
  print("data frames do NOT match")
}

resdataLFC <- resLFC # # # not merging with normalized counts: resdata = merge(as.data.frame(resLFC), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
resdataLFC <- resdataLFC[order(resdataLFC$padj),]
# resdata$log2FoldChange = -resdata$log2FoldChange #inverting fold-chage
head(resdataLFC)

# save raw results for python plotting
if (all(genenames$symbol == resdataLFC$symbol)) {
  print("data frames DO match")
  resdataLFC$genename = genenames$genename
  resdataLFC <- resdataLFC[c(ncol(resdataLFC),1:ncol(resdataLFC)-1)] #not sure why it's adding symbol again
  head(resdataLFC)
  write.csv(resdataLFC, file = "00_bySubtype/USvsML_shrinkFC.csv", row.names=TRUE, quote=FALSE)
  print("saved data frames to csv files")
} else {
  print("data frames do NOT match")
}

resdata[rownames(resdata)=='arr3a',]
# -------------------------------------------------------------------------
# volcano plot
resDF = resdata
resDF = resdataLFC
EnhancedVolcano(resDF,
                lab = rownames(resDF),
                x = 'log2FoldChange',y = 'pvalue',
                xlim = c(-30, 30),
                title = 'UV/S versus M/L',
                pCutoff = 0.05,
                FCcutoff = 0.6,col=c('grey75', 'grey50', 'grey25', 'red','blue'),
                colAlpha = 1)




# ---------------------------------------------------------------------------------
# M vs L -------------------------------------------------------------------
rm(list=ls())

setwd("/Users/angueyraaristjm/Documents/LiMolec/zfRNAseq/20190827/20190827_DESeq2/")
directory <- "/Users/angueyraaristjm/Documents/LiMolec/zfRNAseq/20190827/20190827_DESeq2/"
getwd()

# remove subtrancript info from gCount.csv (sed 's/\.[0-9]//' gCount_ensdart.csv > gCount.csv)
# Read Data and run DESeq2 -------------------------------------------------------------------
countData <- as.matrix(read.csv(paste0(directory,"gCount.csv"), row.names="gene_id"))
#prepDE.py does not observe the order of the provided gtf_list
# countData <- countData[,c(13,14,15,16,17,18,25,26,27,28,29,30,19,20,21,22,23,24,25,7,8,9,10,11,12,1,2,3,4,5,6)]
colData <- read.csv(paste0(directory,"PHENO_DATA_LM.csv"), sep="\t", row.names=1)
# use colData to reorganize countData
countDataSubset <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countDataSubset))
# save total number of mapped reads
# write.csv(colSums(countData), file = "nMappedReads.csv")

dds <- DESeqDataSetFromMatrix(countData = countDataSubset, colData = colData, design = ~subtype)
dds <- DESeq(dds, 'minReplicatesForReplace' = 5)
res <- results(dds)
summary(res)
sprintf('n(DEGenes) = %g (p<0.05) ', sum(res$padj < 0.05, na.rm=TRUE))

# Log fold change shrinkage for visualization and ranking -------------------------------------------------------------------
resLFC <- lfcShrink(dds, coef="subtype_M_vs_L", type="apeglm")
head(resLFC)

# Include Genename (descriptive) -------------------------------------------------
# Run this only if things have changed:
resdata <- res # # # not merging with normalized counts
resdata <- resdata[order(resdata$padj),]
# For some reason this breaks the csv writing and symbol column is replaced by numbers
# resdata$symbol <- tolower(resdata$symbol)
genenames <- mapIds(org.Dr.eg.db, keys=rownames(resdata), column=c("GENENAME"), keytype="SYMBOL", multivals='first')
# write.csv(genenames, file = "genenames.csv", col.names=c("symbol","genename"))
write.csv(genenames, file = "genenames_LvsM.csv", col.names=NA)


genenames <- read.csv("genenames_LvsM.csv", sep=",")
colnames(genenames) <- c("symbol","genename")
genenames$genename <- gsub(",","",genenames$genename)
head(genenames)

# Save Results ------------------------------------------------------------
# results 
resdata <- res # # # not merging with normalized counts: resdata = merge(as.data.frame(resLFC), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
resdata <- resdata[order(resdata$padj),]
# resdata$log2FoldChange = -resdata$log2FoldChange #inverting fold-chage
head(resdata)

# save raw results for python plotting
if (all(genenames$symbol == resdata$symbol)) {
  print("data frames DO match")
  resdata$genename = genenames$genename
  resdata <- resdata[c(ncol(resdata),1:ncol(resdata)-1)] #not sure why it's adding symbol again
  head(resdata)
  write.csv(resdata, file = "00_bySubtype/LvsM_rawFC.csv", row.names=TRUE, quote=FALSE)
  print("saved data frames to csv files")
} else {
  print("data frames do NOT match")
}

resdataLFC <- resLFC # # # not merging with normalized counts: resdata = merge(as.data.frame(resLFC), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
resdataLFC <- resdataLFC[order(resdataLFC$padj),]
# resdata$log2FoldChange = -resdata$log2FoldChange #inverting fold-chage
head(resdataLFC)

# save raw results for python plotting
if (all(genenames$symbol == resdataLFC$symbol)) {
  print("data frames DO match")
  resdataLFC$genename = genenames$genename
  resdataLFC <- resdataLFC[c(ncol(resdataLFC),1:ncol(resdataLFC)-1)] #not sure why it's adding symbol again
  head(resdataLFC)
  write.csv(resdataLFC, file = "00_bySubtype/LvsM_shrinkFC.csv", row.names=TRUE, quote=FALSE)
  print("saved data frames to csv files")
} else {
  print("data frames do NOT match")
}

resdata[rownames(resdata)=='thrb',]
# -------------------------------------------------------------------------
# volcano plot
resDF = resdata
resDF = resdataLFC
EnhancedVolcano(resDF,
                lab = rownames(resDF),
                x = 'log2FoldChange',y = 'pvalue',
                xlim = c(-30, 30),
                title = 'Rods versus Cones',
                pCutoff = 0.05,
                FCcutoff = 0.6,col=c('grey75', 'grey50', 'grey25', 'red','blue'),
                colAlpha = 1)




# ---------------------------------------------------------------------------------
# UV vs S -------------------------------------------------------------------
rm(list=ls())

setwd("/Users/angueyraaristjm/Documents/LiMolec/zfRNAseq/20190827/20190827_DESeq2/")
directory <- "/Users/angueyraaristjm/Documents/LiMolec/zfRNAseq/20190827/20190827_DESeq2/"
getwd()

# remove subtrancript info from gCount.csv (sed 's/\.[0-9]//' gCount_ensdart.csv > gCount.csv)
# Read Data and run DESeq2 -------------------------------------------------------------------
countData <- as.matrix(read.csv(paste0(directory,"gCount.csv"), row.names="gene_id"))
#prepDE.py does not observe the order of the provided gtf_list
# countData <- countData[,c(13,14,15,16,17,18,25,26,27,28,29,30,19,20,21,22,23,24,25,7,8,9,10,11,12,1,2,3,4,5,6)]
colData <- read.csv(paste0(directory,"PHENO_DATA_US.csv"), sep="\t", row.names=1)
# use colData to reorganize countData
countDataSubset <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countDataSubset))
# save total number of mapped reads
# write.csv(colSums(countData), file = "nMappedReads.csv")

dds <- DESeqDataSetFromMatrix(countData = countDataSubset, colData = colData, design = ~subtype)
dds <- DESeq(dds, 'minReplicatesForReplace' = 5)
res <- results(dds)
summary(res)
sprintf('n(DEGenes) = %g (p<0.05) ', sum(res$padj < 0.05, na.rm=TRUE))

# Log fold change shrinkage for visualization and ranking -------------------------------------------------------------------
resLFC <- lfcShrink(dds, coef="subtype_U_vs_S", type="apeglm")
head(resLFC)

# Include Genename (descriptive) -------------------------------------------------
# Run this only if things have changed:
resdata <- res # # # not merging with normalized counts
resdata <- resdata[order(resdata$padj),]
# For some reason this breaks the csv writing and symbol column is replaced by numbers
# resdata$symbol <- tolower(resdata$symbol)
genenames <- mapIds(org.Dr.eg.db, keys=rownames(resdata), column=c("GENENAME"), keytype="SYMBOL", multivals='first')
# write.csv(genenames, file = "genenames.csv", col.names=c("symbol","genename"))
write.csv(genenames, file = "genenames_UvsS.csv", col.names=NA)


genenames <- read.csv("genenames_UvsS.csv", sep=",")
colnames(genenames) <- c("symbol","genename")
genenames$genename <- gsub(",","",genenames$genename)
head(genenames)

# Save Results ------------------------------------------------------------
# results 
resdata <- res # # # not merging with normalized counts: resdata = merge(as.data.frame(resLFC), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
resdata <- resdata[order(resdata$padj),]
# resdata$log2FoldChange = -resdata$log2FoldChange #inverting fold-chage
head(resdata)

# save raw results for python plotting
if (all(genenames$symbol == resdata$symbol)) {
  print("data frames DO match")
  resdata$genename = genenames$genename
  resdata <- resdata[c(ncol(resdata),1:ncol(resdata)-1)] #not sure why it's adding symbol again
  head(resdata)
  write.csv(resdata, file = "00_bySubtype/UvsS_rawFC.csv", row.names=TRUE, quote=FALSE)
  print("saved data frames to csv files")
} else {
  print("data frames do NOT match")
}

resdataLFC <- resLFC # # # not merging with normalized counts: resdata = merge(as.data.frame(resLFC), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
resdataLFC <- resdataLFC[order(resdataLFC$padj),]
# resdata$log2FoldChange = -resdata$log2FoldChange #inverting fold-chage
head(resdataLFC)

# save raw results for python plotting
if (all(genenames$symbol == resdataLFC$symbol)) {
  print("data frames DO match")
  resdataLFC$genename = genenames$genename
  resdataLFC <- resdataLFC[c(ncol(resdataLFC),1:ncol(resdataLFC)-1)] #not sure why it's adding symbol again
  head(resdataLFC)
  write.csv(resdataLFC, file = "00_bySubtype/UvsS_shrinkFC.csv", row.names=TRUE, quote=FALSE)
  print("saved data frames to csv files")
} else {
  print("data frames do NOT match")
}

resdata[rownames(resdata)=='foxq2',]
# -------------------------------------------------------------------------
# volcano plot
resDF = resdata
resDF = resdataLFC
EnhancedVolcano(resDF,
                lab = rownames(resDF),
                x = 'log2FoldChange',y = 'pvalue',
                xlim = c(-30, 30),
                title = 'S versus UV',
                pCutoff = 0.05,
                FCcutoff = 0.6,col=c('grey75', 'grey50', 'grey25', 'red','blue'),
                colAlpha = 1)
