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

# Setup -------------------------------------------------------------------
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
# 2019_09_07: sample S7 is actually an M_cone. Will leave things as they are and just modify PHENO_DATA.
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

# for excel, columns can be rounded using: temp[c("baseMean")]=round_df(temp[c("baseMean")],digits = 2)

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
  write.csv(head(resdata), file = "00_LvsUSM/LvsUSM_test.csv", row.names=FALSE, quote=FALSE)
  write.csv(resdata, file = "00_LvsUSM/LvsUSM_raw.csv", row.names=FALSE, quote=FALSE)
  print("saved data frames to csv files")
} else {
  print("data frames do NOT match")
}
# -------------------------------------------------------------------------

















# Plotters ------------------------------------------------------------

# plot a single gene: counts (normalized by seq depth and +0.5 for log plotting)
test <- plotCounts(dds, gene="rho", intgroup="type", col =c('blue','blue'), fg='white', col.lab ='white', col.main ='white', col.sub ='white', col.axis='white', bg='white')
test <- plotCounts(dds, gene="foxq2", intgroup="subtype", col =c('red','blue'), fg='white', col.lab ='white', col.main ='white', col.sub ='white', col.axis='white', bg='white')

# more customizable plot of a single gene: counts (normalized by seq depth and +0.5 for log plotting)
data <- plotCounts(dds, gene="tbx2a", intgroup=c("subtype"), returnData=TRUE)
data
ggplot(data, aes(x=subtype, y=count, color=subtype)) +
  scale_y_log10() + 
  geom_point(position=position_jitter(width=.1,height=0))

# volcano plot
resDF = resdata
# resDF = resLFC
EnhancedVolcano(resDF,
                lab = rownames(resDF),
                x = 'log2FoldChange',y = 'pvalue',
                xlim = c(-30, 30),
                title = 'L versus UV/S/M',
                pCutoff = 0.05,
                FCcutoff = 0.6,col=c('grey75', 'grey50', 'grey25', 'red','blue'),
                colAlpha = 1)

# visualize as heatmap
n_genes = 100
library("pheatmap")
ntd <- normTransform(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:n_genes]
df <- as.data.frame(colData(dds)[c("subtype")])
assay(ntd)[select,]
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

#built-in
plotPCA(ntd, intgroup=c("subtype"))

##ggplot
pcaData <- plotPCA(ntd, intgroup=c("subtype"), returnData=TRUE)
write.csv(pcaData, file = "00_LvsUSM/pcaData.csv", row.names=FALSE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=subtype)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

plotDispEsts(dds)
#probably means that data should be prefiltered. Maybe remove all lines with very low counts


# launch PCA App
pcaExplorer(dds = dds)



# trying to get PC table
# rld_rods <- rlogTransformation(dds) # using rlog transformation (very slow)
sumExp_log2 <- SummarizedExperiment(log2(counts(dds, normalized=TRUE) + 1),colData=colData(dds)) #using log2+1 pseudocounts
pcaobj <- prcomp(t(assay(sumExp_log2)), rank=3) # calculate the weight of each gene to the first 3 principal components
pcaobj$rotation[1:10,] # this extracts weights for a single gene 
write.csv(pcaobj$rotation, file = "00_LvsUSM/pcaWeights.csv", row.names=TRUE)

# extract the top genes that weigh PC1 the most
PC1_Groups = hi_loadings(pcaobj, whichpc = 1, topN = 20,exprTable=counts(dds))
head(PC1_Groups,40)

# or make plot
hi_loadings(pcaobj, whichpc = 1, topN = 40)
# save plot as 10 x 30 inches pdf