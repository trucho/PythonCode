install.packages("org.Dr.eg.db")
install.packages("org.Dr.eg.db", repos="http://bioconductor.org/packages/3.1/data/annotation")
# Library loading -------------------------------------------------------------------
rm(list=ls())
library("DESeq2")
library("apeglm")
library("ggplot2")
library("org.Dr.eg.db") 
library("ReportingTools")
library("pcaExplorer")

# Use this if loading of Dr.db fails
# options(connectionObserver = NULL)
# Setup -------------------------------------------------------------------
setwd("/Users/angueyraaristjm/Documents/LiMolec/zfRNAseq/20190827/20190827_DESeq2/")
directory <- "/Users/angueyraaristjm/Documents/LiMolec/zfRNAseq/20190827/20190827_DESeq2/"
getwd()

round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  df[,nums] <- round(df[,nums], digits = digits)
  (df)
}

# remove subtrancript info from gCount.csv (sed 's/\.[0-9]//' gCount_ensdart.csv > gCount.csv)
# Read Data and run DESeq2 -------------------------------------------------------------------
countData <- as.matrix(read.csv(paste0(directory,"gCount.csv"), row.names="gene_id"))
# only keep transcription factors
tfList = read.csv(paste0(directory,"00_bySubtype/zfin_TFs.csv"))
countData = subset(countData, rownames(countData) %in% tfList$symbol)
# -------------------------------------------------------------------
# RODS VS CONES
colData <- read.csv(paste0(directory,"PHENO_DATA.csv"), sep="\t", row.names=1)
# 2019_09_07: sample S7 is actually an M_cone. Will leave things as they are and just modify PHENO_DATA.
# use colData to reorganize countData
countData <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countData))
# save total number of mapped reads
write.csv(colSums(countData), file = "nMappedReads_TFs.csv")
# -------------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~type)
dds <- DESeq(dds)
res <- results(dds)
summary(res)
sprintf('n(DEGenes) = %g (p<0.1) ', sum(res$padj < 0.1, na.rm=TRUE))
# -------------------------------------------------------------------


# for excel, columns can be rounded using: temp[c("baseMean")]=round_df(temp[c("baseMean")],digits = 2)

# Log fold change shrinkage for visualization and ranking -------------------------------------------------------------------
resLFC <- lfcShrink(dds, coef="type_Rod_vs_Cone", type="apeglm")
head(resLFC)


# Include Genename (descriptive) -------------------------------------------------
# # Run this only if things have changed:
# resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
# resdata <- resdata[order(resdata$padj),]
# names(resdata)[1] <- 'symbol'
# # For some reason this breaks the csv writing and symbol column is replaced by numbers
# # resdata$symbol <- tolower(resdata$symbol)
# genenames <- mapIds(org.Dr.eg.db, keys=resdata[,c("symbol")], column=c("GENENAME"), keytype="SYMBOL", multivals='first')
# # write.csv(genenames, file = "genenames.csv", col.names=c("symbol","genename"))
# write.csv(genenames, file = "genenamesTFs.csv", col.names=NA)


genenames <- read.csv("genenamesTFs.csv", sep=",")
colnames(genenames) <- c("symbol","genename")
genenames$genename <- gsub(",","",genenames$genename)
head(genenames)

# Save Results ------------------------------------------------------------
# results 
resdata <- res # # # not merging with normalized counts: resdata = merge(as.data.frame(resLFC), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
resdata <- resdata[order(resdata$padj),]
# resdata$log2FoldChange = -resdata$log2FoldChange #inverting fold-change
head(resdata,20)

# save raw results for python plotting
if (all(genenames$symbol == resdata$symbol)) {
  print("data frames DO match")
  resdata$genename = genenames$genename
  resdata <- resdata[c(ncol(resdata),1:ncol(resdata)-1)] #not sure why it's adding symbol again
  head(resdata)
  write.csv(resdata, file = "00_bySubtype/TFs_RvC_rawFC.csv", row.names=TRUE, quote=FALSE)
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
  write.csv(resdataLFC, file = "00_bySubtype/TFs_RvC_shrinkFC.csv", row.names=TRUE, quote=FALSE)
  print("saved data frames to csv files")
} else {
  print("data frames do NOT match")
}


# resdata[rownames(resdata)=='opn1sw1',]
# ------------------------------------------------------------

# plot a single gene: counts (normalized by seq depth and +0.5 for log plotting)
test <- plotCounts(dds, gene="foxq2", intgroup="subtype", col =c('red','green','black','magenta','blue'), fg='white', col.lab ='white', col.main ='white', col.sub ='white', col.axis='white', bg='white')

# more customizable plot of a single gene: counts (normalized by seq depth and +0.5 for log plotting)
data <- plotCounts(dds, gene="nrl", intgroup=c("type"), returnData=TRUE)
data
ggplot(data, aes(x=type, y=count, color=type)) +
  scale_y_log10() + 
  geom_point(position=position_jitter(width=.1,height=0))


# visualize as heatmap
n_genes = 100
library("pheatmap")
ntd <- normTransform(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:n_genes]
df <- as.data.frame(colData(dds)[c("subtype")])
assay(ntd)[select,]
pheatmap(assay(ntd)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df)

#built-in
plotPCA(ntd, intgroup=c("subtype"))

##ggplot
pcaData <- plotPCA(ntd, intgroup=c("subtype"), returnData=TRUE)
write.csv(pcaData, file = "00_bySubtype/TFs_pcaData.csv", row.names=FALSE)
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
pcaobj$rotation[1:40,] # this extracts weights for a single gene 
write.csv(pcaobj$rotation, file = "00_bySubtype/TFs_pcaWeights_RvC.csv", row.names=TRUE)

                 
# extract the top genes that weigh PC1 the most
# manual way
head(pcaobj$rotation[order(pcaobj$rotation[,1]),],20)
tail(pcaobj$rotation[order(pcaobj$rotation[,1]),],20)
#use built-in function
PC1_Groups = hi_loadings(pcaobj, whichpc = 1, topN = 20,exprTable=counts(dds))
head(PC1_Groups,40)

# or make plot
PC1genes = hi_loadings(pcaobj, whichpc = 1, topN = 960,exprTable=counts(dds))
PC1genes
pcaobj# save plot as 10 x 30 inches pdf
