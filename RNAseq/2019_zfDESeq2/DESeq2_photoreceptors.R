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


# Library loading -------------------------------------------------------------------

library("DESeq2")
library("apeglm")
library("ggplot2")
library("org.Dr.eg.db") 
library("ReportingTools")
library("pcaExplorer")

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
#prepDE.py does not observe the order of the provided gtf_list, so reorganizing here to ease plotting later
# countData <- countData[,c(13,14,15,16,17,18,25,26,27,28,29,30,19,20,21,22,23,24,25,7,8,9,10,11,12,1,2,3,4,5,6)]
colData <- read.csv(paste0(directory,"PHENO_DATA.csv"), sep="\t", row.names=1)


countData <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countData))
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~type)
dds <- DESeq(dds)
res <- results(dds)
summary(res)
sprintf('n(DEGenes) = %g (p<0.1) ', sum(res$padj < 0.1, na.rm=TRUE))

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
# write.csv(genenames, file = "genenames.csv", col.names=NA)


genenames <- read.csv("genenames.csv", sep=",")
colnames(genenames) <- c("symbol","genename")
head(genenames)

# Save Results ------------------------------------------------------------
# results + normalized counts as log2((counts/average sequencing depth across samples)+0.5)
resdata <- merge(as.data.frame(resLFC), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
resdata <- resdata[order(resdata$padj),]
names(resdata)[1] <- 'symbol'
head(resdata)

# save raw results for python plotting
if (all(genenames$symbol == resdata$symbol)) {
    print("data frames DO match")
    resdata$genename = genenames$genename
    resdata <- resdata[c(1,ncol(resdata),2:ncol(resdata)-1)] #not sure why it's adding symbol again
    resdata <- resdata[c(1,2,4:ncol(resdata))]
    head(resdata)
    write.csv(resdata, file = "00_rodsVcones/rodVCones_raw.csv")
} else {
    print("data frames do NOT match")
}

# save excel friendly version
# columns can be rounded using: temp[c("baseMean")]=round_df(temp[c("baseMean")],digits = 2)
if (all(genenames$symbol == resdata$symbol)) {
  res_excel <- resdata
  res_excel[c("baseMean")]=round_df(res_excel[c("baseMean")],digits = 2)
  res_excel[c("log2FoldChange")]=round_df(res_excel[c("log2FoldChange")],digits = 4)
  res_excel[c("lfcSE")]=round_df(res_excel[c("lfcSE")],digits = 4)
  
  res_excel[c("R1","R2","R3","R4","R5","R6")]=round_df(res_excel[c("R1","R2","R3","R4","R5","R6")],digits = 2)
  res_excel[c("U1","U2","U3","U4","U5")]=round_df(res_excel[c("U1","U2","U3","U4","U5")],digits = 2)
  res_excel[c("S1","S2","S3","S4","S5","S6","S7")]=round_df(res_excel[c("S1","S2","S3","S4","S5","S6","S7")],digits = 2)
  res_excel[c("M1","M2","M3","M4","M5","M6")]=round_df(res_excel[c("M1","M2","M3","M4","M5","M6")],digits = 2)
  res_excel[c("L1","L2","L3","L4","L5","L6")]=round_df(res_excel[c("L1","L2","L3","L4","L5","L6")],digits = 2)
  
  write.csv(res_excel, file = "00_rodsVcones/rodVCones01_psorted.csv")
  write.csv(res_excel[order(res_excel$symbol),], file = "00_rodsVcones/rodVCones02_abc.csv")
  
  resexcel_pvalue <- subset(res_excel, padj<0.1)
  resexcel_Rods <- subset(resexcel_pvalue, log2FoldChange>0)
  resexcel_Cones <- subset(resexcel_pvalue, log2FoldChange<0)
  
  write.csv(resexcel_pvalue, file = "00_rodsVcones/rodVCones03_pvalue.csv")
  write.csv(resexcel_Rods[order(resexcel_Rods$baseMean),], file = "00_rodsVcones/rodVCones04_rods.csv")
  write.csv(resexcel_Cones[order(resexcel_Cones$baseMean),], file = "00_rodsVcones/rodVCones05_cones.csv")
  
} else {
  print("data frames do NOT match")
}



