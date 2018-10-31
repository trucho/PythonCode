#!/bin/bash

module load trimmomatic || exit 1
java -jar $TRIMMOJAR PE -phred33 /data/angueyraaristjm/20181018_RNAseq/05_zfUV_R1.fastq.gz /data/angueyraaristjm/20181018_RNAseq/05_zfUV_R2.fastq.gz /data/angueyraaristjm/20181018_Trimmed/05_zfUV_fP.fq.gz /data/angueyraaristjm/20181018_Trimmed/05_zfUV_fU.fq.gz /data/angueyraaristjm/20181018_Trimmed/05_zfUV_rP.fq.gz /data/angueyraaristjm/20181018_Trimmed/05_zfUV_rU.fq.gz ILLUMINACLIP:/usr/local/apps/trimmomatic/Trimmomatic-0.36/adapters/TruSeq2-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:70