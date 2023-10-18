#######DESeq analysis of hybrids versus pure forms for misexpression analysis#######
library("DESeq2")
library(edgeR)
library(tidyverse)

sampleFiles <- grep("STAR",list.files(directory),value=TRUE)
metaData <- read.csv('meta.txt', header = TRUE, sep = "\t")
metaData1 <- metaData[order(metaData$ID),]

individual <- dplyr::pull(metaData1, Individual)
conditions <- dplyr::pull(metaData1, Condition)
Population <- dplyr::pull(metaData1, Population)
Batch <- dplyr::pull(metaData1, Batch)

sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=conditions, population=Population, individual=individual, Batch=Batch)
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, design=~condition + Batch+ Population)
dds<-DESeq(ddsHTSeq)

matrix(resultsNames(dds))

dds <- estimateSizeFactors(dds)
idx <- rowSums(counts(dds, normalized=TRUE) >= 5 ) >= 4
dds <- dds[idx,]
dds <- DESeq(dds)

#####output for differential expression between hybrids vs inland pure form samples######
res = results(dds, list( c("population_Hybrid_vs_Inland") ))
res <- res[order(res$padj),]
write.table(res, file='Mixed_CNH1_population_effect_DeSeq_outpuT.csv', sep=',', quote=FALSE)

#####output for differential expression between hybrids vs coastal pure form samples######
res = results(dds, list( c("population_Hybrid_vs_Coastal") ))
res <- res[order(res$padj),]
write.table(res, file='Mixed_CNH1_population_effect_misexpress_DeSeq_outpuT.csv', sep=',', quote=FALSE)

###considered misexpressed if padj < 0.1 in both files above####
