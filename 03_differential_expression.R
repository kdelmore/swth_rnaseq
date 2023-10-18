######DESeq code below for calling differantial expression within a tissue for pure forms samples -- example below from CNH ####
library("DESeq2")
library(edgeR)
library(tidyverse)

###read in individual count files from STAR alignement and HTCount#####
sampleFiles <- grep("STAR",list.files(directory),value=TRUE)
metaData <- read.csv('meta.txt', header = TRUE, sep = "\t")
metaData1 <- metaData[order(metaData$ID),]

individual <- dplyr::pull(metaData1, Individual)
conditions <- dplyr::pull(metaData1, Condition)
Population <- dplyr::pull(metaData1, Population)

sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=conditions, population=Population, individual=individual)
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, design=~condition + population + condition:population)
dds<-DESeq(ddsHTSeq)

matrix(resultsNames(dds))

#####removing genes with low counts#####
dds <- estimateSizeFactors(dds)
idx <- rowSums(counts(dds, normalized=TRUE) >= 5 ) >= 4
dds <- dds[idx,]
dds <- DESeq(dds)

####DE genes for migratory state#####
res <- results(dds, list( c("condition_Winter_vs_Spring") ))
ix = which.min(res$padj) # most significant
ix
res <- res[order(res$padj),]
write.table(res, file='Pure_CNH1_treatment_effect_DeSeq_outpuT.csv', sep=',', quote=FALSE)

####DE genes for subspecies#####
res = results(dds, list( c("population_Inland_vs_Coastal") ))
res <- res[order(res$padj),]
write.table(res, file='Pure_CNH1_population_effect_DeSeq_outpuT.csv', sep=',', quote=FALSE)

####DE genes for GxE#####
res = results(dds, name="conditionWinter.populationInland")
res <- res[order(res$padj),]
write.table(res, file='Pure_CNH1_interaction_effect_DeSeq_outpuT.csv', sep=',', quote=FALSE)
