##### 2023 DEseq2 analysis ###
library(tidyverse)
library(data.table)
library("DESeq2")
args = commandArgs(trailingOnly=TRUE)
#options(scipen = 999)

### File with the counts from subread
countFile = "countReadPairs_using_subread.tab"
if(!file.exists(countFile))
{
    countFile = "countReads_using_subread.tab"
}
print(paste0("Using counts in ",countFile))
countData <- read.table(countFile, header = TRUE)
head(countData)

### Data with the conditions of untreated and treated conditions
samplesFile = "samplesFiles.txt"
print(paste0("Using metaData in ",samplesFile))
metaData <- read.table(samplesFile, header = T)
head(metaData)

# Load the data for the differential expression analysis
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData   = metaData,
                              design= ~condition,
			      tidy = TRUE)
dds

# Compute the FPKMs and write them in a file
fpkmInput <- dds
fpkmInput <- estimateSizeFactors(fpkmInput)
print(mcols(fpkmInput))
print(head(fpkmInput))

genesFile="countReadPairs_using_subread_with_geneLength.tab"
if(!file.exists(genesFile))
{
    genesFile="countReads_using_subread_with_geneLength.tab"
}
genes <- read.table(genesFile, header = TRUE)
print(paste0("Using genes in ",genesFile))
print(head(genes))

mcols(fpkmInput)$basepairs = as.integer(genes$Length)
print(head(fpkmInput))

FPKM <- fpkm(fpkmInput, robust = FALSE)
print(head(FPKM))

FPKMfile="FPKM.tsv"
write.table(FPKM, file = FPKMfile, sep="\t", row.names=TRUE, quote=FALSE)

### We may add here analysis on the FPKM ###