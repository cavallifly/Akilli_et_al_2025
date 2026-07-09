#Do 1 sample_table file per comparison!!!
library(DiffBind)
library(profileplyr)
library(GenomicRanges)
library(ChIPseeker)
BiocManager::install("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
BiocManager::install("BRGenomics")
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(clusterProfiler)
library(ggplot2)
library("GenomicFeatures")
library("BRGenomics")
library(ReactomePA)
library(AnnotationDbi)
library(GenomicTools.fileHandler)
library(enrichplot)
library(biomaRt)
library(tidyverse)
args = commandArgs(trailingOnly=TRUE)
warnings()

comparisons <- c("WD_SUMOnub_vs_WD_NOnub")

#target <- "H3K27ac"
target <- args[1]

for(comparison in comparisons)
{
    inDir <- paste0(comparison)
    condition1 <- unlist(strsplit(comparison,'_vs_'))[1]
    condition2 <- unlist(strsplit(comparison,'_vs_'))[2]
    print(paste0("Performing comparison ",comparison))    
    print(paste0("Condition1 ",condition1))
    print(paste0("Condition2 ",condition2))

    samples <- read.delim("sample_table.tsv")
    print(head(samples))
    K27me3 <- dba(sampleSheet = samples)
    #quit()

    #Counting reads in peaks, summits option is in case you want to center the counting in a range from the summit
    K27me3.counted <- dba.count(K27me3, summits=FALSE, filter=0) #no filter for peaks with low overlapping read counts, summit is for count in a given range from the summit
    print(head(K27me3.counted)) #FRiP: proportion of reads for that sample that overlap a peak in the consensus peakset
    #quit()
    K27me3.counted <- dba.contrast(K27me3.counted, contrast=c("Condition",condition1,condition2))
    K27me3.counted
    dba.show(K27me3.counted,bContrasts=TRUE)
    K27me3.analyzed <- dba.analyze(K27me3.counted, method=DBA_DESEQ2, bBlacklist= FALSE, bGreylist = FALSE) #  DESEq2 without blacklist
    dba.show(K27me3.analyzed,bContrasts=TRUE) # see the amount of differential peaks

    #p1 <- plot(K27me3.analyzed, contrast=1)

    # PCA plot (uses normalized binding matrix)
    #XXX
    # Count reads and normalize
    #samples <- dba.count(samples)
    #samples <- dba.normalize(samples, method=DBA_ALL_METHODS)
    #XXX
    #dba.plotPCA(samples, 
    #        attributes=DBA_CONDITION,   # color/grouping
    #        label=DBA_ID)              # labels = sample IDs

    #dba.plotPCA(K27me3.analyzed,DBA_TISSUE,label=DBA_CONDITION)

    #write report
    report_K27me3 <- dba.report(K27me3.analyzed, contrast = 1, method=DBA_DESEQ2, th=1, bCounts=TRUE)
    report_K27me3 <- as.data.frame(report_K27me3)  
    write.table(report_K27me3, file=paste0("report_",target,"_",condition1,"_vs_",condition2,".tsv"), sep="\t", quote=F, row.names=F)        
    print(head(report_K27me3))
    #quit()
    #scatter plot normalized counts
    colors <- c("blue")

    # Doing the scatter plot
    #p3 <- ggplot(report_K27me3, aes(x=paste0("Conc_",condition2), y=paste0("Conc_",condition1))) +
    p3 <- ggplot(report_K27me3, aes(x=Conc_WD_NOnub, y=Conc_WD_SUMOnub)) +
       geom_point(size=2) +
       geom_point(color="blue") +
       geom_abline(intercept = 0, slope = 1) +
       ylab(paste0("log2 Conc. - ",condition1)) +
       xlab(paste0("log2 Conc. - ",condition2))
    pdf(file = paste0("scatterPlot_Conc_",target,"_",comparison,".pdf"))        
    print(p3)
    dev.off()

    ## Doing the MAplot
    p3_MA <- ggplot(report_K27me3, aes(x=Conc, y= Fold)) +
        geom_point(size=2,alpha=0.5,na.rm=T) +
	geom_point(color="blue") +
	geom_abline(intercept = 0, slope = 0) +
	geom_abline(intercept = -0.58, slope = 0, linetype = "dashed", color= "gray") +
	geom_abline(intercept = 0.58, slope = 0, linetype = "dashed", color= "gray") +
	ylim(-3,3) +
	ggtitle(paste0("Mark - ",condition1," vs ",condition2)) +
	ylab("log2 Fold change") +
	xlab("log2 Conc")
    pdf(file = paste0("MAplot_Conc_",target,"_",comparison,".pdf"))
    print(p3_MA)
    dev.off()
}
