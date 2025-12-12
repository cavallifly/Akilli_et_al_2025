##### 2023 DEseq2 analysis ###
library(tidyverse)
library(data.table)
library("DESeq2")
library("EnhancedVolcano")
args = commandArgs(trailingOnly=TRUE)
#options(scipen = 999)

#toMatch <- c("WT","NO","Gene")
#toMatch <- c("WT","SUMO","Gene")

#toMatch <- c("SUMOnub_Rep2","SUMOnub_Rep3","NOnub_Rep1","NOnub_Rep2","Gene")
#toMatch <- c("SUMO","NO","Gene")

#toMatch <- c("2801_Rep1","2801_Rep3","3601_Rep1","3601_Rep2","Gene")
#toMatch <- c("2801","3601","Gene")

if(is.na(args[1]) && is.na(args[2]))
{
    print(paste0("You didn't provide any in-line values for the p-adj cutoff and the Log2-Fold-Change cutoff, so I will use a default values"))   
}

if(is.na(args[1]))
{
    print(paste0("You didn't provide any in-line value for the p-adj cutoff, so I will use a default value of 0.05"))
    pCutoff  = 0.05
} else {
    pCutoff  <- as.numeric(args[1])
}

if(is.na(args[2]))
{
    print(paste0("You didn't provide any in-line value for the Fold-Change cutoff, so I will use a default value of 1.5,"))
    print(paste0("#NOTE: The Fold-Change you provide will be converted to a Log2-Fold-Change cutoff automatically"))
    FCcutoff     = 1.5
    Log2FCcutoff = log(FCcutoff)/log(2)    
} else {
    FCcutoff     <- as.numeric(args[2])
    Log2FCcutoff <- log(FCcutoff)/log(2)
}
print(paste0("You provided an p-adj cutoff of ",pCutoff," and a Fold-Change cutoff of ",FCcutoff," that "),quote=F,row.names=F)
print(paste0("has been converted to a Log2-Fold-Change cutoff of ",Log2FCcutoff),quote=F,row.names=F)

### File with the counts from subread
countFile = "countReadPairs_using_subread.tab"
if(!file.exists(countFile))
{
    countFile = "countReads_using_subread.tab"
}
print(paste0("Using counts in ",countFile))
countData <- read.table(countFile, header = TRUE)
#countData <- countData[,grep(paste(toMatch,collapse="|"),names(countData))]
#print(head(countData))
#print(nrow(countData))
#row_indices <- apply(countData, 1, function(row) any(as.numeric(row) > 10))
#print(length(row_indices))
#print(row_indices)
#countData <- countData[row_indices,]
#countData <- countData[-grep("NA",),]
#print(head(countData))
#quit()

### Data with the conditions of untreated and treated conditions
samplesFile = "samplesFiles.txt"
print(paste0("Using metaData in ",samplesFile))
metaData <- read.table(samplesFile, header = TRUE)
#metaData <- metaData[grep(paste(toMatch,collapse="|"),metaData$id),]
conditions <- unique(metaData$condition)
head(metaData)

# Load the data for the differential expression analysis
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData   = metaData,
                              design= ~condition,
			      tidy = TRUE)
dds

print(paste0("Pre-filtering of genes with reads sum >= 10"))
head(counts(dds))
keep <- which(rowSums(counts(dds)) >= 10)
dds <- dds[keep,]
head(counts(dds))

#quit()

print(paste0("Running the DEseq2 analysis"))
resDEseq2 <- DESeq(dds)
res <- results(resDEseq2, tidy=TRUE)
print(head(res))
summary(res)

	print(paste0("Doing PCA clustering with vst transformed counts"))
	vsd <- vst(dds, blind=FALSE)
	print(vsd)
	PCAdata <- plotPCA(vsd,returnData=TRUE,intgroup=c("condition", "type", "id"))

	percentVar <- round(100 * attr(PCAdata, "percentVar"))
	p <- ggplot(PCAdata, aes(PC1, PC2, color=condition, shape=type)) +
	  geom_point(size=3) +
          geom_text(aes(label = id), size=2) +
	  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
	  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
	  coord_fixed()
	pdf(paste0("PCAplot_allSample.pdf"))	  
	print(p)
	dev.off()
	quit()

# Sets of genes of interest
#PcGgenes <- read.table("PcG_targets.csv", header=F)
#PcGgenes <- unlist(PcGgenes$V1)
#interestGenes <- list(Hox          = c("Antp","Ubx","lab","Scr","Abd-B","abd-A","pb","Dfd"),
#                      Wing_related = c("cv-2","cv","fuss","cat","notch","ci","spalt","cv-d"),
#                      PcGgenes     = PcGgenes,
#                      none         = c("")
#                      )
#print(interestGenes)
#print(interestGenes$Hox)
#print(intersect(interestGenes$Hox,interestGenes$PcGgenes))
#quit()

#if(is.null(interestGenes))
#{
#    interestGenes <- list(NA = c())
#}

#for(set in names(interestGenes))
#{

set = NA

for(mutant1 in conditions)
{
    for(mutant2 in conditions)
    {
        if(mutant1 == mutant2){next}

	# Write normalized counts
	dds <- estimateSizeFactors(dds)
	normalized_counts <- counts(dds, normalized=TRUE)
	normCountsFile <- paste0("normalizedCounts_using_DEseq2_",mutant2,"_vs_",mutant1,".tsv")
	write.table(normalized_counts, file=normCountsFile, sep="\t", row.names=T, quote=FALSE)

        print(paste0("Doing the DEseq2 analysis of ",mutant2," vs ",mutant1))

	print(paste0("Doing PCA clustering with vst transformed counts"))
	vsd <- vst(dds, blind=FALSE)
	print(vsd)
	PCAdata <- plotPCA(vsd,returnData=TRUE,intgroup=c("condition", "type", "id"))

	percentVar <- round(100 * attr(PCAdata, "percentVar"))
	p <- ggplot(PCAdata, aes(PC1, PC2, color=condition, shape=type)) +
	  geom_point(size=3) +
          geom_text(aes(label = id), size=2) +
	  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
	  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
	  coord_fixed()
	pdf(paste0("PCAplot_",mutant2,"_vs_",mutant1,".pdf"))	  
	print(p)
	dev.off()
	quit()

	res_comparison <- results(resDEseq2, contrast = c("condition", mutant2, mutant1))
        res_comparison <- res_comparison[order(res_comparison$padj),]
	print(head(res_comparison))
        #print(row.names(res_comparison))
        final <- data.frame(
    	   GENEID           = row.names(res_comparison),
	   log2BaseMean     = log2(res_comparison$baseMean),
	   log2Ratio        = res_comparison$log2FoldChange,
	   Stderr_log2Ratio = res_comparison$lfcSE,
	   pvalue           = res_comparison$pvalue,
	   padjust          = res_comparison$padj
	   )

    	outFile <- paste0("MAPlot_",mutant2,"_vs_",mutant1,"_genesOfInterest_",set,".pdf")
	if(!file.exists(outFile))
        {
	   #	    MAplot <- plotMA(res_comparison)
       	    pdf(outFile)
	    plotMA(res_comparison)
	    dev.off()
	    #	    quit()
	}


	# Test different shrinkage procedure
	head(resDEseq2)	
	# because we are interested in treated vs untreated, we set 'coef=2'
	print(resultsNames(resDEseq2))
	resNorm <- lfcShrink(resDEseq2, coef=2, type="normal")
	resAsh  <- lfcShrink(resDEseq2, coef=2, type="ashr")
	resLFC  <- lfcShrink(resDEseq2, coef=2, type="apeglm")
	outFile <- paste0("MAPlot_apeglm_",mutant2,"_vs_",mutant1,"_genesOfInterest_",set,".pdf")
        if(!file.exists(outFile))
        {
            pdf(outFile)
            plotMA(resLFC)
            dev.off()
        }
	outFile <- paste0("MAPlot_Normal_",mutant2,"_vs_",mutant1,"_genesOfInterest_",set,".pdf")
        if(!file.exists(outFile))
        {
            pdf(outFile)
            plotMA(resNorm)
            dev.off()
        }
	outFile <- paste0("MAPlot_ashr_",mutant2,"_vs_",mutant1,"_genesOfInterest_",set,".pdf")
        if(!file.exists(outFile))
        {
            pdf(outFile)
            plotMA(resAsh)
            dev.off()
        }	

    	outfile <- paste0("vulcanoPlot_",mutant2,"_vs_",mutant1,"_genesOfInterest_",set,".pdf")
	if(!file.exists(outfile))
	{
            vulcanoPlot <- EnhancedVolcano(res_comparison,
               lab = NA, #rownames(res_comparison),
  	       x = 'log2FoldChange',
               y = 'padj',
               title = paste0(mutant2," vs ",mutant1),
	       pointSize = 1.5,
	       labSize = 6.0,
	       #	       selectLab = interestGenes[[set]],
	       labCol = 'black',
	       labFace = 'italic',
               pCutoff = pCutoff,
               FCcutoff = Log2FCcutoff,
	       colAlpha = 0.40,
 	       drawConnectors = TRUE,
	       widthConnectors = 0.2
		   )
	    	pdf(paste0("vulcanoPlot_",mutant2,"_vs_",mutant1,"_",set,".pdf"))
            	print(vulcanoPlot)
            	dev.off()
        }   
        write.table(final, file=paste0("DEseq2_results_",mutant2,"_vs_",mutant1,".tsv"), sep="\t", row.names=FALSE, quote=FALSE)

        final_padj_FC <- final[which(final$padjust < pCutoff & abs(final$log2Ratio) > Log2FCcutoff), ]

        final_padj_FC_U <- filter(final_padj_FC, log2Ratio > 0)
        final_padj_FC_D <- filter(final_padj_FC, log2Ratio < 0)

        final_padj_FC_U <- final_padj_FC_U[order(final_padj_FC_U$log2Ratio, decreasing = F),]
        final_padj_FC_D <- final_padj_FC_D[order(final_padj_FC_D$log2Ratio, decreasing = F),]
        print(head(final_padj_FC_U))
        print(head(final_padj_FC_D))	
        write.table(final_padj_FC_U, file=paste0("Up_DEgenes_",mutant2,"_vs_",mutant1,"_padj",pCutoff,"_FC",FCcutoff,".tsv")  , sep="\t", row.names=FALSE, quote=FALSE)
        write.table(final_padj_FC_D, file=paste0("Down_DEgenes_",mutant2,"_vs_",mutant1,"_padj",pCutoff,"_FC",FCcutoff,".tsv"), sep="\t", row.names=FALSE, quote=FALSE)

    }     
}
#}
quit()