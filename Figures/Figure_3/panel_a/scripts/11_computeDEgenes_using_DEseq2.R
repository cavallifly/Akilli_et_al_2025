
##### 2023 DEseq2 analysis ###
library(tidyverse)
library(data.table)
library("DESeq2")
library("EnhancedVolcano")
args = commandArgs(trailingOnly=TRUE)
#options(scipen = 999)


#samples <- c("flybase_only_NOnub_SUMOnub","flybase_only_NOnub_SUMOnub_WT","flybase_Original")
#samples <- c("flybase_only_NOnub_SUMOnub_WT","flybase_Original")
samples <- c("flybase_only_NOnub_SUMOnub")

for(sample in samples)
{

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
countFile = paste0("countReadPairs_using_subread_",sample,".tab")
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
#row_indices[is.na(row_indices)] <- FALSE
#print(row_indices)
#countData <- countData[row_indices,]
#countData <- countData[-grep("NA",),]
print(head(countData))
#quit()

### Data with the conditions of untreated and treated conditions
samplesFile = paste0("samplesFiles_",sample,".txt")
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

print(paste0("Doing PCA clustering with vst transformed counts"))
vsd <- vst(dds, blind=FALSE)
print(vsd)
vsd <- vst(dds, blind = TRUE)

vsd_mat <- assay(vsd)

write.table(vsd_mat,
            file = paste0("vsd_values_",sample,".txt"),
            sep = "\t",
            quote = FALSE,
            col.names = NA)


PCAdata <- plotPCA(vsd,returnData=TRUE,intgroup=c("condition", "type", "id"))

pdf(paste0("PCAplot_",sample,".pdf"))
percentVar <- round(100 * attr(PCAdata, "percentVar"))
ggplot(PCAdata, aes(PC1, PC2, color=condition, shape=type)) +
  geom_point(size=3) +
  geom_text(aes(label = id), size=2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
#print(PCAdata)
dev.off()
#quit()

print(paste0("Running the DEseq2 analysis"))
resDEseq2 <- DESeq(dds)
res <- results(resDEseq2, tidy=TRUE)
print(head(res))
summary(res)

# Write normalized counts
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
normCountsFile <- paste0("normalizedCounts_using_DEseq2_",sample,".tsv")
write.table(normalized_counts, file=normCountsFile, sep="\t", row.names=T, quote=FALSE)

# Sets of genes of interest
PcGgenes <- read.table("PcG_targets_Parreno2024_selected.txt", header=F)
colnames(PcGgenes) <- c("flybase","genes")
genes_of_interest <- PcGgenes
#PcGgenes <- unlist(PcGgenes$genes)
#genes_of_interest <- data.frame(gene = PcGgenes)

PcGgenesAll <- read.table("PcG_targets_Parreno2024.txt", header=F)
colnames(PcGgenesAll) <- c("flybase","genes")
genes_of_interest_all <- PcGgenesAll
#PcGgenesAll <- unlist(PcGgenesAll$genes)
#genes_of_interest_all <- data.frame(gene = PcGgenesAll)

#interestGenes <- list(Hox          = c("Antp","Ubx","lab","Scr","Abd-B","abd-A","pb","Dfd"),
#                      #Wing_related = c("cv-2","cv","fuss","cat","notch","ci","spalt","cv-d"),
#                      PcGgenes     = PcGgenes
#                      #none         = c("")
#                      )
interestGenes <- list(PcGgenes     = PcGgenes)		      
print(interestGenes)
print(interestGenes$PcGgenes)
#print(intersect(interestGenes$Hox,interestGenes$PcGgenes))
#quit()

#if(is.null(interestGenes))
#{
#    interestGenes <- list(NA = c())
#}

for(set in names(interestGenes))
{

#set = NA

#for(mutant1 in conditions)
for(mutant1 in c("WD_NOnub"))
{
    #for(mutant2 in conditions)
    for(mutant2 in c("WD_SUMOnub"))    
    {
        if(mutant1 == mutant2){next}

        print(paste0("Doing the DEseq2 analysis of ",mutant2," vs ",mutant1))

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

    	outFile <- paste0("MAPlot_",mutant2,"_vs_",mutant1,"_genesOfInterest_",set,"_",sample,".pdf")
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
	outFile <- paste0("MAPlot_apeglm_",mutant2,"_vs_",mutant1,"_genesOfInterest_",set,"_",sample,".pdf")
        if(!file.exists(outFile))
        {
            pdf(outFile)
            plotMA(resLFC)
            dev.off()
        }
	outFile <- paste0("MAPlot_Normal_",mutant2,"_vs_",mutant1,"_genesOfInterest_",set,"_",sample,".pdf")
        if(!file.exists(outFile))
        {
            pdf(outFile)
            plotMA(resNorm)
            dev.off()
        }
	outFile <- paste0("MAPlot_ashr_",mutant2,"_vs_",mutant1,"_genesOfInterest_",set,"_",sample,".pdf")
        if(!file.exists(outFile))
        {
            pdf(outFile)
            plotMA(resAsh)
            dev.off()
        }	

    	outfile <- paste0("vulcanoPlot_",mutant2,"_vs_",mutant1,"_genesOfInterest_",set,"_",sample,".pdf")
        #if(file.exists(outFile))
        #{
	#    next;
	#}
	print(head(res_comparison))
	#quit()
	
	if(!file.exists(outfile))
	{

	    df <- res_comparison
	    # Convert DESeq2 results to data frame
	    df <- as.data.frame(df)

	    # Preserve gene names (rownames → column)
	    df$gene <- rownames(df)
	    
	    df$category <- "NS"
	    df$category[df$log2FoldChange >=  FCcutoff & df$padj <= pCutoff] <- "Up"
	    df$category[df$log2FoldChange <= -FCcutoff & df$padj <= pCutoff] <- "Down"

	    df$category[df$gene %in% genes_of_interest_all$flybase & df$category == "NS"] <- "PcG_targets_NS"
	    df$category[df$gene %in% genes_of_interest_all$flybase & df$category == "Up"] <- "PcG_targets_Up"
	    df$category[df$gene %in% genes_of_interest_all$flybase & df$category == "Down"] <- "PcG_targets_Down"	    	    
	    
	    print(head(df))
	    #quit()

	    # Clean data
	    df <- df %>%
	      filter(!is.na(log2FoldChange), !is.na(padj), padj > 0) %>%
	        mutate(
		    neg_log10_padj = -log10(padj),
	      )

	    # Subset for labeling
	    highlight_df <- df[df$gene %in% genes_of_interest$flybase,]
	    highlight_df <- highlight_df[order(highlight_df$gene),]
	    genes_of_interest <- genes_of_interest[genes_of_interest$flybase %in% highlight_df$gene,]
	    genes_of_interest <- genes_of_interest[order(genes_of_interest$flybase),]
	    print(nrow(highlight_df))
	    print(nrow(genes_of_interest))
	    #print(highlight_df[genes_of_interest$flybase],)
	    highlight_df$geneName <- genes_of_interest$gene[match(genes_of_interest$flybase,highlight_df$gene)]
	    print(highlight_df)
	    #quit()

	    print(unique(df$category))
	    #df$category <- factor(df$category, levels = c("Up", "Down", "NS", "PcG_targets_Up", "PcG_targets_Down", "PcG_targets_NS"))
	    #df$category <- factor(df$category, levels = c("Up", "PcG_targets_Up", "Down", "PcG_targets_Down", "NS", "PcG_targets_NS"))
	    df$category <- factor(df$category, levels = c("Down", "PcG_targets_Down", "NS", "PcG_targets_NS", "Up", "PcG_targets_Up"))	    	    
	    counts <- df %>%
	               group_by(category) %>%
		       summarize(n = n(), category=unique(category))
	    print(counts)
	    #quit()

	    y_cap <- 35

	    print(paste0("Number of outliers ",nrow(df[df$neg_log10_padj > y_cap,])," for y_cap ",y_cap))

	    df$neg_log10_padj[df$neg_log10_padj > y_cap] <- y_cap

	    # Plot
	    vulcanoPlot <- ggplot(df, aes(x = log2FoldChange, y = neg_log10_padj)) +
	      geom_point(aes(color = category, fill = category), size = 1.5, alpha = 0.6, shape = 21, stroke = 0.75) +

	      geom_point(
		data = df[df$category == "PcG_targets_Up",],
	        color = "red",
		fill = "yellow",
		shape = 21,
		size = 1.5,
		alpha = 0.6,
		stroke = 0.75
	      ) +
	      geom_point(
		data = df[df$category == "PcG_targets_Down",],
	        color = "blue",
		fill = "yellow",
		shape = 21,
		size = 1.5,
		alpha = 0.6,
		stroke = 0.75
	      ) +
	      geom_point(
		data = df[df$category == "PcG_targets_NS",],
	        color = "gray",
		fill = "yellow",
		shape = 21,
		size = 1.5,
                alpha = 0.6,
		stroke = 0.75
	      ) +	      

	      # Threshold lines
	      geom_vline(xintercept = c(FCcutoff), linetype = "dashed", color = "black", linewidth = 0.4) +
	      geom_vline(xintercept = c(-FCcutoff), linetype = "dashed", color = "black", linewidth = 0.4) +	      
	      geom_hline(yintercept = -log10(pCutoff), linetype = "dashed", color = "black", linewidth = 0.4) +

	      # Labels for top genes
	      geom_label_repel(
	          data = highlight_df,
		  aes(label = geneName),
		  size = 2.5,
		  max.overlaps = Inf,
		  box.padding = 1.2,
		  segment.size = 0.3,
		  point.padding = .10,
		  fill = "white"
	     ) +
  
	     # Colors
	     scale_fill_manual(values = c(
	     			       Up = "red",
				       Down = "blue",
				       NS = "grey",
	     			       PcG_targets_Up = "yellow",
				       PcG_targets_Down = "yellow",
				       PcG_targets_NS = "yellow"				       
				       ), guide = "none") +

	     # Colors
	     scale_color_manual(values = c(
	     			       Up = "red",
				       Down = "blue",
				       NS = "grey",
	     			       PcG_targets_Up = "red",
				       PcG_targets_Down = "blue",
				       PcG_targets_NS = "gray"				       
				       ),
	      			labels = c(
              Up   = paste0("All Up (n=", counts$n[counts$category=="Up"]+counts$n[counts$category=="PcG_targets_Up"], ")"),
              Down = paste0("All Down (n=", counts$n[counts$category=="Down"]+counts$n[counts$category=="PcG_targets_Down"], ")"),
              NS   = paste0("All NS (n=", counts$n[counts$category=="NS"]+counts$n[counts$category=="PcG_targets_NS"], ")"),
              PcG_targets_Up = paste0("PcG Up (n=", counts$n[counts$category=="PcG_targets_Up"], ")"),
              PcG_targets_Down = paste0("PcG Down (n=", counts$n[counts$category=="PcG_targets_Down"], ")"),
              PcG_targets_NS = paste0("PcG NS (n=", counts$n[counts$category=="PcG_targets_NS"], ")")
              )) +   

	      guides(
	        color = guide_legend(
		    override.aes = list(
		    shape = 21,
		    fill = c("blue","yellow","gray","yellow","red","yellow"),
		    stroke = 0.75
	            )
		    )
	       ) +

	     # assi NON distruttivi
  	     #coord_cartesian(
	     #    xlim = range(df$log2FoldChange, na.rm = TRUE),
		# ylim = c(0, y_cap),
		# clip = "off"
		# ) +

               # asse y con simbolo ≥
  	       #scale_y_continuous(
    	       #labels = function(x) ifelse(x == y_cap, paste0("≥", y_cap), x)
	       #) +

             #xlim(min(df$log2FoldChange), max(df$log2FoldChange)) + ylim(-25, max(df$neg_log10_padj)) +
             xlim(min(df$log2FoldChange), max(df$log2FoldChange)) + ylim(-3, y_cap) +	     
				
	     # Labels
	     labs(
	         x = expression(log[2]~"Fold Change"),
		 y = expression(-log[10]~"adj. p-value"),
		 color = NULL
		 ) +
  
	     # Theme (Nature-style minimalism)
	     theme_classic(base_size = 10) +
	       theme(
	            legend.position = 'top',
		    legend.key = element_blank(),
		    legend.key.size = unit(0.5, 'cm'),
		    legend.text = element_text(
		    		  size = 14),
		    title = element_text(
		                 size = 14),
		    axis.text = element_text(size = 14),
		    axis.title = element_text(size = 18),
		    aspect.ratio = 0.5
	     )

	    #vulcanoPlot <- EnhancedVolcano(res_comparison,
            #   lab = rownames(res_comparison),
  	    #   x = 'log2FoldChange',
            #   y = 'padj',
	    #   #col = c("Up" = "red", "Down" = "blue", "NS" = "gray"),
	    #   col = c("NS" = "gray", "Down" = "blue", "Up" = "red"),
	    #   #col = c("gray","gray","gray","red"),	       
	    #   #selectLab = c("Up", "Down"),
	    #   #legendLab = c("Up","Down","NS"),	       
            #   title = "", #paste0(mutant2," vs ",mutant1),
	    #   pointSize = 1.5,
	    #   labSize = 3.0,
	    #   selectLab = interestGenes[[set]],
	    #   labCol = 'black',
	    #   labFace = 'italic',
            #   pCutoff = pCutoff,
            #   FCcutoff = Log2FCcutoff,
	    #   colAlpha = 0.20,
 	    #   drawConnectors = TRUE,
	    #   widthConnectors = 0.2,
	    #   xlim = c(-10,15),
	    #   legendLabSize = 14,
	    #   legendIconSize = 4.0
	    #	   ) + coord_flip()
	    pdf(paste0("vulcanoPlot_",mutant2,"_vs_",mutant1,"_",set,"_",sample,".pdf"), width = 8)
            print(vulcanoPlot)
            dev.off()
        }   
        write.table(final, file=paste0("DEseq2_results_",mutant2,"_vs_",mutant1,"_",sample,".tsv"), sep="\t", row.names=FALSE, quote=FALSE)

        final_padj_FC <- final[which(final$padjust < pCutoff & abs(final$log2Ratio) > FCcutoff), ]

        final_padj_FC_NS <- final[-which(final$padjust < pCutoff & abs(final$log2Ratio) > FCcutoff), ]
        final_padj_FC_U  <- filter(final_padj_FC, log2Ratio > 0)
        final_padj_FC_D  <- filter(final_padj_FC, log2Ratio < 0)

	final_padj_FC_NS <- final_padj_FC_NS[order(final_padj_FC_NS$log2Ratio, decreasing = F),]
        final_padj_FC_U  <- final_padj_FC_U[order(final_padj_FC_U$log2Ratio, decreasing = F),]
        final_padj_FC_D  <- final_padj_FC_D[order(final_padj_FC_D$log2Ratio, decreasing = F),]
        print(head(final_padj_FC_NS))
	print(head(final_padj_FC_U))
        print(head(final_padj_FC_D))	
        write.table(final_padj_FC_NS, file=paste0("NS_DEgenes_",mutant2,"_vs_",mutant1,"_padj",pCutoff,"_FC",FCcutoff,"_",sample,".tsv")  , sep="\t", row.names=FALSE, quote=FALSE)
	write.table(final_padj_FC_U, file=paste0("Up_DEgenes_",mutant2,"_vs_",mutant1,"_padj",pCutoff,"_FC",FCcutoff,"_",sample,".tsv")  , sep="\t", row.names=FALSE, quote=FALSE)
        write.table(final_padj_FC_D, file=paste0("Down_DEgenes_",mutant2,"_vs_",mutant1,"_padj",pCutoff,"_FC",FCcutoff,"_",sample,".tsv"), sep="\t", row.names=FALSE, quote=FALSE)

    }     
}
}
}
quit()