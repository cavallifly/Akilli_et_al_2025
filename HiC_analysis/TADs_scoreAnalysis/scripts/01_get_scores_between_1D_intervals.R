require(reshape2) # cat tabulate to matrix
library(misha)
library(shaman)
library(dplyr)
library(ggplot2)
#library(ggpubr)
#library("ggrepel")

args = commandArgs(trailingOnly=TRUE)
source('/zdata/data/auxFunctions/auxFunctions.R')
options("scipen"=999)

mDBloc <- '/zdata/data/mishaDB/trackdb/'
db <- 'dm6'
dbDir <- paste0(mDBloc,db,'/')
gdb.init(dbDir)
gdb.reload()

obsData <- list(
	WD_2801    = c("hic.hicScores_WD_2801_merge_dm6_NA_1.hicScores_WD_2801_merge_dm6_NA_1_k250_kexp500_step1000kb"),
	WD_3601    = c("hic.hicScores_WD_3601_merge_dm6_NA_1.hicScores_WD_3601_merge_dm6_NA_1_k250_kexp500_step1000kb"),
	WD_NOnub   = c("hic.hicScores_WD_NOnub_merge_dm6_NA_1.hicScores_WD_NOnub_merge_dm6_NA_1_k250_kexp500_step1000kb"),
	WD_SUMOnub = c("hic.hicScores_WD_SUMOnub_merge_dm6_NA_1.hicScores_WD_SUMOnub_merge_dm6_NA_1_k250_kexp500_step1000kb")
		     )
refSample <- "WD_NOnub"
refChrom <- "genomewide"
chrom <- args[1]
print(obsData)

colors <- rainbow(10)

intervals1File <- paste0("./all_gene_fragments_",chrom,"_dmel-all-r6.36_1.gtf")
intervals1     <- read.table(intervals1File, header=F)
intervals1     <- intervals1[,c(1,2,3,4)]
colnames(intervals1) <- c("chrom","start","end","interval")
print(head(intervals1))

intervals2File <- paste0("./all_gene_fragments_",chrom,"_dmel-all-r6.36_2.gtf")
intervals2     <- read.table(intervals2File, header=F)
intervals2     <- intervals2[,c(1,2,3,4)]
colnames(intervals2) <- c("chrom","start","end","interval")
print(head(intervals2))

samples = names(obsData)
print(samples)

outfile <- paste0("scores_trans1Dinterval_",chrom,"_gene_vs_domain.tab")

if(!file.exists(outfile))
{
    header <- paste(c("chrom1","start1","end1","chrom2","start2","end2","score","sample"),collapse="\t")
    write(header,file=outfile,sep="\t",append=T)

    for(sample in samples)
    {
        print(paste0("Analysing ",sample))

        # Get tracks and create virtual tracks: this depends only on the specific sample
        track2D <- obsData[[sample]]

	tmpIntervals1 <- intervals1
        print(head(tmpIntervals1))
	tmpIntervals2 <- intervals2
        print(head(tmpIntervals2))

        chrIntrv2D <- build2D(tmpIntervals1,tmpIntervals2,names1=tmpIntervals1$interval,names2=tmpIntervals2$interval)
	data     <- gextract(track2D,intervals=chrIntrv2D) #,iterator=chrIntrv2D)

        df <- data.frame(data[,1:length(data)-1],sample)			    
        colnames(df) <- c("chrom1","start1","end1","chrom2","start2","end2","score","sample")
        write.table(df,file=outfile,sep="\t",quote=F,row.names=F,col.names=F,append=T)

    } # Close cycle over sample
}
