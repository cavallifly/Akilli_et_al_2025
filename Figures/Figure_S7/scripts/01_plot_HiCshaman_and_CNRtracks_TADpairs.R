
args = commandArgs(trailingOnly=TRUE)

#devtools::load_all("/zssd/scratch/vincent.loubiere/vlite/")
devtools::load_all("./scripts_clean/vlite/")
#hicScreenshot

### Input parameters ###
mainDir	       <- "/work/user/mdistefano/nakilli/CNR/pipeline_analysis/bigWigs/"
HiCresolutions <- c(as.integer(args[1]))
handle <- 200000

chroms = c("chr3L", "chr3R" , "chr2L" , "chr2R" , "chrX" , "chrY" , "chrM" , "chr4")
lengths = c(28110227, 32079331, 23513712, 25286936, 23542271, 3667352, 19524, 1348131)
chromLengths = data.frame(chroms,lengths)
colnames(chromLengths) <- c("chrom","length") 
print(chromLengths)
#quit()

conditions     <- c("NOnub","SUMOnub")
conditionNames <- list("NOnub" = "Control", "SUMOnub" = "SUMO RNAi")

### START For plotting TADs ###
TADs <- read.table("TopDom_domains_hic_WD_NOnub_merge_dm6_NA_window_5_at_5000bp_filtered_refined_at_1000bp_with_state.bed", header=F)
colnames(TADs) <- c("chrom","start","end","name")
#TADs <- TADs[TADs$start > handle,]
#TADs$region <- paste0(TADs$chrom,":",as.integer(as.integer((TADs$start-handle)/HiCresolutions[[1]])*HiCresolutions[[1]]),"-",as.integer(as.integer((TADs$end+handle)/HiCresolutions[[1]])*HiCresolutions[[1]]))
#TADs$span   <- abs(TADs$start-TADs$end)
#print(head(TADs))
#beds       <- setNames(TADs$region, TADs$name)
#entensions <- setNames(TADs$span, TADs$name)
### END For plotting TADs ###

### START For plotting regions of interest ###
interesting_regions <- read.table("regions_to_plot_byDifference.bed", header=T)
interesting_regions$start <- ifelse(interesting_regions$start > handle, interesting_regions$start, handle)
interesting_regions$end <- pmin(
  interesting_regions$end,
  chromLengths$length[match(interesting_regions$chrom, chromLengths$chrom)] - handle
)
#interesting_regions$start <- as.numeric(interesting_regions$start)
#interesting_regions$end <- as.numeric(interesting_regions$end)
print(head(interesting_regions))
interesting_regions$region <- paste0(interesting_regions$chrom,":",as.integer(as.integer((interesting_regions$start-handle)/HiCresolutions[[1]])*HiCresolutions[[1]]),"-",as.integer(as.integer((interesting_regions$end+handle)/HiCresolutions[[1]])*HiCresolutions[[1]]))
interesting_regions$span   <- abs(interesting_regions$start-interesting_regions$end)
print(head(interesting_regions))
#quit()

loops <- read.table("interesting_loops_byDifference.bedpe", header=T)
print(loops)
beds <- setNames(interesting_regions$region, interesting_regions$name)
entensions <- setNames(interesting_regions$span, interesting_regions$name)

### END For plotting regions of interest ###

print(beds)
print(entensions)
#quit()

# For Figure 5
bwMaxs <- list(
    "NOnub"   = c(180, 370, 150, 740),
    "SUMOnub" = c(180, 150, 150, 740)
)
bwHeight <- 50


### Input parameters ###
for(geneCategory in c("UP","DOWN"))
{
# Selected genes
if(geneCategory == "UP")
{
  geneSetFile <- "UPgenes.txt"
}
if(geneCategory == "DOWN")
{
  geneSetFile <- "DOWNgenes.txt"
}
geneSet <- read.table(geneSetFile)
colnames(geneSet) <- c("Flybase","GeneID","chrom","TSS")
print(head(geneSet))


for(regionTag in names(beds))
{
    bed <- beds[[regionTag]]
    print(paste0("Region to visualize ",regionTag," ",bed), quote = FALSE)   

    coor <- unlist(tstrsplit(bed, ":|-"))
    chrom <- coor[1]
    start <- as.integer(coor[2])
    end <- as.integer(coor[3])
    print(paste0(chrom," ",start," ",end))

    for(HiCresolution in HiCresolutions)
    {

        CNRresolution  <- HiCresolution / 10

        for(condition in conditions)
    	{
            conditionName <- conditionNames[condition]
            bwMax         <- bwMaxs[[condition]]
	    print(paste0("Condition ",condition," ",conditionName), quote = FALSE)
	    print(paste0("Tracks Maxima ",paste(bwMax,collapse = " ")), quote = FALSE)

            outFile   <- paste0("HiCshaman_",HiCresolution,"bp_and_CNRtracks_plot_WD_",condition,"_",regionTag,"_",geneCategory,"genes.pdf")
	    if(file.exists(outFile)){next;}

	    mishaTrack <- paste0("hicScores_WD_",condition,"_merge_dm6_NA_1_k250_kexp500_step1000kb")    

	    nBins = as.integer((entensions[[regionTag]])/CNRresolution)
	    print(paste0("Number of bins to plot ",nBins))

	    TAD1 <- unlist(tstrsplit(regionTag, "_"))[1]
	    TAD2 <- unlist(tstrsplit(regionTag, "_"))[2]
	    print(paste0(TAD1," ",TAD2))
	    print(head(loops))
	    loopsTmp <- loops[loops$name == regionTag,c(2,3,4,5,6,7)]
	    colnames(loopsTmp) <- c("chrom1","start1","end1","chrom2","start2","end2")
	    print(head(loopsTmp))
	    write.table(loopsTmp,file="interesting_loops_tmp.bedpe",sep="\t", row.names=FALSE, quote=FALSE, col.names=T, append=F)	    

	    # Selected genes
	    geneSetTmp <- geneSet[geneSet$chrom==chrom,]
	    geneSetTmp1 <- geneSetTmp[loopsTmp$start1 <= geneSetTmp$TSS & geneSetTmp$TSS <= loopsTmp$end1,]$GeneID
	    geneSetTmp2 <- geneSetTmp[loopsTmp$start2 <= geneSetTmp$TSS & geneSetTmp$TSS <= loopsTmp$end2,]$GeneID
	    geneSetTmp <- c(geneSetTmp1,geneSetTmp2)
	    print(geneSetTmp)
	    selGenes <- geneSetTmp
	    if(length(geneSetTmp)==0)
	    {
	      next
	    }
	    print(paste0("Selected genes",selGenes))

	    for(state in c("PcG","Active","Het","Null"))
	    {
	      outFileTADsTmp <- paste0("TopDom_domains_",state,"state_to_plot.bed")
	      if (file.exists(outFileTADsTmp))
	      {
	        file.remove(outFileTADsTmp)
	      }
	      file.copy("Empty.bed", outFileTADsTmp)
	      TADsTmp <- TADs[grep(state,TADs$name),]
	      TADsTmp <- TADsTmp[TADsTmp$name == TAD1,c(1,2,3,4)]
	      print(paste0("Select TADs"))
	      write.table(TADsTmp,file=outFileTADsTmp,sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE, append=T)
	      TADsTmp <- TADs[grep(state,TADs$name),]
	      TADsTmp <- TADsTmp[TADsTmp$name == TAD2,c(1,2,3,4)]
	      print(paste0("Select TADs"))
	      write.table(TADsTmp,file=outFileTADsTmp,sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE, append=T)	      
	    }

    	    hicScreenshot(
	        map.name = paste0("Hi-C map ",conditionName),
                #mcool.file = mcoolFile,
		misha.track = mishaTrack,		
  	    	region     = bed,
		loops.file = "/work/user/mdistefano/nakilli/HiC/03_TADplots_cHiC_and_ChIPseq/interesting_loops_tmp.bedpe",
  	    	resolution = HiCresolution,
  	    	pdf.file   = outFile,		
  	    	tracks = c(
	    	    "Empty.bed",
		    "TopDom_domains_hic_WD_NOnub_merge_dm6_NA_window_5_at_5000bp_filtered_refined_at_1000bp_with_PcGstate.bed",
                    "TopDom_domains_hic_WD_NOnub_merge_dm6_NA_window_5_at_5000bp_filtered_refined_at_1000bp_with_Activestate.bed",
                    "TopDom_domains_hic_WD_NOnub_merge_dm6_NA_window_5_at_5000bp_filtered_refined_at_1000bp_with_Hetstate.bed",
                    "TopDom_domains_hic_WD_NOnub_merge_dm6_NA_window_5_at_5000bp_filtered_refined_at_1000bp_with_Nullstate.bed",		    
		    #"TopDom_domains_PcGstate_to_plot.bed",	    
		    #"TopDom_domains_Activestate_to_plot.bed",
		    #"TopDom_domains_Hetstate_to_plot.bed",
	    	    #"TopDom_domains_Nullstate_to_plot.bed",
	    	    "Empty.bed",		    
    	            paste0(mainDir,"cnr_H3K27me3_WD_",condition,"_merge_all_dm6_NA_normRPKM.bw"),
    	  	    paste0(mainDir,"cnr_PcXlinked_WD_",condition,"_merge_all_dm6_NALF2_normRPKM.bw"),
    	    	    paste0(mainDir,"cnr_H3K27ac_WD_",condition,"_merge_all_dm6_NA_normRPKM.bw"),
    	    	    paste0(mainDir,"cnr_H3K9me3_WD_D907C_merge_all_dm6_GS_normRPKM.bw")
  	        ),
  	    	gtf = "./dmel-all-r6.36.gtf",
		bw.max = bwMax,
	    	col = c("white","blue","red","darkgreen","black","white","blue","black","red","darkgreen"),
  	    	track.names= c(paste0(""),paste0("PcG"), paste0("Active"), paste0("Het"), paste0("Null"), paste0(""), paste0("H3K27me3"), paste0("Pc"), paste0("H3K27ac"), paste0("H3K9me3")),
		sel.gene.symbols = selGenes,
            	#sel.gene.symbols = c("dac","inv","Antp","Ubx","upd1", "upd2", "upd3", "zfh1", "lab", "Abd-B", "en", "Dronc", "h"),
	    	gene.height = 3,
	    	cex.gene.symbol = .9,
		space.height= 12,
		border.lwd = 2,
		border.col = c("white","white","white","white","white","white",NA,NA,NA,NA),
		nbins = nBins
            )
	    #quit()
        }
    }
}
}
for(state in c("PcG","Active","Het","Null"))
{
  outFileTADsTmp <- paste0("TopDom_domains_",state,"state_to_plot.bed")
  if (file.exists(outFileTADsTmp))
   {
     file.remove(outFileTADsTmp)
   }
}
file = "/work/user/mdistefano/nakilli/HiC/03_TADplots_cHiC_and_ChIPseq/interesting_loops_tmp.bedpe"
if (file.exists(file))
{
  file.remove(file)
}
