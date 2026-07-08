#setwd("/zssd/scratch/vincent.loubiere/projects/marco/")
devtools::load_all("./scripts_clean/vlite/")

args = commandArgs(trailingOnly=TRUE)

### Input parameters ###
mainDir	       <- "/work/user/mdistefano//nakilli/CNR/pipeline_analysis/bigWigs/"

HiCresolutions <- c(as.integer(args[1]))
if(is.na(HiCresolutions))
{
  print(paste0("Hi-C resolution not provided. Exiting!"))
  quit()
}
conditions     <- c("NOnub","SUMOnub")
conditionNames <- list("NOnub" = "Control", "SUMOnub" = "SUMO RNAi")
beds <- list(
    "PcG_zfh1_region" = "chr3R:30,330,000-31,400,000"
)
entensions <- list(
    "PcG_zfh1_region" = abs(30330000-31400000)
)

insMax <- list(
    "PcG_zfh1_region" = 1.2
)

insMin <- list(
    "PcG_zfh1_region" = -0.76
)

loops <- read.table("interesting_loops_for_ROI.bedpe", header=T)
print(loops)

bwMaxs <- list(
    "NOnub"   = c(NA,370, 370, 150, 740),
    "SUMOnub" = c(NA,370, 150, 150, 740)
)
bwMins <- list(
    "NOnub"   = c(NA,0, 0, 0, 0),
    "SUMOnub" = c(NA,0, 0, 0, 0)
)

bwHeight <- 50
### Input parameters ###

#for(vmin.hic in c(0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001))
for(vmin.hic in c(0.001, 0.00075, 0.00050, 0.00025, 0.00010))
{
#for(vmax.hic in c(0.5, 0.1, 0.05, 0.01, 0.005, 0.001))
for(vmax.hic in c(0.125, 0.1, 0.075))
{
if(vmin.hic >= vmax.hic)
{
  next
}
print(paste0(vmin.hic," ",vmax.hic))
for(regionTag in names(beds))
{
    bed <- beds[[regionTag]]
    print(paste0("Region to visualize ",bed), quote = FALSE)

    for(HiCresolution in HiCresolutions)
    {

        CNRresolution  <- HiCresolution / 10

        for(condition in conditions)
    	{
            conditionName <- conditionNames[condition]
            bwMax         <- bwMaxs[[condition]]
            bwMin         <- bwMins[[condition]]
	    bwMax[[1]]    <- insMax[[regionTag]]
	    bwMin[[1]]    <- insMin[[regionTag]]	    
	    
	    print(paste0("Condition ",condition," ",conditionName), quote = FALSE)
	    print(paste0("Tracks Minima ",paste(bwMin,collapse = " ")), quote = FALSE)
	    print(paste0("Tracks Maxima ",paste(bwMax,collapse = " ")), quote = FALSE)	    

	    gamma= 1 
	    #outFile   <- paste0("HiChiglass_",HiCresolution,"bp_and_CNRtracks_plot_WD_",condition,"_",regionTag,"_VLscale_vmin_",vmin.hic,"_vmax_",vmax.hic,"_gamma1.pdf")	    
	    outFile   <- paste0("HiChiglass_",HiCresolution,"bp_and_CNRtracks_plot_WD_",condition,"_",regionTag,"_vmin_",vmin.hic,"_vmax_",vmax.hic,"_gamma",gamma,".pdf")	    
	    if(file.exists(outFile)){next;}
    	    mcoolFile <- paste0("/work/user/mdistefano/nakilli/HiC/01_cool_files/hic_contacts_WD_",condition,"_merge_all_dm6_NA.mcool")

	    nBins = as.integer((entensions[[regionTag]])/CNRresolution)
	    print(paste0("Number of bins to plot ",nBins))

	    print(head(loops))
	    loopsTmp <- loops[loops$name == regionTag,c(2,3,4,5,6,7)]
	    colnames(loopsTmp) <- c("chrom1","start1","end1","chrom2","start2","end2")
	    print(head(loopsTmp))
	    write.table(loopsTmp,file="interesting_loops_tmp.bedpe",sep="\t", row.names=FALSE, quote=FALSE, col.names=T, append=F)

    	    hicScreenshot(
	        map.name = paste0("Hi-C map ",conditionName),
                mcool.file = mcoolFile,
  	    	region     = bed,
  	    	resolution = HiCresolution,
  	    	pdf.file   = outFile,
		vmin.hic   = vmin.hic,
		vmax.hic   = vmax.hic,
		#Cc.hic     = c("black","darkblue", "blue", "cornflowerblue", "grey90", "tomato", "red","darkorange","yellow"), # Color-scale from Vincent
		gamma.hic  = gamma,
		useRaster= F,
		loops.file = "/work/user/mdistefano/nakilli/HiC/03_TADplots_cHiC_and_ChIPseq/interesting_loops_tmp.bedpe",
  	    	tracks = c(
	    	    "Empty.bed",
	    	    #"TopDom_domains_hic_WD_NOnub_merge_dm6_NA_window_5_at_5000bp_filtered_refined_at_1000bp_with_PcGstate.bed",			    
	    	    #"TopDom_domains_hic_WD_NOnub_merge_dm6_NA_window_5_at_5000bp_filtered_refined_at_1000bp_with_Activestate.bed",
		    #"TopDom_domains_hic_WD_NOnub_merge_dm6_NA_window_5_at_5000bp_filtered_refined_at_1000bp_with_Hetstate.bed",
	    	    #"TopDom_domains_hic_WD_NOnub_merge_dm6_NA_window_5_at_5000bp_filtered_refined_at_1000bp_with_Nullstate.bed",
	    	    #"Empty.bed",
		    paste0("insulationProfiles_raw_WD_",condition,"_merge_at_1000bp.bw"),		    
    	            paste0(mainDir,"cnr_H3K27me3_WD_",condition,"_merge_all_dm6_NA_normRPKM.bw"),
    	  	    paste0(mainDir,"cnr_PcXlinked_WD_",condition,"_merge_all_dm6_NALF2_normRPKM.bw"),
    	    	    paste0(mainDir,"cnr_H3K27ac_WD_",condition,"_merge_all_dm6_NA_normRPKM.bw"),
    	    	    paste0(mainDir,"cnr_H3K9me3_WD_D907C_merge_all_dm6_GS_normRPKM.bw")
  	        ),
  	    	gtf = "./dmel-all-r6.36.gtf",
		bw.max = bwMax,
		bw.min = bwMin,		
	    	#col = c("white","blue","red","darkgreen","black","white",NA,"blue","black","red","darkgreen"),
	    	col = c("white",NA,"blue","black","red","darkgreen"),		
  	    	#track.names= c(paste0(""),paste0("PcG"), paste0("Active"), paste0("Het"), paste0("Null"), paste0(""), paste0("Insulation"), paste0("H3K27me3"), paste0("Pc"), paste0("H3K27ac"), paste0("H3K9me3")),
  	    	track.names= c(paste0(""), paste0("Insulation"), paste0("H3K27me3"), paste0("Pc"), paste0("H3K27ac"), paste0("H3K9me3")),				
            	sel.gene.symbols = c("dac","inv","Antp","Ubx","upd1", "upd2", "upd3", "zfh1", "lab", "Abd-B", "en", "Dronc", "h"),
	    	gene.height = 3,
	    	cex.gene.symbol = .9,
		space.height= 12,
		border.lwd = 2,
		#border.col = c("white","white","white","white","white","white","gray60",NA,NA,NA,NA),
		border.col = c("white","gray60",NA,NA,NA,NA),		
		nbins = nBins
            )
	    #quit()
        }
    }
}
}
}