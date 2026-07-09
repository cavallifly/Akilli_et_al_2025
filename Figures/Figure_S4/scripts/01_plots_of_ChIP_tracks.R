setwd("./")

devtools::load_all("./scripts_clean/vlite/")
mainDir="../../CNR/pipeline_analysis/bigWigs/"

for(regionTag in c("differentialPcG_regions"))
{

  if(regionTag == "differentialPcG_regions")
  {
    bed <- c("chrX:340,000-370,000","chr2L:7,595,000-7,625,000","chr2R:5,095,000-5,125,000")
  }
  print(bed)

  outFile = paste0("./plots_of_ChIP_tracks_",regionTag,"_FigS4b.pdf")

  if(file.exists(outFile)){next;}
  pdf(outFile, width = 9, height = 5)
  vl_par(mai= c(.9, 2, .9, .9))

  vlite::bwScreenshot(
    bed = bed,
    tracks = c(
      paste0("/work/user/mdistefano/nakilli/HiC/03_TADplots_cHiC_and_ChIPseq/Empty.bed"),    
      paste0("/work/user/mdistefano/nakilli/HiC/03_TADplots_cHiC_and_ChIPseq/differentialH3K27me3peaks.bed"),
      paste0("/work/user/mdistefano/nakilli/HiC/03_TADplots_cHiC_and_ChIPseq/Empty.bed"),          
      paste0(mainDir,"cnr_H3K27me3_WD_NOnub_merge_all_dm6_NA_normRPKM.bw"),
      paste0(mainDir,"cnr_H3K27me3_WD_SUMOnub_merge_all_dm6_NA_normRPKM.bw")
    ),
    gtf = "./dmel-all-r6.36_differentialH3K27me3peaks.gtf",
    #bw.max = 24,
    #bw.n.breaks = NULL,
    # nbins = NULL,
    col = c("white","blue","white","blue","blue"),
    track.names= c(paste0(""),paste0("DP"),paste0(""),paste0("H3K27me3 - Control"), paste0("H3K27me3 - SUMO RNAi")),
    sel.gene.symbols = c("yellow","CG3107","gus","CG45370","Cyp4d21","CG6739","CG3777","CG43181","CG43182","lncRNA:yar"),
    gene.height = 3,
    bw.max = c(300, 300),
    cex.gene.symbol = .6
  )
  dev.off()
}