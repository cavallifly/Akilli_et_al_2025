setwd("./")

devtools::load_all("./scripts_clean/vlite/")
mainDir="/work/user/mdistefano/nakilli/CNR/pipeline_analysis/bigWigs/"

PcColor = "black"

for(regionTag in c("PcG_regions"))
{

#udp1-3
#Antp chr3R:6896253-6999228
#Ubx chr3R:16656623-16734426
#inv chr2R:11474465-11509695
#zfh1 chr3R:30765926-30789202.
if(regionTag == "PcG_regions")
{
    #regionTag <- "PcG_regions"
    bed <- c("chr3R:6,640,000-7,090,000", "chr3R:16,630,000-17,010,000", "chr2R:11,440,000-11,610,000", "chr3R:30,700,000-30,970,000")
}

outFile = paste0("./plots_of_ChIP_tracks_",regionTag,".pdf")
if(file.exists(outFile)){next;}
pdf(outFile, width = 9, height = 5)
vl_par(mai= c(.9, 2, .9, .9))
print(bed)
vlite::bwScreenshot(
  bed = bed,
  tracks = c(
      
    paste0(mainDir,"cnr_H3K27me3_WD_NOnub_merge_all_dm6_NA_normRPKM.bw"),
    paste0(mainDir,"cnr_H3K27me3_WD_SUMOnub_merge_all_dm6_NA_normRPKM.bw"),
    paste0(mainDir,"cnr_PcXlinked_WD_NOnub_merge_all_dm6_NALF2_normRPKM.bw"),
    paste0(mainDir,"cnr_PcXlinked_WD_SUMOnub_merge_all_dm6_NALF2_normRPKM.bw"),
    paste0(mainDir,"cnr_H3K27ac_WD_NOnub_merge_all_dm6_NA_normRPKM.bw"),
    paste0(mainDir,"cnr_H3K27ac_WD_SUMOnub_merge_all_dm6_NA_normRPKM.bw"),
    paste0(mainDir,"cnr_H3K9me3_WD_D907C_merge_all_dm6_GS_normRPKM.bw")
  ),
  gtf = "./dmel-all-r6.36.gtf",
  # bw.max = 24,
#  bw.n.breaks = NULL,
  # nbins = NULL,
  col = c("blue","blue",PcColor,PcColor,"red","red","darkgreen"),
  track.names= c(paste0("H3K27me3 - Control"), paste0("H3K27me3 - SUMO RNAi"), paste0("Pc - Control"), paste0("Pc - SUMO RNAi"), paste0("H3K27ac - Control"), paste0("H3K27ac - SUMO RNAi"), paste0("H3K9me3 - Control D907C")),
  #sel.gene.symbols = c("dac","inv","Antp","Ubx","upd1", "upd2", "upd3", "zfh1", "lab", "Abd-B", "en", "Dronc", "h"),
  sel.gene.symbols = c("dac","inv","Antp","Ubx","upd1", "upd2", "upd3", "zfh1", "lab", "Abd-B", "en", "Dronc", "h"),  
  gene.height = 3,
  bw.max = c(460, 460, 500, 210, 180, 180, 740),  
  cex.gene.symbol = .6 #,
#  col.gene.ps = "grey60" ,
#  col.gene.ns = "grey60" 
)
dev.off()
}
