#install.packages("ggstatsplot")
#install.packages("ggpubr",force=)

#library(ggstatsplot)
library(ggplot2)
library(ggpubr)
library(rtracklayer)
library(tidyr)

scoresSUMO_NO_3L <- read.delim("./avg_scores_trans1Dinterval_chr3L_SUMOnub_NOnub_all_domains_vs_all_domains.tab",header = FALSE)
scoresSUMO_NO_2L <- read.delim("./avg_scores_trans1Dinterval_chr2L_SUMOnub_NOnub_all_domains_vs_all_domains.tab",header = FALSE)
scoresSUMO_NO_2R <- read.delim("./avg_scores_trans1Dinterval_chr2R_SUMOnub_NOnub_all_domains_vs_all_domains.tab",header = FALSE)
scoresSUMO_NO_3R <- read.delim("./avg_scores_trans1Dinterval_chr3R_SUMOnub_NOnub_all_domains_vs_all_domains.tab",header = FALSE)
scoresSUMO_NO_X  <- read.delim("./avg_scores_trans1Dinterval_chrX_SUMOnub_NOnub_all_domains_vs_all_domains.tab",header = FALSE)

colnames(scoresSUMO_NO_3L)<-c("Chr1","start1","end1","Chr2","start2","end2","interval1","interval2","av_score","sample")
colnames(scoresSUMO_NO_3R)<-c("Chr1","start1","end1","Chr2","start2","end2","interval1","interval2","av_score","sample")
colnames(scoresSUMO_NO_2L)<-c("Chr1","start1","end1","Chr2","start2","end2","interval1","interval2","av_score","sample")
colnames(scoresSUMO_NO_2R)<-c("Chr1","start1","end1","Chr2","start2","end2","interval1","interval2","av_score","sample")
colnames(scoresSUMO_NO_X)<-c("Chr1","start1","end1","Chr2","start2","end2","interval1","interval2","av_score","sample")
#Domains for ChiaSig

###ChiaSig is the all NOnub and all SUMOnub domains

ChiaSig<-rbind(scoresSUMO_NO_3L,scoresSUMO_NO_3R,scoresSUMO_NO_2R,scoresSUMO_NO_2L,scoresSUMO_NO_X)
########
max(ChiaSig$av_score)
min(ChiaSig$av_score)

ChiaSig$av_score<-ChiaSig$av_score+100
head(ChiaSig)

head(ChiaSig)
ChiaSig2<-ChiaSig[,c(1:6,9,10)]
ChiaSig_NO<-ChiaSig2[ChiaSig2$sample=="WD_NOnub",]
ChiaSig_SUMO<-ChiaSig2[ChiaSig2$sample=="WD_SUMOnub",]
ChiaSig_NO<-ChiaSig_NO[,c(1:7)]
ChiaSig_SUMO<-ChiaSig_SUMO[,c(1:7)]

write.table(ChiaSig_NO, file = "ChiaSig_NO.bedpe", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(ChiaSig_SUMO, file = "ChiaSig_SUMO.bedpe", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

nrow(ChiaSig_NO)
nrow(ChiaSig_SUMO)

getwd()


scoresSUMO_NO_3R_f<-scoresSUMO_NO_3R[,7:10]
head(scoresSUMO_NO_3R_f)
scoresSUMO_NO_3L_f<-scoresSUMO_NO_3L[,7:10]
head(scoresSUMO_NO_3L_f)
scoresSUMO_NO_2L_f<-scoresSUMO_NO_2L[,7:10]
head(scoresSUMO_NO_2L_f)
scoresSUMO_NO_2R_f<-scoresSUMO_NO_2R[,7:10]
head(scoresSUMO_NO_2R_f)
scoresSUMO_NO_X_f<-scoresSUMO_NO_X[,7:10]
head(scoresSUMO_NO_X_f)

domains1_X<-scoresSUMO_NO_X[,c(1,2,3,7)]
domains2_X<-scoresSUMO_NO_X[,c(4,5,6,8)]
colnames(domains1_X)<-c("Chr","Start","End","Type")
colnames(domains2_X)<-c("Chr","Start","End","Type")
DomainsX<-rbind(domains1_X,domains2_X)
DomainsX<-DomainsX[!duplicated(DomainsX),]
DomainsX_GR<-GenomicRanges::makeGRangesFromDataFrame(DomainsX,seqnames.field="Chr",
                                                      start.field='Start',end.field="End",keep.extra.columns = TRUE,na.rm=TRUE)



domains1_3L<-scoresSUMO_NO_3L[,c(1,2,3,7)]
domains2_3L<-scoresSUMO_NO_3L[,c(4,5,6,8)]
colnames(domains1_3L)<-c("Chr","Start","End","Type")
colnames(domains2_3L)<-c("Chr","Start","End","Type")
Domains3L<-rbind(domains1_3L,domains2_3L)
Domains3L<-Domains3L[!duplicated(Domains3L),]
Domains3L_GR<-GenomicRanges::makeGRangesFromDataFrame(Domains3L,seqnames.field="Chr",
                                                      start.field='Start',end.field="End",keep.extra.columns = TRUE,na.rm=TRUE)


domains1_2L<-scoresSUMO_NO_2L[,c(1,2,3,7)]
domains2_2L<-scoresSUMO_NO_2L[,c(4,5,6,8)]
colnames(domains1_2L)<-c("Chr","Start","End","Type")
colnames(domains2_2L)<-c("Chr","Start","End","Type")
Domains2L<-rbind(domains1_2L,domains2_2L)
Domains2L<-Domains2L[!duplicated(Domains2L),]
Domains2L_GR<-GenomicRanges::makeGRangesFromDataFrame(Domains2L,seqnames.field="Chr",
                                                      start.field='Start',end.field="End",keep.extra.columns = TRUE,na.rm=TRUE)

domains1_2R<-scoresSUMO_NO_2R[,c(1,2,3,7)]
domains2_2R<-scoresSUMO_NO_2R[,c(4,5,6,8)]
colnames(domains1_2R)<-c("Chr","Start","End","Type")
colnames(domains2_2R)<-c("Chr","Start","End","Type")
Domains2R<-rbind(domains1_2R,domains2_2R)
Domains2R<-Domains2R[!duplicated(Domains2R),]
Domains2R_GR<-GenomicRanges::makeGRangesFromDataFrame(Domains2R,seqnames.field="Chr",
                                                      start.field='Start',end.field="End",keep.extra.columns = TRUE,na.rm=TRUE)


domains1_3R<-scoresSUMO_NO_3R[,c(1,2,3,7)]
domains2_3R<-scoresSUMO_NO_3R[,c(4,5,6,8)]
colnames(domains1_3R)<-c("Chr","Start","End","Type")
colnames(domains2_3R)<-c("Chr","Start","End","Type")
Domains3R<-rbind(domains1_3R,domains2_3R)
Domains3R<-Domains3R[!duplicated(Domains3R),]
Domains3R_GR<-GenomicRanges::makeGRangesFromDataFrame(Domains3R,seqnames.field="Chr",
                                                      start.field='Start',end.field="End",keep.extra.columns = TRUE,na.rm=TRUE)

AllDomainsSN<-c(Domains3R_GR,Domains3L_GR,Domains2R_GR,Domains2L_GR,DomainsX_GR)

SUMO<-read.table("./DEseq2_results_WD_SUMOnub_vs_WD_NOnub_table.tsv",header=T)

RNAseq_GR<-GenomicRanges::makeGRangesFromDataFrame(SUMO,seqnames.field="chrom",
                                                 start.field='start',end.field="end",keep.extra.columns = TRUE,na.rm=TRUE)
RNA_Domains<-mergeByOverlaps(AllDomainsSN,RNAseq_GR)


colnames(RNA_Domains)<-c("Domains3R_GR","Interval","SUMO_sig_GR","GeneID","FBID" ,
                         "TSS","FPKM_WD_NOnub_Rep1","FPKM_WD_NOnub_Rep2","FPKM_WD_NOnub_Rep3",
                         "FPKM_WD_SUMOnub_Rep1","FPKM_WD_SUMOnub_Rep2","FPKM_WD_SUMOnub_Rep3",
                         "NormCounts_WD_NOnub_Rep1","NormCounts_WD_NOnub_Rep2","NormCounts_WD_NOnub_Rep3",
                         "NormCounts_WD_SUMOnub_Rep1","NormCounts_WD_SUMOnub_Rep2","NormCounts_WD_SUMOnub_Rep3",
                         "log2BaseMean","log2Ratio","Stderr_log2Ratio","pvalue","padjust")

#Assign interaction classes
All_S_N_avS<-rbind(scoresSUMO_NO_3R_f,scoresSUMO_NO_3L_f,scoresSUMO_NO_2L_f,scoresSUMO_NO_2R_f,scoresSUMO_NO_X_f)


#chrom color
cat<-c('PcG','Het','Null','Active')
All_S_N_avS$domains<-'x'
All_S_N_avS$domains2<-'y'

#Assign if inter or intra domain contact

i<-1

for (i in 1:4) {
  
  All_S_N_avS[grep(cat[i],All_S_N_avS$interval1),]$domains<-cat[i]
  All_S_N_avS[grep(cat[i],All_S_N_avS$interval2),]$domains2<-cat[i]   
  i<-i+1
}
All_S_N_avS$Int_er_ra<-'x'
All_S_N_avS[All_S_N_avS$interval1 == All_S_N_avS$interval2,]$Int_er_ra<-'Intra'
All_S_N_avS[All_S_N_avS$interval1 != All_S_N_avS$interval2,]$Int_er_ra<-'Inter'

head(All_S_N_avS)

#use columns needed only
All_S_N_avS2 <- All_S_N_avS %>% unite("Domains", 5:6, na.rm = TRUE, remove = TRUE)
head(All_S_N_avS2)
All_S_N_avS3 <- All_S_N_avS2 %>% unite("Interval", 1:2, na.rm = TRUE, remove = FALSE)
head(All_S_N_avS3)

##Separate NOnub and SUMOnub in the beginning
All_S_N_avS3_NO<-All_S_N_avS3[All_S_N_avS3$sample=="WD_NOnub",]
All_S_N_avS3_SUMO<-All_S_N_avS3[All_S_N_avS3$sample=="WD_SUMOnub",]
head(All_S_N_avS3_SUMO)
head(All_S_N_avS3_NO)
#Cbind for future need
All_S_N_avS4<-merge(All_S_N_avS3_NO,All_S_N_avS3_SUMO,by="Interval")
All_S_N_avS5<-All_S_N_avS4[,c(1:7,10,11)]
All_S_N_avS5$dif<-All_S_N_avS5$av_score.y-All_S_N_avS5$av_score.x
All_S_N_avS6<-All_S_N_avS5[!duplicated(All_S_N_avS5),]
head(All_S_N_avS6)
nrow(All_S_N_avS6)/2

head(All_S_N_avS6)
print(unique(All_S_N_avS6$Domains))

states <- c("Active","Het","Null","PcG")

for(s1 in seq(from=1,by=1,to=length(states)))
{
    state1 = states[[s1]]
    
    for(s2 in seq(from=1,by=1,to=length(states)))
    {
        state2 = states[[s2]]

	if(s2 > s1)
	{
   	    print(paste0("Substitute ",paste0(state1,"_",state2)," with ",paste0(state2,"_",state1)))
	    All_S_N_avS6$Domains <- gsub(paste0(state1,"_",state2), paste0(state2,"_",state1), All_S_N_avS6$Domains)
	    print(unique(All_S_N_avS6$Domains))	    
	}
    }
}

data <- data.frame()
allData <- data.frame()
All_S_N_avS6Inter <- All_S_N_avS6[All_S_N_avS6$Int_er_ra == "Inter",]
for(statePair in unique(All_S_N_avS6Inter$Domains))
{
    print(paste0("NOnub ",statePair," ",mean(All_S_N_avS6Inter[All_S_N_avS6Inter$Domains == statePair,]$av_score.x)))

    data <- cbind(paste0(statePair,"_",unique(All_S_N_avS6Inter$sample.x)),All_S_N_avS6Inter[All_S_N_avS6Inter$Domains == statePair,]$av_score.x)	
    if(nrow(allData) == 0)
    {
        allData <- data
    } else {
        allData <- rbind(allData,data)
    }   
    #print(head(data))

    print(paste0("SUMOnub ",statePair," ",mean(All_S_N_avS6Inter[All_S_N_avS6Inter$Domains == statePair,]$av_score.y)))
    data <- cbind(paste0(statePair,"_",unique(All_S_N_avS6Inter$sample.y)),All_S_N_avS6Inter[All_S_N_avS6Inter$Domains == statePair,]$av_score.y)
    allData <- rbind(allData,data)    
    print(nrow(allData))
}
outFile <- "avgScores_trans1Dinterval_allChrom_allStatePairs_data.tsv"
write.table(allData , file=outFile,  row.names=F, sep="\t", append=F, quote=F, col.names=F)
quit()
