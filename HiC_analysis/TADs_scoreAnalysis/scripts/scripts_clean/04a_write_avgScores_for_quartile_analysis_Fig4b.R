#install.packages("ggstatsplot")
#install.packages("ggpubr",force=)

#library(ggstatsplot)
library(ggplot2)
library(ggpubr)
#library(rtracklayer)
library(tidyr)

data <- read.table("avgScores_trans1Dinterval_allChrom_all_domains_vs_all_domains.tab")
colnames(data) <- c("chr1","start1","end1","chr2","start2","end2","interval1","interval2","av_score","sample")
data$sample <- gsub("WD_","",data$sample)


data$statePair <- paste0(gsub("[0-9]+", "", data$interval1),"_",gsub("[0-9]+", "", data$interval2))
data$statePair <- gsub("Active_PcG" ,"PcG_Active" ,data$statePair)
data$statePair <- gsub("Null_PcG"   ,"PcG_Null"   ,data$statePair)
data$statePair <- gsub("Het_PcG"   ,"PcG_Het"   ,data$statePair)
data$statePair <- gsub("Null_Active","Active_Null",data$statePair)
data$statePair <- gsub("Het_Active","Active_Het",data$statePair)
data$statePair <- gsub("Null_Het","Het_Null",data$statePair)
print(head(data))

states <- c("PcG","Active","Het","Null")
samples <- c("NOnub","SUMOnub")
refSample <- "NOnub"	

for(s1 in states)
{
    state1 = s1
    
    #for(s2 in seq(from=s1,by=1,to=length(states)-1))
    for(s2 in states)
    {
        state2 = s2
	statePair <- paste0(state1,"_",state2)

	outFile <- paste0("avgScores_trans1Dinterval_all_domains_vs_all_domains_",statePair,"_data.tab")
	if(file.exists(outFile)){file.remove(outFile);}

	dataStates <- data[grep(statePair,data$statePair),c("sample","interval1","interval2","av_score")]
	if(nrow(dataStates) == 0){next;}
	print(statePair)

	dataStates$intervals <- paste0(dataStates$interval1,"_",dataStates$interval2)

	print(paste0("Keep only the trans-interval scores that are present in both the conditions"))
	dataStates <- dataStates[dataStates$interval1 != dataStates$interval2,]
	counts <- table(dataStates$intervals)
	counts <- counts[counts == 2]
	dataStates <- dataStates[dataStates$intervals %in% names(counts),]
	###DONE###

	print(paste0("Write the table as interval1 interval2 av_score_in_sample1 av_score_in_sample2"))
	dataStatesNew <- data.frame()
	for(sample in samples)
	{
	    dataStates1 <- dataStates[grep(sample,dataStates$sample),]
	    dataStates1 <- data.frame(dataStates1$interval1, dataStates1$interval2, dataStates1$av_score)
	    colnames(dataStates1) <- c("interval1","interval2",paste0("av_score_",sample))
	    if(nrow(dataStatesNew) == 0)
	    {
	        dataStatesNew <- dataStates1
		colnames(dataStatesNew) <- c("interval1","interval2",paste0("av_score_",sample))
	    } else {
	        print(names(dataStatesNew))
	        print(names(dataStates1))		
	        dataStatesNew <- merge(dataStatesNew,dataStates1,by = c("interval1","interval2"))
#		colnames(dataStatesNew) <- gsub("WD_","",colnames(dataStatesNew))
	    }
   	    print(head(dataStatesNew))	    
	}
	dataStates <- dataStatesNew
	###DONE###

	#newPcG <- c("PcG27","PcG147","PcG205","PcG106","PcG191","PcG24","PcG76","PcG160")
	#print(paste0("New PcG-PcG avg. scores ",nrow(dataStates[dataStates$interval1 %in% newPcG | dataStates$interval2 %in% newPcG,])))
	# Check worked!
	###DONE###

	print(paste0("Writing all the average scores (All)"))
	if(nrow(dataStates) > 0)
	{
	    print(statePair)
 	    print(head(dataStates))

	    for(sample in samples)
	    {
	        print(paste0("av_score_",sample))
		df <- dataStates[,c(paste0("av_score_",sample))]
		df <- data.frame(paste0(sample,"All"),df)
	        write.table(df, file = outFile, row.names = F, sep = "\t", append = T, quote = F, col.names = F)
	    }	    
        }

	print(paste0("Dividing the dataset in quartiles depending on the avgScore in ",refSample))
	qControl <- quantile(dataStates[,c(paste0("av_score_",refSample))], na.rm = TRUE)
	print(qControl)

	for(i in seq(from=1,to=4,by=1))
	{
	      quartile <- dataStates[dataStates[,c(paste0("av_score_",refSample))] >= qControl[[i]] & dataStates[,c(paste0("av_score_",refSample))] <= qControl[[i+1]],]
	      print(head(quartile))

	      for(sample in samples)
	      {
                  print(paste0("av_score_",sample))
                  df <- quartile[,c(paste0("av_score_",sample))]
                  df <- data.frame(paste0(sample,"Q",i),df)
                  write.table(df, file = outFile, row.names = F, sep = "\t", append = T, quote = F, col.names = F)
              }
	}
	allPoints <- read.table(outFile)
	colnames(allPoints) <- c("category","av_score")
	print(table(allPoints$category))
    }
}

quit()

# Now, the corrected loop:
for (pair_name in names(all_inter_dfs))
{
    current_df <- all_inter_dfs[[pair_name]]
  
    print(paste0("Processing domain pair: ", pair_name))
  
  # Skip if the data frame is empty (e.g., if a particular pair had no 'Inter' interactions)
  if (nrow(current_df) == 0) {
    print(paste0("No 'Inter' data found for ", pair_name, ". Skipping."))
    next
  }
  
  # Define output file name dynamically based on the pair_name
  outFile <- paste0("avg_scores_", pair_name, "_inter_data.tsv")
  
  # Write initial data to file (NOnubAll and SUMOnubAll)
  # Using column names `av_score.x` and `av_score.y` from your `head()` output
  write.table(cbind("NOnubAll", current_df$av_score.x), file = outFile, row.names = F, sep = "\t", append = F, quote = F, col.names = F)
  write.table(cbind("SUMOnubAll", current_df$av_score.y), file = outFile, row.names = F, sep = "\t", append = T, quote = F, col.names = F)
  
  # Check if there's enough unique data for quantile calculation
  if (length(unique(current_df$av_score.x)) < 2) {
    print(paste0("Not enough unique values for quantile calculation for NOnub in ", pair_name, ". Skipping quantiles for this pair."))
    next
  }
  
  qNOnub <- quantile(current_df$av_score.x, na.rm = TRUE)
  print(paste0("Quantiles for ", pair_name, ":"))
  print(qNOnub)
  
  # Filter for Q1 (bottom quantile)
  interQ1 <- current_df[current_df$av_score.x >= qNOnub[[1]] & current_df$av_score.x <= qNOnub[[2]],]
  print(paste0("Head of Q1 for ", pair_name))
  print(head(interQ1))
  
  # Filter for Q4 (top quantile)
  interQ4 <- current_df[current_df$av_score.x >= qNOnub[[4]] & current_df$av_score.x <= qNOnub[[5]],]
  print(paste0("Head of Q4 for ", pair_name))
  print(head(interQ4))
  
  # Write quantile data to file
  write.table(cbind("NOnubQ1", interQ1$av_score.x), file = outFile, row.names = F, sep = "\t", append = T, quote = F, col.names = F)
  write.table(cbind("SUMOnubQ1", interQ1$av_score.y), file = outFile, row.names = F, sep = "\t", append = T, quote = F, col.names = F)
  write.table(cbind("NOnubQ4", interQ4$av_score.x), file = outFile, row.names = F, sep = "\t", append = T, quote = F, col.names = F)
  write.table(cbind("SUMOnubQ4", interQ4$av_score.y), file = outFile, row.names = F, sep = "\t", append = T, quote = F, col.names = F)
  
  print(paste0("Finished processing ", pair_name))
}
quit()

quit()

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

#S_N_avS6[All_S_N_avS6$Domains == "PcG_PcG",]

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
for(statePair in unique(All_S_N_avS6$Domains))
{
    print(paste0("NOnub ",statePair," ",mean(All_S_N_avS6[All_S_N_avS6$Domains == statePair,]$av_score.x)))

    data <- cbind(paste0(statePair,"_",unique(All_S_N_avS6$sample.x)),All_S_N_avS6[All_S_N_avS6$Domains == statePair,]$av_score.x)	
    if(nrow(allData) == 0)
    {
        allData <- data
    } else {
        allData <- rbind(allData,data)
    }   
    #print(head(data))

    print(paste0("SUMOnub ",statePair," ",mean(All_S_N_avS6[All_S_N_avS6$Domains == statePair,]$av_score.y)))
    data <- cbind(paste0(statePair,"_",unique(All_S_N_avS6$sample.y)),All_S_N_avS6[All_S_N_avS6$Domains == statePair,]$av_score.y)
    allData <- rbind(allData,data)    
    print(nrow(allData))
}
quit()
outFile <- "avgScores_trans1Dinterval_allChrom_allStatePairs_data.tsv"
write.table(allData , file=outFile,  row.names=F, sep="\t", append=T, quote=F, col.names=F)
quit()


quit()

##DomainPairs##
#PcG_PcG inter or intra
SUMO_P_P       <- All_S_N_avS6[All_S_N_avS6$Domains == "PcG_PcG",]
#PcG_PcG inter
SUMO_P_P_inter <- SUMO_P_P[SUMO_P_P$Int_er_ra.x=="Inter",]
#PcG_PcG intra
SUMO_P_P_intra <- SUMO_P_P[SUMO_P_P$Int_er_ra.x=="Intra",]

outFile = "avg_scores_PcG_PcG_inter_data.tsv"
write.table(cbind("NOnubAll",SUMO_P_P_inter[,c(4)])  , file=outFile,  row.names=F, sep="\t", append=F, quote=F, col.names=F)
write.table(cbind("SUMOnubAll",SUMO_P_P_inter[,c(8)]), file=outFile,  row.names=F, sep="\t", append=T, quote=F, col.names=F)
print(head(SUMO_P_P_inter))
qNOnub <- quantile(SUMO_P_P_inter[,c(4)])
print(paste0("Bottom quantile"))
P_P_interQ1 <- SUMO_P_P_inter[qNOnub[[1]] <= SUMO_P_P_inter[,c(4)] & SUMO_P_P_inter[,c(4)] <= qNOnub[[2]],]
print(head(P_P_interQ1))
#SUMO_P_P_inter[qNOnub[[1]] <= SUMO_P_P_inter[,c(4)] & SUMO_P_P_inter[,c(4)] <= qNOnub[[2]],]$qNOnub <- "Q1"
#SUMO_P_P_inter[qNOnub[[2]] <= SUMO_P_P_inter[,c(4)] & SUMO_P_P_inter[,c(4)] <= qNOnub[[3]],]$qNOnub <- "Q2"
#SUMO_P_P_inter[qNOnub[[3]] <= SUMO_P_P_inter[,c(4)] & SUMO_P_P_inter[,c(4)] <= qNOnub[[4]],]$qNOnub <- "Q3"
#SUMO_P_P_inter[qNOnub[[4]] <= SUMO_P_P_inter[,c(4)] & SUMO_P_P_inter[,c(4)] <= qNOnub[[5]],]$qNOnub <- "Q4"
print(paste0("Top quantile"))
P_P_interQ4 <- SUMO_P_P_inter[qNOnub[[4]] <= SUMO_P_P_inter[,c(4)] & SUMO_P_P_inter[,c(4)] <= qNOnub[[5]],]
print(head(SUMO_P_P_inter))

write.table(cbind("NOnubQ1"  ,P_P_interQ1[,c(4)]), file=outFile,  row.names=F, sep="\t", append=T, quote=F, col.names=F)
write.table(cbind("SUMOnubQ1",P_P_interQ1[,c(8)]), file=outFile,  row.names=F, sep="\t", append=T, quote=F, col.names=F)
write.table(cbind("NOnubQ4"  ,P_P_interQ4[,c(4)]), file=outFile,  row.names=F, sep="\t", append=T, quote=F, col.names=F)
write.table(cbind("SUMOnubQ4",P_P_interQ4[,c(8)]), file=outFile,  row.names=F, sep="\t", append=T, quote=F, col.names=F)

SUMO_P_A       <- All_S_N_avS6[All_S_N_avS6$Domains == "PcG_Active" | All_S_N_avS6$Domains == "Active_PcG",]
SUMO_P_A_inter <- SUMO_P_A[SUMO_P_A$Int_er_ra.x=="Inter",]

SUMO_P_N       <- All_S_N_avS6[All_S_N_avS6$Domains == "PcG_Null"   | All_S_N_avS6$Domains == "Null_PcG",]
SUMO_P_N_inter <- SUMO_P_N[SUMO_P_N$Int_er_ra.x=="Inter",]

SUMO_P_H       <- All_S_N_avS6[All_S_N_avS6$Domains == "PcG_Het"    | All_S_N_avS6$Domains == "Het_PcG",]
SUMO_P_H_inter <- SUMO_P_H[SUMO_P_H$Int_er_ra.x=="Inter",]

# Create a NAMED list of the data frames you want to loop through
# The names of the list elements will be used for file naming
all_inter_dfs <- list(
  "PcG_PcG" = SUMO_P_P_inter,
  "PcG_Active" = SUMO_P_A_inter, # Naming this consistently
  "PcG_Null" = SUMO_P_N_inter,   # Naming this consistently
  "PcG_Het" = SUMO_P_H_inter     # Naming this consistently
)



























  
  
  









print(head(P_P_interQ4))
print(paste0("Quantiles"))
print(qNOnub)

quit()

SUMO_P_P_inter2<-SUMO_P_P_inter[,c(4,8)]
colnames(SUMO_P_P_inter2)<-c("Control","SUMO")  
#pheatmap(SUMO_P_P_inter2)
head(SUMO_P_P_inter2)

#library()


l_SUMO_P_P_inter2 <- SUMO_P_P_inter2 %>%
  rownames_to_column("Interval") %>%  # If row names are Intervals, convert to a column
  pivot_longer(cols = c(Control, SUMO), names_to = "Condition", values_to = "Value")

#ggbetweenstats(
#  data  = l_SUMO_P_P_inter2,
#  x     = Condition,
#  y     = Value ,
#  title = "Inter PcG in Control vs SUMO all "
#)
int<-median(SUMO_P_P_inter2$SUMO)
allPP<-ggplot(l_SUMO_P_P_inter2, aes(x = Condition, y = Value, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.5,show.legend = FALSE) + # Violin plot with transparency
  theme_minimal() + # Clean theme
  labs(title = "",
       x="",
       y = "PcG-PcG") +
  scale_fill_manual(values = c("Control" = "darkgreen", "SUMO" = "orange"))+ # Custom colors
  stat_compare_means(method = "t.test", label = "p.format",aes(x=1.5),vjust = 0,size=3)+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+theme(plot.title = element_text(size = 10,hjust = 0.5)) +
  geom_hline(yintercept = int, linetype = "dashed", color = "lightgrey")+
  geom_boxplot(width = 0.1, outlier.shape = NA,show.legend = FALSE)  # Small box plot
medianControlP<-median(SUMO_P_P_inter2$Control)
quantile(SUMO_P_P_inter2$Control)

intLow<-median(lowC_P_P$SUMO)

lowC_P_P<-SUMO_P_P_inter2[SUMO_P_P_inter2$Control <0,]


l_lowC_P_P <- lowC_P_P %>%
  rownames_to_column("Interval") %>%  # If row names are Intervals, convert to a column
  pivot_longer(cols = c(Control, SUMO), names_to = "Condition", values_to = "Value")

HighC_P_P<-SUMO_P_P_inter2[SUMO_P_P_inter2$Control >20,]

intHigh<-median(HighC_P_P$SUMO)

l_HighC_P_P <- HighC_P_P %>%
  rownames_to_column("Interval") %>%  # If row names are Intervals, convert to a column
  pivot_longer(cols = c(Control, SUMO), names_to = "Condition", values_to = "Value")

Low_PP<-ggplot(l_lowC_P_P, aes(x = Condition, y = Value, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.5,show.legend = FALSE) + # Violin plot with transparency
  theme_minimal() + # Clean theme
  labs(title = "",
       x="",
       y = "PcG-PcG") +
  scale_fill_manual(values = c("Control" = "darkgreen", "SUMO" = "orange"))+ # Custom colors
  stat_compare_means(method = "t.test", label = "p.format",aes(x=1.5),vjust = 0,size=3)+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+theme(plot.title = element_text(size = 10,hjust = 0.5)) +
  geom_hline(yintercept = intLow, linetype = "dashed", color = "lightgrey")+
  geom_boxplot(width = 0.1, outlier.shape = NA,show.legend = FALSE)  # Small box plot


High_PP<-ggplot(l_HighC_P_P, aes(x = Condition, y = Value, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.5,show.legend = FALSE) + # Violin plot with transparency
  theme_minimal() + # Clean theme
  labs(title = "",
       x="",
       y = "PcG-PcG") +
  scale_fill_manual(values = c("Control" = "darkgreen", "SUMO" = "orange"))+ # Custom colors
  stat_compare_means(method = "t.test", label = "p.format",aes(x=1.5),vjust = 0,size=3)+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+theme(plot.title = element_text(size = 10,hjust = 0.5)) +
  geom_hline(yintercept = intHigh, linetype = "dashed", color = "lightgrey")+
  geom_boxplot(width = 0.1, outlier.shape = NA,show.legend = FALSE)  # Small box plot


combinesPP<-allPP + Low_PP + High_PP

#Active_PcG #PcG_Active

SUMO_P_A<-All_S_N_avS6[All_S_N_avS6$Domains=="PcG_Active"| All_S_N_avS6$Domains=="Active_PcG",]
#Active_PcG inter
SUMO_P_A_inter<-SUMO_P_A[SUMO_P_A$Int_er_ra.x=="Inter",]
rownames(SUMO_P_A_inter)<-SUMO_P_A_inter$Interval
#PcG_PcG intra

SUMO_P_A_inter2<-SUMO_P_A_inter[,c(4,8)]
colnames(SUMO_P_A_inter2)<-c("Control","SUMO")  
#pheatmap(SUMO_P_A_inter2)
head(SUMO_P_A_inter2)

l_SUMO_P_A_inter2 <- SUMO_P_A_inter2 %>%
  rownames_to_column("Interval") %>%  # If row names are Intervals, convert to a column
  pivot_longer(cols = c(Control, SUMO), names_to = "Condition", values_to = "Value")


intAllA<-median(SUMO_P_A_inter2$SUMO)
print(l_SUMO_P_A_inter2)
#quit()

#All_PA<-ggplot(l_SUMO_P_A_inter2, aes(x = Condition, y = Value, fill = Condition)) +
#  geom_violin(trim = FALSE, alpha = 0.5,show.legend = FALSE) + # Violin plot with transparency
#  theme_minimal() + # Clean theme
#  labs(title = "",
#       x="",
#       y = "PcG-Active") +
#  scale_fill_manual(values = c("Control" = "darkgreen", "SUMO" = "orange"))+ # Custom colors
#  #stat_compare_means(method = "t.test", label = "p.format",aes(x=1.5),vjust = 0,size=3)+
#  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+theme(plot.title = element_text(size = 10,hjust = 0.5)) +
#  geom_hline(yintercept = intAllA, linetype = "dashed", color = "lightgrey")+
#  geom_boxplot(width = 0.1, outlier.shape = NA,show.legend = FALSE)  # Small box plot


#ggbetweenstats(
#  data  = l_SUMO_P_A_inter2,
#  x     = Condition,
#  y     = Value ,
#  title = "All PcG to Active in Control vs SUMO all "
#)


q <- quantile(SUMO_P_A_inter2$Control)
print(q)
quit()

lowC_P_A<-SUMO_P_A_inter2[SUMO_P_A_inter2$Control < -6,]


intLowA<-median(lowC_P_A$SUMO)

l_intLowA<- lowC_P_A %>%
  rownames_to_column("Interval") %>%  # If row names are Intervals, convert to a column
  pivot_longer(cols = c(Control, SUMO), names_to = "Condition", values_to = "Value")

Low_PA<-ggplot(l_intLowA, aes(x = Condition, y = Value, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.5,show.legend = FALSE) + # Violin plot with transparency
  theme_minimal() + # Clean theme
  labs(title = "",
       x="",
       y = "PcG-Active") +
  scale_fill_manual(values = c("Control" = "darkgreen", "SUMO" = "orange"))+ # Custom colors
  stat_compare_means(method = "t.test", label = "p.format",aes(x=1.5),vjust = 0,size=3)+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+theme(plot.title = element_text(size = 10,hjust = 0.5)) +
  geom_hline(yintercept = intLowA, linetype = "dashed", color = "lightgrey")+
  geom_boxplot(width = 0.1, outlier.shape = NA,show.legend = FALSE)  # Small box plot



HighC_P_A<-SUMO_P_A_inter2[SUMO_P_A_inter2$Control >10,]

intHighA<-median(HighC_P_A$SUMO)

l_HighC_P_A <- HighC_P_A %>%
  rownames_to_column("Interval") %>%  # If row names are Intervals, convert to a column
  pivot_longer(cols = c(Control, SUMO), names_to = "Condition", values_to = "Value")

High_PA<-ggplot(l_HighC_P_A, aes(x = Condition, y = Value, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.5,show.legend = FALSE) + # Violin plot with transparency
  theme_minimal() + # Clean theme
  labs(title = "",
       x="",
       y = "PcG-Active") +
  scale_fill_manual(values = c("Control" = "darkgreen", "SUMO" = "orange"))+ # Custom colors
  stat_compare_means(method = "t.test", label = "p.format",aes(x=1.5),vjust = 0,size=3)+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+theme(plot.title = element_text(size = 10,hjust = 0.5)) +
  geom_hline(yintercept = intHighA, linetype = "dashed", color = "lightgrey")+
  geom_boxplot(width = 0.1, outlier.shape = NA,show.legend = FALSE)  # Small box plot




#Null_PcG #PcG_Null

SUMO_P_N<-All_S_N_avS6[All_S_N_avS6$Domains=="PcG_Null"| All_S_N_avS6$Domains=="Null_PcG",]
#Active_PcG inter
SUMO_P_N_inter<-SUMO_P_N[SUMO_P_N$Int_er_ra.x=="Inter",]
rownames(SUMO_P_N_inter)<-SUMO_P_N_inter$Interval
#PcG_PcG intra

SUMO_P_N_inter2<-SUMO_P_N_inter[,c(4,8)]
colnames(SUMO_P_N_inter2)<-c("Control","SUMO")  
#pheatmap(SUMO_P_A_inter2)
head(SUMO_P_N_inter2)

l_SUMO_P_N_inter2 <- SUMO_P_N_inter2 %>%
  rownames_to_column("Interval") %>%  # If row names are Intervals, convert to a column
  pivot_longer(cols = c(Control, SUMO), names_to = "Condition", values_to = "Value")

intAllN<-median(SUMO_P_N_inter2$SUMO)



All_PN<-ggplot(l_SUMO_P_N_inter2, aes(x = Condition, y = Value, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.5,show.legend = FALSE) + # Violin plot with transparency
  theme_minimal() + # Clean theme
  labs(title = "",
       x="",
       y = "PcG-Null") +
  scale_fill_manual(values = c("Control" = "darkgreen", "SUMO" = "orange"))+ # Custom colors
  stat_compare_means(method = "t.test", label = "p.format",aes(x=1.5),vjust = 0,size=3)+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+theme(plot.title = element_text(size = 10,hjust = 0.5)) +
  geom_hline(yintercept = intAllN, linetype = "dashed", color = "lightgrey")+
  geom_boxplot(width = 0.1, outlier.shape = NA,show.legend = FALSE)  # Small box plot






quantile(SUMO_P_N_inter2$Control)


lowC_P_N<-SUMO_P_N_inter2[SUMO_P_N_inter2$Control < -8,]


intLowN<-median(lowC_P_N$SUMO)

l_intLowN<- lowC_P_N %>%
  rownames_to_column("Interval") %>%  # If row names are Intervals, convert to a column
  pivot_longer(cols = c(Control, SUMO), names_to = "Condition", values_to = "Value")

Low_PN<-ggplot(l_intLowN, aes(x = Condition, y = Value, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.5,show.legend = FALSE) + # Violin plot with transparency
  theme_minimal() + # Clean theme
  labs(title = "",
       x="",
       y = "PcG-Null") +
  scale_fill_manual(values = c("Control" = "darkgreen", "SUMO" = "orange"))+ # Custom colors
  stat_compare_means(method = "t.test", label = "p.format",aes(x=1.5),vjust = 0,size=3)+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+theme(plot.title = element_text(size = 10,hjust = 0.5)) +
  geom_hline(yintercept = intLowN, linetype = "dashed", color = "lightgrey")+
  geom_boxplot(width = 0.1, outlier.shape = NA,show.legend = FALSE)  # Small box plot



HighC_P_N<-SUMO_P_N_inter2[SUMO_P_N_inter2$Control >10,]

intHighN<-median(HighC_P_N$SUMO)

l_HighC_P_N <- HighC_P_N %>%
  rownames_to_column("Interval") %>%  # If row names are Intervals, convert to a column
  pivot_longer(cols = c(Control, SUMO), names_to = "Condition", values_to = "Value")

High_PN<-ggplot(l_HighC_P_N, aes(x = Condition, y = Value, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.5,show.legend = FALSE) + # Violin plot with transparency
  theme_minimal() + # Clean theme
  labs(title = "",
       x="",
       y = "PcG-Null") +
  scale_fill_manual(values = c("Control" = "darkgreen", "SUMO" = "orange"))+ # Custom colors
  stat_compare_means(method = "t.test", label = "p.format",aes(x=1.5),vjust = 0,size=3)+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+theme(plot.title = element_text(size = 10,hjust = 0.5)) +
  geom_hline(yintercept = intHighN, linetype = "dashed", color = "lightgrey")+
  geom_boxplot(width = 0.1, outlier.shape = NA,show.legend = FALSE)  # Small box plot






#Het_PcG #PcG_Het

SUMO_P_H<-All_S_N_avS6[All_S_N_avS6$Domains=="PcG_Het"| All_S_N_avS6$Domains=="Het_PcG",]
#Active_PcG inter
SUMO_P_H_inter<-SUMO_P_H[SUMO_P_H$Int_er_ra.x=="Inter",]
rownames(SUMO_P_H_inter)<-SUMO_P_H_inter$Interval
#PcG_PcG intra

SUMO_P_H_inter2<-SUMO_P_H_inter[,c(4,8)]
colnames(SUMO_P_H_inter2)<-c("Control","SUMO")  
#pheatmap(SUMO_P_A_inter2)
head(SUMO_P_H_inter2)
SUMO_P_H_inter2<-na.omit(SUMO_P_H_inter2)

l_SUMO_P_H_inter2 <- SUMO_P_H_inter2 %>%
  rownames_to_column("Interval") %>%  # If row names are Intervals, convert to a column
  pivot_longer(cols = c(Control, SUMO), names_to = "Condition", values_to = "Value")

intAllH<-median(SUMO_P_H_inter2$SUMO)
is.numeric(SUMO_P_H_inter2$SUMO)


All_PH<-ggplot(l_SUMO_P_H_inter2, aes(x = Condition, y = Value, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.5,show.legend = FALSE) + # Violin plot with transparency
  theme_minimal() + # Clean theme
  labs(title = "",
       x="",
       y = "PcG-Het") +
  scale_fill_manual(values = c("Control" = "darkgreen", "SUMO" = "orange"))+ # Custom colors
  stat_compare_means(method = "t.test", label = "p.format",aes(x=1.5),vjust = 0,size=3)+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+theme(plot.title = element_text(size = 10,hjust = 0.5)) +
  geom_hline(yintercept = intAllH, linetype = "dashed", color = "lightgrey")+
  geom_boxplot(width = 0.1, outlier.shape = NA,show.legend = FALSE)  # Small box plot






quantile(SUMO_P_H_inter2$Control)


lowC_P_H<-SUMO_P_H_inter2[SUMO_P_H_inter2$Control < -7,]


intLowH<-median(lowC_P_H$SUMO)

l_intLowH<- lowC_P_H %>%
  rownames_to_column("Interval") %>%  # If row names are Intervals, convert to a column
  pivot_longer(cols = c(Control, SUMO), names_to = "Condition", values_to = "Value")

Low_PH<-ggplot(l_intLowH, aes(x = Condition, y = Value, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.5,show.legend = FALSE) + # Violin plot with transparency
  theme_minimal() + # Clean theme
  labs(title = "",
       x="",
       y = "PcG-Het") +
  scale_fill_manual(values = c("Control" = "darkgreen", "SUMO" = "orange"))+ # Custom colors
  stat_compare_means(method = "t.test", label = "p.format",aes(x=1.5),vjust = 0,size=3)+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+theme(plot.title = element_text(size = 10,hjust = 0.5)) +
  geom_hline(yintercept = intLowH, linetype = "dashed", color = "lightgrey")+
  geom_boxplot(width = 0.1, outlier.shape = NA,show.legend = FALSE)  # Small box plot



HighC_P_H<-SUMO_P_H_inter2[SUMO_P_H_inter2$Control >10,]

intHighH<-median(HighC_P_H$SUMO)

l_HighC_P_H<- HighC_P_H %>%
  rownames_to_column("Interval") %>%  # If row names are Intervals, convert to a column
  pivot_longer(cols = c(Control, SUMO), names_to = "Condition", values_to = "Value")

High_PH<-ggplot(l_HighC_P_H, aes(x = Condition, y = Value, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.5,show.legend = FALSE) + # Violin plot with transparency
  theme_minimal() + # Clean theme
  labs(title = "",
       x="",
       y = "PcG-Het") +
  scale_fill_manual(values = c("Control" = "darkgreen", "SUMO" = "orange"))+ # Custom colors
  stat_compare_means(method = "t.test", label = "p.format",aes(x=1.5),vjust = 0,size=3)+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+theme(plot.title = element_text(size = 10,hjust = 0.5)) +
  geom_hline(yintercept = intHighH, linetype = "dashed", color = "lightgrey")+
  geom_boxplot(width = 0.1, outlier.shape = NA,show.legend = FALSE)  # Small box plot



library(cowplot)

combined <- plot_grid(allPP ,Low_PP,High_PP,All_PA,Low_PA,High_PA,All_PN,Low_PN,High_PN,All_PH,Low_PH,High_PH
                      , ncol = 3)

# Save
ggsave("PcG_to_All.pdf", combined, width = 10, height = 7, dpi = 350)

# Print in R
#combined


row1 <- plot_grid(allPP, Low_PP, High_PP, ncol = 3)
row2 <- plot_grid(All_PA, Low_PA, High_PA, ncol = 3)
row3 <- plot_grid(All_PN, Low_PN, High_PN, ncol = 3)
row4 <- plot_grid(All_PH, Low_PH, High_PH, ncol = 3)

library(cowplot)
library(grid)


color_panel <- function(plot_row, fill = "grey90") {
  ggdraw() +
    # Draw a full background rectangle
    draw_plot(
      ggplot() + 
        theme_void() + 
        theme(plot.background = element_rect(fill = fill, color = NA))
    ) +
    # Overlay the actual row of plots
    draw_plot(plot_row)
}



frame_row <- function(plot_row, color = "black", size = 2) {
  ggdraw() +
    draw_grob(
      rectGrob(gp = gpar(fill = NA, col = color, lwd = size))
    ) +
    draw_plot(plot_row)
}

row1 <- plot_grid(allPP, Low_PP, High_PP, ncol = 3)
row2 <- plot_grid(All_PA, Low_PA, High_PA, ncol = 3)
row3 <- plot_grid(All_PN, Low_PN, High_PN, ncol = 3)
row4 <- plot_grid(All_PH, Low_PH, High_PH, ncol = 3)




row1_framed <- frame_row(row1, color = "royalblue", size = 3)
row2_framed <- frame_row(row2, color = "red", size = 3)
row3_framed <- frame_row(row3, color = "lightgreen", size = 3)
row4_framed <- frame_row(row4, color = "black", size = 3)


#row1_colored <- plot_grid(row1, fill = "lightblue")
#row2_colored <- plot_grid(row2, fill = "pink")
#row3_colored <- plot_grid(row3, fill = "grey")
#row4_colored <- plot_grid(row4, fill = "lightgreen")



combined <- plot_grid(
  row1_framed,
  row2_framed,
  row3_framed,
  row4_framed,
  ncol = 1,
  rel_heights = c(1, 1, 1, 1)  # adjust row heights if needed
)



ggsave("PcG_to_All.pdf", combined, width = 10, height = 7.5, dpi = 350)
combined






#####Check PcG intras

#PcG_PcG intra
SUMO_P_P_intra<-SUMO_P_P[SUMO_P_P$Int_er_ra.x=="Intra",]

SUMO_P_P_intra2<-SUMO_P_P_intra[,c(4,8)]


colnames(SUMO_P_P_intra2)<-c("Control","SUMO") 

l_SUMO_P_P_intra2 <- SUMO_P_P_intra2 %>%
  rownames_to_column("Interval") %>%  # If row names are Intervals, convert to a column
  pivot_longer(cols = c(Control, SUMO), names_to = "Condition", values_to = "Value")

intIntraPP<-median(SUMO_P_P_intra2$Control)

P_P_intra<-ggplot(l_SUMO_P_P_intra2, aes(x = Condition, y = Value, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.5,show.legend = FALSE) + # Violin plot with transparency
  theme_minimal() + # Clean theme
  labs(title = "",
       x="",
       y = "PcG_intra") +
  scale_fill_manual(values = c("Control" = "darkgreen", "SUMO" = "orange"))+ # Custom colors
  stat_compare_means(method = "t.test", label = "p.format",aes(x=1.5),vjust = 0,size=3)+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+theme(plot.title = element_text(size = 10,hjust = 0.5)) +
  geom_hline(yintercept = intIntraPP, linetype = "dashed", color = "lightgrey")+
  geom_boxplot(width = 0.1, outlier.shape = NA,show.legend = FALSE)  # Small box plot
medianControlP<-median(SUMO_P_P_inter2$Control)
quantile(SUMO_P_P_inter2$Control)

lowIntraPP<-SUMO_P_P_intra2[SUMO_P_P_intra2$Control < 3.5,]

l_lowIntraPP <- lowIntraPP %>%
  rownames_to_column("Interval") %>%  # If row names are Intervals, convert to a column
  pivot_longer(cols = c(Control, SUMO), names_to = "Condition", values_to = "Value")

lowintIntraPP<-median(lowIntraPP$Control)

low_P_P_intra<-ggplot(l_lowIntraPP, aes(x = Condition, y = Value, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.5,show.legend = FALSE) + # Violin plot with transparency
  theme_minimal() + # Clean theme
  labs(title = "",
       x="",
       y = "PcG_intra") +
  scale_fill_manual(values = c("Control" = "darkgreen", "SUMO" = "orange"))+ # Custom colors
  stat_compare_means(method = "t.test", label = "p.format",aes(x=1.5),vjust = 0,size=3)+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+theme(plot.title = element_text(size = 10,hjust = 0.5)) +
  geom_hline(yintercept = lowintIntraPP, linetype = "dashed", color = "lightgrey")+
  geom_boxplot(width = 0.1, outlier.shape = NA,show.legend = FALSE)  # Small box plot




HighIntraPP<-SUMO_P_P_intra2[SUMO_P_P_intra2$Control > 20,]

l_highIntraPP <- HighIntraPP %>%
  rownames_to_column("Interval") %>%  # If row names are Intervals, convert to a column
  pivot_longer(cols = c(Control, SUMO), names_to = "Condition", values_to = "Value")

highintIntraPP<-median(HighIntraPP$Control)

high_P_P_intra<-ggplot(l_highIntraPP, aes(x = Condition, y = Value, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.5,show.legend = FALSE) + # Violin plot with transparency
  theme_minimal() + # Clean theme
  labs(title = "",
       x="",
       y = "PcG_intra") +
  scale_fill_manual(values = c("Control" = "darkgreen", "SUMO" = "orange"))+ # Custom colors
  stat_compare_means(method = "t.test", label = "p.format",aes(x=1.5),vjust = 0,size=3)+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+theme(plot.title = element_text(size = 10,hjust = 0.5)) +
  geom_hline(yintercept = highintIntraPP, linetype = "dashed", color = "lightgrey")+
  geom_boxplot(width = 0.1, outlier.shape = NA,show.legend = FALSE)  # Small box plot




library(cowplot)

combined_Intra_PP <- plot_grid(P_P_intra,low_P_P_intra,high_P_P_intra, ncol = 3)
getwd()
# Save
ggsave("combined_Intra_PP.pdf", combined_Intra_PP, width = 10, height = 4, dpi = 350)

# Print in R
combined



########

RNA_Domains[RNA_Domains$GeneID == "Antp",]
RNA_Domains[RNA_Domains$GeneID == "Ubx",]
RNA_Domains[RNA_Domains$GeneID == "vg",]


#Where are the downregulated genes
AllRNAdoms<-RNA_Domains$Interval
RNA_Domains_PcG_all<-AllRNAdoms[grep("PcG",AllRNAdoms)]


RNA_Domains2<-RNA_Domains[RNA_Domains$log2Ratio< -1.5,]
RNA_Domains2 <- RNA_Domains2[!is.na(RNA_Domains2$padjust) & RNA_Domains2$padjust <= 0.05, ]
head(RNA_Domains2)
nrow(RNA_Domains2)

RNA_Domains4<-RNA_Domains[RNA_Domains$log2Ratio > 1.5,]
RNA_Domains4 <- RNA_Domains4[!is.na(RNA_Domains4$padjust) & RNA_Domains4$padjust <= 0.05, ]


#Below are the domains in which down regulated genes are found

RNA_Domains3<-RNA_Domains2$Interval
RNA_Domains3_Actives_down<-RNA_Domains3[grep("Active",RNA_Domains3)]
RNA_Domains3_Nul_down<-RNA_Domains3[grep("Null",RNA_Domains3)]
RNA_Domains3_Het_down<-RNA_Domains3[grep("Het",RNA_Domains3)]
RNA_Domains3_PcG_downs<-RNA_Domains3[grep("PcG",RNA_Domains3)]



RNA_Domains3_Actives_down<-RNA_Domains3_Actives_down[!duplicated(RNA_Domains3_Actives_down)]
RNA_Domains3_Nul_down<-RNA_Domains3_Nul_down[!duplicated(RNA_Domains3_Nul_down)]
RNA_Domains3_Het_down<-RNA_Domains3_Het_down[!duplicated(RNA_Domains3_Het_down)]
RNA_Domains3_PcG_downs<-RNA_Domains3_PcG_downs[!duplicated(RNA_Domains3_PcG_downs)]

#Below are the domains in which up regulated genes are found

RNA_Domains_ups<-RNA_Domains4$Interval
RNA_Domains3_Actives_ups<-RNA_Domains_ups[grep("Active",RNA_Domains_ups)]
RNA_Domains3_Null_ups<-RNA_Domains_ups[grep("Null",RNA_Domains_ups)]
RNA_Domains3_Het_ups<-RNA_Domains_ups[grep("Het",RNA_Domains_ups)]
RNA_Domains3_PcG_ups<-RNA_Domains_ups[grep("PcG",RNA_Domains_ups)]



RNA_Domains3_Actives_ups<-RNA_Domains3_Actives_ups[!duplicated(RNA_Domains3_Actives_ups)]
RNA_Domains3_Null_ups<-RNA_Domains3_Null_ups[!duplicated(RNA_Domains3_Null_ups)]
RNA_Domains3_Het_ups<-RNA_Domains3_Het_ups[!duplicated(RNA_Domains3_Het_ups)]
RNA_Domains3_PcG_ups<-RNA_Domains3_PcG_ups[!duplicated(RNA_Domains3_PcG_ups)]


common <- intersect(RNA_Domains3_PcG_ups, RNA_Domains3_PcG_downs)
common


onlyDownPcG <- setdiff(RNA_Domains3_PcG_downs, common)
onlyUpPcG<-setdiff(RNA_Domains3_PcG_ups, common)




onlyDownPcG_scores<-All_S_N_avS6[All_S_N_avS6$interval2.x %in% onlyDownPcG,]
onlyDownPcG_scores.1<-onlyDownPcG_scores[,c(2,3,10)]

onlyDownPcG_scores2<-All_S_N_avS6[All_S_N_avS6$interval1.x %in% onlyDownPcG,]
onlyDownPcG_scores2.1<-onlyDownPcG_scores2[,c(2,3,10)]

colnames(onlyDownPcG_scores.1)<-c("Active","PcG","Dif")
colnames(onlyDownPcG_scores2.1)<-c("Active","PcG","Dif")

#Downs_A_H<-rbind(Active_Het_Down1.1,Active_Het_Down2.1)

ggplot(Downs_A_H, aes(x = Het, y = Active, fill = Dif)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Heatmap of Dif Values", x = "Active", y = "Het", fill = "Dif") +
  theme_minimal()










highActive<-SUMO_P_A_inter[SUMO_P_A_inter$dif>30,]
genes_in_high_A_P<-RNA_Domains[RNA_Domains$Interval %in% highActive$interval1.x | RNA_Domains$Interval %in% highActive$interval2.x ,]
head(genes_in_high_A_P)
genes_in_high_A_P<-genes_in_high_A_P[,c(2:ncol(genes_in_high_A_P))]

EnhancedVolcano(genes_in_high_A_P,lab = genes_in_high_A_P$GeneID,
                x = 'log2Ratio',
                y = 'pvalue',pCutoff = 5*10e-2,
                FCcutoff = 1.5)
Genes_Actives<-genes_in_high_A_P[grep("Active",genes_in_high_A_P$Interval),]

EnhancedVolcano(Genes_Actives,lab = Genes_Actives$GeneID,
                x = 'log2Ratio',
                y = 'pvalue',pCutoff = 5*10e-2,
                FCcutoff = 1.5)





#All genes under all Active domains
RNA_Domains11<-RNA_Domains$Interval
RNA_Domains12_Actives<-RNA_Domains11[grep("Active",RNA_Domains11)]

#All inter Actives
SUMO_P_A<-All_S_N_avS6[All_S_N_avS6$Domains=="PcG_Active"| All_S_N_avS6$Domains=="Active_PcG",]
SUMO_H_A<-All_S_N_avS6[All_S_N_avS6$Domains=="Het_Active"| All_S_N_avS6$Domains=="Active_Het",]
SUMO_N_A<-All_S_N_avS6[All_S_N_avS6$Domains=="Null_Active"| All_S_N_avS6$Domains=="Active_Null",]
Active_to_nonA<-rbind(SUMO_P_A,SUMO_H_A,SUMO_N_A)

head(Active_to_nonA)
#All non Actives gaining active
SUMO_N_A_gain<- Active_to_nonA[Active_to_nonA$dif> 0,]
head(SUMO_N_A_gain)
#All genes under Active which gain interaction with others
x<-RNA_Domains12_Actives[RNA_Domains12_Actives %in% SUMO_N_A_gain$interval1.x |RNA_Domains12_Actives %in% SUMO_N_A_gain$interval2.x]
t<-RNA_Domains[RNA_Domains$Interval %in% x,]
t2<-t[,2:ncol(t)]

EnhancedVolcano(t2,lab = t2$GeneID,
                x = 'log2Ratio',
                y = 'pvalue',pCutoff = 5*10e-2,
                FCcutoff = 1.5)



RNA_Domains3_PcG_ups<-RNA_Domains4.2[grep("PcG",RNA_Domains4.2)]

RNA_Domains3_PcG_ups<-unique(RNA_Domains3_PcG_ups)
#RNA_Domains3_PcG_ups, RNA_Domains3_PcG_downs


RNA_Domains3_PcG_downs<-unique(RNA_Domains3_PcG_downs)
RNA_Domains3_PcG_ups<-unique(RNA_Domains3_PcG_ups)



RNA_Domains3_PcG_check_d<-SUMO_P_P_inter[SUMO_P_P_inter$interval2.x %in% RNA_Domains3_PcG_downs,]

RNA_Domains3_PcG_check_d<-RNA_Domains3_PcG_check_d[,c(4,8)]
colnames(RNA_Domains3_PcG_check_d)<-c("Control","SUMO")  

RNA_Domains3_PcG_londowns <- RNA_Domains3_PcG_check_d %>%
  rownames_to_column("Interval") %>%  # If row names are Intervals, convert to a column
  pivot_longer(cols = c(Control, SUMO), names_to = "Condition", values_to = "Value")



ggplot(RNA_Domains3_PcG_londowns, aes(x = Condition, y = Value, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.5,show.legend = FALSE) + # Violin plot with transparency
  theme_minimal() + # Clean theme
  labs(title = "",
       x="",
       y = "PcG-PcG") +
  scale_fill_manual(values = c("Control" = "darkgreen", "SUMO" = "orange"))+ # Custom colors
  stat_compare_means(method = "t.test", label = "p.format",aes(x=1.5),vjust = 0,size=3)







###Compare the PcG-PcG contac change of down PcG genes

RNA_Domains3_PcG_check<-SUMO_P_P_inter[SUMO_P_P_inter$interval2.x %in% RNA_Domains3_PcG_downs,]

RNA_Domains3_PcG_check<-RNA_Domains3_PcG_check[,c(4,8)]
colnames(RNA_Domains3_PcG_check)<-c("Control","SUMO")  

RNA_Domains3_PcG_long <- RNA_Domains3_PcG_check %>%
  rownames_to_column("Interval") %>%  # If row names are Intervals, convert to a column
  pivot_longer(cols = c(Control, SUMO), names_to = "Condition", values_to = "Value")



ggplot(RNA_Domains3_PcG_long, aes(x = Condition, y = Value, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.5,show.legend = FALSE) + # Violin plot with transparency
  theme_minimal() + # Clean theme
  labs(title = "",
       x="",
       y = "PcG-PcG") +
  scale_fill_manual(values = c("Control" = "darkgreen", "SUMO" = "orange"))+ # Custom colors
  stat_compare_means(method = "t.test", label = "p.format",aes(x=1.5),vjust = 0,size=3)+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+theme(plot.title = element_text(size = 10,hjust = 0.5)) +
  geom_hline(yintercept = int, linetype = "dashed", color = "lightgrey")+
  geom_boxplot(width = 0.1, outlier.shape = NA,show.legend = FALSE)  # Small box plot

####Compare upPcGdifs to downPcGdifs

#RNA_Domains3_PcG_ups, RNA_Domains3_PcG_downs
RNA_Domains3_PcG_ups_dif<-SUMO_P_P_inter[SUMO_P_P_interRNA_Domains3_PcG_ups
SUMO_P_P_inter

#How many down genes under Active
nrow(RNA_Domains2[RNA_Domains2$Interval %in% RNA_Domains3_Actives,])
nrow(RNA_Domains2[RNA_Domains2$Interval %in% RNA_Domains3_Null,])
nrow(RNA_Domains2[RNA_Domains2$Interval %in% RNA_Domains3_Het,])
nrow(RNA_Domains2[RNA_Domains2$Interval %in% RNA_Domains3_PcG,])











#DownGeneDomainsAll<-c(RNA_Domains3_Actives,RNA_Domains3_Null,RNA_Domains3_Het,RNA_Domains3_PcG)

length(RNA_Domains3_Actives)
length(RNA_Domains3_Null)
length(RNA_Domains3_Het)
length(RNA_Domains3_PcG)

##Are the Null and Active Domains in which the Down genes found are the ones with ectopic contacts with the 
#repressıve domains?
head(All_S_N_avS6)
All_S_N_avS7<-All_S_N_avS6[All_S_N_avS6$dif>0,]#Any domain with incr. contact in SUMO
#Increased PcG Actives
#SUMO_P_A_inc<-All_S_N_avS7[All_S_N_avS7$Domains=="Active_PcG"|All_S_N_avS7$Domains=="PcG_Active",]
#Increased Het Actives
#SUMO_H_A_inc<-All_S_N_avS7[All_S_N_avS7$Domains=="Active_Het"|All_S_N_avS7$Domains=="Het_Active",]
#All Active_Het
Active_Het<-All_S_N_avS6[All_S_N_avS6$Domains=="Active_Het"|All_S_N_avS6$Domains=="Het_Active",]
Active_Het_Down1<-Active_Het[Active_Het$interval1.x %in% RNA_Domains3_Actives ,]
Active_Het_Down1.1<-Active_Het_Down1[,c(2,3,10)]
Active_Het_Down2<-Active_Het[Active_Het$interval2.x %in% RNA_Domains3_Actives,]
Active_Het_Down2.1<-Active_Het_Down2[,c(2,3,10)]

colnames(Active_Het_Down2.1)<-c("Het","Active","Dif")
colnames(Active_Het_Down1.1)<-c("Active","Het","Dif")
Downs_A_H<-rbind(Active_Het_Down1.1,Active_Het_Down2.1)

ggplot(Downs_A_H, aes(x = Het, y = Active, fill = Dif)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Heatmap of Dif Values", x = "Active", y = "Het", fill = "Dif") +
  theme_minimal()

mean(Downs_A_H$Dif)

Downs_A_H_summed <- Downs_A_H %>%
  group_by(Active) %>%
  summarise(Total_Dif = sum(Dif, na.rm = TRUE))

# Print the new dataframe
head(Downs_A_H_summed)







#All PcG_Active
PcG_Active<-All_S_N_avS6[All_S_N_avS6$Domains=="Active_PcG"|All_S_N_avS6$Domains=="PcG_Active",]
PcG_Active_Down1<-PcG_Active[PcG_Active$interval1.x %in% RNA_Domains3_PcG ,]
PcG_Active_Down1.1<-PcG_Active_Down1[,c(2,3,10)]
PcG_Active_Down2<-PcG_Active[PcG_Active$interval2.x %in% RNA_Domains3_PcG,]
PcG_Active_Down2.1<-PcG_Active_Down2[,c(2,3,10)]

colnames(PcG_Active_Down2.1)<-c("Active","PcG","Dif")
colnames(PcG_Active_Down1.1)<-c("PcG","Active","Dif")
Downs_P_A<-rbind(PcG_Active_Down1.1,PcG_Active_Down2.1)

ggplot(Downs_P_A, aes(x = PcG, y = Active, fill = Dif)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Heatmap of Dif Values", x = "Active", y = "Het", fill = "Dif") +
  theme_minimal()


sum(Downs_P_A$Dif)
min(Downs_P_A$Dif)
max(Downs_P_A$Dif)
sum(Downs_P_A[Downs_P_A$Dif<0,]$Dif)
sum(Downs_P_A[Downs_P_A$Dif>0,]$Dif)

#All PcG_Het
PcG_Het<-All_S_N_avS6[All_S_N_avS6$Domains=="Het_PcG"|All_S_N_avS6$Domains=="PcG_Het",]
PcG_Het_Down1<-PcG_Het[PcG_Het$interval1.x %in% RNA_Domains3_PcG ,]
PcG_Het_Down1.1<-PcG_Het_Down1[,c(2,3,10)]
PcG_Het_Down2<-PcG_Het[PcG_Het$interval2.x %in% RNA_Domains3_PcG,]
PcG_Het_Down2.1<-PcG_Het_Down2[,c(2,3,10)]

colnames(PcG_Het_Down2.1)<-c("Het","PcG","Dif")
colnames(PcG_Het_Down1.1)<-c("PcG","Het","Dif")
Downs_P_H<-rbind(PcG_Het_Down1.1,PcG_Het_Down2.1)

ggplot(Downs_P_A, aes(x = PcG, y = Active, fill = Dif)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Heatmap of Dif Values", x = "Active", y = "Het", fill = "Dif") +
  theme_minimal()


sum(Downs_P_H$Dif)
min(Downs_P_H$Dif)
max(Downs_P_H$Dif)
sum(Downs_P_H[Downs_P_H$Dif<0,]$Dif)
sum(Downs_P_H[Downs_P_H$Dif>0,]$Dif)


#All PcG_Null
PcG_Null<-All_S_N_avS6[All_S_N_avS6$Domains=="Null_PcG"|All_S_N_avS6$Domains=="PcG_Null",]
PcG_Null_Down1<-PcG_Null[PcG_Null$interval1.x %in% RNA_Domains3_PcG ,]
PcG_Null_Down1.1<-PcG_Null_Down1[,c(2,3,10)]
PcG_Null_Down2<-PcG_Null[PcG_Null$interval2.x %in% RNA_Domains3_PcG,]
PcG_Null_Down2.1<-PcG_Null_Down2[,c(2,3,10)]

colnames(PcG_Null_Down1.1)<-c("PcG","Null","Dif")
colnames(PcG_Null_Down2.1)<-c("Null","PcG","Dif")
Downs_P_N<-rbind(PcG_Null_Down1.1,PcG_Null_Down2.1)

ggplot(Downs_P_N, aes(x = PcG, y = Null, fill = Dif)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Heatmap of Dif Values", x = "Active", y = "Het", fill = "Dif") +
  theme_minimal()

sum(Downs_P_N$Dif)
min(Downs_P_N$Dif)
max(Downs_P_N$Dif)
sum(Downs_P_N[Downs_P_N$Dif<0,]$Dif)
sum(Downs_P_N[Downs_P_N$Dif>0,]$Dif)


## Inter PcG-PcG: SUMO_P_P_inter

PcG_int_Down1<-SUMO_P_P_inter[SUMO_P_P_inter$interval1.x %in% RNA_Domains3_PcG ,]
PcG_int_Down1.1<-PcG_int_Down1[,c(2,3,10)]
PcG_int_Down2<-SUMO_P_P_inter[SUMO_P_P_inter$interval2.x %in% RNA_Domains3_PcG,]
PcG_int_Down2.1<-PcG_int_Down2[,c(2,3,10)]

colnames(PcG_int_Down1.1)<-c("PcG","PcG2","Dif")
colnames(PcG_int_Down2.1)<-c("PcG2","PcG","Dif")
Downs_P_P<-rbind(PcG_int_Down1.1,PcG_int_Down2.1)

ggplot(Downs_P_P, aes(x = PcG, y = PcG2, fill = Dif)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Heatmap of Dif Values", x = "PcG2", y = "PcG", fill = "Dif") +
  theme_minimal()

sum(Downs_P_P$Dif)
min(Downs_P_P$Dif)
max(Downs_P_P$Dif)
sum(Downs_P_P[Downs_P_P$Dif<0,]$Dif)
sum(Downs_P_P[Downs_P_P$Dif>0,]$Dif)
nrow(Downs_P_P)
nrow(Downs_P_P[Downs_P_P$PcG2 %in% RNA_Domains3_PcG,])
sum(Downs_P_P[Downs_P_P$PcG2 %in% RNA_Domains3_PcG,]$Dif)
sum(Downs_P_P[!Downs_P_P$PcG2 %in% RNA_Domains3_PcG,]$Dif)


nrow(Downs_P_P[!duplicated(Downs_P_P$PcG),])
nrow(Downs_P_P[Downs_P_P$PcG2 %in% Downs_P_P$PcG,])

length(Downs_P_P)
length(Downs_P_P[Downs_P_P$Dif<0,]$Dif)
length(Downs_P_P[Downs_P_P$Dif>0,]$Dif)




sum(Downs_P_A$Dif)


# Print the new dataframe
head(Downs_N_P_summed)
colnames(Downs_N_P_summed)<-c("DownGenesUnderNull","PcG_Sum")
colnames(Downs_N_H_summed)<-c("DownGenesUnderNull","Het_Sum")


Downs_Null_To_Repressive<-merge(Downs_N_P_summed,Downs_N_H_summed,by="DownGenesUnderNull")
rownames(Downs_Null_To_Repressive)<-Downs_Null_To_Repressive$DownGenesUnderNull
Downs_Null_To_Repressive<-Downs_Null_To_Repressive[,-1]

pheatmap(Downs_Null_To_Repressive,fontsize_row=5)






SUMO_P_P_inter2<-SUMO_P_P_inter[,c(4,8)]
colnames(SUMO_P_P_inter2)<-c("Control","SUMO")  
pheatmap(SUMO_P_P_inter2)
head(SUMO_P_P_inter2)






#All inter Nulls
SUMO_P_N<-All_S_N_avS6[All_S_N_avS6$Domains=="PcG_Null"| All_S_N_avS6$Domains=="Null_PcG",]
SUMO_H_N<-All_S_N_avS6[All_S_N_avS6$Domains=="Het_Null"| All_S_N_avS6$Domains=="Null_Het",]
#SUMO_N_A<-All_S_N_avS6[All_S_N_avS6$Domains=="Null_Active"| All_S_N_avS6$Domains=="Active_Null",]
Nullto_HP<-rbind(SUMO_P_N,SUMO_H_N)

Nullto_HP_gain<-Nullto_HP[Nullto_HP$dif>35,]


#All down genes under Null domainsm these null domains have increased contacts with Het and PcG
#Down genes under null domains : RNA_Domains3_Null
#Null domains gaining contacts with Het or PcG: Nullto_HP_gain

Nulls_DownGenes_HPgain<-RNA_Domains3_Null[RNA_Domains3_Null %in% Nullto_HP_gain$interval1.x |RNA_Domains3_Null %in% Nullto_HP_gain$interval2.x]
RNAs_Under_Null<-RNA_Domains[grep('Null',RNA_Domains$Interval),]
RNAs_Under_Null<-RNAs_Under_Null[,-1]
EnhancedVolcano(RNAs_Under_Null,lab = RNAs_Under_Null$GeneID,
                x = 'log2Ratio',
                y = 'pvalue',pCutoff = 5*10e-2,
                FCcutoff = 1.5)


memory.limit(size = 16000)
RNAs_increasesHP_underN <- RNAs_Under_Null[
  !is.na(match(RNAs_Under_Null$Interval, Nullto_HP_gain$interval1.x)) |
    !is.na(match(RNAs_Under_Null$Interval, Nullto_HP_gain$interval2.x)), ]


EnhancedVolcano(RNAs_increasesHP_underN,lab = RNAs_increasesHP_underN$GeneID,
                x = 'log2Ratio',
                y = 'pvalue',pCutoff = 5*10e-2,
                FCcutoff = 1.5)





#Active_Active #Null_Null #Het_Het
#Active_PcG #PcG_Active
#PcG_Null #Null_PcG
#Het_PcG #PcG_Het
#Het_Null #Null_Het
#Act_Het #Het_Active
#Null_Act #Act_Null

#In which domains arethe down-regulated genes


#All inter Nulls:AllNulltoAll
SUMO_P_N<-All_S_N_avS6[All_S_N_avS6$Domains=="PcG_Null"| All_S_N_avS6$Domains=="Null_PcG",]
SUMO_H_N<-All_S_N_avS6[All_S_N_avS6$Domains=="Het_Null"| All_S_N_avS6$Domains=="Null_Het",]
SUMO_A_N<-All_S_N_avS6[All_S_N_avS6$Domains=="Null_Active"| All_S_N_avS6$Domains=="Active_Null",]
Nullto_HP<-rbind(SUMO_P_N,SUMO_H_N)
AllNulltoAll<-rbind(SUMO_P_N,SUMO_H_N,SUMO_A_N)
#All inter Actives:Active_to_nonA
#All inter Hets:AllNulltoAll
SUMO_H_N<-All_S_N_avS6[All_S_N_avS6$Domains=="Het_Null"| All_S_N_avS6$Domains=="Null_Het",]
SUMO_H_P<-All_S_N_avS6[All_S_N_avS6$Domains=="PcG_Het"| All_S_N_avS6$Domains=="PcG_Het",]
SUMO_H_A<-All_S_N_avS6[All_S_N_avS6$Domains=="Het_Active"| All_S_N_avS6$Domains=="Active_Het",]
AllHetoAll<-rbind(SUMO_H_N,SUMO_H_P,SUMO_H_A)


###########################################################################################################
#####Check for the significant interactions according to ChiaSig############

Sig_NO<-read.delim2("/Users/nazliakilli/Desktop/HiC_scores/ChiaSig_py3_NOnub.sig.bedpe",header = FALSE)
head(Sig_NO)
nrow(Sig_NO)

Sig_SUMO<-read.delim2("/Users/nazliakilli/Desktop/HiC_scores/ChiaSig_py3_SUMOnub.sig.bedpe",header = FALSE)
head(Sig_SUMO)
nrow(Sig_SUMO)

###this is the list of all domains: AllDomainsSN
library(GenomicRanges)
library(dplyr)
###Start with SUMO
gr1 <- GRanges(seqnames = Sig_SUMO$V1,
               ranges = IRanges(start = Sig_SUMO$V2, end = Sig_SUMO$V3))

# Second genomic region
gr2 <- GRanges(seqnames = Sig_SUMO$V4,
               ranges = IRanges(start = Sig_SUMO$V5, end = Sig_SUMO$V6))


colnames(mcols(AllDomainsSN))
colnames(mcols(AllDomainsSN)) <- "domain"
hits2 <- findOverlaps(gr2, AllDomainsSN)
domain2 <- rep(NA, length(gr2))
domain2[queryHits(hits2)] <- mcols(AllDomainsSN)$domain[subjectHits(hits2)]
hits1 <- findOverlaps(gr1, AllDomainsSN)
domain1 <- rep(NA, length(gr1))
domain1[queryHits(hits1)] <- mcols(AllDomainsSN)$domain[subjectHits(hits1)]
Sig_SUMO$Domain1 <- domain1
Sig_SUMO$Domain2 <- domain2


##Now with Control

gr3 <- GRanges(seqnames = Sig_NO$V1,
               ranges = IRanges(start = Sig_NO$V2, end = Sig_NO$V3))

# Second genomic region
gr4 <- GRanges(seqnames = Sig_NO$V4,
               ranges = IRanges(start = Sig_NO$V5, end = Sig_NO$V6))


colnames(mcols(AllDomainsSN))
colnames(mcols(AllDomainsSN)) <- "domain"
hits2 <- findOverlaps(gr4, AllDomainsSN)
domain2 <- rep(NA, length(gr4))
domain2[queryHits(hits2)] <- mcols(AllDomainsSN)$domain[subjectHits(hits2)]
hits1 <- findOverlaps(gr3, AllDomainsSN)
domain1 <- rep(NA, length(gr3))
domain1[queryHits(hits1)] <- mcols(AllDomainsSN)$domain[subjectHits(hits1)]
Sig_NO$Domain1 <- domain1
Sig_NO$Domain2 <- domain2


###Take significant NOnub domains


All_S_N_avS3_NO$pair_key <- apply(All_S_N_avS3_NO[, c("interval1", "interval2")], 1, function(x) {
  paste(sort(x), collapse = "_")
})

# For Sig_NO: create a consistent key
Sig_NO$pair_key <- apply(Sig_NO[, c("Domain1", "Domain2")], 1, function(x) {
  paste(sort(x), collapse = "_")
})


filtered_NO <- All_S_N_avS3_NO[All_S_N_avS3_NO$pair_key %in% Sig_NO$pair_key, ]
nrow(filtered_NO)
###Take significant NOnub domains


All_S_N_avS3_SUMO$pair_key <- apply(All_S_N_avS3_SUMO[, c("interval1", "interval2")], 1, function(x) {
  paste(sort(x), collapse = "_")
})

# For Sig_SUMO: create a consistent key
Sig_SUMO$pair_key <- apply(Sig_SUMO[, c("Domain1", "Domain2")], 1, function(x) {
  paste(sort(x), collapse = "_")
})

nrow(Sig_SUMO)

filtered_SUMO <- All_S_N_avS3_SUMO[All_S_N_avS3_SUMO$pair_key %in% Sig_SUMO$pair_key, ]
nrow(filtered_SUMO)
#Cbind for future need
filtered_both<-merge(filtered_NO,filtered_SUMO,by="pair_key")
filtered_both2<-filtered_both[,c(2:8,12,13)]
filtered_both2$dif<-filtered_both2$av_score.y-filtered_both2$av_score.x
SigDoms<-filtered_both2[!duplicated(filtered_both2),]
head(SigDoms)

nrow(SigDoms)



SigDoms[SigDoms$Int_er_ra=="Intra",]

#PcG_PcG inter or intra
Sig_P_P<-SigDoms[SigDoms$Domains=="PcG_PcG",]
nrow(Sig_P_P)
#PcG_PcG inter
Sig_P_P_inter<-Sig_P_P[Sig_P_P$Int_er_ra=="Inter",]
Sig_P_P_inter_filt<-Sig_P_P_inter[abs(Sig_P_P_inter$dif) > 5,]
#PcG_PcG intra
Sig_P_P_intra<-Sig_P_P[Sig_P_P$Int_er_ra=="Intra",]
nrow(Sig_P_P_intra)


Sig_P_P_inter2<-Sig_P_P_inter[,c(4,8)]
colnames(Sig_P_P_inter2)<-c("Control","SUMO")  
#pheatmap(Sig_P_P_inter2)
head(Sig_P_P_inter2)

l_Sig_P_P_inter2<- Sig_P_P_inter2 %>%
  rownames_to_column("Interval") %>%  # If row names are Intervals, convert to a column
  pivot_longer(cols = c(Control, SUMO), names_to = "Condition", values_to = "Value")


#library(ggstatsplot)
ggbetweenstats(
  data  = l_Sig_P_P_inter2,
  x     = Condition,
  y     = Value ,
  title = "Inter PcG in Control vs SUMO sig scores"
)

intersept<-median(Sig_P_P_inter2$SUMO)

ggplot(l_Sig_P_P_inter2, aes(x = Condition, y = Value, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.5,show.legend = FALSE) + # Violin plot with transparency
  theme_minimal() + # Clean theme
  labs(title = "",
       x="",
       y = "PcG-PcG") +
  scale_fill_manual(values = c("Control" = "darkgreen", "SUMO" = "orange"))+ # Custom colors
  stat_compare_means(method = "t.test", label = "p.format",aes(x=1.5),vjust = 0,size=3)+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+theme(plot.title = element_text(size = 10,hjust = 0.5)) +
  geom_hline(yintercept = intersept, linetype = "dashed", color = "lightgrey")+
  geom_boxplot(width = 0.1, outlier.shape = NA,show.legend = FALSE)  # Small box plot









##PcG Active

#PcG_PcG inter or intra
Sig_P_A<-SigDoms[SigDoms$Domains=="PcG_Active" |SigDoms$Domains=="Active_PcG" ,]
#PcG_PcG inter
Sig_P_A_inter<-Sig_P_A[Sig_P_A$Int_er_ra=="Inter",]

Sig_P_A_inter2<-Sig_P_A_inter[,c(4,8)]
colnames(Sig_P_A_inter2)<-c("Control","SUMO")  




l_Sig_P_A_inter2<- Sig_P_A_inter2 %>%
  rownames_to_column("Interval") %>%  # If row names are Intervals, convert to a column
  pivot_longer(cols = c(Control, SUMO), names_to = "Condition", values_to = "Value")

ggbetweenstats(
  data  = l_Sig_P_A_inter2,
  x     = Condition,
  y     = Value ,
  title = " PcG Active in Control vs SUMO sigs "
)


#RNA_Domains3_PcG_ups, RNA_Domains3_PcG_downs


## Inter PcG-PcG: SUMO_P_P_inter
#Downs
PcG_int_Down1<-Sig_P_P_inter_filt[Sig_P_P_inter_filt$interval1.x %in% RNA_Domains3_PcG_downs ,]
PcG_int_Down1.1<-PcG_int_Down1[,c(2,3,10)]
PcG_int_Down2<-Sig_P_P_inter_filt[Sig_P_P_inter_filt$interval2.x %in% RNA_Domains3_PcG_downs,]
PcG_int_Down2.1<-PcG_int_Down2[,c(2,3,10)]

colnames(PcG_int_Down1.1)<-c("PcG","PcG2","Dif")
colnames(PcG_int_Down2.1)<-c("PcG2","PcG","Dif")
Downs_P_P<-rbind(PcG_int_Down1.1,PcG_int_Down2.1)

ggplot(Downs_P_P, aes(x = PcG, y = PcG2, fill = Dif)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Heatmap of Dif Values", x = "PcG2", y = "PcG", fill = "Dif") +
  theme_minimal()

sum(Downs_P_P$Dif)



#Ups
PcG_int_Up1<-Sig_P_P_inter_filt[Sig_P_P_inter_filt$interval1.x %in% RNA_Domains3_PcG_ups ,]
PcG_int_Up1.1<-PcG_int_Up1[,c(2,3,10)]
PcG_int_Up2<-Sig_P_P_inter_filt[Sig_P_P_inter_filt$interval2.x %in% RNA_Domains3_PcG_ups,]
PcG_int_Up2.1<-PcG_int_Up2[,c(2,3,10)]

colnames(PcG_int_Up1.1)<-c("PcG","PcG2","Dif")
colnames(PcG_int_Up2.1)<-c("PcG2","PcG","Dif")
Ups_P_P<-rbind(PcG_int_Up1.1,PcG_int_Up2.1)

ggplot(Ups_P_P, aes(x = PcG, y = PcG2, fill = Dif)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Heatmap of Dif Values", x = "PcG2", y = "PcG", fill = "Dif") +
  theme_minimal()

sum(Ups_P_P$Dif)

#Get the non overlapping domains in Up and Down

Ups_P_P

Ups_P_P$Domains <- apply(Ups_P_P[, c("PcG", "PcG2")], 1, function(x) {
  paste(sort(x), collapse = "_")
})

Downs_P_P$Domains <- apply(Downs_P_P[, c("PcG", "PcG2")], 1, function(x) {
  paste(sort(x), collapse = "_")
})


Ups_P_P_only <- anti_join(Ups_P_P, Downs_P_P, by = "Domains")
Downs_P_P_only <- anti_join(Downs_P_P, Ups_P_P, by = "Domains")


ggplot(Downs_P_P_only, aes(x = PcG, y = PcG2, fill = Dif)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Heatmap of Dif Values", x = "PcG2", y = "PcG", fill = "Dif") +
  theme_minimal()


ggplot(Ups_P_P_only, aes(x = PcG, y = PcG2, fill = Dif)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Heatmap of Dif Values", x = "PcG2", y = "PcG", fill = "Dif") +
  theme_minimal()

sum(Ups_P_P_only$Dif)
mean(Ups_P_P_only$Dif)
mean(Downs_P_P_only$Dif)
median(Downs_P_P_only$Dif)
median(Ups_P_P_only$Dif)




############################################
Sig_P_P_inter2.2<-Sig_P_P_inter2[Sig_P_P_inter2$Control< 0,]
Sig_P_P_inter2.3<-Sig_P_P_inter2[Sig_P_P_inter2$Control > 0,]

Sig_P_P_inter2.2<- Sig_P_P_inter2.2 %>%
  rownames_to_column("Interval") %>%  # If row names are Intervals, convert to a column
  pivot_longer(cols = c(Control, SUMO), names_to = "Condition", values_to = "Value")

Sig_P_P_inter2.3<- Sig_P_P_inter2.3 %>%
  rownames_to_column("Interval") %>%  # If row names are Intervals, convert to a column
  pivot_longer(cols = c(Control, SUMO), names_to = "Condition", values_to = "Value")





nrow(l_Sig_P_P_inter2)
# Convert Interval to a factor (if categorical)
l_Null_PcG3$Interval <- as.factor(l_Null_PcG3$Interval)
l_Null_PcG2$Interval <- as.factor(l_Null_PcG2$Interval)
l_Null_PcG4$Interval <- as.factor(l_Null_PcG4$Interval)

library(ggpubr)
# Create the violin plot
h<-ggplot(l_Sig_P_P_inter2, aes(x = Condition, y = Value, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.5,show.legend = FALSE) + # Violin plot with transparency
  theme_minimal() + # Clean theme
  labs(title = "",
       x="",
       y = "PcG-PcG") +
  scale_fill_manual(values = c("Control" = "darkgreen", "SUMO" = "orange"))+ # Custom colors
  stat_compare_means(method = "t.test", label = "p.format",aes(x=1.5),vjust = 0,size=3)+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+theme(plot.title = element_text(size = 10,hjust = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "lightgrey")+
  geom_boxplot(width = 0.1, outlier.shape = NA,show.legend = FALSE)  # Small box plot







