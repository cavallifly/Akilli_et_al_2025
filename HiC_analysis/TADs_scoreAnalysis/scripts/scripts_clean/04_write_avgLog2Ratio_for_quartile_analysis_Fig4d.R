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

data <- data[data$chr1 == data$chr2,]

print(paste0("Log2ratio per TAD definition"))
Log2ratio <- read.table("avgLog2ratio_per_TAD_basedOnGeneTSS.tsv",header=F)
colnames(Log2ratio) <- c("interval","chr","start","end","sumLog2ratio","sum2Log2ratio","nGenes","avgLog2ratio","stddevLog2ratio")
Log2ratio <- Log2ratio[,c("interval","interval","avgLog2ratio")]
colnames(Log2ratio) <- c("interval1","interval2","avgLog2ratio")
print(head(Log2ratio))

print(paste0("Associate avgLog2ratio per TAD"))
data <- merge(data, Log2ratio, by = "interval1")
data <- data[,c("chr1","start1","end1","chr2","start2","end2","interval1","avgLog2ratio","interval2.x","av_score","sample")]
colnames(data) <- c("chr1","start1","end1","chr2","start2","end2","interval1","avgLog2ratio1","interval2","av_score","sample")
print(head(data))

data <- merge(data, Log2ratio, by = "interval2")
data <- data[,c("chr1","start1","end1","chr2","start2","end2","interval1.x","avgLog2ratio1","interval2","avgLog2ratio","av_score","sample")]
colnames(data) <- c("chr1","start1","end1","chr2","start2","end2","interval1","avgLog2ratio1","interval2","avgLog2ratio2","av_score","sample")

print(paste0("Compute avgLog2ratio per TAD pair as the average of the avgLog2ratio of each of the two TADs"))
data$avgLog2ratio <- (data$avgLog2ratio1+data$avgLog2ratio2)*0.5
data <- data[,c("chr1","start1","end1","chr2","start2","end2","interval1","avgLog2ratio1","interval2","avgLog2ratio2","avgLog2ratio","av_score","sample")]
print(head(data))

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
    
    for(s2 in states)
    {
        state2 = s2
	statePair <- paste0(state1,"_",state2)

	outFile <- paste0("avgLog2ratio_trans1Dinterval_all_domains_vs_all_domains_",statePair,"_data.tab")
	if(file.exists(outFile)){file.remove(outFile);}

	dataStates <- data[grep(statePair,data$statePair),c("sample","interval1","interval2","avgLog2ratio","av_score")]
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
	    dataStates1 <- data.frame(dataStates1$interval1, dataStates1$interval2, dataStates1$avgLog2ratio, dataStates1$av_score)
	    colnames(dataStates1) <- c("interval1","interval2","avgLog2ratio",paste0("av_score_",sample))
	    if(nrow(dataStatesNew) == 0)
	    {
	        dataStatesNew <- dataStates1
		colnames(dataStatesNew) <- c("interval1","interval2","avgLog2ratio",paste0("av_score_",sample))
	    } else {
	        print(names(dataStatesNew))
	        print(names(dataStates1))		
	        dataStatesNew <- merge(dataStatesNew,dataStates1,by = c("interval1","interval2","avgLog2ratio"))
	    }
   	    print(head(dataStatesNew))	    
	}
	dataStates <- dataStatesNew
	###DONE###

	print(paste0("Writing all the average scores (All)"))
	if(nrow(dataStates) > 0)
	{
	    print(statePair)
 	    print(head(dataStates))

	    for(sample in samples)
	    {
	        print(paste0("av_score_",sample))
		#df <- dataStates[,c(paste0("av_score_",sample))]
		df <- dataStates[,c("avgLog2ratio")]		
		df <- data.frame(paste0(sample,"All"),df)
	        write.table(df, file = outFile, row.names = F, sep = "\t", append = T, quote = F, col.names = F)
	    }	    
        }

	print(paste0("Dividing the dataset in quartiles depending on the 1Dav_score in ",refSample))
	qControl <- quantile(dataStates[,c(paste0("av_score_",refSample))], na.rm = TRUE)
	print(qControl)

	for(i in seq(from=1,to=4,by=1))
	{
	      quartile <- dataStates[dataStates[,c(paste0("av_score_",refSample))] >= qControl[[i]] & dataStates[,c(paste0("av_score_",refSample))] <= qControl[[i+1]],]
	      print(head(quartile))

	      for(sample in samples)
	      {
                  print(paste0("av_score_",sample))
                  #df <- quartile[,c(paste0("av_score_",sample))]
                  df <- quartile[,c("avgLog2ratio")]		  
                  df <- data.frame(paste0(sample,"Q",i),df)
                  write.table(df, file = outFile, row.names = F, sep = "\t", append = T, quote = F, col.names = F)
              }
	}
	allPoints <- read.table(outFile)
	colnames(allPoints) <- c("category","avgLog2ratio")
	print(table(allPoints$category))
    }
}

quit()
