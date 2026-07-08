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
print(paste0("1Ddistance definition"))
data$distance <- abs((data$start1+data$end1)*0.5-(data$start2+data$end2)*0.5)

data$statePair <- paste0(gsub("[0-9]+", "", data$interval1),"_",gsub("[0-9]+", "", data$interval2))
data$statePair <- gsub("Active_PcG" ,"PcG_Active" ,data$statePair)
data$statePair <- gsub("Null_PcG"   ,"PcG_Null"   ,data$statePair)
data$statePair <- gsub("Het_PcG"   ,"PcG_Het"   ,data$statePair)
data$statePair <- gsub("Null_Active","Active_Null",data$statePair)
data$statePair <- gsub("Het_Active","Active_Het",data$statePair)
data$statePair <- gsub("Null_Het","Het_Null",data$statePair)
print(head(data))
#outFileAll <- paste0("1Ddistances_cisAndTrans1Dinterval_all_domains_vs_all_domains_complete.tab")
#write.table(data, file = outFileAll, row.names = F, sep = "\t", append = F, quote = F, col.names = T)
outFileAll <- paste0("1Ddistances_cisAndTrans1Dinterval_all_domains_vs_all_domains.tab")
if(file.exists(outFileAll)){file.remove(outFileAll);}
#quit()

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

	outFile <- paste0("1Ddistances_cisAndTrans1Dinterval_all_domains_vs_all_domains_",statePair,"_data.tab")
	if(file.exists(outFile)){file.remove(outFile);}

	dataStates <- data[grep(statePair,data$statePair),c("sample","interval1","interval2","distance","av_score")]
	if(nrow(dataStates) == 0){next;}
	print(statePair)

	dataStates$intervals <- paste0(dataStates$interval1,"_",dataStates$interval2)

	print(paste0("Keep only the trans-interval scores that are present in both the conditions"))
	#dataStates <- dataStates[dataStates$interval1 != dataStates$interval2,]
	counts <- table(dataStates$intervals)
	counts <- counts[counts == 2]
	dataStates <- dataStates[dataStates$intervals %in% names(counts),]
	###DONE###

	print(paste0("Write the table as interval1 interval2 av_score_in_sample1 av_score_in_sample2"))
	dataStatesNew <- data.frame()
	for(sample in samples)
	{
	    dataStates1 <- dataStates[grep(sample,dataStates$sample),]
	    dataStates1 <- data.frame(dataStates1$interval1, dataStates1$interval2, dataStates1$distance, dataStates1$av_score)
	    colnames(dataStates1) <- c("interval1","interval2","distance",paste0("av_score_",sample))
	    if(nrow(dataStatesNew) == 0)
	    {
	        dataStatesNew <- dataStates1
		colnames(dataStatesNew) <- c("interval1","interval2","distance",paste0("av_score_",sample))
	    } else {
	        print(names(dataStatesNew))
	        print(names(dataStates1))		
	        dataStatesNew <- merge(dataStatesNew,dataStates1,by = c("interval1","interval2","distance"))
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
		df <- dataStates[,c("distance")]		
		df <- data.frame(paste0(sample,"All"),df)
	        #write.table(df, file = outFile, row.names = F, sep = "\t", append = T, quote = F, col.names = F)
	    }	    
        }

	print(paste0("Dividing the dataset in quartiles depending on the 1Dav_score in ",refSample))
	qControl <- quantile(dataStates[,c(paste0("av_score_",refSample))], na.rm = TRUE)
	print(qControl)

	for(i in seq(from=1,to=4,by=1))
	{
	      quartile <- dataStates[dataStates[,c(paste0("av_score_",refSample))] >= qControl[[i]] & dataStates[,c(paste0("av_score_",refSample))] <= qControl[[i+1]],]
	      print(head(quartile))
	      toPrint <- cbind(quartile,paste0(statePair,"_Q",i))
	      colnames(toPrint) <- c(colnames(quartile),"quartile")
	      write.table(toPrint, file = outFileAll, row.names = F, sep = "\t", append = T, quote = F, col.names = T)	     
	      #quit()

	      for(sample in samples)
	      {
                  print(paste0("av_score_",sample))
                  #df <- quartile[,c(paste0("av_score_",sample))]
                  df <- quartile[,c("distance")]		  
                  df <- data.frame(paste0(sample,"Q",i),df)
                  #write.table(df, file = outFile, row.names = F, sep = "\t", append = T, quote = F, col.names = F)
              }
	}
	#allPoints <- read.table(outFile)
	#colnames(allPoints) <- c("category","distance")
	#print(table(allPoints$category))
    }
}

quit()
