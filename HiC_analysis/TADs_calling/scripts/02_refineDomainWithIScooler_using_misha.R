require(misha)
require(shaman)

require(doParallel)

### misha working DB
mDBloc <-  '/zdata/data/mishaDB/trackdb/'
db <- 'dm6'
dbDir <- paste0(mDBloc,db,'/')
gdb.init(dbDir)
gdb.reload()

source("/zdata/data/auxFunctions/auxFunctions.R")
options(scipen=20,gmax.data.size=0.5e8,shaman.sge_support=1)

domainRes <-  5000
isRes     <-  1000
window    <- 25000

insTracks <- list(
	WD_2801_merge    = paste0("insulation.hic_WD_2801_merge_dm6_NA_w",window,"bp_r",isRes,"bp_cooler"),	
	WD_3601_merge    = paste0("insulation.hic_WD_3601_merge_dm6_NA_w",window,"bp_r",isRes,"bp_cooler"),	
	WD_NOnub_merge   = paste0("insulation.hic_WD_NOnub_merge_dm6_NA_w",window,"bp_r",isRes,"bp_cooler"),		
	WD_SUMOnub_merge = paste0("insulation.hic_WD_SUMOnub_merge_dm6_NA_w",window,"bp_r",isRes,"bp_cooler")	
	)

insIntrv <- gintervals.all()

# Retrieve the insulation tracks from mishaDB
print("Get the min, mean, stddev, and max of the insulation tracks")
insData <- list()
for(set in names(insTracks))
{
    print(set, row.names=F, quote=F)
    track <- insTracks[[set]]
    data <- gextract(track,insIntrv,iterator=1e3,colnames="INS")	
    data$INS  <- clipQuantile(-data$INS,0.999)
    data$INS  <- scaleData(data$INS,0,1,1e-9)
    
    trackName  <- set
    #trackName1 <- unlist(strsplit(trackName,"_"))[1:2]
    #trackName1 <- paste(trackName1,collapse="_")
    print(trackName)
    window     <- trackName #paste0("w",unlist(strsplit(trackName,"_"))[3])
    if(length(insData[[trackName]]) == 0)
    {
        insData[[trackName]] <- data[c(1,2,3,4)]
	colnames(insData[[trackName]]) <- c("chrom","start","end",window)
    } else {
        cNames <- c(colnames(insData[[trackName]]),window)
        insData[[trackName]] <- cbind(insData[[trackName]],data$INS)
	colnames(insData[[trackName]]) <- cNames
    }
    
}

print("")
for(trackName in names(insData))
{
    print(trackName)
    print(head(insData[[trackName]]))

    #insData[[trackName]]$min    <- apply(insData[[trackName]][4:length(insData[[trackName]])],1,min)    
    insData[[trackName]]$mean   <- apply(insData[[trackName]][4:length(insData[[trackName]])],1,mean)
    #insData[[trackName]]$stddev <- apply(insData[[trackName]][4:length(insData[[trackName]])],1,sd)
    #insData[[trackName]]$max    <- apply(insData[[trackName]][4:length(insData[[trackName]])],1,max)

    outProfile <- paste0("insulationProfiles_",trackName,".tsv")
    write.table(insData[[trackName]], file = outProfile, sep="\t", row.names=FALSE, quote=FALSE)

    print("", row.names=F, quote=F)
}
#print(head(insData))

# Retrieve the TopDom domains
for(trackName in names(insData))
{
    print(trackName)
    TopDomFile <- list.files(path="./",patter=paste0("Top.*",trackName,".*","_",5,"_.*",domainRes,"bp.*filtered.*"))
    #TopDomFile <- list.files(path="./",patter=paste0("Top.*",trackName,".*","_",5,"_.*",domainRes,"bp.t.*"))
    #TopDomFile <- TopDomFile[-grep("refined",TopDomFile)]    
    print(TopDomFile)

    #TopDomDomains <- read.table(TopDomFile,header=F,fill=TRUE)
    #TopDomDomains <- TopDomDomains[,c(2,3,4)]
    #TopDomDomains <- TopDomDomains[-grep("start",TopDomDomains),]
    #outFile <- gsub(".tsv",paste0("_refined_at_",isRes,"bp.tsv"),TopDomFile)    

    TopDomDomains <- read.table(TopDomFile,header=F,fill=TRUE)
    TopDomDomains <- TopDomDomains[,c(1,2,3)]
    outFile <- gsub(".bed",paste0("_refined_at_",isRes,"bp.tsv"),TopDomFile)    

    colnames(TopDomDomains) <- c("chrom","start","end")
    print(head(TopDomDomains))

    currentTrack <- insData[[trackName]]
    print(head(currentTrack))

    refinedTopDomDomains <- data.frame()

    if(file.exists(outFile)){next;}
    print(outFile)

    for(d in 1:nrow(TopDomDomains))
    {
	domain <- TopDomDomains[d,]

        boundary1 <- data.frame(domain$chrom,as.numeric(domain$start)-domainRes,as.numeric(domain$start)+domainRes)
	colnames(boundary1) <- c("chrom","start","end")
	print(boundary1)
        boundary2 <- data.frame(domain$chrom,as.numeric(domain$end)-domainRes,as.numeric(domain$end)+domainRes)
	colnames(boundary2) <- c("chrom","start","end")
	print(boundary2)

	print(currentTrack[currentTrack$chrom == boundary1$chrom & (boundary1$start <= currentTrack$start & currentTrack$end <= boundary1$end),])
	print(currentTrack[currentTrack$chrom == boundary2$chrom & (boundary2$start <= currentTrack$start & currentTrack$end <= boundary2$end),])

	df1 <- currentTrack[currentTrack$chrom == boundary1$chrom & (boundary1$start <= currentTrack$start & currentTrack$end <= boundary1$end),]
	df2 <- currentTrack[currentTrack$chrom == boundary2$chrom & (boundary2$start <= currentTrack$start & currentTrack$end <= boundary2$end),]

	if(nrow(df1)!=0)
	{
	    print(df1[df1$mean == max(df1$mean),])

	    #refStart <- (df1[df1$mean == max(df1$mean),]$start + df1[df1$mean == max(df1$mean),]$end) / 2
	    refStart <- df1[df1$mean == max(df1$mean),]$start	    
	    if(length(refStart) > 1)
	    {
	        refStart <- as.numeric(domain$start)
	    }
	} else {
	    refStart <- as.numeric(domain$start)
	}
	if(nrow(df2)!=0)
	{
	    print(df2[df2$mean == max(df2$mean),])
	    #refEnd   <- (df2[df2$mean == max(df2$mean),]$start + df2[df2$mean == max(df2$mean),]$end) / 2
	    refEnd   <- max(df2[df2$mean == max(df2$mean),]$start)
	    if(length(refEnd) > 1)
            {
              	refEnd <- as.numeric(domain$end)
            }
	} else {
	    refEnd <- as.numeric(domain$end)
	}

	refinedTopDomDomain <- data.frame(boundary1$chrom,refStart,refEnd)
	colnames(refinedTopDomDomain) <- c("chrom","start","end")
	print(refinedTopDomDomain)

	if(is.null(refinedTopDomDomains))
	{
	    refinedTopDomDomains <- refinedTopDomDomain
	} else {
	    refinedTopDomDomains <- rbind(refinedTopDomDomains,refinedTopDomDomain)
	}
    }
    colnames(refinedTopDomDomains) <- c("chrom","start","end")
    write.table(refinedTopDomDomains, file = outFile, sep="\t", row.names=FALSE, quote=FALSE)    
}
