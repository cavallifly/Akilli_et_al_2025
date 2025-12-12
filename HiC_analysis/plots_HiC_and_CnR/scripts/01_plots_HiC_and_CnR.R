options("scipen"=999, max.print=999999)

library(misha)
library(shaman)
library(reshape2)
library(plyr)
source('./scripts/auxFunctions.R')

mDBloc <- '/zdata/data/mishaDB/trackdb/'
db <- 'dm6'
dbDir <- paste0(mDBloc,db,'/')
gdb.init(dbDir)
gdb.reload()

###Region to visualize
#Visualization: chr18:53,859,000-56,456,400
#Modelling region -> chr18:53700000-56700000
domains <- read.table('TopDom_domains_hic_WD_NOnub_merge_dm6_NA_window_5_at_5000bp_filtered_refined_at_1000bp_with_state.bed')
colnames(domains) <- c("chrom1","start1","end1","chrom2","start2","end2","name","value")
print(head(domains))

regions <- read.table('region_for_dac.bed')
colnames(regions) <- c("chrom1","start1","end1","chrom2","start2","end2","name","value")
print(head(regions))
extension = 20000
zmax <-  100
zmin <- -100

Capture    = 0 # 1 is FALSE 0 is TRUE
Genes      = 0
ViewPoints = 1 # 1 is FALSE 0 is TRUE
Insulation = 0 # 1 is FALSE 0 is TRUE
ChIPseqs   = 0 # 1 is FALSE 0 is TRUE

###Resolutions
chicRes    = 5000
virt4CRes  = 5000
insRes     = 1000

###Colors
scoreCol <- colorRampPalette(c("darkblue","white","darkred"))(200)
# Plot the colorbar once for all
pdf(paste0('colorbar_scores.pdf'),width=5,height=8)
labels <- c(as.character(as.integer(zmin)),as.character(as.integer(zmin/2.)),"0",as.character(as.integer(zmax/2.)),as.character(as.integer(zmax)))
x <- length(scoreCol)
plot(y=1:x,x=rep(0,x),col=scoreCol, pch=15, xlab='', ylab='', yaxt='n', xaxt='n',frame=F,ylim=c(-x,x),xlim=c(-5,5), main='');
sapply(1:length(labels), function(y){text(y=seq(0,1,length.out=length(labels))[y]*x,x=0,labels=labels[y],pos=2)});
dev.off()
#quit()

insCols     =  c("#b22222")
#ChIPseq tracks      R45     G45      B45    A255
colChip = rgb(45 /255, 45/255, 45/255, alpha=255/255)
#RNAseq tracks       R194   G37      B92    A255
colRNA  = rgb(194/255, 37/255, 92/255, alpha=255/255)

###Initialization of variables
allViewPoints = list()

# Plots per condition
celltypes  <- list('WD')

for(r in 1:nrow(regions))
{
   region <- regions[r,]  

   chrom1 <- regions[r,]$chrom1
   start1 <- regions[r,]$start1
   end1   <- regions[r,]$end1

   chrom2 <- regions[r,]$chrom2
   start2 <- regions[r,]$start2
   end2   <- regions[r,]$end2
   name   <- regions[r,]$name
   value  <- regions[r,]$value

   print(region)

   if(is.na(end2) || is.na(end1) || is.na(start2) || is.na(start1))
   {
       next
   }      

   if(start1 < 100001)
   {
       next
   }
   if(start2 < 100001)
   {
       next
   }
   if(start1 >= end2)
   {
       next
   }      

   if((end1+extension) > gintervals.all()[gintervals.all()$chrom == chrom1,]$end)
   {
       next
   }
   if((end2+extension) > gintervals.all()[gintervals.all()$chrom == chrom2,]$end)
   {
       next
   }   
   print(paste0(chrom1,' ',start1,' ',end1,' ',chrom2,' ',start2,' ',end2))

   interval1D  <- gintervals(chrom1,start1-extension,end2+extension)
   print(paste0("Interval1D ",paste(interval1D,collapse=" ")))
   interval1D1 <- gintervals(chrom1,start1,end1)
   interval1D2 <- gintervals(chrom2,start2,end2)
   
for(ChIPres in c(2000))
{
    for(celltype in celltypes)
    {
        if(celltype == "WD")
	{
	    conditions <- list('NOnub','SUMOnub')
	    refCondition <- 'NOnub'

	    chipDatasets <- c('H3K27me3','Pc','H3K27ac','H3K9me3')
    	    plotOrder    <- chipDatasets
	    plotTitles   <- chipDatasets

	    #Ivana's suggestion
	    plotMin      <- list(0,0,0,0)
	    if(ChIPres == 2000)
	    {
		plotMax      <- list(1000, 150, 1000, 1000)
	    }
	    #plotCols     = c(rep(colChip,length(plotTitles)-1),colRNA)
	    plotCols     = c("blue","cyan","red","green")
	}

        for(condition in conditions)
        {
	    ChIPseqCondition = condition
	    print(ChIPseqCondition)
	    if(condition == "NOnub")
	    {
	        plotAnn = NA
	    }
	    if(condition == "SUMOnub")
	    {
	        plotAnn = read.table("./scripts/insulation_changes_SUMOnub.bed",header=T)
	    }	    
	    if(condition == "VELOnub")
	    {
	        plotAnn = read.table("./scripts/insulation_changes_VELOnub.bed",header=T)
	    }	    

            if(ViewPoints == 0)
            {
	        #### List of all viepoints of the virtual 4C tracks
	        #viewpoint at the following points
	        #proximal promoter (P)     - deletion coordinates -    chr18:54,986,882-54,991,228
	        #alternative promoter (aP)                             chr18:54,991,244-54,993,455
	        #enhancer A (A)            - deletion coordinates -    chr18:55,544,055-55,552,304
	        #enhancer B (B)            - deletion coordinates -    chr18:55,811,387-55,818,779
	        allViewPoints <- list(
			        proximalPromoter    = list('chr18', 54986882, 54991228),
					  )
	        plotOrdervirt4C  <- c("proximalPromoter","alternativePromoter","enhancerA","enhancerB")	
	        plotTitlesvirt4C <- c("Proximal promoter","Alternative promoter","A","B")
	        #v4C tracks               R57      G106   B156  A255
	        plotColsvirt4C   <- c(rep(rgb(0.2235294, 0.4156863, 0.6117647, alpha=1.0),length(allViewPoints)))
	        plotMaxvirt4C    <- c(rep(100,length(allViewPoints)))
            }
            insTracks <- c()
            if(Insulation == 0)
            {
	        ### Reference Insulation ###
		#hic_WD_NOnub_merge_dm6_NA_w25000bp_r1000bp_cooler.track
		insTracks <- c(gtrack.ls('insulation',celltype,refCondition,paste0('_w',25000,'bp'),paste0('_r',insRes,'bp')))				
		#insTracks <- insTracks[-grep('500kb|450kb|400kb|350kb',insTracks)]
		print(insTracks)

		#trackNames <- gsub("kb_r", " ", gsub("_w"," w",insTracks))
		trackNames <- gsub("insulation.|merge_dm6_NA_|_cooler|hic_","",insTracks)
		print(trackNames)
		insRefData <- gextract(insTracks,intervals=interval1D,iterator=insRes, colnames=trackNames)
		print(head(insRefData))
		insRefData[is.na(insRefData)] = 0
		#quit()
		for(t in c(trackNames,trackNames))
        	{
		    insRefData[,t] <- clipQuantile(-insRefData[,t],0.999)
		    insRefData[,t] <- scaleData(insRefData[,t],0,1,1e-9)
		    print(max(insRefData[,t]))
		    print(min(insRefData[,t]))
		}

		#insRefData$min    <- apply(insRefData[1:length(insRefData)],1,min)    
		#insRefData$mean   <- apply(insRefData[1:length(insRefData)],1,mean)
		#insRefData$stddev <- apply(insRefData[1:length(insRefData)],1,sd)
		#insRefData$max    <- apply(insRefData[1:length(insRefData)],1,max)
		#insRefData$min    <- apply(insRefData,1,min)    
		#insRefData$mean   <- apply(insRefData,1,mean)
		#insRefData$stddev <- apply(insRefData,1,sd)
		#insRefData$max    <- apply(insRefData[,5:length(insRefData)-1],1,max)    		
		#print(head(insRefData))		
		#quit()
		print(head(insRefData))		
		insRefData = insRefData[(start1 <= insRefData$start & insRefData$start <= end1) & (start1 <= insRefData$end & insRefData$end <= end1),]
		#insRefData = insRefData[,c("chrom","start","end","mean","stddev")]
		#insRefData = insRefData[,c("chrom","start","end","mean")]
		insRefData = insRefData[,c(1,2,3,4)]
		insRefData$stddev <- 0
		colnames(insRefData) < c("chrom","start","end","mean","stddev")
		print(head(insRefData))				
		print(insTracks)

                ### Insulation ###
	        insDatasets <- c('25000bp')
    	        insOrder    <- c('25000bp')
    	        insTitles   <- c(paste0('Insulation score 25000bp window ',insRes,'bp'))
    	        insMax      <- c(1.0)

		#insTracks <- c(gtrack.ls('insulation',celltype,condition,'w400kb','_r',insRes/1e3,'kb'))
		insTracks <- c(gtrack.ls('insulation',celltype,condition,'w25000bp','_r',insRes,'bp'))				
		#trackNames <- gsub("kb_r", " ", gsub("_w"," w",insTracks))
		trackNames <- gsub("insulation.|merge_dm6_NA_|_cooler|hic_","",insTracks)
		print(trackNames)
		insData <- gextract(insTracks,intervals=interval1D,iterator=insRes, colnames=trackNames)
		print(head(insData))
		insData[is.na(insData)] = 0

		for(t in trackNames)
        	{
		    insData[,t] <- clipQuantile(-insData[,t],0.999)
		    insData[,t] <- scaleData(insData[,t],0,1,1e-9)
		    print(max(insData[,t]))
		    print(min(insData[,t]))
		}

		#insData$min    <- apply(insData[4:length(insData)],1,min)    
		#insData$mean   <- apply(insData[4:length(insData)],1,mean)
		#insData$stddev <- apply(insData[4:length(insData)],1,sd)
		#insData$max    <- apply(insData[4:length(insData)],1,max)    
		#print(head(insData))		

		#print(head(insData))		
		#insData = insData[(start1 <= insData$start & insData$start <= end1) & (start1 <= insData$end & insData$end <= end1),]
		#insData = insData[,c("chrom","start","end","mean","stddev")]		
		print(head(insData))				
		print(insTracks)

            }
	    if(ChIPseqs == 0)
	    {
    	        ### ChIPseq ###
      	        print(chipDatasets)
		chipTracks <- c()
		#print(gtrack.ls('chipseq','NA',celltype,condition))
		#		quit()
		for(chip in chipDatasets)
		{
		    chipTracks <- c(chipTracks,gtrack.ls('chipseq','NA',celltype,condition,chip))
		}
		for(chip in c("Pc"))
		{
		    chipTracks <- c(chipTracks,gtrack.ls('cnr','NA',celltype,condition,chip))
		}		
	        print(chipTracks)
		#quit()
	    }

       	    hicDatasets <- c(celltype)[1]	
    	    tag <- paste0(celltype,'_',condition)

	    print(tag)

	    #outfile <- paste0("./plots_shamanScores_",chrom1,"_",start1,"_",end2,"_",celltype,"_",gsub("_","",condition),"_chicRes",chicRes,"bp_insRes",insRes,"bp_",name,"_",value,".png")
	    outfile <- paste0("./plots_shamanScores_",chrom1,"_",start1,"_",end2,"_",celltype,"_",gsub("_","",condition),"_chicRes",chicRes,"bp_insRes",insRes,"bp_",name,"_",value,".pdf")	    
	    print(outfile)
	    if (file.exists(outfile)){next}

	    ##############################################################################
	    #plotRatios <- list(hic=1.5,genes=0.3,ins=0.5,virt4C=0.5,chip=0.2,clusters=0.075,scale=0.1)
	    plotRatios <- list(hic=1.,genes=0.1,ins=0.15,virt4C=0.3,chip=0.15,clusters=0.075,scale=0.1)	    
	    band <- 1
	    ##############################################################################

	    nPlots <- length(hicDatasets) + length(chipDatasets) + length(insTracks) + 1 ### + 1 is for genes!
	    print(nPlots)
	    # PNG
	    #lM <- 0.4 # 3.5
	    lM <- 3. # 3.5	    
	    #bM <- (0.2*band)/2
	    #tM <- (0.2*band)/2
	    bM <- lM
	    tM <- lM
	    rM <- lM

	    width=550
	    # par("cra") = (width,height) size of default character in "rasters" (pixels) https://stackoverflow.com/questions/10277292/how-to-determine-symbol-size-in-x-and-y-units
	    # For pch in 0:25 the default size is about 75% of the character height (see par("cin")) https://stat.ethz.ch/R-manual/R-patched/library/graphics/html/points.html
	    # https://stackoverflow.com/questions/17213293/how-to-get-r-plot-window-size
	    #		    print()

	    # Plot squared matrix		        
	    print(width-96*(lM+rM))
	    print(par("pin")*96)

	    print(paste0("Width ",width))
	    #mai		    A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.
	    #png(outfile, width=width,height=(length(hicDatasets)*band*plotRatios[['hic']]+length(insTracks)*plotRatios[['ins']]+plotRatios[['genes']]+length(chipDatasets)*plotRatios[['chip']])*width, type = "cairo")
	    pdf(outfile, width=width,height=(length(hicDatasets)*band*plotRatios[['hic']]+length(insTracks)*plotRatios[['ins']]+plotRatios[['genes']]+length(chipDatasets)*plotRatios[['chip']])*width)	    

	    par(mai=c(bM,lM,tM,rM), xaxs="i", yaxs="i", family='sans', cex=1)	    

	    layout(matrix(1:nPlots, ncol=1), height=c(rep(plotRatios[['hic']],length(hicDatasets))*band, rep(plotRatios[['genes']],1), rep(plotRatios[['virt4C']],length(allViewPoints)), rep(plotRatios[['ins']],length(insTracks)), rep(plotRatios[['chip']],length(chipDatasets)), width=2, respect=T))


            if(Capture == 0)
            {



		interval2D <- gintervals.2d(chrom1,start1-extension,end2+extension,chrom1,start1-extension,end2+extension)
		print(paste0("Interval2D ",paste(interval1D,collapse=" ")))
		#quit()
		binnedIterator <- giterator.intervals(intervals=interval2D, iterator=c(chicRes,chicRes))

		for(chicTrack in gtrack.ls('hic.','hicScores_',paste0(celltype,"_",condition),chrom1,"k250"))
		{
        	    print(chicTrack)
	    	    if(length(chicTrack) == 0)
	    	    {
		        print(chicTrack)
		    	next
		    }

		    #		    XXX

		    if(!is.null(chicTrack))
		    {
	
			print("##### Plot squared matrix #####")
			data <- gextract(chicTrack, binnedIterator, iterator=binnedIterator,colnames = c('score'))

			plotData <- data
			print(head(plotData))
			print(tail(plotData))
			plotData$score[is.na(plotData$score)] = 0

			mtx    <- acast(data,start1~start2,value.var="score")

			#counts <- gextract(track,g2d(intrv),colnames="score")
			#regReads <- nrow(counts)/2

			image(mtx, main="", useRaster=TRUE, xlab="", ylab="", axes=FALSE, cex.main=4, col=scoreCol, zlim=c(zmin,zmax))

			# Annotation on the diagonal
			#plotRectangles2D(gintervals.2d(chrom1,start1,end1,chrom1,start1,end1),interval2D,lwd=1,image=TRUE,col='black')
			#plotRectangles2D(gintervals.2d(chrom2,start2,end2,chrom2,start2,end2),interval2D,lwd=1,image=TRUE,col='black')

			# TADs annotation
			for(d in 1:nrow(domains))
			{
			    #print(domains[d,])
			    chrom1D <- domains[d,]$chrom1
			    start1D <- domains[d,]$start1
			    end1D   <- domains[d,]$end1
			    chrom2D <- domains[d,]$chrom2
			    start2D <- domains[d,]$start2
			    end2D   <- domains[d,]$end2

			    if(chrom1D != chrom1)
			    {
			        next
			    }
			    if((start1 <= start1D && start1D <= end1) || (start1 <= end1D && end1D <= end1))
			    {
			    	domainAnnotation <- gintervals.2d(chrom1D,start1D,end1D,chrom2D,start2D,end2D)
			    	print(domainAnnotation)
			    	print(interval2D)
				print(domains[d,]$value)

				colorD <- "red"
				if(length(grep("PcG",domains[d,]$value))  != 0)
				{
				    colorD <- "blue"
				}
				if(length(grep("Null",domains[d,]$value)) != 0)
				{
				    colorD <- "black"
				}
				if(length(grep("Het",domains[d,]$value)) != 0)
				{
				    colorD <- "green"
				}								
				print(colorD)
			    	plotRectangles2D(domainAnnotation,interval2D,lwd=3,image=TRUE,col=colorD)
				
			    }
			}
		    }
		    ### Plot loops on the map ###
		    
		    #data <- read.table(paste0("./scripts/loops_",celltype,"_",condition,"_10000bp.tsv"),header=T)
		    #data <- as.data.frame(data)
		    #print(head(data))
		    #col <- "yellow"
		    #lwd = 4
		    #lty = 1
		    #if(nrow(data) == 0)
		    #{
		    #    for(g in 1:nrow(data))
		    #    {			
			#    row <- data[g,]
			 #   if(row$start1 != row$start2)
			 #   {
		#		if ( interval1D$chrom == row$chrom1 & (interval1D$start < row$start1 & row$start1 < interval1D$end) & (interval1D$start < row$start2 & row$start2 < interval1D$end))
		#		{				
	    	 #                   print(row)
		#		    x1 <- row$start1 + (row$start2-row$start1)/2
		#		    x2 <- row$start1 + (row$end2-row$start1)/2
		#		    x3 <- row$end1 + (row$end2-row$end1)/2
		#		    x4 <- row$end1 + (row$start2-row$end1)/2

		#		    y1 <- row$start2-row$start1
		#		    y2 <- row$end2-row$start1
		#		    y3 <- row$end2-row$end1
		#		    y4 <- row$start2-row$end1


		#		    segments(x1,y1,x2,y2,col=col,lwd=lwd,lty=lty)
		#		    segments(x2,y2,x3,y3,col=col,lwd=lwd,lty=lty)

		#		    segments(x3,y3,x4,y4,col=col,lwd=lwd,lty=lty)
		#		    segments(x4,y4,x1,y1,col=col,lwd=lwd,lty=lty)
		#		}
		#	    }
		#	}
		#    }

		}
	    }
		    
	    if(Genes == 0)
	    {
	        print("# Plot gene annotation")
	        ### Genes
	        ### Load TSS and GENE coordinates...
	        tssCoordinates <- gintervals.load('intervals.ucscCanTSS') # Ritorna un dataframe
	        rownames(tssCoordinates) <- tssCoordinates$geneName

		#TADcoordinates <- gintervals.load('domains.TAD_hic_LE_WTref_merge_dm6_YOVL_r10kb')
		#print(TADcoordinates)

		geneCoordinates <- gintervals.load('intervals.ucscCanGenes')
#		geneCoordinates <- geneCoordinates[geneCoordinates$geneName == geneName,]		
#		print(head(geneCoordinates))

		rownames(geneCoordinates) <- geneCoordinates$geneName
	        genes <- gintervals.neighbors(geneCoordinates,interval1D) # Non e' 100% sicuro che funzioni bene. E' un intersect tra intervalli e geni presenti.
	        TSSs  <- gintervals.neighbors(tssCoordinates,interval1D) # Non e' 100% sicuro che funzioni bene. E' un intersect tra intervalli e geni presenti.
                # Puoi impostare il numero di primi vicini...o usare massimo distance minimum distance: centre+/-2kb
         	                     # Usa un altro exact_neighbors che e' basato su stringhe!

		genesPlot(genes,plotLim=c(interval1D$start,interval1D$end),cex=1.0,rHeight=90)
	    }

	    if(Insulation == 0)
	    {
		if(length(insTracks) > 0)
	        {
		    print(head(insData))
	            print("# Plot insulation tracks")

		    plotIns(insData,plotOrder=insOrder,plotCols=insCols,plotTitles=insTitles,plotStart=interval1D$start,plotEnd=interval1D$end,chipRes=insRes,plotMax=c(1.0),chrom=chrom1,main=FALSE)

		    #		    plotIns(insData, chipStart=NA, plotOrder=NA, plotCols=NA, plotTitles=NA, vertical=FALSE, chipRes=insRes, plotAnn=F, plotStart=F, plotEnd=F, inputType='avg', cex=1, plotMax=1.0, main=TRUE, adj = 0.5, plotAxis=TRUE, CTCFDir=NA, peaksDir=NA, chrom=NA, condition=NA, peakSet=NA)

		    #plot(insRefData$mean, cex.main=1.2, type='l', xlab='', ylab='', yaxt='n', xaxt='n', lwd=1.5, col="black", border="black", axes=F, main="", xpd = F, space=0, ylim=c(0,1),xlim=c(start1/1000,end1/1000))
		    #polygon(c(1:length(insRefData$mean), length(insRefData$mean):1), c(insRefData$mean-insRefData$stddev,rev(insRefData$mean+insRefData$stddev)), col = rgb(0,0,0, max=255, alpha=50), border=F)

		    #axis(1, at=c(1,length(insRefData$mean)), labels=c("",""), cex.axis=1, las=1)
		    #axis(1, at=c(0,1), labels=c(0,signif(1.0,2)), cex.axis=1, las=3)
		    #axis(2, at=c(0,0.5,1), labels=c(0,0.5,signif(1.0,2)), cex.axis=1, las=1)


		    #if(condition != refCondition)
		    #{
   		    #    lines(insData$mean, type='l', col="blue", lwd=3)
			#polygon(c(1:length(insData$mean), length(insData$mean):1), c(insData$mean-insData$stddev,rev(insData$mean+insData$stddev)), col = rgb(0,0,1, max=255, alpha=50), border=F)
		    #}



 	       	}
	    }

	    if(ViewPoints == 0)
            {
	        ##### Plot virtual4C track
	        virt4CData <- c()
	        for(viewPoint in names(allViewPoints))
	        {
	            print(viewPoint)
	    	    locus <- allViewPoints[[viewPoint]]
	 	    locus <- unlist(locus)
		    print(locus)

	  	    interval2Dvirt4C <- gintervals.2d(chrom, start, end, locus[1], as.numeric(locus[2]), as.numeric(locus[3]))
		    print(interval2Dvirt4C)
		    binnedIterator <- data.frame(matrix(ncol = 6, nrow = as.integer((end-start)/virt4CRes)+2))
		    colnames(binnedIterator) <- c("chrom1","start1","end1","chrom2","start2","end2")
		    binnedIterator$chrom1	 <- chrom
		    binnedIterator$chrom2	 <- chrom
		    print(colnames(binnedIterator))
		    nbin = 0
		    s1 = as.integer(start/virt4CRes+1)*virt4CRes
		    s2 = as.integer(end/virt4CRes)*virt4CRes
		    resolution1 = virt4CRes
		    resolution2 = as.numeric(locus[3])-as.numeric(locus[2])
		    print(s1)
		    print(s2)		
		    print(resolution1)
		    print(resolution2)		
		
		if(s1>start)
		{
		    nbin = nbin+1
		    binnedIterator$start1[nbin] <- start
		    binnedIterator$end1[nbin]   <- s1
		    binnedIterator$start2[nbin] <- as.numeric(locus[2])
		    binnedIterator$end2[nbin] <- as.numeric(locus[3])			    		    
		}

		for(binStart1 in seq(from=s1, to=s2-resolution1, by=resolution1))
		{
		    for(binStart2 in seq(from=as.numeric(locus[2]), to=as.numeric(locus[3])-1, by=resolution2))
		    {
			nbin = nbin+1
			#print(paste(nbin,binStart1,binStart1+resolution1,binStart2,binStart2+resolution2))
			binnedIterator$start1[nbin] <- binStart1
			binnedIterator$end1[nbin]   <- binStart1+resolution1
			binnedIterator$start2[nbin] <- binStart2
			binnedIterator$end2[nbin]   <- binStart2+resolution2-1
		    }
		}
		if(s2<end)
		{
		    nbin = nbin+1
		    binnedIterator$start1[nbin] <- s2
		    binnedIterator$end1[nbin]   <- end
		    binnedIterator$start2[nbin] <- as.numeric(locus[2])
		    binnedIterator$end2[nbin]   <- as.numeric(locus[3])			    		    
		}
		#print(binnedIterator)
		#binnedIterator <- giterator.intervals(intervals=interval2Dvirt4C, iterator=c(virt4CRes,as.numeric(locus[3])-as.numeric(locus[2])))
		#quit()

		virt4CTrack <- gextract(chicTrack, binnedIterator, iterator=binnedIterator, colnames= c(viewPoint))
		#print(virt4CTrack)
		print(length(virt4CTrack))
		print(nrow(virt4CTrack))

		if(is.null(virt4CData))
		{
		    virt4CData <- virt4CTrack[,c("chrom1","start1","end1",viewPoint)]
		    } else {
		    virt4CData <- cbind(virt4CData,virt4CTrack[c(viewPoint)])
		}

		    virt4CData[is.na(virt4CData)] <- 0
		}

		outfile <- paste0("./virtual4C_shamanScores_Zfp608_",celltype,"_",gsub("_","",condition),"_chic_",chicRes,"bp.tsv")	    	    
		virt4CData[is.na(virt4CData)] <- 0
		print(virt4CData)	    

		print(head(virt4CData))	

		plotAnn=allViewPoints
		print(plotAnn)

		if(condition == "NOnub")
		{
		    plotVirt4C(virt4CData,plotOrder=plotOrdervirt4C,plotCols=plotColsvirt4C,plotTitles=plotTitlesvirt4C,plotStart=interval1D$start,plotEnd=interval1D$end,chipRes=virt4CRes,plotMax=plotMaxvirt4C,chrom=chrom,main=FALSE)
		} else {
		    plotVirt4C(virt4CData,plotOrder=plotOrdervirt4C,plotCols=plotColsvirt4C,plotTitles=plotTitlesvirt4C,plotStart=interval1D$start,plotEnd=interval1D$end,chipRes=virt4CRes,plotMax=plotMaxvirt4C,chrom=chrom,main=FALSE,plotAnn=plotAnn)
		}
	    }

            if(ChIPseqs == 0)
            {
                if(length(chipTracks) > 0)
                {
                    chipData <- gextract(chipTracks,intervals=interval1D, iterator=ChIPres, colnames = gsub('_merge_|dm6_|NA|chipseq.','',chipTracks))
		    print("# Plot ChIPseq tracks")
                    print(chipTracks)
  	            print(head(chipData))
		    for(i in 4:length(colnames(chipData)))
		    {
			print(colnames(chipData)[i])
		        print(paste0("Min ",min(chipData[,i])))
		        print(paste0("Max ",max(chipData[,i])))
		    }
		    #quit()
		    #if(condition == "NOnub")
		    #{
		        plotChIP(chipData,plotOrder=plotOrder,plotCols=plotCols,plotTitles=plotTitles,plotStart=interval1D$start,plotEnd=interval1D$end,chipRes=ChIPres,plotMax=plotMax,chrom=chrom,main=FALSE,plotMin=plotMin)			    
		    #} else {			    
		    #    plotAnn = plotAnn[plotAnn$chrom == chrom,]
		    #    print(plotAnn)

                    #    CTCFDir <- "/work/cavalli/mdistefano/2021_04_13_Project_with_Ivana/CTCF_peaks_directionality/pwmscan_mm10_19159_19671_all_possible_CTCF_motifs.bed"
		    #    plotChIP(chipData,plotOrder=plotOrder,plotCols=plotCols,plotTitles=plotTitles,plotStart=interval1D$start,plotEnd=interval1D$end,chipRes=ChIPres,plotMax=plotMax,chrom=chrom,main=FALSE,plotMin=plotMin,plotAnn=plotAnn)
		    #}
            }
	}
        dev.off()

    }
}
#	quit()
    }
}

warnings()
quit()
