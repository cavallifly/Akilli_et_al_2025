library(InteractionSet)
library(GenomicRanges)
library(HiCExperiment)
library(HiContactsData)
library(TopDom)
library(dplyr)
library(tidyr)
library(purrr)
options(scipen=999)

args = commandArgs(trailingOnly=TRUE)

# Code elaborated from https://bioconductor.org/books/devel/OHCA/pages/interoperability.html
HiCExperiment2TopDom <- function(hic, chr) {
    data <- list()
    cm <- as(hic[chr], 'ContactMatrix')
    data$counts <- as.matrix(cm) |> base::as.matrix()
    data$counts[is.na(data$counts)] <- 0
    data$bins <- regions(cm) |> 
    as.data.frame() |> 
    select(seqnames, start, end) |>
    mutate(seqnames = as.character(seqnames)) |>
    mutate(id = 1:n(), start = start - 1) |> 
    relocate(id) |> 
    dplyr::rename(chr = seqnames, from.coord = start, to.coord = end)
    class(data) <- 'TopDomData'
    return(data)
}

condition  <- args[1]
resolution <- as.numeric(args[2])


for(mcoolFile in list.files(path="../../01_cool_files/",pattern=paste0(condition,".*merge.*.mcool$"),full.names = TRUE))
{
    print(mcoolFile)
    name <- gsub(".mcool","",basename(mcoolFile))

    for(chrom in c("chr2L","chr2R","chr3L","chr3R","chrX"))
    {

        hic <- zoom(import(mcoolFile, format = 'mcool', focus=chrom), resolution)
        print(paste0("Working on chromosome ",chrom))
        hic_topdom <- HiCExperiment2TopDom(hic, chrom)

	for(window in c(3,4,5,6,7,9,12,15,20,25,30))
        {
	    print(paste0("Working on window ",window))
	    checkFile <- paste0(name,"_",chrom,"_window_",window,"_at_",resolution,"bp.bed")
	    if(file.exists(checkFile))
	    {
	        print(paste0("File ",checkFile," exists!"))
	        next
	    }
            outFile <- paste0("TopDom_domains_",name,"_window_",window,"_at_",resolution,"bp.tsv")	    
            print(outFile)

	    domains <- TopDom::TopDom(hic_topdom, window.size = window, outFile=paste0(name,"_",chrom,"_window_",window,"_at_",resolution,"bp"))

	    topologicalFeatures(hic, 'domain') <- domains$bed |> 
                mutate(chromStart = chromStart + 1) |> 
    	        filter(name == 'domain') |> 
    	        makeGRangesFromDataFrame()
    	        tf <- topologicalFeatures(hic, 'domain')
            write.table(tf,outFile,sep="\t",quote=F,append=T)
        }
    }
}
# Check dac boundary #grep chr2L TopDom_domains_hic_WD_??01_Rep* | sort -k 3,3n | grep 16350001