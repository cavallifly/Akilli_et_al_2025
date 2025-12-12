library(Gviz)
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)

options(ucscChromosomeNames = FALSE)

# ---- Load your local annotation ----
# (Example: a downloaded GTF from Ensembl or FlyBase)
#gtf_file <- "dm6.refGene.gtf"
gtf_file <- "dmel-all-r6.36.gtf"

# Import into R
gtf <- rtracklayer::import(gtf_file)
gtf <- gtf[gtf$type %in% c("gene"), ]

data(cpgIslands)
class(cpgIslands)
print(head(cpgIslands))

# ---- Create the Gviz tracks ----
genome <- "dm6"
chr <- "chr2L"
regionStart <- 16330000
regionEnd   <- 16620000

# Filter for the dac locus
region <- subset(
  gtf,
  seqnames == chr & start >= regionStart & end <= regionEnd & (end-start) > 8000
)
class(gtf)
print(head(gtf))
#range(region)

ideoTrack <- IdeogramTrack(genome = genome, chromosome = chr)
axisTrack <- GenomeAxisTrack()

atrack <- AnnotationTrack(region, name = "Genes", genome = genome, chromosome = gsub("chr","",chr))

geneTrack <- GeneRegionTrack(
  region,
  genome = genome,
  chromosome = chr,
  start      = regionStart,
  end        = regionEnd,
  name = "dm6 Genes",
  exon = TRUE,
  strand = TRUE,
  showId = TRUE,  
  geneSymbol = TRUE,
  symbol = region$gene_symbol,
  shape = "arrow",
  arrowHeadMaxWidth=10,
  ucscChromosomeNames=FALSE,
  fontsize.group=24,
  stackHeight=0.3,
  stacking = "squish"
)


# ---- Plot to PDF ----
pdf("local_annotation.pdf", width = 10, height = 6)

p <- plotTracks(list(ideoTrack, axisTrack, geneTrack),
           from = regionStart, to = regionEnd,
           main = "dm6 gene annotation")
print(p)
dev.off()

quit()

data(cpgIslands)
class(cpgIslands)
## [1] "GRanges"
## attr(,"package")
## [1] "GenomicRanges"
chr <- as.character(unique(seqnames(cpgIslands)))
gen <- genome(cpgIslands)
atrack <- AnnotationTrack(cpgIslands, name = "CpG")

# ---- Plot to PDF ----
pdf("local_annotation.pdf", width = 10, height = 6)
plotTracks(atrack)
dev.off()
quit()




quit()
# ---- Define region ----
genome <- "dm6"
chr <- "2L"
start <- regionStart
end <- regionEnd

# Create the genome axis
axisTrack <- GenomeAxisTrack()

# Gene annotation track from UCSC
refseqTrack <- UcscTrack(
  genome = genome,
  chromosome = chr,
  track = "NCBI RefSeq",
  table = "ncbiRefSeq",
  trackType = "GeneRegionTrack",
  rstarts = start,
  rends = end,
  from = start,
  to = end,
  name = "RefSeq Genes",
  showId = TRUE,
  geneSymbol = TRUE
)

# Ensembl transcript-level track
transcriptTrack <- UcscTrack(
  genome = genome,
  chromosome = chr,
  track = "Ensembl Genes",
  table = "ensGene",
  trackType = "GeneRegionTrack",
  from = start, to = end,
  name = "Ensembl Transcripts",
  showId = TRUE,
  geneSymbol = TRUE
)

# ✅ FlyBase annotations (alternative, often richer for dm6)
#flybaseTrack <- UcscTrack(
#  genome = genome,
#  chromosome = chr,
#  track = "FlyBase Genes",
#  table = "flyBaseGene",
#  trackType = "GeneRegionTrack",
#  from = start, to = end,
#  name = "FlyBase Genes",
#  showId = TRUE,
#  geneSymbol = TRUE
#)

# Optional: Add an ideogram for chromosome context
ideoTrack <- IdeogramTrack(genome = genome, chromosome = chr)

pdf("test.pdf")
# Combine and plot
plotTracks(
  list(ideoTrack, axisTrack, refseqTrack, transcriptTrack),
  from = start, to = end,
  main = paste0(chr,":",start,"-",end," (Drosophila melanogaster, dm6)")
)

#print(p)
dev.off()
