# For Venn-diagram
library(dplyr)
library(ggvenn)
library(IRanges)
library(GenomicRanges)

require(doParallel)

### misha working DB
mDBloc <-  '/zdata/data/mishaDB/trackdb/'
db <- 'dm6'
dbDir <- paste0(mDBloc,db,'/')

source("/zdata/data/auxFunctions/auxFunctions.R")
options(scipen=20,gmax.data.size=0.5e8,shaman.sge_support=1)

domainRes  <-  5000
isRes      <-  400
window     <- 25000
tolerance <-  2500

#chrom	start	end	insnoNub1	insnoNub2	insSUMOnub1	insSUMOnub2
print(paste0("First set of TopDom domains"))
TopDomDomainsFile1 <- "TopDom_domains_hic_WD_NOnub_merge_dm6_NA_window_5_at_5000bp_filtered_refined_at_1000bp_with_WD_NOnub_and_SUMOnub_merge_insulation.tsv"
TopDomDomains1 	   <- read.table(TopDomDomainsFile1,header=T,fill=TRUE)
TopDomBorders11    <- TopDomDomains1[,c(1,2,4,6)]
colnames(TopDomBorders11) <- c("chrom","start","insNub","insSUMO")
TopDomBorders12    <- TopDomDomains1[,c(1,3,5,7)]
colnames(TopDomBorders12) <- c("chrom","start","insNub","insSUMO")
TopDomBorders1    <- rbind(TopDomBorders11,TopDomBorders12)
TopDomBorders1    <- bind_cols(as.character(TopDomBorders1$chrom),as.integer(TopDomBorders1$start-tolerance),as.integer(TopDomBorders1$start+tolerance),as.numeric(TopDomBorders1$insNub),as.numeric(TopDomBorders1$insSUMO))
colnames(TopDomBorders1)  <- c("chrom","start","end","insNub","insSUMO")
TopDomBorders1    <- unique(TopDomBorders1)
TopDomBorders1    <- as.data.frame(TopDomBorders1)
head(TopDomDomains1)
print(TopDomBorders1)
tail(TopDomBorders1)
print(paste0("Total number of borders in set 1 ", nrow(TopDomBorders1)))

print(paste0("Second set of TopDom domains"))
TopDomDomainsFile2 <- "TopDom_domains_hic_WD_SUMOnub_merge_dm6_NA_window_5_at_5000bp_filtered_refined_at_1000bp_with_WD_NOnub_and_SUMOnub_merge_insulation.tsv"
TopDomDomains2 	  <- read.table(TopDomDomainsFile2,header=T,fill=TRUE)	
TopDomBorders21    <- TopDomDomains2[,c(1,2,4,6)]
colnames(TopDomBorders21) <- c("chrom","start","insNub","insSUMO")
TopDomBorders22    <- TopDomDomains2[,c(1,3,5,7)]
colnames(TopDomBorders22) <- c("chrom","start","insNub","insSUMO")
TopDomBorders2    <- rbind(TopDomBorders21,TopDomBorders22)
TopDomBorders2    <- bind_cols(as.character(TopDomBorders2$chrom),as.integer(TopDomBorders2$start-tolerance),as.integer(TopDomBorders2$start+tolerance),as.numeric(TopDomBorders2$insNub),as.numeric(TopDomBorders2$insSUMO))
colnames(TopDomBorders2)  <- c("chrom","start","end","insNub","insSUMO")
TopDomBorders2    <- unique(TopDomBorders2)
TopDomBorders2    <- as.data.frame(TopDomBorders2)
head(TopDomDomains2)
head(TopDomBorders2)
tail(TopDomBorders2)
print(paste0("Total number of borders in set 2 ", nrow(TopDomBorders2)))

print(paste0("1. Venn-diagram to show the overlap of the borders."))
# Convert to GenomicRanges objects (for overlap calculations)
gr1 <- makeGRangesFromDataFrame(TopDomBorders1[,c(1,2,3)], keep.extra.columns = TRUE)
gr2 <- makeGRangesFromDataFrame(TopDomBorders2[,c(1,2,3)], keep.extra.columns = TRUE)
gr1 <- reduce(gr1)
gr2 <- reduce(gr2)
print(paste0("Total number of borders in set 1 ", length(gr1)))
print(paste0("Total number of borders in set 2 ", length(gr2)))

# Find overlaps (intersection) and unique regions
overlap <- findOverlaps(gr1, gr2)
print(length(overlap))
# Extract overlapping indices
query_hits   <- queryHits(overlap)   # Indices in gr1
subject_hits <- subjectHits(overlap) # Indices in gr2

print(head(gr1[query_hits]))
print(length(gr1[query_hits]))

print(head(gr2[subject_hits]))
print(length(gr2[subject_hits]))

# Create Venn diagram
df1 <- c(paste0("common_",seq(1,length(overlap),1)), paste0("gr1_",seq(length(overlap)+1,length(gr1),1)))
df2 <- c(paste0("common_",seq(1,length(overlap),1)), paste0("gr2_",seq(length(overlap)+1,length(gr2),1)))
print(paste0("Total number of borders in set 1 ", length(df1)))
print(paste0("Total number of borders in set 2 ", length(df2)))
venn_list <- list(
  "Control"   = df1,
  "SUMO RNAi" = df2
)

# Plot with ggvenn
pdf("Venn_diagram.pdf")
p <- ggvenn(
    venn_list,
    fill_color = c("#1f77b4", "#ff7f0e", "#2ca02c"),
    stroke_size = 0.5,
    set_name_size = 10,   # Increase set name size
    text_size = 7,         # Increase number size
    auto_scale = TRUE
    )
print(p)
dev.off()
print(paste0("End Venn-diagram analysis"))

print(paste0("1.1 Write a list of the condition-specific borders."))
specific_to_condition1 <- gr1[-query_hits]
specific_to_condition2 <- gr2[-subject_hits]

# Convert to data frames
df_condition1 <- data.frame(
  chrom = as.character(seqnames(specific_to_condition1)),
  start = start(specific_to_condition1),
  end = end(specific_to_condition1),
  stringsAsFactors = FALSE
)

df_condition2 <- data.frame(
  chrom = as.character(seqnames(specific_to_condition2)),
  start = start(specific_to_condition2),
  end = end(specific_to_condition2),
  stringsAsFactors = FALSE
)

# Write to files
write.table(
  df_condition1,
  file = "TopDom_domains_hic_WD_NOnub_merge_dm6_NA_window_5_at_5000bp_filtered_refined_at_1000bp_specificBorders.bed",
  sep = "\t",  # Tab-separated
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

write.table(
  df_condition2,
  file = "TopDom_domains_hic_WD_SUMOnub_merge_dm6_NA_window_5_at_5000bp_filtered_refined_at_1000bp_specificBorders.bed",
  sep = "\t",  # Tab-separated
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

print(paste0("2. Violin plots of the insulation values in Control (WD_NOnub) vs. SUMO RNAi (SUMOnub)."))

print(paste0("2.1 Violin plots of the insulation values in Control (WD_NOnub) vs. SUMO RNAi (SUMOnub) with lines connecting corresponding borders."))

print(paste0("2.2 Violin plots of the insulation values in Control (WD_NOnub) vs. SUMO RNAi (SUMOnub) with lines connecting corresponding borders with the most changing borders in red/blue."))

