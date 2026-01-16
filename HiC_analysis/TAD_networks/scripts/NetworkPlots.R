library(tidyr)
library(GenomicRanges)
library(dplyr)
library(igraph)
library(randomForest)
library(ranger)
library(caret)
library(ggplot2)
library(pROC)
library(rtracklayer)
library(MLmetrics)
#library(smotefamily) # For handling class imbalance for ML


###Get HiC Scores

scoresSUMO_NO_3L<-read.delim("/Users/nazliakilli/Desktop/HiC_scores/avg_scores_trans1Dinterval_chr3L_SUMOnub_NOnub_all_domains_vs_all_domains.tab",header = FALSE)
scoresSUMO_NO_2L<-read.delim("/Users/nazliakilli/Desktop/HiC_scores/avg_scores_trans1Dinterval_chr2L_SUMOnub_NOnub_all_domains_vs_all_domains.tab",header = FALSE)
scoresSUMO_NO_2R<-read.delim("/Users/nazliakilli/Desktop/HiC_scores/avg_scores_trans1Dinterval_chr2R_SUMOnub_NOnub_all_domains_vs_all_domains.tab",header = FALSE)
scoresSUMO_NO_3R<-read.delim("/Users/nazliakilli/Desktop/HiC_scores/avg_scores_trans1Dinterval_chr3R_SUMOnub_NOnub_all_domains_vs_all_domains.tab",header = FALSE)
scoresSUMO_NO_X<-read.delim("/Users/nazliakilli/Desktop/HiC_scores/avg_scores_trans1Dinterval_chrX_SUMOnub_NOnub_all_domains_vs_all_domains.tab",header = FALSE)

colnames(scoresSUMO_NO_3L)<-c("Chr1","start1","end1","Chr2","start2","end2","interval1","interval2","av_score","sample")
colnames(scoresSUMO_NO_3R)<-c("Chr1","start1","end1","Chr2","start2","end2","interval1","interval2","av_score","sample")
colnames(scoresSUMO_NO_2L)<-c("Chr1","start1","end1","Chr2","start2","end2","interval1","interval2","av_score","sample")
colnames(scoresSUMO_NO_2R)<-c("Chr1","start1","end1","Chr2","start2","end2","interval1","interval2","av_score","sample")
colnames(scoresSUMO_NO_X)<-c("Chr1","start1","end1","Chr2","start2","end2","interval1","interval2","av_score","sample")
#Domains for ChiaSig

###ChiaSig is the all NOnub and all SUMOnub domains

ChiaSig<-rbind(scoresSUMO_NO_3L,scoresSUMO_NO_3R,scoresSUMO_NO_2R,scoresSUMO_NO_2L,scoresSUMO_NO_X)
ChiaSig$av_score<-ChiaSig$av_score+100
head(ChiaSig)
ChiaSig_Control<-ChiaSig[ChiaSig$sample=="WD_NOnub",]
ChiaSig_SUMO<-ChiaSig[ChiaSig$sample=="WD_SUMOnub",]

###Gett All domain info

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

# Load CNR

Control_H3K27me3<-import.bw("/Users/nazliakilli/Desktop/DiffBind/cnr_H3K27me3_WD_NOnub_merge_all_dm6_NA/cnr_H3K27me3_WD_NOnub_merge_all_dm6_NA_normRPKM.bw")
head(Control_H3K27me3)
SUMO_H3K27me3<-import.bw("/Users/nazliakilli/Desktop/DiffBind/cnr_H3K27me3_WD_SUMOnub_merge_all_dm6_NA/cnr_H3K27me3_WD_SUMOnub_merge_all_dm6_NA_normRPKM.bw")
head(SUMO_H3K27me3)
head(Control_H3K27me3)

Control_Pc<-import.bw("/Users/nazliakilli/Desktop/DiffBind/cnr_PcXlinked_WD_NOnub_merge_all_dm6_NALF2/cnr_PcXlinked_WD_NOnub_merge_all_dm6_NALF2_normRPKM.bw")
head(Control_Pc)
SUMO_Pc<-import.bw("/Users/nazliakilli/Desktop/DiffBind/cnr_PcXlinked_WD_SUMOnub_merge_all_dm6_NALF2/cnr_PcXlinked_WD_SUMOnub_merge_all_dm6_NALF2_normRPKM.bw")
head(SUMO_Pc)

Control_H2Aub118<-import.bw("/Users/nazliakilli/Desktop/DiffBind/cnr_H2Aub118_WD_NOnub_merge_all_dm6_NALF2/cnr_H2Aub118_WD_NOnub_merge_all_dm6_NALF2_normRPKM.bw")
head(Control_H2Aub118)
SUMO_H2Aub118<-import.bw("/Users/nazliakilli/Desktop/DiffBind/cnr_H2Aub118_WD_SUMOnub_merge_all_dm6_NALF2/cnr_H2Aub118_WD_SUMOnub_merge_all_dm6_NALF2_normRPKM.bw")
head(SUMO_H2Aub118)

Control_H3K27ac<-import.bw("/Users/nazliakilli/Desktop/DiffBind/cnr_H3K27ac_WD_NOnub_merge_all_dm6_NA/cnr_H3K27ac_WD_NOnub_merge_all_dm6_NA_normRPKM.bw")
head(Control_H3K27ac)
SUMO_H3K27ac<-import.bw("/Users/nazliakilli/Desktop/DiffBind/cnr_H3K27ac_WD_SUMOnub_merge_all_dm6_NA/cnr_H3K27ac_WD_SUMOnub_merge_all_dm6_NA_normRPKM.bw")
head(SUMO_H3K27ac)

# ==================================================
# Enhanced TAD Contact Pipeline with ML + Rigor
# This version adds new features and handles class imbalance
# ==================================================


# NOTE: This script assumes you have already loaded the following data:
# AllDomainsSN
# Control_H3K27me3, SUMO_H3K27me3
# Control_H3K27ac, SUMO_H3K27ac
# Control_H2Aub118, SUMO_H2Aub118
# Control_Pc, SUMO_Pc
# ChiaSig_Control, ChiaSig_SUMO
head(ChiaSig_Control)


# Rename the existing RNAseq_GR object to RNA_GR for consistency
# This is a good practice to match the variable name in the main script.
RNA<-read.table("/Users/nazliakilli/Desktop/RNAseq/DEseq2_results_WD_SUMOnub_vs_WD_NOnub_table.tsv",header=T)
RNAseq_GR<-GenomicRanges::makeGRangesFromDataFrame(RNA,seqnames.field="chrom",
                                                   start.field='start',end.field="end",keep.extra.columns = TRUE,na.rm=TRUE)
RNA_GR <- RNAseq_GR
head(RNA_GR)
RNAseq_GR[is.na(RNAseq_GR$log2Ratio),]
# Add a new column named 'ExpressionValue' to the RNA_GR object,
# using the 'log2Ratio' values as the expression measure.
mcols(RNA_GR)$ExpressionValue <- mcols(RNA_GR)$log2Ratio
RNA_GR<-na.omit(RNA_GR)
# Check the new column to confirm it's been created correctly
head(mcols(RNA_GR)$ExpressionValue)
# RNA-seq data as a GenomicRanges object.
# Assumed to be created from a data frame like this:
# RNA_GR <- GenomicRanges::makeGRangesFromDataFrame(RNA_data, seqnames.field="chrom", start.field='start', end.field="end", keep.extra.columns = TRUE, na.rm=TRUE)
# Make sure the 'RNA_data' data frame has a column with expression values.
# For example, a column named 'log2FoldChange' or similar.
# For this script, we'll assume the expression value column is named 'ExpressionValue'.
# YOU MAY NEED TO CHANGE THE 'signal_col' ARGUMENT BELOW.

# ==================================================
# 1. Function to compute average ChIP per TAD
# ==================================================
compute_avg_score <- function(tads_gr, chip_gr, signal_col = NULL) {
  if (is.null(signal_col)) {
    num_cols <- which(sapply(mcols(chip_gr), is.numeric))
    if (length(num_cols) == 0) stop("No numeric columns in chip_gr for signal!")
    signal_col <- num_cols[1]
  }
  signal_values <- mcols(chip_gr)[[signal_col]]
  hits <- findOverlaps(tads_gr, chip_gr)
  avg_scores <- rep(0, length(tads_gr))
  if (length(hits) > 0) {
    avg_values <- tapply(signal_values[subjectHits(hits)], queryHits(hits), mean, na.rm = TRUE)
    avg_scores[as.integer(names(avg_values))] <- avg_values
  }
  return(avg_scores)
}

# ==================================================
# 2. Annotate TADs with histone and RNA-seq scores
# ==================================================
AllDomainsSN$avgH3K27me3Control <- compute_avg_score(AllDomainsSN, Control_H3K27me3)
AllDomainsSN$avgH3K27me3SUMO    <- compute_avg_score(AllDomainsSN, SUMO_H3K27me3)
AllDomainsSN$avgH3K27acControl  <- compute_avg_score(AllDomainsSN, Control_H3K27ac)
AllDomainsSN$avgH3K27acSUMO     <- compute_avg_score(AllDomainsSN, SUMO_H3K27ac)
AllDomainsSN$avgH2Aub118Control <- compute_avg_score(AllDomainsSN, Control_H2Aub118)
AllDomainsSN$avgH2Aub118SUMO    <- compute_avg_score(AllDomainsSN, SUMO_H2Aub118)
AllDomainsSN$avgPcControl       <- compute_avg_score(AllDomainsSN, Control_Pc)
AllDomainsSN$avgPcSUMO          <- compute_avg_score(AllDomainsSN, SUMO_Pc)




# NEW: Annotate TADs with average gene expression (on the AllDomainsSN GRanges object).
AllDomainsSN$avg_gene_expression <- compute_avg_score(AllDomainsSN, RNA_GR, signal_col = "ExpressionValue")

# NEW: Annotate TADs with gene count (on the AllDomainsSN GRanges object).
AllDomainsSN$num_genes <- countOverlaps(AllDomainsSN, RNA_GR)
head(AllDomainsSN)



# ==================================================
# 3. CONVERT TO DATAFRAME & VALIDATE UNIQUENESS (FIXED)
# ==================================================
# Convert the fully annotated GRanges object to a data frame.
df_AllDomainsSN <- as.data.frame(AllDomainsSN, row.names = NULL)


# --- Apply the FIXED logic to create IDs and Types ---

# 1. FIXED: Create DomainID from genomic coordinates (the unique identifier).
df_AllDomainsSN$DomainID <- paste(
  df_AllDomainsSN$seqnames, 
  df_AllDomainsSN$start, 
  df_AllDomainsSN$end, 
  sep = ":"
)

# 2. FIXED: Proper TAD classification: Create TADType from the 'Type' column
cat("\n=== Checking original Type values ===\n")
type_summary <- table(df_AllDomainsSN$Type)
cat("Unique Type values and their counts:\n")
print(type_summary)

# CRITICAL FIX 2: Extract TADType (e.g., "Null545" -> "Null")
df_AllDomainsSN$TADType <- gsub("[0-9]+$", "", as.character(df_AllDomainsSN$Type))

# Verify the classification worked correctly
cat("\n=== TADType classification results ===\n")
cat("TADType distribution:\n")
print(table(df_AllDomainsSN$TADType))

# Check distribution across chromosomes for ALL TAD types
cat("\nTAD types by chromosome:\n")
tadtype_chr_dist <- df_AllDomainsSN %>%
  count(seqnames, TADType) %>%
  # NOTE: library(tidyr) must be loaded for pivot_wider to work.
  pivot_wider(names_from = TADType, values_from = n, values_fill = 0)
print(tadtype_chr_dist)

# Specifically check PcG distribution
cat("\nPcG TADs by chromosome:\n")
pcg_chr_dist <- df_AllDomainsSN %>%
  filter(TADType == "PcG") %>%
  count(seqnames)
print(pcg_chr_dist)

# --- Redundancy Check (Should now show Redundancy_Factor = 1 for all) ---
cat("\n=== REDUNDANCY CHECK ===\n")

# 1. Check the total number of rows per chromosome BEFORE filtering
before_distinct <- table(df_AllDomainsSN$seqnames)

# 2. Check the number of unique Domain IDs per chromosome
unique_ids_per_chr <- df_AllDomainsSN %>%
  group_by(seqnames) %>%
  summarise(Unique_Domains = n_distinct(DomainID)) %>%
  as.data.frame()

# 3. Calculate the redundancy factor
redundancy_check <- merge(unique_ids_per_chr, as.data.frame(before_distinct),
                          by.x = "seqnames", by.y = "Var1") %>%
  rename(Total_Rows = Freq) %>%
  mutate(Redundancy_Factor = Total_Rows / Unique_Domains) %>%
  arrange(desc(Redundancy_Factor))

print(redundancy_check)
cat("==============================================================\n")

# ==================================================
# 4. ASSIGN CANONICAL TAD IDs & ADD NEW FEATURES (FINAL CLEANING) - FIXED
# ==================================================

df_AllDomainsSN <- df_AllDomainsSN %>%
  # Ensure only unique DomainIDs (now based on coords) are kept
  distinct(DomainID, .keep_all = TRUE) %>%
  mutate(
    # Create the TAD index (0 to N-1) for graph construction
    TAD_index = 0:(n()-1), 
    TAD_index_char = as.character(TAD_index)
  ) %>%
  mutate(
    # Recalculate size and density based on final unique domains
    TAD_size = end - start + 1,
    gene_density = num_genes / TAD_size
  ) %>%
  # Ensure TADType is a factor with all types preserved
  mutate(TADType = factor(TADType, levels = unique(TADType)))

# Verify we still have all TAD types
cat("\n=== Final TAD Type Verification ===\n")
cat("TAD types after processing:\n")
print(table(df_AllDomainsSN$TADType))
cat("Total domains:", nrow(df_AllDomainsSN), "\n")

# --- CRITICAL POINT FOR DEBUGGING: Map Creation - FIXED ---
# The domain_index_map is built using the original full 'Type' name (e.g., "Null545")
# as the key, and the new canonical 'TAD_index_char' as the value.
domain_index_map <- setNames(df_AllDomainsSN$TAD_index_char, df_AllDomainsSN$Type)

cat("\n--- Map Diagnostics (COMPREHENSIVE) ---\n")
cat(paste0("Number of unique domains used to create map: ", length(domain_index_map), "\n"))


# --- CRITICAL POINT FOR DEBUGGING: Map Creation - FIXED ---

# Create the domain index map using the original Type column
domain_index_map <- setNames(df_AllDomainsSN$TAD_index_char, df_AllDomainsSN$Type)

cat("\n--- Map Diagnostics (COMPREHENSIVE) ---\n")
cat(paste0("Number of unique domains used to create map: ", length(domain_index_map), "\n"))

# Check distribution of all TAD types in the map
cat("\nTAD types in domain index map:\n")
map_tadtypes <- sapply(strsplit(names(domain_index_map), "[0-9]"), function(x) x[1])
print(table(map_tadtypes))

# Specifically check PcG domains in the map
pcg_domains_in_map <- names(domain_index_map)[grepl("^PcG", names(domain_index_map))]
cat("\nPcG domains in map:", length(pcg_domains_in_map), "\n")

# Check PcG domain distribution by chromosome in the map
if (length(pcg_domains_in_map) > 0) {
  pcg_indices <- domain_index_map[pcg_domains_in_map]
  pcg_chromosomes <- df_AllDomainsSN$seqnames[match(pcg_indices, df_AllDomainsSN$TAD_index_char)]
  cat("PcG domains in map by chromosome:\n")
  print(table(pcg_chromosomes))
  
  # Show some examples
  cat("Example PcG domains in map: ", paste(head(pcg_domains_in_map), collapse = ", "), "\n")
} else {
  cat("WARNING: No PcG domains found in the map!\n")
}

# Check other TAD types too
for (tadtype in c("Active", "Null", "Het")) {
  domains_in_map <- names(domain_index_map)[grepl(paste0("^", tadtype), names(domain_index_map))]
  cat(paste0(tadtype, " domains in map: ", length(domains_in_map), "\n"))
}

cat("Example Keys in Map: ", paste(head(names(domain_index_map)), collapse = ", "), "\n")
cat("-------------------------------\n")

# ==================================================
# 5. BUILD EDGE LISTS FOR CONTROL & SUMO (ChIA-PET) - WITH ENHANCED DIAGNOSTICS
# ==================================================

# --- COMPREHENSIVE DIAGNOSTICS FOR EDGE INDEXING ---
cat("\n=== COMPREHENSIVE EDGE INDEXING DIAGNOSTICS ===\n")

# First, let's check what PcG domains exist in ChiaSig data
cat("\n1. PcG domains in ChiaSig data:\n")
chiasig_pcg_domains_control <- unique(c(
  ChiaSig_Control$interval1[grepl("^PcG", ChiaSig_Control$interval1)],
  ChiaSig_Control$interval2[grepl("^PcG", ChiaSig_Control$interval2)]
))

chiasig_pcg_domains_sumo <- unique(c(
  ChiaSig_SUMO$interval1[grepl("^PcG", ChiaSig_SUMO$interval1)],
  ChiaSig_SUMO$interval2[grepl("^PcG", ChiaSig_SUMO$interval2)]
))

chiasig_pcg_domains_all <- unique(c(chiasig_pcg_domains_control, chiasig_pcg_domains_sumo))

cat("Unique PcG domains in ChiaSig Control:", length(chiasig_pcg_domains_control), "\n")
cat("Unique PcG domains in ChiaSig SUMO:", length(chiasig_pcg_domains_sumo), "\n")
cat("Total unique PcG domains in ChiaSig:", length(chiasig_pcg_domains_all), "\n")

# Check which of these exist in our domain_index_map
existing_chiasig_pcg <- chiasig_pcg_domains_all[chiasig_pcg_domains_all %in% names(domain_index_map)]
missing_chiasig_pcg <- chiasig_pcg_domains_all[!chiasig_pcg_domains_all %in% names(domain_index_map)]

cat("PcG domains from ChiaSig that exist in our map:", length(existing_chiasig_pcg), "\n")
cat("PcG domains from ChiaSig that are MISSING from our map:", length(missing_chiasig_pcg), "\n")

if (length(missing_chiasig_pcg) > 0) {
  cat("Sample missing PcG domains: ", paste(head(missing_chiasig_pcg), collapse = ", "), "\n")
}

# Check chromosome distribution of existing PcG domains
if (length(existing_chiasig_pcg) > 0) {
  existing_indices <- domain_index_map[existing_chiasig_pcg]
  existing_chromosomes <- df_AllDomainsSN$seqnames[match(existing_indices, df_AllDomainsSN$TAD_index_char)]
  cat("Chromosome distribution of PcG domains that exist in both datasets:\n")
  print(table(existing_chromosomes))
}

# --- DEBUGGING CHIASIG CONTROL INDEXING ---
cat("\n--- Control Edge List Indexing Diagnostic ---\n")
cat(paste0("Input ChiaSig_Control rows: ", nrow(ChiaSig_Control), "\n"))
quantile(ChiaSig_Control$av_score)
# Count PcG-PcG edges in the original ChiaSig data BEFORE mapping
pcg_pcg_edges_original <- ChiaSig_Control %>%
  filter(grepl("^PcG", interval1) & grepl("^PcG", interval2))

cat("PcG-PcG edges in original ChiaSig_Control:", nrow(pcg_pcg_edges_original), "\n")

if (nrow(pcg_pcg_edges_original) > 0) {
  cat("Chromosome distribution of original PcG-PcG edges:\n")
  # Extract chromosome from domain names (assuming format like "PcG133")
  from_chr <- sapply(strsplit(pcg_pcg_edges_original$interval1, "G"), function(x) ifelse(length(x) > 1, substr(x[2], 1, 1), "Unknown"))
  to_chr <- sapply(strsplit(pcg_pcg_edges_original$interval2, "G"), function(x) ifelse(length(x) > 1, substr(x[2], 1, 1), "Unknown"))
  cat("From chromosomes: ", table(from_chr), "\n")
  cat("To chromosomes: ", table(to_chr), "\n")
}

# The core indexing block
edge_list_indexed_Control <- ChiaSig_Control %>%
  select(Domain1 = interval1, Domain2 = interval2, ContactScore = av_score) %>%
  mutate(
    from = domain_index_map[Domain1],
    to   = domain_index_map[Domain2]
  )

# Enhanced diagnostics for mapping failures
num_na_from <- sum(is.na(edge_list_indexed_Control$from))
num_na_to <- sum(is.na(edge_list_indexed_Control$to))

cat(paste0("Rows where 'from' node is NA (map failed): ", num_na_from, " / ", nrow(ChiaSig_Control), "\n"))
cat(paste0("Rows where 'to' node is NA (map failed): ", num_na_to, " / ", nrow(ChiaSig_Control), "\n"))

# Check what types of domains are failing to map
if (num_na_from > 0) {
  failed_domains_from <- unique(edge_list_indexed_Control$Domain1[is.na(edge_list_indexed_Control$from)])
  cat("Sample 'from' domains that failed to map: ", paste(head(failed_domains_from), collapse = ", "), "\n")
  cat("Types of failed 'from' domains: ", table(grepl("^PcG", failed_domains_from)), "\n")
}

if (num_na_to > 0) {
  failed_domains_to <- unique(edge_list_indexed_Control$Domain2[is.na(edge_list_indexed_Control$to)])
  cat("Sample 'to' domains that failed to map: ", paste(head(failed_domains_to), collapse = ", "), "\n")
  cat("Types of failed 'to' domains: ", table(grepl("^PcG", failed_domains_to)), "\n")
}

# The actual filter that drops the rows
edge_list_indexed_Control <- edge_list_indexed_Control %>%
  filter(!is.na(from) & !is.na(to), from != to) %>%
  distinct()

cat(paste0("Final indexed Control edges: ", nrow(edge_list_indexed_Control), "\n"))

# Check PcG-PcG edges after mapping
if (nrow(edge_list_indexed_Control) > 0) {
  from_types <- df_AllDomainsSN$TADType[match(edge_list_indexed_Control$from, df_AllDomainsSN$TAD_index_char)]
  to_types <- df_AllDomainsSN$TADType[match(edge_list_indexed_Control$to, df_AllDomainsSN$TAD_index_char)]
  
  pcg_pcg_edges_mapped <- sum(from_types == "PcG" & to_types == "PcG")
  cat("PcG-PcG edges after mapping: ", pcg_pcg_edges_mapped, "/", nrow(edge_list_indexed_Control), "\n")
  
  if (pcg_pcg_edges_mapped > 0) {
    # Get chromosome distribution of mapped PcG-PcG edges
    pcg_edges <- edge_list_indexed_Control[from_types == "PcG" & to_types == "PcG", ]
    from_chr <- df_AllDomainsSN$seqnames[match(pcg_edges$from, df_AllDomainsSN$TAD_index_char)]
    to_chr <- df_AllDomainsSN$seqnames[match(pcg_edges$to, df_AllDomainsSN$TAD_index_char)]
    
    cat("Chromosome distribution of mapped PcG-PcG edges:\n")
    print(table(c(from_chr, to_chr)))
  }
}

cat("--------------------------------------------\n")

# --- EDGE LIST FOR SUMO ---
edge_list_indexed_SUMO <- ChiaSig_SUMO %>%
  select(Domain1 = interval1, Domain2 = interval2, ContactScore = av_score) %>%
  mutate(
    from = domain_index_map[Domain1],
    to   = domain_index_map[Domain2]
  ) %>%
  filter(!is.na(from) & !is.na(to), from != to) %>%
  distinct()

# ==================================================
# 6. MERGE EDGES ACROSS CONDITIONS & CLASSIFY CHANGES - WITH DIAGNOSTICS
# ==================================================
all_edges <- full_join(
  edge_list_indexed_Control %>% rename(ContactScore_C = ContactScore),
  edge_list_indexed_SUMO      %>% rename(ContactScore_S = ContactScore),
  by = c("from", "to")
) %>%
  mutate(ContactScore_C = ifelse(is.na(ContactScore_C), 0, ContactScore_C),
         ContactScore_S = ifelse(is.na(ContactScore_S), 0, ContactScore_S),
         Change = ContactScore_S - ContactScore_C,
         Label = case_when(
           ContactScore_C == 0 & ContactScore_S > 10 ~ "Gained",
           ContactScore_C > 10 & ContactScore_S == 0 ~ "Lost",
           ContactScore_C > 0 & ContactScore_S > 0 & Change > 10 ~ "Increased",
           ContactScore_C > 0 & ContactScore_S > 0 & Change < -10 ~ "Decreased",
           TRUE ~ "Stable"
         )) %>%
  filter(Label != "Stable")

# Final check for PcG-PcG contacts in differential edges
cat("\n=== FINAL PcG-PcG CONTACT CHECK ===\n")
if (nrow(all_edges) > 0) {
  # Get TAD types for both ends of each edge
  from_types <- df_AllDomainsSN$TADType[match(all_edges$from, df_AllDomainsSN$TAD_index_char)]
  to_types <- df_AllDomainsSN$TADType[match(all_edges$to, df_AllDomainsSN$TAD_index_char)]
  
  # Count PcG-PcG edges
  pcg_pcg_edges <- sum(from_types == "PcG" & to_types == "PcG")
  cat("PcG-PcG edges in final differential data:", pcg_pcg_edges, "/", nrow(all_edges), "\n")
  
  if (pcg_pcg_edges > 0) {
    # Get chromosome distribution of PcG-PcG edges
    pcg_edges <- all_edges[from_types == "PcG" & to_types == "PcG", ]
    from_chr <- df_AllDomainsSN$seqnames[match(pcg_edges$from, df_AllDomainsSN$TAD_index_char)]
    to_chr <- df_AllDomainsSN$seqnames[match(pcg_edges$to, df_AllDomainsSN$TAD_index_char)]
    
    cat("PcG-PcG edges by chromosome in final data:\n")
    print(table(c(from_chr, to_chr)))
  } else {
    cat("WARNING: No PcG-PcG edges found in final differential data!\n")
    cat("This suggests the mapping between ChiaSig domain names and AllDomainsSN Type names is failing.\n")
  }
}

# ==================================================
# 7. FEATURE ENGINEERING: JOIN TAD FEATURES - INTER-TAD ONLY (SAME CHROMOSOME)
# ==================================================
tad_feats <- df_AllDomainsSN %>%
  select(TAD_index_char, start, end, seqnames, TADType, 
         avgH3K27me3Control, avgH3K27me3SUMO,
         avgH3K27acControl, avgH3K27acSUMO,
         avgH2Aub118Control, avgH2Aub118SUMO,
         avgPcControl, avgPcSUMO,
         TAD_size, gene_density, avg_gene_expression)

ml_edges <- all_edges %>%
  # Join TAD A features (using 'from')
  left_join(tad_feats, by = c("from" = "TAD_index_char")) %>%
  rename_with(~ paste0("a_", .x), c(
    "start", "end", "seqnames", "TADType", 
    starts_with("avg", vars = names(.)), 
    starts_with("TAD_", vars = names(.)), 
    starts_with("gene_", vars = names(.))
  )) %>%
  # Join TAD B features (using 'to')
  left_join(tad_feats, by = c("to" = "TAD_index_char")) %>%
  rename_with(~ paste0("b_", .x), c(
    "start", "end", "seqnames", "TADType", 
    starts_with("avg", vars = names(.)), 
    starts_with("TAD_", vars = names(.)), 
    starts_with("gene_", vars = names(.))
  )) %>%
  # NEW FEATURE: Genomic distance between TADs
  mutate(distance = abs(a_end - b_start)) %>%
  # CRITICAL FILTER: Keep only INTER-TAD contacts on SAME CHROMOSOME for plotting
  # (we keep inter-chromosomal in the data but don't plot them now)
  filter(a_seqnames == b_seqnames) %>%
  # Create canonical TAD Pair Type (Symmetry enforced by sort())
  rowwise() %>%
  mutate(TAD_pair_type = paste(
    sort(c(as.character(a_TADType), as.character(b_TADType))), 
    collapse = "-")) %>%
  ungroup()

# ==================================================
# R Setup (Ensure these libraries are loaded at the top of your full script)
# ==================================================
# library(tidyverse) # Includes dplyr and tidyr
# library(igraph)

# ==================================================
# R Setup and Global Variables
# ==================================================
library(tidyverse) # Includes dplyr and tidyr
library(igraph)

# Define the color palette (CRITICAL FIX for 'invalid color name' error)
tad_color_palette <- c(
  "PcG" = "royalblue",      # Purple for PcG
  "Active" = "#FF0000",   # Red for Active
  "Null" = "#A9A9A9",     # Dark Gray for Null/Inactive
  "Het" = "#008000"       # Green for Heterochromatin (Het)
)

# Note: The data frames (ml_edges, ChiaSig_Control, ChiaSig_SUMO,
# df_AllDomainsSN, domain_index_map) are assumed to be correctly defined
# and populated by your earlier script sections (1 through 9).

# ==================================================
# 8. ENHANCED ML FEATURE TABLE (FINAL FEATURE CALCULATION) - INTER-TAD ONLY (SAME CHROMOSOME)
# ==================================================
ml_tbl <- ml_edges %>%
  mutate(
    # Differences between TAD ends
    d_H3K27me3_C = a_avgH3K27me3Control - b_avgH3K27me3Control,
    d_H3K27me3_S = a_avgH3K27me3SUMO - b_avgH3K27me3SUMO,
    d_H3K27ac_C = a_avgH3K27acControl- b_avgH3K27acControl,
    d_H3K27ac_S= a_avgH3K27acSUMO - b_avgH3K27acSUMO,
    d_H2Aub_C = a_avgH2Aub118Control - b_avgH2Aub118Control,
    d_H2Aub_S = a_avgH2Aub118SUMO - b_avgH2Aub118SUMO,
    d_Pc_C= a_avgPcControl - b_avgPcControl,
    d_Pc_S= a_avgPcSUMO - b_avgPcSUMO,
    
    # Mean values for context
    m_H3K27me3_C = (a_avgH3K27me3Control + b_avgH3K27me3Control) / 2,
    m_H3K27me3_S = (a_avgH3K27me3SUMO + b_avgH3K27me3SUMO) / 2,
    m_H3K27ac_C = (a_avgH3K27acControl + b_avgH3K27acControl) / 2,
    m_H3K27ac_S = (a_avgH3K27acSUMO + b_avgH3K27acSUMO) / 2,
    m_H2Aub_C = (a_avgH2Aub118Control + b_avgH2Aub118Control) / 2,
    m_H2Aub_S = (a_avgH2Aub118SUMO + b_avgH2Aub118SUMO) / 2,
    m_Pc_C= (a_avgPcControl+ b_avgPcControl) / 2,
    m_Pc_S= (a_avgPcSUMO + b_avgPcSUMO) / 2,
    
    # Ratio features
    ratio_H3K27me3_a = (a_avgH3K27me3SUMO + 1) / (a_avgH3K27me3Control + 1),
    ratio_H3K27me3_b = (b_avgH3K27me3SUMO + 1) / (b_avgH3K27me3Control + 1),
    ratio_H3K27ac_a = (a_avgH3K27acSUMO + 1) / (a_avgH3K27acControl + 1),
    ratio_H3K27ac_b = (b_avgH3K27acSUMO + 1) / (b_avgH3K27acControl + 1),
    
    # Change magnitude features
    total_epigenetic_change = abs(a_avgH3K27me3SUMO - a_avgH3K27me3Control) +
      abs(b_avgH3K27me3SUMO - b_avgH3K27me3Control) +
      abs(a_avgH3K27acSUMO - a_avgH3K27acControl) +
      abs(b_avgH3K27acSUMO - b_avgH3K27acControl),
    
    # Contact strength context
    contact_strength_context = (ContactScore_C + ContactScore_S) / 2,
    
    # Epigenetic similarity between TADs
    epigenetic_similarity_C = 1 / (1 + abs(a_avgH3K27me3Control - b_avgH3K27me3Control)),
    epigenetic_similarity_S = 1 / (1 + abs(a_avgH3K27me3SUMO - b_avgH3K27me3SUMO)),
    
    # Distance/Size/Gene features
    log_distance = log10(distance + 1),
    total_TAD_size = a_TAD_size + b_TAD_size,
    total_gene_density = a_gene_density + b_gene_density,
    diff_gene_expression = a_avg_gene_expression - b_avg_gene_expression,
    avg_gene_expression = (a_avg_gene_expression + b_avg_gene_expression) / 2
  ) %>%
  select(
    # Select all calculated numeric features
    starts_with("d_"), starts_with("m_"), starts_with("ratio_"),
    total_epigenetic_change, contact_strength_context,
    epigenetic_similarity_C, epigenetic_similarity_S,
    log_distance, total_TAD_size, total_gene_density,
    diff_gene_expression, avg_gene_expression,
    
    # Select the symmetric categorical feature (TAD_pair_type) and the Label
    TAD_pair_type,
    Label
  ) %>%
  mutate(Label = factor(Label, levels = unique(Label)))




# ==================================================
# 10. PcG-PcG NETWORK ANALYSIS - CLEAN AND EFFICIENT
# ==================================================

cat("\n\n--- Starting PcG-PcG Network Analysis ---\n")

# FIXED: Create a clean edge list without duplicates
create_clean_pcg_edges <- function(chiasig_data, condition_name) {
  cat(paste("Creating clean PcG edges for", condition_name, "...\n"))
  
  # Step 1: Filter PcG-PcG contacts
  pcg_edges <- chiasig_data %>%
    dplyr::filter(grepl("PcG", interval1) & grepl("PcG", interval2)) %>%
    dplyr::select(Domain1 = interval1, Domain2 = interval2, ContactScore = av_score)
  
  cat(paste("Raw PcG-PcG edges:", nrow(pcg_edges), "\n"))
  
  # Step 2: Create canonical edge keys to identify duplicates
  pcg_edges_clean <- pcg_edges %>%
    dplyr::mutate(
      # Create canonical key (alphabetically sorted domain names)
      edge_key = ifelse(Domain1 < Domain2,
                        paste(Domain1, Domain2, sep = "|"),
                        paste(Domain2, Domain1, sep = "|")),
      # Ensure consistent ordering
      canonical_from = ifelse(Domain1 < Domain2, Domain1, Domain2),
      canonical_to = ifelse(Domain1 < Domain2, Domain2, Domain1)
    )
  
  # Step 3: Remove duplicates by keeping the highest contact score
  pcg_edges_clean <- pcg_edges_clean %>%
    dplyr::group_by(edge_key, canonical_from, canonical_to) %>%
    dplyr::summarise(
      ContactScore = max(ContactScore), # Keep strongest contact
      .groups = 'drop'
    ) %>%
    dplyr::select(Domain1 = canonical_from, Domain2 = canonical_to, ContactScore)
  
  cat(paste("After duplicate removal:", nrow(pcg_edges_clean), "\n"))
  
  # Step 4: Map to TAD indices
  pcg_edges_clean <- pcg_edges_clean %>%
    dplyr::mutate(
      from = domain_index_map[Domain1],
      to = domain_index_map[Domain2]
    ) %>%
    dplyr::filter(!is.na(from) & !is.na(to)) %>%
    dplyr::select(from, to, ContactScore)
  
  cat(paste(" After mapping to TAD indices:", nrow(pcg_edges_clean), "\n"))
  
  return(pcg_edges_clean)
}

# Rebuild clean edge lists
pcg_edges_Control_clean <- create_clean_pcg_edges(ChiaSig_Control, "Control")
pcg_edges_SUMO_clean <- create_clean_pcg_edges(ChiaSig_SUMO, "SUMO")

# Verify no duplicates in clean data
cat("\n=== VERIFYING CLEAN DATA ===\n")
control_dups <- pcg_edges_Control_clean %>%
  dplyr::count(from, to) %>%
  dplyr::filter(n > 1)
cat("Duplicates in clean Control data:", nrow(control_dups), "\n")

sumo_dups <- pcg_edges_SUMO_clean %>%
  dplyr::count(from, to) %>%
  dplyr::filter(n > 1)
cat("Duplicates in clean SUMO data:", nrow(sumo_dups), "\n")



# ==================================================
# R Setup (Ensure these libraries are loaded at the top of your full script)
# ==================================================
library(tidyverse) 
library(igraph)

# Define the color palette for TAD types
tad_color_palette <- c(
  "PcG" = "royalblue",        # Purple for PcG
  "Active" = "#FF0000",        # Red for Active
  "Null" = "#A9A9A9",          # Dark Gray for Null/Inactive
  "Het" = "#008000"           # Green for Heterochromatin (Het)
)

# Define a color palette for edge strength
strength_color_palette <- colorRampPalette(c("grey", "gray50", "red"))

# NOTE: The data frames (pcg_edges_Control_clean, pcg_edges_SUMO_clean,
# df_AllDomainsSN, etc.) must be defined and populated by prior steps.
# For this script to be runnable, these variables must exist in the environment.

#
# ==================================================
# 11. BUILD GRAPHS WITH CLEAN DATA
# ==================================================

cat("\n=== BUILDING GRAPHS WITH CLEAN DATA ===\n")

# Function to create the set of unique TAD vertices
get_pcg_vertices <- function(pcg_edges_control, pcg_edges_sumo) {
  all_from <- unique(c(pcg_edges_control$from, pcg_edges_sumo$from))
  all_to <- unique(c(pcg_edges_control$to, pcg_edges_sumo$to))
  all_vertices <- unique(c(all_from, all_to))
  
  # Get vertex information from df_AllDomainsSN (assumed external data frame)
  vertices <- df_AllDomainsSN %>%
    dplyr::filter(TAD_index_char %in% all_vertices) %>%
    # CRITICAL: Include 'seqnames' for chromosome filtering later
    dplyr::select(name = TAD_index_char, TADType, seqnames, width, start, end, TAD_size)
  
  return(vertices)
}

# --- Mock Data for Executability (REPLACE WITH YOUR REAL DATA LOADING) ---
# NOTE: In a real environment, you must load these dataframes correctly.
# This section is for structure only.
if (!exists("pcg_edges_Control_clean")) {
  df_AllDomainsSN <- tibble(
    TAD_index_char = c("TAD_A1", "TAD_A2", "TAD_B1", "TAD_B2", "TAD_C1", "TAD_C2"),
    TADType = c("PcG", "PcG", "Active", "PcG", "PcG", "Active"),
    seqnames = c("chr2R", "chr2R", "chr3R", "chr3R", "chr4", "chr4")
  )
  pcg_edges_Control_clean <- tibble(
    from = c("TAD_A1", "TAD_A2", "TAD_B1", "TAD_B2", "TAD_C1"), 
    to = c("TAD_A2", "TAD_B2", "TAD_C1", "TAD_A1", "TAD_C2"), 
    ContactScore = c(15, 8, 20, 5, 45) # Variety of scores
  )
  pcg_edges_SUMO_clean <- tibble(
    from = c("TAD_A1", "TAD_B1", "TAD_B2", "TAD_C1", "TAD_C2"), 
    to = c("TAD_A2", "TAD_C1", "TAD_A2", "TAD_B2", "TAD_C1"), 
    ContactScore = c(5, 30, 18, 55, 30) # Variety of scores
  )
}
# --- END Mock Data ---

pcg_vertices_clean <- get_pcg_vertices(pcg_edges_Control_clean, pcg_edges_SUMO_clean)

# CRITICAL FIX: Add the color attribute to the vertex data frame before graph building
pcg_vertices_clean$color <- tad_color_palette[as.character(pcg_vertices_clean$TADType)]

# Function to build the comprehensive igraph object
build_pcg_graph_comprehensive <- function(edges, vertices, cond = "Control") {
  cat(paste0("Building graph for ", cond, "...\n"))
  
  edges_for_graph <- edges %>%
    dplyr::mutate(
      from = as.character(from),
      to = as.character(to),
      ContactScore = as.numeric(ContactScore)
    ) %>%
    dplyr::select(from, to, ContactScore)
  
  # CRITICAL: Explicitly Remove self-loops (Intra-TAD contacts)
  edges_for_graph <- edges_for_graph %>% dplyr::filter(from != to)
  
  g <- graph_from_data_frame(
    edges_for_graph,
    directed = FALSE,
    vertices = vertices
  )
  
  edge_attrs <- edge_attr_names(g)
  if ("ContactScore" %in% edge_attrs) {
    scores <- E(g)$ContactScore
    cat(paste0(" ContactScore present - Range: ", round(min(scores), 2), " to ", round(max(scores), 2), "\n"))
  }
  
  g$condition <- cond
  return(g)
}

# Build graphs
g_Control_PcG <- build_pcg_graph_comprehensive(pcg_edges_Control_clean, pcg_vertices_clean, "Control")
g_SUMO_PcG <- build_pcg_graph_comprehensive(pcg_edges_SUMO_clean, pcg_vertices_clean, "SUMO")

# Build differential graph (Edges)
control_edges_diff_clean <- pcg_edges_Control_clean %>%
  dplyr::select(from, to, ContactScore_C = ContactScore)

sumo_edges_diff_clean <- pcg_edges_SUMO_clean %>%
  dplyr::select(from, to, ContactScore_S = ContactScore)

pcg_diff_edges_clean <- dplyr::full_join(
  control_edges_diff_clean,
  sumo_edges_diff_clean,
  by = c("from", "to")
) %>%
  dplyr::filter(from != to) %>% 
  dplyr::mutate(
    ContactScore_C = ifelse(is.na(ContactScore_C), 0, ContactScore_C),
    ContactScore_S = ifelse(is.na(ContactScore_S), 0, ContactScore_S),
    Change = ContactScore_S - ContactScore_C,
    Label = dplyr::case_when(
      ContactScore_C == 0 & ContactScore_S > 10 ~ "Gained",
      ContactScore_C > 10 & ContactScore_S == 0 ~ "Lost",
      ContactScore_C > 0 & ContactScore_S > 0 & Change > 10 ~ "Increased",
      ContactScore_C > 0 & ContactScore_S > 0 & Change < -10 ~ "Decreased",
      TRUE ~ "Stable"
    )
  ) %>%
  dplyr::filter(Label != "Stable")

cat("Clean differential PcG-PcG edges:", nrow(pcg_diff_edges_clean), "\n")

# Build differential graph (Graph object)
diff_edges_for_graph_clean <- pcg_diff_edges_clean %>%
  dplyr::select(from, to, Label, Change)

g_Diff_PcG <- graph_from_data_frame(
  diff_edges_for_graph_clean,
  directed = FALSE,
  vertices = pcg_vertices_clean
)

# Apply edge colors for differential graph
E(g_Diff_PcG)$color <- dplyr::recode(E(g_Diff_PcG)$Label,
                                     "Gained" = "red",
                                     "Lost" = "blue",
                                     "Increased" = "orange",
                                     "Decreased" = "lightblue")
g_Diff_PcG$condition <- "Diff"

# ==================================================
# 11A. QUARTILE CALCULATION AND FILTERING
# ==================================================

cat("\n\n=== CALCULATING CONTACT SCORE QUARTILES ===\n")

all_pcg_scores <- c(pcg_edges_Control_clean$ContactScore, pcg_edges_SUMO_clean$ContactScore)
score_quantiles <- quantile(all_pcg_scores, probs = c(0.25, 0.75), na.rm = TRUE)

Q1_threshold_max <- score_quantiles["25%"] 
Q4_threshold_min <- score_quantiles["75%"] 

# 2. Filter Edge Lists for Q1 and Q4 Contacts
pcg_edges_Q1_Control <- pcg_edges_Control_clean %>% dplyr::filter(ContactScore <= Q1_threshold_max)
pcg_edges_Q1_SUMO <- pcg_edges_SUMO_clean %>% dplyr::filter(ContactScore <= Q1_threshold_max)
pcg_edges_Q4_Control <- pcg_edges_Control_clean %>% dplyr::filter(ContactScore >= Q4_threshold_min)
pcg_edges_Q4_SUMO <- pcg_edges_SUMO_clean %>% dplyr::filter(ContactScore >= Q4_threshold_min)

# 3. Build Q1 and Q4 Graphs
g_Control_PcG_Q1 <- build_pcg_graph_comprehensive(pcg_edges_Q1_Control, pcg_vertices_clean, "Control_Q1")
g_SUMO_PcG_Q1 <- build_pcg_graph_comprehensive(pcg_edges_Q1_SUMO, pcg_vertices_clean, "SUMO_Q1")
g_Control_PcG_Q4 <- build_pcg_graph_comprehensive(pcg_edges_Q4_Control, pcg_vertices_clean, "Control_Q4")
g_SUMO_PcG_Q4 <- build_pcg_graph_comprehensive(pcg_edges_Q4_SUMO, pcg_vertices_clean, "SUMO_Q4")

# 4. Build Differential Graphs for Q1 and Q4
build_diff_graph_q_subset <- function(edges_C, edges_S, q_label) {
  cat(paste0("Building Differential Graph for ", q_label, "...\n"))
  
  diff_edges <- dplyr::full_join(
    edges_C %>% dplyr::select(from, to, ContactScore_C = ContactScore),
    edges_S %>% dplyr::select(from, to, ContactScore_S = ContactScore),
    by = c("from", "to")
  ) %>%
    dplyr::filter(from != to) %>% 
    dplyr::mutate(
      ContactScore_C = ifelse(is.na(ContactScore_C), 0, ContactScore_C),
      ContactScore_S = ifelse(is.na(ContactScore_S), 0, ContactScore_S),
      Change = ContactScore_S - ContactScore_C,
      Label = dplyr::case_when(
        ContactScore_C == 0 & ContactScore_S > 0 ~ "Gained", 
        ContactScore_C > 0 & ContactScore_S == 0 ~ "Lost",
        ContactScore_C > 0 & ContactScore_S > 0 & Change > 10 ~ "Increased",
        ContactScore_C > 0 & ContactScore_S > 0 & Change < -10 ~ "Decreased",
        TRUE ~ "Stable"
      )
    ) %>%
    dplyr::filter(Label != "Stable")
  
  if (nrow(diff_edges) == 0) {
    g_Diff <- make_empty_graph(n = 0, directed = FALSE)
    V(g_Diff)$name <- character(0)
    g_Diff$condition <- paste0("Diff_", q_label)
    return(g_Diff)
  }
  
  g_Diff <- graph_from_data_frame(
    diff_edges %>% dplyr::select(from, to, Label, Change),
    directed = FALSE,
    vertices = pcg_vertices_clean
  )
  
  E(g_Diff)$color <- dplyr::recode(E(g_Diff)$Label,
                                   "Gained" = "red",
                                   "Lost" = "blue",
                                   "Increased" = "orange",
                                   "Decreased" = "lightblue")
  g_Diff$condition <- paste0("Diff_", q_label)
  return(g_Diff)
}

g_Diff_PcG_Q1 <- build_diff_graph_q_subset(pcg_edges_Q1_Control, pcg_edges_Q1_SUMO, "Q1")
g_Diff_PcG_Q4 <- build_diff_graph_q_subset(pcg_edges_Q4_Control, pcg_edges_Q4_SUMO, "Q4")

# ==================================================
# 12. PLOT AND SAVE GRAPHS (FINALIZED PLOTTING FUNCTION)
# ==================================================

# Define the output directory name
output_dir <- "pcg_network_plots_by_chromosome_2" 

# Create the directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
  cat(paste0("Created output directory: ", output_dir, "\n"))
} else {
  cat(paste0("Output directory already exists: ", output_dir, "\n"))
}

# FINAL PLOTTING FUNCTION
plot_and_save_graph_by_chr <- function(graph, filename, plot_title, fixed_layout = NULL, chr_name) {
  
  # Define edge properties based on graph condition
  if (graph$condition %in% c("Control", "SUMO", "Control_Q1", "SUMO_Q1", "Control_Q4", "SUMO_Q4")) {
    # --- Non-differential graphs (use gradient) ---
    max_score <- max(E(graph)$ContactScore)
    min_score <- min(E(graph)$ContactScore) 
    
    # Scale scores from 0-1 for uniform color mapping
    if ((max_score - min_score) == 0) {
      scaled_scores <- rep(0.5, ecount(graph)) 
    } else {
      scaled_scores <- (E(graph)$ContactScore - min_score) / (max_score - min_score)
    }
    
    # Map scaled scores to the white-to-dark-grey strength_color_palette (100 shades)
    edge_color_to_use <- strength_color_palette(100)[round(scaled_scores * 99) + 1] 
    
    edge_widths <- E(graph)$ContactScore / ifelse(max_score > 0, max_score, 1) * 3
    
  } else if (grepl("^Diff", graph$condition)) {
    # --- Differential graphs (use discrete colors) ---
    max_change <- max(abs(E(graph)$Change))
    edge_widths <- abs(E(graph)$Change) / ifelse(max_change > 0, max_change, 1) * 3
    edge_color_to_use <- E(graph)$color
    
  } else {
    # Fallback
    edge_widths <- 1
    edge_color_to_use <- "gray50"
  }
  
  # Start plotting to a PNG file
  png_path <- file.path(output_dir, filename)
  png(filename = png_path, width = 1200, height = 1200, res = 100)
  
  plot(
    graph,
    layout = fixed_layout,
    vertex.size = 10,
    vertex.label = V(graph)$name,
    vertex.label.cex = 0.8,
    vertex.color = V(graph)$color,
    vertex.frame.color = "gray",
    edge.width = edge_widths,
    edge.color = edge_color_to_use, 
    main = paste(plot_title, "PcG-PcG Network on", chr_name),
    sub = paste("Vertices:", vcount(graph), ", Edges:", ecount(graph)),
    margin = c(0, 0, 0, 0)
  )
  
  # Add legend(s) based on graph type
  if (grepl("^Diff", graph$condition)) {
    # Differential Legend (Change Types)
    legend(
      "bottomleft",
      legend = c("Gained", "Lost", "Increased", "Decreased"), 
      col = c("red", "blue", "orange", "lightblue"),
      lwd = 5, bty = "n", cex = 1.2, title = "Edge Change"
    )
  } else {
    # Non-Differential Legends (TAD Type and Strength Gradient)
    
    # 1. TAD Type Legend
    legend(
      "topright",
      legend = names(tad_color_palette),
      col = tad_color_palette,
      pch = 19, bty = "n", cex = 1.2, title = "TAD Type"
    )
    
    # 2. Contact Strength Legend (Gradient representation)
    legend_scores <- seq(min_score, max_score, length.out = 5)
    legend_colors <- strength_color_palette(5) # Get 5 representative colors
    
    # Place this below the TAD type legend
    legend(
      "bottomright",
      legend = round(legend_scores, 1),
      col = legend_colors,
      lwd = 5, # Make lines thick to represent color blocks
      bty = "n", cex = 1.2, title = "Contact Strength"
    )
  }
  
  dev.off()
  cat(paste0("Successfully saved: ", png_path, "\n"))
}

# --- Execution Loop for Chromosome-Specific Plots (FINAL) ---

# Define the list of graph versions to loop through
graph_versions <- list(
  list(ctrl = g_Control_PcG, sumo = g_SUMO_PcG, diff = g_Diff_PcG, suffix = "All"),
  list(ctrl = g_Control_PcG_Q1, sumo = g_SUMO_PcG_Q1, diff = g_Diff_PcG_Q1, suffix = "Q1"),
  list(ctrl = g_Control_PcG_Q4, sumo = g_SUMO_PcG_Q4, diff = g_Diff_PcG_Q4, suffix = "Q4")
)

# 1. Identify all chromosomes present in the universal network
all_chromosomes <- unique(V(g_Control_PcG)$seqnames)

cat("\nStarting chromosome-specific plotting for ALL, Q1, and Q4 network versions.\n")

for (chr in all_chromosomes) {
  
  chr_vertices <- V(g_Control_PcG)[V(g_Control_PcG)$seqnames == chr]$name
  
  if (length(chr_vertices) < 2) {
    cat(paste0("  Skipping ", chr, ": Only ", length(chr_vertices), " vertex(es).\n"))
    next
  }
  
  cat(paste0("\nProcessing ", chr, " (", length(chr_vertices), " vertices)...\n"))
  
  # Generate the layout once based on the ALL-Contacts Control graph for consistency
  g_Control_All_chr <- induced_subgraph(g_Control_PcG, vids = chr_vertices)
  
  if (vcount(g_Control_All_chr) > 1) {
    chr_layout <- layout_nicely(g_Control_All_chr) 
  } else {
    chr_layout <- layout_in_circle(make_full_graph(length(chr_vertices)))
  }
  
  
  for (version in graph_versions) {
    
    ctrl_g <- version$ctrl
    sumo_g <- version$sumo
    diff_g <- version$diff
    suffix <- version$suffix
    
    # --- A. Control Network ---
    g_Control_chr <- induced_subgraph(ctrl_g, vids = chr_vertices)
    if (vcount(g_Control_chr) > 1 && ecount(g_Control_chr) > 0) {
      plot_and_save_graph_by_chr(
        g_Control_chr,
        filename = paste0("Control_", chr, "_PcG_Network_", suffix, ".png"),
        plot_title = paste("Control Condition", suffix, "Contacts"),
        fixed_layout = chr_layout, 
        chr_name = chr
      )
    } else {
      cat(paste0("  Skipping Control ", suffix, " plot for ", chr, ": No edges found.\n"))
    }
    
    # --- B. SUMO Network ---
    g_SUMO_chr <- induced_subgraph(sumo_g, vids = chr_vertices)
    if (vcount(g_SUMO_chr) > 1 && ecount(g_SUMO_chr) > 0) {
      plot_and_save_graph_by_chr(
        g_SUMO_chr,
        filename = paste0("SUMO_", chr, "_PcG_Network_", suffix, ".png"),
        plot_title = paste("SUMO Condition", suffix, "Contacts"),
        fixed_layout = chr_layout, 
        chr_name = chr
      )
    } else {
      cat(paste0("  Skipping SUMO ", suffix, " plot for ", chr, ": No edges found.\n"))
    }
    
    # --- C. Differential Network ---
    g_Diff_chr <- induced_subgraph(diff_g, vids = chr_vertices)
    if (vcount(g_Diff_chr) > 1 && ecount(g_Diff_chr) > 0) {
      plot_and_save_graph_by_chr(
        g_Diff_chr,
        filename = paste0("Differential_", chr, "_PcG_Network_", suffix, ".png"),
        plot_title = paste("SUMO vs Control Differential", suffix, "Contacts"),
        fixed_layout = chr_layout, 
        chr_name = chr
      )
    } else {
      cat(paste0("  Skipping Differential ", suffix, " plot for ", chr, ": No significant change edges found.\n"))
    }
  }
}

cat("\nChromosome-specific plotting complete for All, Q1, and Q4 network versions.\n")
