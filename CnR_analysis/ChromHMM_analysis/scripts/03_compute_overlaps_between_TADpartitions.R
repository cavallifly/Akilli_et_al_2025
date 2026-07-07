
# -----------------------------
# 0. LIBRERIE
# -----------------------------
library(data.table)
library(pheatmap)
library(reshape2)
library(grid)

# -----------------------------
# 1. FILES - gene lists
# -----------------------------
base_files <- list(
  PcG1    = "cluster_PcG1_TADs.txt",
  Active1 = "cluster_Active1_TADs.txt",
  Null1   = "cluster_Null1_TADs.txt",
  Het1   = "cluster_Het1_TADs.txt"
)
#base_files_order = c("PcG1","Active1","Null1","Het1")
base_files_order = c(1,2,3,4)

rec_files <- list(
  PcG2    = "cluster_PcG2_TADs.txt",
  Active2 = "cluster_Active2_TADs.txt",
  Null2   = "cluster_Null2_TADs.txt",
  Het2   = "cluster_Het2_TADs.txt"  
  )
rec_files_order = c(1,2,3,4)  

all_files <- c(base_files, rec_files)
all_files <- unname(unlist(all_files))
print(all_files)

colorScaleMax <- 1.5
colorScaleMin <- -colorScaleMax

# -----------------------------
# 2. CLEAN FUNCTION (FIX DEFINITIVO)
# -----------------------------
clean_genes <- function(x){
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\\s+", "", x)
  x <- toupper(x)
  x <- gsub("[^A-Z0-9]", "", x)
  x <- x[!is.na(x) & x != ""]
  unique(x)
}

# -----------------------------
# 3. LOAD + CLEAN
# -----------------------------
gene_sets <- lapply(all_files, function(f) {
  clean_genes(scan(f, what = character(), quiet = TRUE))
})
names(gene_sets) <- c(names(base_files),names(rec_files))
print(names(gene_sets))

# -----------------------------
# 4. GLOBAL GENE SET
# -----------------------------
all_genes <- unique(unlist(gene_sets))

# -----------------------------
# 5. BINARY MATRIX
# -----------------------------
genes_dt <- data.table(symbol = all_genes)

for (set_name in names(gene_sets)) {
  genes_dt[[set_name]] <- genes_dt$symbol %in% gene_sets[[set_name]]
}

# -----------------------------
# 6. FISHER TEST
# -----------------------------
results <- data.table()

for (g1 in names(rec_files)[rec_files_order]) {
  for (g2 in names(base_files)[base_files_order]) {
    
    test <- fisher.test(
      genes_dt[[g1]],
      genes_dt[[g2]],
      alternative = "greater"
    )
    
    overlap <- sum(genes_dt[[g1]] & genes_dt[[g2]])
    
    n_set1 <- sum(genes_dt[[g1]])
    n_set2 <- sum(genes_dt[[g2]])
    
    results <- rbind(results, data.table(
      set1 = g1,
      set2 = g2,
      n_set1 = n_set1,
      n_set2 = n_set2,
      overlap = overlap,
      odds_ratio = unname(test$estimate),
      pvalue = test$p.value
    ))
  }
}
print(results)

# -----------------------------
# 7. MULTIPLE TEST CORRECTION
# -----------------------------
results[, padj := p.adjust(pvalue, method = "BH")]
results[, odds_ratio := as.numeric(odds_ratio)]
results[, log2OR := log2(odds_ratio)]
print(results)

# fix Inf
finite_vals <- results[is.finite(log2OR), log2OR]

if (length(finite_vals) > 0) {
  min_finite <- min(finite_vals)
  results[is.infinite(log2OR) & log2OR < 0, log2OR := min_finite]
}

# -----------------------------
# 8. HEATMAP MATRICES
# -----------------------------
results$set1 <- factor(results$set1, levels = names(rec_files))
results$set2 <- factor(results$set2, levels = names(base_files))

heat_mat <- dcast(results, set1 ~ set2, value.var = "log2OR")
rownames(heat_mat) <- heat_mat[[1]]
heat_mat <- as.matrix(heat_mat[, -1])
print(heat_mat)

padj_mat <- dcast(results, set1 ~ set2, value.var = "padj")
rownames(padj_mat) <- padj_mat[[1]]
padj_mat <- as.matrix(padj_mat[, -1])

ov_mat <- dcast(results, set1 ~ set2, value.var = "overlap")
rownames(ov_mat) <- ov_mat[[1]]
ov_mat <- as.matrix(ov_mat[, -1])

# -----------------------------
# 9. CLEAN MATRICES
# -----------------------------

heat_mat[heat_mat == Inf]  <- colorScaleMax
heat_mat[heat_mat == -Inf] <- colorScaleMin
#heat_mat[!is.finite(heat_mat)] <- 0
padj_mat[is.na(padj_mat)] <- 1
ov_mat[is.na(ov_mat)] <- 0


# -----------------------------
# 10. CLUSTER SIZE LABELS (FIX RICHIESTO)
# -----------------------------
cluster_sizes <- sapply(gene_sets, length)

rownames(heat_mat) <- paste0(
  rownames(heat_mat),
  " (n=", cluster_sizes[rownames(heat_mat)], ")"
)

colnames(heat_mat) <- paste0(
  colnames(heat_mat),
  " (n=", cluster_sizes[colnames(heat_mat)], ")"
)

rownames(padj_mat) <- rownames(heat_mat)
colnames(padj_mat) <- colnames(heat_mat)

rownames(ov_mat) <- rownames(heat_mat)
colnames(ov_mat) <- colnames(heat_mat)

# -----------------------------
# 11. STARS MATRIX
# -----------------------------
stars <- padj_mat
stars[padj_mat <= 0.05] <- "*"
stars[padj_mat < 0.01] <- "**"
stars[padj_mat < 0.001] <- "***"
stars[padj_mat < 0.0001] <- "****"
stars[padj_mat > 0.05] <- "N.S."

# -----------------------------
# 12. HEATMAP OUTPUT
# -----------------------------
pdf("heatmap_withNAH3K9me3vswithGSH3K9me3.pdf", height = 4, width = 5)

par(mar = c(5, 5, 4, 2))

p <- pheatmap::pheatmap(
  heat_mat,
  breaks = seq(colorScaleMin, colorScaleMax, length.out = 21),
  col = colorRampPalette(c("cornflowerblue", "white", "tomato"))(21),
  display_numbers = matrix(
    paste0(stars, "\n", ov_mat),
    nrow = nrow(stars)
  ),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main = "Log2 Odd-Ratio of TAD states overlaps"
)

#grid.newpage()
#grid.draw(p$gtable)

# Y-axis label
#grid.text(
#  "Rows (Set1)",
#  x = 0.02, y = 0.5,
#  rot = 90,
#  gp = gpar(fontsize = 12)
#)

# X-axis label
#grid.text(
#  "Columns (Set2)",
#  x = 0.5, y = 0.02,
#  gp = gpar(fontsize = 12)
#)



dev.off()
