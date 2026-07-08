# Load required libraries
library(ggplot2)
library(ggpubr)

# Read your data file (replace with your file path)
# Expected format: 3 columns (e.g., "region", "condition1", "condition2")
#inFile <- "TopDom_domains_hic_WD_NOnub_merge_dm6_NA_window_5_at_5000bp_filtered_refined_at_1000bp_with_WD_NOnub_and_SUMOnub_merge_raw_insulation_data.tsv"
inFile <- "TopDom_domains_hic_WD_NOnub_merge_dm6_NA_window_5_at_5000bp_filtered_refined_at_1000bp_with_WD_NOnub_and_SUMOnub_merge_raw_insulation_scatterData.tsv"
#inFile <- "TopDom_domains_hic_WD_SUMOnub_merge_dm6_NA_window_5_at_5000bp_filtered_refined_at_1000bp_with_WD_NOnub_and_SUMOnub_merge_insulation_data.tsv"
data <- read.table(inFile, header = T, sep = "\t")
print(head(data))

name <- gsub(".tsv","",inFile)
outFile <- paste0("scatterPlot_",name,"_with_stats.pdf")

# Extract the values for the two conditions
dataScatter <- data.frame(Control = data$insNOnub, SUMORNAi = data$insSUMOnub)

# Clean the data: Remove NA/Inf
dataScatter_clean <- na.omit(dataScatter)
dataScatter_clean <- dataScatter_clean[
  is.finite(dataScatter_clean$Control) &
  is.finite(dataScatter_clean$SUMORNAi), ]
print(head(dataScatter))
#quit()

# Calculate correlation coefficients
pearson_r    <- cor(as.numeric(dataScatter_clean$Control), as.numeric(dataScatter_clean$SUMORNAi), method = "pearson")
spearman_rho <- cor(as.numeric(dataScatter_clean$Control), as.numeric(dataScatter_clean$SUMORNAi), method = "spearman")
r_squared    <- pearson_r^2

# Print correlation coefficients to the console
#cat("Pearson's r:", round(pearson_r, 3), "\n")
#cat("R²:", round(r_squared, 3), "\n")
#cat("Spearman's ρ:", round(spearman_rho, 3), "\n")
cat("Pearson's r:", pearson_r, "\n")
cat("R²:", r_squared, "\n")
cat("Spearman's ρ:", spearman_rho, "\n")
#quit()
data <- data.frame(dataScatter$Control, dataScatter$SUMORNAi)
colnames(data) <- c("condition1","condition2")
print(head(data))

# Create the scatter plot with regression line
p <- ggplot(data, aes(x = condition1, y = condition2)) +
  geom_point(alpha = 0.6, color = "darkgreen", size=1) +  # Scatter points
  #geom_smooth(method = "lm", se = FALSE, color = "red") +  # Linear regression line
  labs(
    title = "Scatter Plot with Correlation Coefficients",
    x = "TADs Insulation Control",
    y = "TADs Insulation SUMO RNAi"
  ) +
  theme_classic()

# Add correlation annotations to the plot
p <- p +
  annotate(
    "text",
    x = mean(data$condition1), y = mean(data$condition2),
    label = paste0(
      "Pearson r = ", round(pearson_r, 3), "\n",
      "R² = ", round(r_squared, 3), "\n",
      "Spearman ", expression(rho), " = ", round(spearman_rho, 3)
    ),
    #hjust = 1.1, vjust = 1.1,
    size = 8,        
    color = "black"
  ) +
  theme(legend.position = "none",
        plot.title = element_text(size = 18, color = "black", face= "bold"),
        plot.subtitle = element_text(size = 16, face = "bold.italic"),
	axis.text  = element_text(size = 16),
	axis.title = element_text(size = 16, face = "bold"),	
	axis.ticks = element_line(color = "black", size = 0.5)
        )
	     

#quit()

# Print the plot
pdf(outFile)
print(p)
dev.off()

quit()

