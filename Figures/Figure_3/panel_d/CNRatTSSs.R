library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggpubr)

H2Aub118<-read.delim2("./report_H2Aub118_WD_SUMOnub_vs_WD_NOnub_withGeneNames.tsv",header= TRUE)
nrow(H2Aub118)
H2Aub118<-H2Aub118[-1,]

H2Aub118<- H2Aub118 %>%
  rownames_to_column("id") %>%                                 # bring rownames into a column
  extract(
    id,
    into = c("gene", "FBgn", "status"),
    regex = "^(.*)_(FBgn\\d+)_(.*)$",                          # gene_FBgn#####_status
    remove = TRUE
  )

head(H2Aub118)
H2Aub118_clean <- na.omit(H2Aub118)

H2Aub118_clean <- H2Aub118 %>%
  filter(complete.cases(.))


# Reshape to long format
H2Aub118_long <- H2Aub118_clean %>%
  select(gene, status, Conc_WD_NOnub, Conc_WD_SUMOnub) %>%
  pivot_longer(cols = starts_with("Conc"),
               names_to = "Condition",
               values_to = "Concentration")

# Make condition labels nicer
H2Aub118_long$Condition <- recode(H2Aub118_long$Condition,
                                  Conc_WD_NOnub = "Control",
                                  Conc_WD_SUMOnub = "SUMO")


# Make a combined condition variable
H2Aub118_long <- H2Aub118_long %>%
  mutate(Group = case_when(
    status == "notChanged" & Condition == "Control" ~ "NS_C",
    status == "notChanged" & Condition == "SUMO"    ~ "NS_S",
    status == "Up"         & Condition == "Control" ~ "Up_C",
    status == "Up"         & Condition == "SUMO"    ~ "Up_S",
    status == "Down"       & Condition == "Control" ~ "Down_C",
    status == "Down"       & Condition == "SUMO"    ~ "Down_S"
  ))
print(H2Aub118_long)

# Make sure Group is an ordered factor (for plotting order)
H2Aub118_long$Group <- factor(H2Aub118_long$Group,
                              levels = c("NS_C","NS_S",
                                         "Up_C","Up_S",
                                         "Down_C","Down_S"))
print(H2Aub118_long)
print(table(H2Aub118_long$Group))


# Define comparisons
comparisons <- list(
  c("NS_C", "NS_S"),
  c("Up_C", "Up_S"),
  c("Down_C", "Down_S")
)

# Convert Concentration to numeric
H2Aub118_long$Concentration <- as.numeric(H2Aub118_long$Concentration)






comparisons <- list(c("NS_C", "NS_S"), c("Up_C", "Up_S"), c("Down_C", "Down_S"))
print(head(H2Aub118_long))

#maxmedian <- H2Aub118_long %>%
#  group_by(Group) %>%
#  summarise(med = median(Concentration, na.rm = TRUE)) %>%
#  pull(med) %>%
#  max()

#Select PcG targets:PcGtargets2

PcG_H2AubNS <- read.table("PcG_targets_Parreno2024_NS.txt")
PcG_H2AubNS$status <- "NS"
colnames(PcG_H2AubNS) <- c("gene_name","status")
PcG_H2AubUp <- read.table("PcG_targets_Parreno2024_Up.txt")
PcG_H2AubUp$status <- "Up"
colnames(PcG_H2AubUp) <- c("gene_name","status")
PcG_H2AubDown <- read.table("PcG_targets_Parreno2024_Down.txt")
PcG_H2AubDown$status <- "Down"
colnames(PcG_H2AubDown) <- c("gene_name","status")
PcGtargets2 <- rbind(PcG_H2AubNS,PcG_H2AubUp,PcG_H2AubDown)
PcG_H2Aub <- H2Aub118_long[H2Aub118_long$gene %in% PcGtargets2$gene_name,]

print(paste0(length(unique(PcG_H2Aub[PcG_H2Aub$status == "Up",]$gene))))
print(paste0(length(unique(PcG_H2Aub[PcG_H2Aub$status == "Down",]$gene))))
print(paste0(length(unique(PcG_H2Aub[PcG_H2Aub$status == "NS",]$gene))))
#nrow(PcG_H2Aub[PcG_H2Aub$status == "Down",])
#nrow(PcG_H2Aub[PcG_H2Aub$status == "NS",])
print(PcG_H2Aub)

###Make a combined plot for all

H3K27me3<-read.delim2("./report_H3K27me3_WD_SUMOnub_vs_WD_NOnub_withGeneNames.tsv",header= TRUE)
nrow(H3K27me3)
H3K27me3<-H3K27me3[-1,]

H3K27ac<-read.delim2("./report_H3K27ac_WD_SUMOnub_vs_WD_NOnub_withGeneNames.tsv",header= TRUE)
nrow(H3K27ac)
H3K27ac<-H3K27ac[-1,]

Pc<-read.delim2("./report_PcXlinked_WD_SUMOnub_vs_WD_NOnub_withGeneNames.tsv",header= TRUE)
nrow(Pc)
Pc<-Pc[-1,]

H2Aub118<-read.delim2("./report_H2Aub118_WD_SUMOnub_vs_WD_NOnub_withGeneNames.tsv",header= TRUE)
nrow(H2Aub118)
H2Aub118<-H2Aub118[-1,]

process <- function(df) {
  df %>%
    # bring rownames into a column
    tibble::rownames_to_column("id") %>%
    
    # extract gene, FBgn, status from rownames
    tidyr::extract(
      id,
      into = c("gene", "FBgn", "status"),
      regex = "^(.*)_(FBgn\\d+)_(.*)$",  
      remove = TRUE
    ) %>%
    
    
    # remove rows with NA
    filter(complete.cases(.)) %>%
    
    # **filter for specific genes**
    #filter(gene %in% PcGtargets2$gene_name) %>%
    
    # reshape to long
    select(gene, status, Conc_WD_NOnub, Conc_WD_SUMOnub) %>%
    pivot_longer(
      cols = starts_with("Conc"),
      names_to = "Condition",
      values_to = "Concentration"
    ) %>%
    
    # relabel conditions
    mutate(Condition = recode(Condition,
                              Conc_WD_NOnub = "Control",
                              Conc_WD_SUMOnub = "SUMO")) %>%
    
    # build group labels
    mutate(Group = case_when(
      status == "notChanged" & Condition == "Control" ~ "NS_C",
      status == "notChanged" & Condition == "SUMO"    ~ "NS_S",
      status == "Up"         & Condition == "Control" ~ "Up_C",
      status == "Up"         & Condition == "SUMO"    ~ "Up_S",
      status == "Down"       & Condition == "Control" ~ "Down_C",
      status == "Down"       & Condition == "SUMO"    ~ "Down_S"
    )) %>%
    
    # order groups
    mutate(Group = factor(Group,
                          levels = c("NS_C","NS_S",
                                     "Up_C","Up_S",
                                     "Down_C","Down_S")),
           # force numeric
           Concentration = as.numeric(Concentration))
        
}

process_PcG <- function(df) {
  df %>%
    # bring rownames into a column
    tibble::rownames_to_column("id") %>%
    
    # extract gene, FBgn, status from rownames
    tidyr::extract(
      id,
      into = c("gene", "FBgn", "status"),
      regex = "^(.*)_(FBgn\\d+)_(.*)$",  
      remove = TRUE
    ) %>%
    
    
    # remove rows with NA
    filter(complete.cases(.)) %>%
    
    # **filter for specific genes**
    filter(gene %in% PcGtargets2$gene_name) %>%
    
    # reshape to long
    select(gene, status, Conc_WD_NOnub, Conc_WD_SUMOnub) %>%
    pivot_longer(
      cols = starts_with("Conc"),
      names_to = "Condition",
      values_to = "Concentration"
    ) %>%
    
    # relabel conditions
    mutate(Condition = recode(Condition,
                              Conc_WD_NOnub = "Control",
                              Conc_WD_SUMOnub = "SUMO")) %>%
    
    # build group labels
    mutate(Group = case_when(
      status == "notChanged" & Condition == "Control" ~ "NS_C",
      status == "notChanged" & Condition == "SUMO"    ~ "NS_S",
      status == "Up"         & Condition == "Control" ~ "Up_C",
      status == "Up"         & Condition == "SUMO"    ~ "Up_S",
      status == "Down"       & Condition == "Control" ~ "Down_C",
      status == "Down"       & Condition == "SUMO"    ~ "Down_S"
    )) %>%
    
    # order groups
    mutate(Group = factor(Group,
                          levels = c("NS_C","NS_S",
                                     "Up_C","Up_S",
                                     "Down_C","Down_S")),
           # force numeric
           Concentration = as.numeric(Concentration))
        
}

# ---- Example usage ----
# say you have several datasets: H3K27me3, Pc, H2Aub118 df4
datasets <- list(H3K27me3 = H3K27me3, Pc = Pc, H2Aub118 = H2Aub118, H3K27ac= H3K27ac)

processed <- lapply(datasets, process)
processed_PcG <- lapply(datasets, process_PcG)

# access results
processed$H3K27me3 %>% head()
processed$Pc %>% head()
processed$H2Aub118 %>% head()
processed$H3K27ac %>% head()

# access results
processed_PcG$H3K27me3 %>% head()
processed_PcG$Pc %>% head()
processed_PcG$H2Aub118 %>% head()
processed_PcG$H3K27ac %>% head()

plot_violin <- function(df, title = "All TSSs") {
  
  # define comparisons for stats
  comparisons <- list(
    c("NS_C", "NS_S"),
    c("Up_C", "Up_S"),
    c("Down_C", "Down_S")
  )
  
  # compute dashed line at max median across groups
  #maxmedian <- df %>%
  #  group_by(Group) %>%
  #  summarise(med = median(Concentration, na.rm = TRUE)) %>%
  #  pull(med) %>%
  #  max(na.rm = TRUE)

  summ <- df %>%
    group_by(Group) %>%
    summarize(n = n(), y=min(df$Concentration)-0.5)
summ <- unique(summ)
print(as.data.frame(summ), quote=F)   
summ <- summ[-grep("_S",summ$Group),]


ggplot(df, aes(x = Group, y = Concentration, fill = Group)) +
    geom_violin(trim = T, alpha = 0.6, width=0.9) +
    geom_boxplot(width = 0.05, outlier.shape = NA, show.legend = FALSE) +
    stat_compare_means(comparisons = comparisons,
                       method = "wilcox.test",
                       exact = FALSE,
                       label.x.npc = "center",
                       #label="p.signif",
		       label="p.adj",		       
                       size = 3, step.increase = 0,bracket.size = 0.4, tip.length = 0) +
    #geom_hline(yintercept = maxmedian, linetype = "dashed", color = "grey", size = 0.7) +
    geom_text(data=summ, aes(x=Group, y=y, label = paste0("n = ",n)), size=3., color="black",position = position_nudge(x = 0.5)) +
    scale_fill_manual(values = c(
      "NS_C" = "gray30", "NS_S" = "gray30",
      "Up_C" = "skyblue4","Up_S" = "skyblue4",
      "Down_C" = "firebrick","Down_S" = "firebrick"
    )) +
    labs(title = title,
         x = "", y = "Enrichment at TSSs") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.4,size = 10),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 10)
    )
}

plot_violin_PcG <- function(df, title = "All TSSs") {
  
  # define comparisons for stats
  comparisons <- list(
    c("NS_C", "NS_S"),
    c("Up_C", "Up_S"),
    c("Down_C", "Down_S")
  )
  
  # compute dashed line at max median across groups
  #maxmedian <- df %>%
  #  group_by(Group) %>%
  #  summarise(med = median(Concentration, na.rm = TRUE)) %>%
  #  pull(med) %>%
  #  max(na.rm = TRUE)

  summ <- df %>%
    group_by(Group) %>%
    summarize(n = n(), y=min(df$Concentration)-0.5)
summ <- unique(summ)
print(as.data.frame(summ), quote=F)   
summ <- summ[-grep("_S",summ$Group),]


ggplot(df, aes(x = Group, y = Concentration, fill = Group)) +
    geom_violin(trim = T, alpha = 0.6, width= 0.9) +
    geom_boxplot(width = 0.05, outlier.shape = NA, show.legend = FALSE) +
    stat_compare_means(comparisons = comparisons,
                       method = "wilcox.test",
                       exact = FALSE,
                       label.x.npc = "center",
                       #label="p.signif",
                       label="p.adj",		       
                       size = 3, step.increase = 0, bracket.size = 0.4, tip.length = 0) +
    #geom_hline(yintercept = maxmedian, linetype = "dashed", color = "grey", size = 0.7) +
    geom_text(data=summ, aes(x=Group, y=y, label = paste0("n = ",n)), size=3., color="black", position = position_nudge(x = 0.5)) +
    scale_fill_manual(values = c(
      "NS_C" = "gray30", "NS_S" = "gray30",
      "Up_C" = "skyblue4","Up_S" = "skyblue4",
      "Down_C" = "firebrick","Down_S" = "firebrick"
    )) +
    labs(title = title,
         x = "", y = "Enrichment at TSS of PcG target genes") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.4,size = 10),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 10)
    )
}

# Generate plots for each dataset
p1 <- plot_violin(processed$H3K27me3, title = "H3K27me3")
p2 <- plot_violin(processed$Pc, title = "Pc")
p3 <- plot_violin(processed$H2Aub118, title = "H2Aub118")
p4 <- plot_violin(processed$H3K27ac, title = "H3K27ac")

combined <- ggarrange(p1, p2, p3, p4,
                      ncol = 2, nrow = 2,   # 2x2 grid
                      common.legend = TRUE, # share legend if you want
                      legend = "none")


outFile <- "Enrichment_at_TSS_of_all_genes.tab"
write.table(processed,file=outFile,sep="\t",quote=F,row.names=F,)
print(processed)

#H3K27me3.gene	H3K27me3.status	H3K27me3.Condition	H3K27me3.Concentration	H3K27me3.Group	Pc.gene	Pc.status	Pc.Condition	Pc.Concentration	Pc.Group	H2Aub118.geneH2Aub118.status	H2Aub118.Condition	H2Aub118.Concentration	H2Aub118.Group	H3K27ac.gene	H3K27ac.status	H3K27ac.Condition	H3K27ac.Concentration	H3K27ac.Group


# Generate plots for each dataset
p1_PcG <- plot_violin_PcG(processed_PcG$H3K27me3, title = "H3K27me3")
p2_PcG <- plot_violin_PcG(processed_PcG$Pc, title = "Pc")
p3_PcG <- plot_violin_PcG(processed_PcG$H2Aub118, title = "H2Aub118")
p4_PcG <- plot_violin_PcG(processed_PcG$H3K27ac, title = "H3K27ac")

outFile <- "Enrichment_at_TSS_of_PcGtargets_genes.tab"
write.table(processed_PcG,file=outFile,sep="\t",quote=F,row.names=F,)

combined_PcG <- ggarrange(p1_PcG, p2_PcG, p3_PcG, p4_PcG,
                      ncol = 2, nrow = 2,   # 2x2 grid
                      common.legend = TRUE, # share legend if you want
                      legend = "none")

pdf(paste0("violinPlots_for_allGenes.pdf"), width=8)
print(combined)
print(combined_PcG)
dev.off