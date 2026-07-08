# Script inspired by https://www.datanovia.com/en/lessons/anova-in-r/#two-way-independent-anova

library(tidyverse)
library(ggpubr)
library(rstatix)
library(dunn.test)
library(dplyr)
library(FSA)

inFiles = list.files("./",pattern="1Ddistances_trans1Dinterval_all_domains_vs_all_domains_.*_data.tab")
#inFiles = list.files("./",pattern="1Ddistances_trans1Dinterval_all_domains_vs_all_domains_PcG_PcG_data.tab")
print(inFiles)

#subsetting <- "FALSE"
subsetting <- "TRUE"

#levels obtained from the input file
#rgb(red, green, blue, alpha, names = NULL, maxColorValue = 1)
# Using Generic RGB
#colors made below from the input file

for(inFile in inFiles)
{

    print(inFile, quote=F)

    allData <- read.table(inFile, header=F)
    colnames(allData) <- c("condition","values")
    allData <- allData[grep("NOnub",allData$condition),]
    allData$condition <- gsub("NOnub","",allData$condition)
    print(head(allData), quote=F)

    name <- gsub(".tab","",inFile)    
    outFile <- paste0("violinPlot_",name,"_with_stats.pdf")
    #if(file.exists(outFile)){next;}    
    data <- allData

    # List of specific comparisons (custom pairs)
    #    c("control" = "#0072B2", "treatment" = "#D55E00")
    levels <- unique(data$condition)
    print(levels)

    #colors = rep(c("#0072B2","#D55E00"),as.integer(length(levels)/2))
    colors = rep("gray80",as.integer(length(levels)))     
    custom_pairs <- list()
    n = 1
    for(i in 2:(length(levels)))
    {
        for(j in i:(length(levels)))
        {
            if(i==j){next;}
            custom_pairs[[n]] <- c(levels[[i]],levels[[j]])
	    n = n + 1
        }
    }
    print(levels)
    print(colors)
    print(custom_pairs)
    data$condition = factor(data$condition, levels=levels)
    print(head(data))

    print(paste0("### Get summary statistics ###"), quote=F)
    dataStats <- data %>%
	  group_by(condition) %>%
	  get_summary_stats(values, type = "mean_sd")
    print(as.data.frame(dataStats), quote=F)
    #quit()

    summ <- data %>%
           group_by(condition) %>%
	   summarize(n = n(), condition=condition, score = min(min(data$values),-min(data$values))*(50))
    summ <- unique(summ)
    print(as.data.frame(summ), quote=F)   
    #quit()

    bxp <- ggviolin(data, x="condition", y="values", color="black", fill="condition", alpha=.50, scale = "width", trim = TRUE) + scale_fill_manual(values=colors) +
           geom_boxplot(position=position_dodge(1), width=0.1, color="black", outlier.shape = NA) +
           #geom_point(position = position_jitter(width=0.4), size = 0.5, color="black", alpha=1.0) +	   	   
	   labs(y="Genomic domain-domain distance", x = "") +	   
	   geom_text(data=summ, aes(x=condition, y=score, label = paste0("n = ",n)), size=4., color="black") +    
	   theme_classic() +
   	   theme(axis.text.x = element_text(angle = 60, hjust=1), legend.position = "none")
	
    outFile <- paste0("violinPlot_",name,".pdf")
    #pdf(outFile)
    #print(bxp)
    #dev.off()
    #quit()
    
    print(paste0("### Check assumptions ###"), quote=F)
    print(paste0("1) Checking the presence of outliers"), quote=F)
    checkOutliers <- data %>%
	  group_by(condition) %>%
	  identify_outliers(values)
    if(nrow(checkOutliers) == 0)
    {
        print(paste0("No outliers found: This assumption is verified"), quote=F)
    } else {
        print(paste0("### WARNING: The following data-points are outiers!"), quote=F)
        print(checkOutliers, quote=F)
        print(paste0("Note that, in the situation where you have extreme outliers, this can be due to:"), quote=F)
        print(paste0("1) data entry errors, measurement errors or unusual values."), quote=F)
        print(paste0("You can include the outlier in the analysis anyway if you do not believe the result will be substantially affected."), quote=F)
        print(paste0("This can be evaluated by comparing the result of the ANOVA test with and without the outlier."), quote=F)
        print(paste0("It’s also possible to keep the outliers in the data and perform robust ANOVA test using the WRS2 package."), quote=F)
        print("", quote=F)     	    
    }
    print("", quote=F)


    print(paste0("2) Checking normality assumption"), quote=F)
    print(paste0("# Building the linear model"), quote=F)
    model  <- lm(values ~ condition,
    	         data = data)
    print(paste0("# Creating a QQ plot of residuals"), quote=F)
    p <- ggqqplot(residuals(model))
    pdf(paste0("qqplot_",name,".pdf")) 
    print(p)
    dev.off()	
    print(paste0("### Compute Shapiro-Wilk test of normality ###"), quote=F)
    #checkNormality <- shapiro_test(residuals(model)[1:5])
    #print(checkNormality, quote=F)
    #if(checkNormality$p.value > 0.05)
    #{
    #    print(paste0("In the QQ plot, as all the points fall approximately along the reference line, we can assume normality."), quote=F)
#	print(paste0("This conclusion is supported by the Shapiro-Wilk test. The p-value is not significant (p = ",checkNormality$p.value,"),"), quote=F)
#	print(paste0("so we can assume normality."), quote=F)
#    } else {
#        print(paste0("The normality assumption is not supported by the Shapiro-Wilk test. The p-value is significant (p = ",checkNormality$p.value,"),"), quote=F)
#        print(paste0("so we cannot assume normality for this dataset."), quote=F)
        #    quit()
#    }

    print("", quote=F)
    print(paste0("3) Checking homogeneity of variance assumption"), quote=F)
    print(paste0("This can be checked using the Levene’s test:"), quote=F)
    checkHomogenityOfVariance <- data %>% levene_test(values ~ condition)
    print(checkHomogenityOfVariance, quote=F)
    if(checkHomogenityOfVariance$p > -0.05)
    {
        print(paste0("The Levene’s test is not significant (p > ",checkHomogenityOfVariance$p,"). Therefore, we can assume the homogeneity of variances in the different groups."), quote=F)
	print(paste0("### Computation of the one-way ANOVA test ###"), quote=F)
	print(paste0("In the R code below, the asterisk represents the interaction effect and the main effect of each variable (and all lower-order interactions)."), quote=F)
	res.aov <- kruskal_test(values ~ condition, data = data)
	#	res.aov <- data %>% anova_test(values ~ condition)
	print(res.aov, quote=F)

	significantInteraction = res.aov[res.aov$p < 0.05,]
	if(nrow(significantInteraction) > 0)
	{
	    print(paste0("There was a statistically significant interaction between locis"), quote=F)
	    print(significantInteraction, quote=F)
	} else {
	    print(paste0("There was no statistically significant interaction between loci"), quote=F)
	    #next
	}

	print("", quote=F)	
	print(paste0("### Post-hoct tests ###"), quote=F)

	print(paste0("### Procedure for significant one-way interaction ###"), quote=F)
	print(paste0("### Compute pairwise comparisons ###"), quote=F)

	print(paste0("To determine which group means are different. We’ll now perform multiple pairwise comparisons between the different conditions."), quote=F)

	print(paste0("You can run and interpret all possible pairwise comparisons using a Bonferroni adjustment."), quote=F)
	print(paste0("This can be easily done using the function emmeans_test() [rstatix package], a wrapper around the emmeans package,"), quote=F)
	print(paste0("which needs to be installed. Emmeans stands for estimated marginal means (aka least square means or adjusted means)."), quote=F)
	

	print(paste0("Compare the values of the different cellType levels by condition levels:"), quote=F)
	pairwiseTestsCondition <- dunn_test(values ~ condition, data = data, p.adjust.method = "BH")
	#pairwiseTestsCondition <- data %>%
	#		       pairwise_t_test(
	#		       values ~ condition, 
	#		       p.adjust.method = "BH"
	#		    )		      

	if(subsetting == "TRUE")
	{
	# Filter only desired comparisons
	if(exists("subset_t_test"))
        {
	    rm(subset_t_test)
	}
	i = 0
	ypos <- c()
	for(pair in custom_pairs)
	{
	    delta <- 2
	    yposStart <- max(data$values)+delta
	    
	    print(paste0(pair[[1]]," ",pair[[2]]))
	    subset <- pairwiseTestsCondition[(pairwiseTestsCondition$group1 == pair[[1]] & pairwiseTestsCondition$group2 == pair[[2]]) | (pairwiseTestsCondition$group1 == pair[[2]] & pairwiseTestsCondition$group2 == pair[[1]]),]
	    print(subset)
	    #padjust = pvaluesFromDEseq2[pvaluesFromDEseq2$condition2 == pair[[2]]]$pvalue
	    subset$p.adj
	    ypos <- c(ypos,yposStart + 0 * delta)
	    if(exists("subset_t_test"))
	    {
	        subset_t_test <- rbind(subset_t_test,subset)
	    } else {
		subset_t_test <- subset	       
	    }	    
	    i = i + 1
	}
	
	# Adjust p-values only for this subset
	subset_t_test$p.adj <- p.adjust(subset_t_test$p, method = "BH")
	pairwiseTestsCondition <- subset_t_test
	}
	print(as.data.frame(pairwiseTestsCondition), quote=F)
	#quit()

	pairwiseTestsCondition <- pairwiseTestsCondition %>%
                       dplyr::mutate(p.adj.formatted = sprintf("%.3g", p.adj))

    } else {
       print(paste0("The Levene’s test is significant (p = ",checkHomogenityOfVariance$p," < 0.05). Therefore, we cannot assume the homogeneity of variances in the different groups."), quote=F)

       print(paste0("The Welch one-way test is an alternative to the standard one-way ANOVA"), quote=F)
       print(paste0("in the situation where the homogeneity of variance can’t be assumed (i.e., Levene test is significant)."), quote=F)

       res.aov <- data %>% welch_anova_test(values ~ condition)
       print(res.aov, quote=F)

       print(paste0("In this case, the Games-Howell post hoc test or pairwise t-tests (with no assumption of equal variances) can be used to compare all possible combinations of group differences."), quote=F)

       pairwiseTestsCondition <- data %>% games_howell_test(values ~ condition)
	    
    }
    print("", quote=F)	

    print(paste0("### Visualization: barplots with p-values"), quote=F)
    #pairwiseTestsCondition <- as.data.frame(pairwiseTestsCondition)    		    
    print(as.data.frame(pairwiseTestsCondition), quote=F)
    outFileText <- paste0("violinPlot_",name,"_with_stats.tsv")
    print(outFileText)
    write.table(as.data.frame(pairwiseTestsCondition), file=outFileText,  row.names=F, sep="\t")    
    pairwiseTestsCondition <- pairwiseTestsCondition %>% 
    			      			       add_xy_position(x = "condition")

    print(pairwiseTestsCondition$y.position)
    print(min(pairwiseTestsCondition$y.position))
    d = pairwiseTestsCondition$y.position[[2]]-pairwiseTestsCondition$y.position[[1]]
    x <- 0:(length(pairwiseTestsCondition$y.position)-1)
    pairwiseTestsCondition$y.position <- min(pairwiseTestsCondition$y.position)+(d*6)*x
    print(pairwiseTestsCondition$y.position)

    print(as.data.frame(pairwiseTestsCondition), quote=F)

    finalBxp <- bxp +
	   #stat_pvalue_manual(data=pairwiseTestsCondition[pairwiseTestsCondition$p.adj < 0.05,], label="p.adj.formatted", label.size = 4, bracket.size = 0.1, tip.length = 0, ) +
	   stat_pvalue_manual(data=pairwiseTestsCondition, label="p.adj.formatted", label.size = 4, bracket.size = 0.4, tip.length = 0, ) +
	   labs(subtitle = get_test_label(res.aov, detailed = TRUE), caption = get_pwc_label(pairwiseTestsCondition)) +
	   theme(legend.position = "none",
               plot.title = element_text(size = 18, color = "black", face= "bold"),
               plot.subtitle = element_text(size = 16, face = "bold.italic")
               )

    outFile <- paste0("violinPlot_",name,"_with_stats.pdf")
    pdf(outFile)
    print(finalBxp)
    dev.off()	
    #quit()
    print("", quote=F)
    #quit()
} # Close cycle over inFile