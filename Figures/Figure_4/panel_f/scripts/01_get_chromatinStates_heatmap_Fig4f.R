# Run clustering
library("ggdendro")
#library("data.table")
#suppressPackageStartupMessages(library(dendextend))
#suppressPackageStartupMessages(library(reshape2))
library(dendextend)
library("ggplot2")
library("gplots")
library("RColorBrewer")
#library("grid")
#library("reshape2")

# Color definition

dat <- read.table("_tmp",header=T,row.names=1)

#dat <- dat[ c(3, 1, 2), c(3, 1, 2)]
dat <- dat[ c(1, 2, 3, 4), c(1, 2, 3, 4)]

#dat <- as.matrix(dat)
#dat <- as.dist(dat, diag=F, upper=TRUE)
print(dat)

par(cex=0.5)

#hc <- hclust(dat, method="ward.D2")
#plot(hc)

#quit()

#dendro <- as.dendrogram(hclust(d = dat, method="ward.D2"))

#labels_colors(dendro) <- colorCodes[groupCodes][order.dendrogram(dendro)]
#dendro <- dendro %>%
#      color_branches(k = 100, col = colorCodes[groupCodes][order.dendrogram(dendro)]) %>%
#      set("branches_lwd", c(2)) %>%
#      set("branches_lty", c(1)) 

# Heatmap
#cols = c("darkred", "red", "white", "blue", "darkblue")
cols = c("darkblue", "blue", "gray", "red", "darkred")
mypalette <- colorRampPalette(cols)(100)
print(length(seq(XXXminXXX,XXXmaxXXX,length.out=101)))

#layout(matrix(1:1, ncol=1), height=500, width=4, respect=T)
#png(file="Rplots.png", width=1300,height=1200)

pdf(file="plots.pdf", width = 14, height = 12, pointsize = 10)
par(mar = c(6, 6, 4, 2))  # bottom, left, top, right
#sM <- 3.5
#tM <- (sM)/2	    
#par(mfrow = c(1, 1), mar = c(1, 1, 1, 1))

#p <- heatmap.2(as.matrix(dat), Rowv=F, Colv=F, breaks = seq(XXXminXXX,XXXmaxXXX,length.out=101), col = mypalette, symmbreaks=T, symkey=T, symm = TRUE, margins=c(7,14), trace="none", density.info="none", keysize=0.75, key.xlab="avg. scores", key.title=NA, cexRow=3, cexCol=3, key.par = list(cex=1.5), cellnote=round(dat, 2), notecex=3.0, notecol="black", srtCol = 90, labRow = rownames(dat))

heatmap.2(as.matrix(dat), Rowv=F, Colv=F, breaks = seq(XXXminXXX,XXXmaxXXX,length.out=101), col = mypalette, symmbreaks=T, symkey=T, symm = TRUE, margins=c(7,14), trace="none", density.info="none", keysize=0.75, key.xlab="avg. scores", key.title=NA, cexRow=3, cexCol=3, key.par = list(cex=1.5), cellnote=round(dat, 2), notecex=3.0, notecol="white", srtCol = 90, labRow = rownames(dat))
#print(p)
#ggsave("Rplots.pdf")
dev.off()

