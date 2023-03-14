file_1 <- "FILE1"
file_2 <- "FILE2" 

## Mouse
#file_1 <- "3868_PAS_ENSEMBL75_overlaps_genes"
#file_2 <- "4210_PAS_ENSEMBL75_overlaps_genes"

## Human
#file_1 <- "3843_PAS_ENSEMBL75_overlaps_genes"
#file_2 <- "4211_PAS_ENSEMBL75_overlaps_genes"

## Process sample 1
process_sample <- function (file) {
	# Read file
	sample_1 <- read.table(file_1, sep="\t", header=F, stringsAsFactors=F)[,c(5,15)]
	# Remove semicolons
	test$V15 <- gsub(";", "", test$V15)
	# 
	
	# Aggregate counts by gene ID
	counts_1 <- aggregate(sample_1$V5, list(sample_1$V15), sum)
	# 
	counts_1[,1] <- gsub("gene_id ", "", counts_1[,1])
	total_1 <- sum(counts_1[,2])
	counts_1 <- cbind(counts_1, counts_1[,2] / total_1 * 1000000)
}

library(plyr)
data.frame(test[,1:14], do.call(rbind, strsplit(as.character(test$V15),' ')))

before <- data.frame(attr = c(1,30,4,6), type=c('foo_and_bar','foo_and_bar_2'))  
out <- strsplit(as.character(before$type),'; ') 
do.call(rbind, out)
data.frame(before$attr, do.call(rbind, out))

[,1]  [,2]   
[1,] "foo" "bar"  
[2,] "foo" "bar_2"
[3,] "foo" "bar"  
[4,] "foo" "bar_2"
And to combine:
		


sample_1 <- read.table(file_1, sep="\t", header=F, stringsAsFactors=F)[,c(5,15)]
counts_1 <- aggregate(sample_1$V5, list(sample_1$V15), sum)
counts_1[,1] <- gsub("gene_id ", "", counts_1[,1])
total_1 <- sum(counts_1[,2])
counts_1 <- cbind(counts_1, counts_1[,2] / total_1 * 1000000)

## Process sample 2
sample_2 <- read.table(file_2, sep="\t", header=F, stringsAsFactors=F)[,c(5,15)]
counts_2 <- aggregate(sample_2$V5, list(sample_2$V15), sum)
counts_2[,1] <- gsub("gene_id ", "", counts_2[,1])
total_2 <- sum(counts_2[,2])
counts_2 <- cbind(counts_2, counts_2[,2] / total_2 * 1000000)

## Merge and produce x and y vectors
merged_counts_all <- merge(counts_1[,-2], counts_2[,-2], by=1, all=TRUE)
merged_counts <- merge(counts_1[,-2], counts_2[,-2], by=1)
x <- log2(merged_counts[,2])
y <- log2(merged_counts[,3])

## Plot
# fit y axis label
# include zeros?
# include gene numbers (repl. 1, repl. 2, plotted, all)?
# make more script-like/customizable (add correlation, generate labels as automatically as possible)
##
suppressPackageStartupMessages(library(geneplotter))
suppressPackageStartupMessages(library(RColorBrewer))


max <- ceiling(max(x,y))

#pdf("Jurkat_replicates_agreement.pdf")
colors  <- densCols(x, y, colramp=colorRampPalette(c("grey80","black")))
order <- rev(order(colors))
plot(0, type="n", xlim=c(0,max), ylim=c(0,max), xlab="Normalized A-seq reads (log2 CPM) per gene\nReplicate 1", ylab="Normalized A-seq reads (log2 CPM) per gene\nReplicate 2", main="Agreement of A-seq counts between replicates\nCell line: Jurkat (Homo sapiens)")
points(x[order], y[order], col=colors[order], pch=16)
text(max*0.125,max*0.95, "R = 0.9759", cex=1.25)
abline(0,1)
dev.off()
~                     