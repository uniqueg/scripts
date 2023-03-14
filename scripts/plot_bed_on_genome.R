### Author: Alexander Kanitz
### Created: 22-JAN-2013
### Modified: 22-JAN-2013
### Description: 
### Arguments: bed file
### Output: 
### Usage: perl 

### A. PRE-REQUISITES
# Initialize and get command line arguments
args <- commandArgs(trailingOnly=TRUE)
###

## INDUCED
# Load data 
i <- read.table("chip_peaks_induced_mid", col.names = c("chromosome", "position"))
# Split by chromosome
ils <- split(i, i$chromosome)
## Traverse through each chromosome
dump <- sapply(attributes(ils)$names, function(chr) {
	## Print to pdf
	pdf(paste0(chr,"_induced.pdf"), width = 24, height = 8)
	# Plot histogram: Counts vs chromosome position; bin size = 20000
	hist(ils[[chr]]$position, breaks=seq(1, max(ils[[chr]]$position) + 50000, by= 50000), main=paste0("ChIP peaks across chromosome ", chr,", induced"), xlab="Chromosome position", ylab="Counts")
	dev.off()			
})

## UNINDUCED
# Load data
u <- read.table("chip_peaks_uninduced_mid", col.names = c("chromosome", "position"))
# Split by chromosome
uls <- split(u, u$chromosome)
## Traverse through each chromosome
dump <- sapply(attributes(uls)$names, function(chr) {
	## Plot
	pdf(paste0(chr,"_uninduced.pdf"), width = 24, height = 8)
	# Plot histogram: Counts vs chromosome position; bin size = 20000
	hist(uls[[chr]]$position, breaks=seq(1, max(uls[[chr]]$position) + 50000, by= 50000), main=paste0("ChIP peaks across chromosome ", chr,", uninduced"), xlab="Chromosome position", ylab="Counts")
	dev.off()			
})