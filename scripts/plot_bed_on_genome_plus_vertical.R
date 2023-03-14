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

## AsiSI positions
v <- read.table("/home/kanitz/Dropbox/Work/BioZ/sync/AsiSI_occurrences/AsiSI_orig_mid", col.names = c("chromosome", "position"))
vls <- split(v, v$chromosome)

## INDUCED
# Load data 
i <- read.table("/home/kanitz/Dropbox/Work/Eclipse/general/test_files/chip_peaks_induced_mid", col.names = c("chromosome", "position"))
# Split by chromosome
ils <- split(i, i$chromosome)
## Traverse through each chromosome
dump <- sapply(attributes(ils)$names, function(chr) {
	## Print to pdf
	pdf(paste0(chr,"_induced.pdf"), width = 48, height = 16)
	# Plot histogram: Counts vs chromosome position; bin size = 20000
	hist(ils[[chr]]$position, breaks=seq(1, max(ils[[chr]]$position) + 50000, by= 50000), main=paste0("ChIP peaks across chromosome ", chr,", induced"), xlab="Chromosome position", ylab="Counts", ylim=c(0,20), axes=FALSE)
	axis(1,at=seq(0,max(ils[[chr]]$position)+50000,500000), labels=F)
	axis(1,at=seq(0,max(ils[[chr]]$position)+50000,10000000), lwd.ticks=3)
	axis(2,at=seq(0,20,2))
	abline(v = vls[[chr]]$position, col="red")
	dev.off()			
})

## UNINDUCED
# Load data
u <- read.table("/home/kanitz/Dropbox/Work/Eclipse/general/test_files/chip_peaks_uninduced_mid", col.names = c("chromosome", "position"))
# Split by chromosome
uls <- split(u, u$chromosome)
## Traverse through each chromosome
dump <- sapply(attributes(uls)$names, function(chr) {
	## Plot
	pdf(paste0(chr,"_uninduced.pdf"), width = 48, height = 16)
	# Plot histogram: Counts vs chromosome position; bin size = 20000
	hist(uls[[chr]]$position, breaks=seq(1, max(uls[[chr]]$position) + 50000, by= 50000), main=paste0("ChIP peaks across chromosome ", chr,", uninduced"), xlab="Chromosome position", ylab="Counts", ylim=c(0,20), axes=FALSE)
	axis(1,at=seq(0,max(uls[[chr]]$position)+50000,500000), labels=F)
	axis(1,at=seq(0,max(uls[[chr]]$position)+50000,10000000), lwd.ticks=3)
	axis(2,at=seq(0,25,2))	
	abline(v = vls[[chr]]$position, col="red")
	dev.off()			
})