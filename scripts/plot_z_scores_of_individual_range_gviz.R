### Author: Alexander Kanitz
### Created: 22-MAR-2013
### Modified: 22-MAR-2013
### Description: Plot Z scores and AsiSI restriction sites for the indicated range
### Arguments: 1. R data object [FILE]; 2. chromosome [CHAR]; 3. start of range [INT]; 4. end of range [INT]; 5. path to output files [PATH]
### Output: PDF and PNG plots of Z scores (induced and uninduced) and AsiSI sites 
### Usage: ./Rscript ./plot_z_scores_of_individual_range_gviz.R ./z_score_plots_essential.R chr1 10000000 20000000 ./

### A. PRE-REQUISITES
# Initialize and get command line arguments
args <- commandArgs(trailingOnly=TRUE)
## Load libraries
library(Gviz)
library(rtracklayer)
# Disable scientific notation
options(scipen=999)
###

### B. LOAD DATA & PASS ARGUMENTS
load(args[1])					## Objects: z_ind_chr_ls, z_unind_chr_ls, asi_chr_ls, chr_ls, gen, max_z_chr_ls, max_z_all -> lists of Z scores (induced and uninduced; data frames), AsiSI sites (GRanges), chromosomes (char vectors), genome (scalar: "hg19") & maximum Z scores per chromosome (numeric vectors) and overall maximum Z score (not needed for script) 
chr <- args[2]
start <- as.integer(args[3])
end <- as.integer(args[4])
dir <- paste(sub("/$", "", args[5]), "/", sep="")
###

### C. SUBSET DATA
z_ind_local <- z_ind_chr_ls[[chr]][z_ind_chr_ls[[chr]][,4] >= start & z_ind_chr_ls[[chr]][,4] <= end, ]
z_unind_local <- z_unind_chr_ls[[chr]][z_unind_chr_ls[[chr]][,4] >= start & z_unind_chr_ls[[chr]][,4] <= end, ]
asi_local <- asi_chr_ls[[chr]][end(asi_chr_ls[[chr]]) >= start & start(asi_chr_ls[[chr]]) <= end, ]
	#names(asi_local) <- asi_local$name
	#asi_local <- ranges(asi_local) # GRanges to IRanges
max_z_local <- max(c(z_ind_local[,7], z_unind_local[,7]))

### D. PREPARE TRACKS TO PLOT
i_track <- IdeogramTrack(genome=gen, chromosome=chr, showId=FALSE, collapse=TRUE)
g_track <- GenomeAxisTrack(range=ranges(asi_local), fill.range="brown", col.range="brown") #, showId=TRUE, cex.id=2, col.id=NA)
z_ind_track <- DataTrack(start=z_ind_local[,4], width=1, data=z_ind_local[,7], chromosome=chr, genome=gen, name="Z score induced", type=c("smooth", "p"), ylim=c(0, ceiling(max_z_local)), lwd=2.5, col.line="blue", span=1/20, evaluation=500, col.title="black", col.axis="black", background.title="transparent", cex.title=0.9, baseline=3, col.baseline="red", lwd.baseline=1.5, lty.baseline="dashed")
z_unind_track <- DataTrack(start=z_unind_local[,4], width=1, data=z_unind_local[,7], chromosome=chr, genome=gen, name="Z score uninduced", type=c("smooth", "p"), ylim=c(ceiling(max_z_local), 0), lwd=2.5, col="orange2", col.line="brown", span=1/20, evaluation=500, col.title="black", col.axis="black", background.title="transparent", cex.title=0.9, baseline=3, col.baseline="red", lwd.baseline=1.5, lty.baseline="dashed")

### E. PLOT
## Prepare file basename & plot title
basename <- paste(dir, paste("z_scores", chr, start, end, sep="_"), sep="") 
title <- paste("Chromosome ", sub("^chr", "", chr, ignore.case=TRUE), ", ", start, "-", end, sep="")
## Plot PNG
png(paste(basename, ".png", sep=""), units="in", res=600, width=12, height=6)
plotTracks(list(i_track, z_ind_track, g_track, z_unind_track), from=start, to=end, chromosome=chr, main=title)
dev.off()
## Plot PDF
pdf(paste(basename, ".pdf", sep=""), width=12, height=6)
plotTracks(list(i_track, z_ind_track, g_track, z_unind_track), from=start, to=end, chromosome=chr, main=title)
dev.off()
# Re-enable scientific notation
options(scipen=0)
###

## PROBLEMS
# - size of regions in GenomeAxisTrack
# - problem of scale line widths in png > no regions!!
# - incorporate probable cut sites, i.e. selected yH2Ax clusters