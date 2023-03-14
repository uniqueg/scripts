### Author: Alexander Kanitz
### Created: 22-MAR-2013
### Modified: 22-MAR-2013
### Description: For each chromosome, plots the Z scores and AsiSI restriction sites
### Arguments: 1. R data object; 2. path to output files
### Output: PDF and PNG plots of Z scores (induced and uninduced) and AsiSI sites; the y (Z score) axis max labeling will be derived from the highest Z score found in each chromosome (induced or uninduced), i.e. different scales will be used for different chromosomes
### Usage: ./Rscript ./plot_z_scores_all_chromosomes_gviz.R ./z_score_plots_essential.R ./

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
dir <- paste(sub("/$", "", args[2]), "/", sep="")
###

### C. PLOT EACH CHROMOSOME INDIVIDUALLY WITH MAPPLY
dump <- lapply(chr_ls, function(chr, asi, ind, unind, max) {

	### C1. PREPARE TRACKS TO PLOT
	i_track <- IdeogramTrack(genome=gen, chromosome=chr, showId=FALSE, fill=FALSE, col=FALSE, collapse=TRUE)
	g_track <- GenomeAxisTrack(range=ranges(asi[[chr]]), fill.range="brown", col.range="brown") #, showId=TRUE, cex.id=2, col.id=NA)
	z_ind_track <- DataTrack(start=ind[[chr]][,4], width=1, data=ind[[chr]][,7], chromosome=chr, genome=gen, name="Z score induced", type=c("smooth", "p"), ylim=c(0, ceiling(max[[chr]])), lwd=2.5, col.line="blue", span=1/20, evaluation=500, col.title="black", col.axis="black", background.title="transparent", cex.title=0.9, baseline=3, col.baseline="red", lwd.baseline=1.5, lty.baseline="dashed")
	z_unind_track <- DataTrack(start=unind[[chr]][,4], width=1, data=unind[[chr]][,7], chromosome=chr, genome=gen, name="Z score uninduced", type=c("smooth", "p"), ylim=c(ceiling(max[[chr]]), 0), lwd=2.5, col="orange2", col.line="brown", span=1/20, evaluation=500, col.title="black", col.axis="black", background.title="transparent", cex.title=0.9, baseline=3, col.baseline="red", lwd.baseline=1.5, lty.baseline="dashed")
	
	### C2. PLOT
	## Prepare file basename & plot title
	basename <- paste(dir, "z_scores_", chr, sep="") 
	title <- paste("Chromosome", sub("^chr", "", chr, ignore.case=TRUE), sep=" ")
	## Plot PNG
	png(paste(basename, ".png", sep=""), units="in", res=600, width=12, height=6)
	plotTracks(list(i_track, z_ind_track, g_track, z_unind_track), from=1, to=seqlengths(asi[[chr]])[chr], chromosome=chr, main=title)
	dev.off()
	## Plot PDF
#	pdf(paste(basename, ".pdf", sep=""), width=12, height=6)
#	plotTracks(list(i_track, z_ind_track, g_track, z_unind_track), from=1, to=seqlengths(asi[[chr]])[chr], chromosome=chr, main=title)
#	dev.off()
			
}, asi_chr_ls, z_ind_chr_ls, z_unind_chr_ls, max_z_chr_ls)
# Re-enable scientific notation
options(scipen=0)
###

## IMPLEMENT / BUGS
# - write script for same scale for all (use max_z_all or optionally allow specification of max)
# - add plot type manually
# - bug with apparent negative widths?!?! checked: not true! got better when using IRanges instead of GRanges, but in chr7 (and possibly chr8/9) there is still a problem, probably for AsiSI sites (all other chromosomes ok)
# - size of regions in GenomeAxisTrack
# - problem of scale line widths in png > no regions!!
# - incorporate probable cut sites, i.e. selected yH2Ax clusters
# - only plots canonical 24 chromosomes >> change in data preparation script if required!