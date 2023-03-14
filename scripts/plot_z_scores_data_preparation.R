## PLOT YH2AX Z-SCORES AND ASISI SITES 

# Read Z score tables
z_ind <- read.table("z_all_induced", sep="\t")
z_unind <- read.table("z_all_uninduced", sep="\t")
# Save/load raw data objects: save(z_ind, z_unind, file="z_all_raw.R") / load("z_all_raw.R")

# Filter positive Z scores
z_ind_pos <- z_ind[z_ind[,7] > 0, ]
z_unind_pos <- z_unind[z_unind[,7] > 0, ]

# Import AsiSI motif ranges as GRanges object
library(rtracklayer)
asi <- import("AsiSI_orig.bed", genome="hg19", asRangedData=FALSE)

# Split data by chromosome
z_ind_chr_ls <- split(z_ind_pos, z_ind_pos[,1], drop=TRUE)
z_unind_chr_ls <- split(z_unind_pos, z_unind_pos[,1], drop=TRUE)
asi_chr_ls <- split(asi, seqnames(asi), drop=TRUE)

# Extract genome and chromosome information
gen <- as.character(unique(genome(asi_chr_ls)))
chr <- unique(c(names(z_ind_chr_ls), names(z_unind_chr_ls), names(asi_chr_ls)))
chr <- chr[!grepl("_|M", chr)]					# Remove mitochondrial and incomplete/unknown chromosomes
chr_ls <- split(chr, chr)

# Determine max Z score in whole dataset (used as common y axis limit in the plotting of all chromosomes)
max_z_ind <- sapply(z_ind_chr_ls, function(x) {max(x[,7])})
max_z_unind <- sapply(z_unind_chr_ls, function(x) {max(x[,7])})
max_z_chr <- pmax(max_z_ind, max_z_unind)
max_z_chr_ls <- split(max_z_chr, names(max_z_chr))
max_z_all <- max(c(max_z_ind, max_z_unind))

# Save essential objects:
save(z_ind_chr_ls, z_unind_chr_ls, asi_chr_ls, chr_ls, max_z_chr_ls, gen, max_z_all, file="z_score_plots_essential.R")
# Load with: load("z_score_plots_essential.R")


# - ? use ALL chromosomes, not just canonical --> remove line: chr <- chr[!grepl("_|M", chr)]
# - add defined yH2Ax cluster data