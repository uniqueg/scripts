#######
### A. PRE-REQUISITES
### B. IMPORT CUT/SUT DATA
### C. UPDATE COORDINATES
### D. GROUP OBJECTS & CLEAN-UP
#######

### A. PRE-REQUISITES
# Remove all objects
rm(list = ls())
# Load libraries
library(gdata)
library(rtracklayer)
###

### B. IMPORT CUT/SUT DATA
# Xu et al., 2009, Nature 457:1033-1037
# Import sheets for CUTs, SUTs and all transcripts
xu_CUTs_xls <- read.xls("http://www.nature.com/nature/journal/v457/n7232/extref/nature07728-s2.xls",sheet=3)
xu_SUTs_xls <- read.xls("http://www.nature.com/nature/journal/v457/n7232/extref/nature07728-s2.xls",sheet=2)
xu_all_transcripts_xls <- read.xls("http://www.nature.com/nature/journal/v457/n7232/extref/nature07728-s2.xls", sheet=1)
# Edit chromosome identifiers to match chain file
xu_CUTs_seq_levels <- paste0("chr", xu_CUTs_xls$chr)
xu_SUTs_seq_levels <- paste0("chr", xu_SUTs_xls$chr)
xu_all_transcripts_seq_levels <- paste0("chr", xu_all_transcripts_xls$chr)
# GRanges coercion
xu_CUTs_SacCer1 <- GRanges(xu_CUTs_seq_levels, IRanges(xu_CUTs_xls$start, xu_CUTs_xls$end), strand=xu_CUTs_xls$strand)
xu_SUTs_SacCer1 <- GRanges(xu_SUTs_seq_levels, IRanges(xu_SUTs_xls$start, xu_SUTs_xls$end), strand=xu_SUTs_xls$strand)
xu_all_transcripts_SacCer1 <- GRanges(xu_all_transcripts_seq_levels, IRanges(xu_all_transcripts_xls$start, xu_all_transcripts_xls$end), strand=xu_all_transcripts_xls$strand)
# Remove unused objects
rm(xu_CUTs_xls, xu_SUTs_xls, xu_all_transcripts_xls, xu_CUTs_seq_levels, xu_SUTs_seq_levels, xu_all_transcripts_seq_levels)
###

### C. MIGRATION TO CURRENT COORDINATES
# Update genomic coordinates from SacCer1 to SacCer3
# Download chain file for rtracklayer::liftOver from UCSC
download.file("http://hgdownload.cse.ucsc.edu/goldenPath/sacCer1/liftOver/sacCer1ToSacCer3.over.chain.gz", destfile="~/riboRNA/genome_files/sacCer1ToSacCer3.over.chain.gz")
# Unpack chain file
system("gunzip ~/riboRNA/genome_files/sacCer1ToSacCer3.over.chain.gz")
# Import chain file
SacCer1_3_chain <- import.chain("/home/akanitz/riboRNA/genome_files/sacCer1ToSacCer3.over.chain")
# Use rtracklayer::liftOver method to migrate genomic coordinates
xu_CUTs <- liftOver(xu_CUTs_SacCer1, SacCer1_3_chain)
xu_SUTs <- liftOver(xu_SUTs_SacCer1, SacCer1_3_chain)
xu_all_transcripts <- liftOver(xu_all_transcripts_SacCer1, SacCer1_3_chain)
# Convert GRangesList to GRanges objects
xu_CUTs <- unlist(xu_CUTs)
xu_SUTs <- unlist(xu_SUTs)
xu_all_transcripts <- unlist(xu_all_transcripts)
# Save objects
save(xu_CUTs, file=file.path("/home/akanitz/riboRNA/R_files/", "xu_CUTs.rda"))
save(xu_SUTs, file=file.path("/home/akanitz/riboRNA/R_files/", "xu_SUTs.rda"))
save(xu_all_transcripts, file=file.path("/home/akanitz/riboRNA/R_files/", "xu_all_transcripts.rda"))
# Remove unused objects
rm(xu_CUTs_SacCer1, xu_SUTs_SacCer1, xu_all_transcripts_SacCer1, SacCer1_3_chain)
###

### D. GROUP OBJECTS & CLEAN-UP
# GRangesList of all objects
list_xu_objects <- GRangesList()
list_xu_objects$xu_CUTs <- xu_CUTs
list_xu_objects$xu_SUTs <- xu_SUTs
list_xu_objects$xu_all <- xu_all_transcripts
save(list_xu_objects, file=file.path("~/riboRNA/R_files/", "list_xu_objects.rda"))
# Image file of all objects
save.image(file = "image_xu_objects.Rdata")
# Remove all objects
rm(list = ls())
# Remove unused packages
detach(package:gdata)
detach(package:rtracklayer)
detach(package:GenomicRanges)
detach(package:IRanges)
detach(package:BiocGenerics)
###
