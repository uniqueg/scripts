#######
### A. PRE-REQUISITES
### B. IMPORT AS GAPPEDALIGNMENTS OBJECTS
### C. IMPORT AS BAMFILELIST
### D. CLEAN-UP
#######

### A. PRE-REQUISITES
# Remove all objects
rm(list = ls())
# Load libraries
library(Rsamtools)
###

### B. IMPORT AS GAPPEDALIGNMENTS OBJECTS
# Read as GappedAlignments objects
reads_A <- readGappedAlignments("/home/akanitz/riboRNA/alignment_files/riboA_s.bam")
reads_A2 <- readGappedAlignments("/home/akanitz/riboRNA/alignment_files/riboA2_s.bam")
reads_B <- readGappedAlignments("/home/akanitz/riboRNA/alignment_files/riboB_s.bam")
# Save
save(reads_A, file=file.path("/home/akanitz/riboRNA/R_files/", "reads_A.rda"))
save(reads_A2, file=file.path("/home/akanitz/riboRNA/R_files/", "reads_A2.rda"))
save(reads_B, file=file.path("/home/akanitz/riboRNA/R_files/", "reads_B.rda"))
###

### C. IMPORT AS BAMFILELIST
# Use base::list.files method to create a char vector or .bam read files
read_files <- list.files("/home/akanitz/riboRNA/alignment_files","_s.bam$",full=TRUE)
# Use Rsamtools::BamFileList method to create .bam file list from char vector
list_read_files <- BamFileList(read_files)
# Save file list
save(list_read_files, file=file.path("~/riboRNA/R_files/", "list_read_files.rda"))
# Remove unused files
rm(read_files)
###

### E. GROUP OBJECTS & CLEAN-UP
# List of all objects
list_read_objects <- as.list.environment(.GlobalEnv)
save(list_read_objects, file=file.path("~/riboRNA/R_files/", "list_read_objects.rda"))
# Image file of all objects
save.image(file = "image_read_objects.Rdata")
# Remove all objects
rm(list = ls())
# Remove unused packages
detach(package:Rsamtools)
detach(package:Biostrings)
detach(package:GenomicRanges)
detach(package:IRanges)
detach(package:BiocGenerics)
###
