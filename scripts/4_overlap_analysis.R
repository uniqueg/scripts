#######
### A. PRE-REQUISITES
### B. LOAD OBJECTS
### C. OVERLAP ANALYSIS GENOME
### D. OVERLAP ANALYSIS CUTS/SUTS
### E. GROUP OBJECTS & CLEAN-UP
#######

### A. PRE-REQUISITES
# Remove all objects in environment
rm(list = ls())
# Load libraries
library(Rsamtools)
###

### B. LOAD OBJECTS
# Genomic objects
load("image_sc68_objects.Rdata")
# CUT/SUT objects
load("image_xu_objects.Rdata")
# Read objects
load("list_read_files.rda")
###

### C. OVERLAP ANALYSIS GENOME
# GenomicRanges::summarizeOverlaps

## C1. MODE: UNION
# GENOMIC FEATURES
overlap_union_reads_genomic_type <- summarizeOverlaps(list_sc68_genomic_type, list_read_files, mode = "Union", ignore.strand = FALSE)
colnames(overlap_union_reads_genomic_type) <- gsub('^/.*(A2|A|B).*','union_\\1', colnames(overlap_union_reads_genomic_type))
save(overlap_union_reads_genomic_type, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_union_reads_genomic_type.rda"))
# GENE/EXON TYPES
overlap_union_reads_gene_type <- summarizeOverlaps(list_sc68_gene_type, list_read_files, mode = "Union", ignore.strand = FALSE)
colnames(overlap_union_reads_gene_type) <- gsub('^/.*(A2|A|B).*','union_\\1', colnames(overlap_union_reads_gene_type))
save(overlap_union_reads_gene_type, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_union_reads_gene_type.rda"))
# GENES
overlap_union_reads_genes <- summarizeOverlaps(sc68_genes, list_read_files, mode = "Union", ignore.strand = FALSE)
colnames(overlap_union_reads_genes) <- gsub('^/.*(A2|A|B).*','union_\\1', colnames(overlap_union_reads_genes))
save(overlap_union_reads_genes, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_union_reads_genes.rda"))
# GENES (ANTISENSE)
overlap_union_reads_genes_as <- summarizeOverlaps(sc68_genes_as, list_read_files, mode = "Union", ignore.strand = FALSE)
colnames(overlap_union_reads_genes_as) <- gsub('^/.*(A2|A|B).*','union_\\1', colnames(overlap_union_reads_genes_as))
save(overlap_union_reads_genes_as, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_union_reads_genes_as.rda"))
# EXONS
overlap_union_reads_exons <- summarizeOverlaps(sc68_exons, list_read_files, mode = "Union", ignore.strand = FALSE)
colnames(overlap_union_reads_exons) <- gsub('^/.*(A2|A|B).*','union_\\1', colnames(overlap_union_reads_exons))
save(overlap_union_reads_exons, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_union_reads_exons.rda"))
# EXONS (ANTISENSE)
overlap_union_reads_exons_as <- summarizeOverlaps(sc68_exons_as, list_read_files, mode = "Union", ignore.strand = FALSE)
colnames(overlap_union_reads_exons_as) <- gsub('^/.*(A2|A|B).*','union_\\1', colnames(overlap_union_reads_exons_as))
save(overlap_union_reads_exons_as, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_union_reads_exons_as.rda"))
# INTRONS
overlap_union_reads_introns <- summarizeOverlaps(sc68_introns, list_read_files, mode = "Union", ignore.strand = FALSE)
colnames(overlap_union_reads_introns) <- gsub('^/.*(A2|A|B).*','union_\\1', colnames(overlap_union_reads_introns))
save(overlap_union_reads_introns, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_union_reads_introns.rda"))
# INTRONS (ANTISENSE)
overlap_union_reads_introns_as <- summarizeOverlaps(sc68_introns_as, list_read_files, mode = "Union", ignore.strand = FALSE)
colnames(overlap_union_reads_introns_as) <- gsub('^/.*(A2|A|B).*','union_\\1', colnames(overlap_union_reads_introns_as))
save(overlap_union_reads_introns_as, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_union_reads_introns_as.rda"))
##

## C2. MODE: INTERSECTION STRICT
# GENOMIC FEATURES
overlap_int_strict_reads_genomic_type <- summarizeOverlaps(list_sc68_genomic_type, list_read_files, mode = "IntersectionStrict", ignore.strand = FALSE)
colnames(overlap_int_strict_reads_genomic_type) <- gsub('^/.*(A2|A|B).*','int_strict_\\1', colnames(overlap_int_strict_reads_genomic_type))
save(overlap_int_strict_reads_genomic_type, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_strict_reads_genomic_type.rda"))
# GENE/EXON TYPES
overlap_int_strict_reads_gene_type <- summarizeOverlaps(list_sc68_gene_type, list_read_files, mode = "IntersectionStrict", ignore.strand = FALSE)
colnames(overlap_int_strict_reads_gene_type) <- gsub('^/.*(A2|A|B).*','int_strict_\\1', colnames(overlap_int_strict_reads_gene_type))
save(overlap_int_strict_reads_gene_type, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_strict_reads_gene_type.rda"))
# GENES
overlap_int_strict_reads_genes <- summarizeOverlaps(sc68_genes, list_read_files, mode = "IntersectionStrict", ignore.strand = FALSE)
colnames(overlap_int_strict_reads_genes) <- gsub('^/.*(A2|A|B).*','int_strict_\\1', colnames(overlap_int_strict_reads_genes))
save(overlap_int_strict_reads_genes, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_strict_reads_genes.rda"))
# GENES (ANTISENSE)
overlap_int_strict_reads_genes_as <- summarizeOverlaps(sc68_genes_as, list_read_files, mode = "IntersectionStrict", ignore.strand = FALSE)
colnames(overlap_int_strict_reads_genes_as) <- gsub('^/.*(A2|A|B).*','int_strict_\\1', colnames(overlap_int_strict_reads_genes_as))
save(overlap_int_strict_reads_genes_as, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_strict_reads_genes_as.rda"))
# EXONS
overlap_int_strict_reads_exons <- summarizeOverlaps(sc68_exons, list_read_files, mode = "IntersectionStrict", ignore.strand = FALSE)
colnames(overlap_int_strict_reads_exons) <- gsub('^/.*(A2|A|B).*','int_strict_\\1', colnames(overlap_int_strict_reads_exons))
save(overlap_int_strict_reads_exons, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_strict_reads_exons.rda"))
# EXONS (ANTISENSE)
overlap_int_strict_reads_exons_as <- summarizeOverlaps(sc68_exons_as, list_read_files, mode = "IntersectionStrict", ignore.strand = FALSE)
colnames(overlap_int_strict_reads_exons_as) <- gsub('^/.*(A2|A|B).*','int_strict_\\1', colnames(overlap_int_strict_reads_exons_as))
save(overlap_int_strict_reads_exons_as, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_strict_reads_exons_as.rda"))
# INTRONS
overlap_int_strict_reads_introns <- summarizeOverlaps(sc68_introns, list_read_files, mode = "IntersectionStrict", ignore.strand = FALSE)
colnames(overlap_int_strict_reads_introns) <- gsub('^/.*(A2|A|B).*','int_strict_\\1', colnames(overlap_int_strict_reads_introns))
save(overlap_int_strict_reads_introns, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_strict_reads_introns.rda"))
# INTRONS (ANTISENSE)
overlap_int_strict_reads_introns_as <- summarizeOverlaps(sc68_introns_as, list_read_files, mode = "IntersectionStrict", ignore.strand = FALSE)
colnames(overlap_int_strict_reads_introns_as) <- gsub('^/.*(A2|A|B).*','int_strict_\\1', colnames(overlap_int_strict_reads_introns_as))
save(overlap_int_strict_reads_introns_as, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_strict_reads_introns_as.rda"))
##

## C3. MODE: INTERSECTION NOT EMPTY
# GENOMIC FEATURES
overlap_int_not_empty_reads_genomic_type <- summarizeOverlaps(list_sc68_genomic_type, list_read_files, mode = "IntersectionNotEmpty", ignore.strand = FALSE)
colnames(overlap_int_not_empty_reads_genomic_type) <- gsub('^/.*(A2|A|B).*','int_not_empty_\\1', colnames(overlap_int_not_empty_reads_genomic_type))
save(overlap_int_not_empty_reads_genomic_type, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_not_empty_reads_genomic_type.rda"))
# GENE/EXON TYPES
overlap_int_not_empty_reads_gene_type <- summarizeOverlaps(list_sc68_gene_type, list_read_files, mode = "IntersectionNotEmpty", ignore.strand = FALSE)
colnames(overlap_int_not_empty_reads_gene_type) <- gsub('^/.*(A2|A|B).*','int_not_empty_\\1', colnames(overlap_int_not_empty_reads_gene_type))
save(overlap_int_not_empty_reads_gene_type, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_not_empty_reads_gene_type.rda"))
# GENES
overlap_int_not_empty_reads_genes <- summarizeOverlaps(sc68_genes, list_read_files, mode = "IntersectionNotEmpty", ignore.strand = FALSE)
colnames(overlap_int_not_empty_reads_genes) <- gsub('^/.*(A2|A|B).*','int_not_empty_\\1', colnames(overlap_int_not_empty_reads_genes))
save(overlap_int_not_empty_reads_genes, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_not_empty_reads_genes.rda"))
# GENES (ANTISENSE)
overlap_int_not_empty_reads_genes_as <- summarizeOverlaps(sc68_genes_as, list_read_files, mode = "IntersectionNotEmpty", ignore.strand = FALSE)
colnames(overlap_int_not_empty_reads_genes_as) <- gsub('^/.*(A2|A|B).*','int_not_empty_\\1', colnames(overlap_int_not_empty_reads_genes_as))
save(overlap_int_not_empty_reads_genes_as, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_not_empty_reads_genes_as.rda"))
# EXONS
overlap_int_not_empty_reads_exons <- summarizeOverlaps(sc68_exons, list_read_files, mode = "IntersectionNotEmpty", ignore.strand = FALSE)
colnames(overlap_int_not_empty_reads_exons) <- gsub('^/.*(A2|A|B).*','int_not_empty_\\1', colnames(overlap_int_not_empty_reads_exons))
save(overlap_int_not_empty_reads_exons, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_not_empty_reads_exons.rda"))
# EXONS (ANTISENSE)
overlap_int_not_empty_reads_exons_as <- summarizeOverlaps(sc68_exons_as, list_read_files, mode = "IntersectionNotEmpty", ignore.strand = FALSE)
colnames(overlap_int_not_empty_reads_exons_as) <- gsub('^/.*(A2|A|B).*','int_not_empty_\\1', colnames(overlap_int_not_empty_reads_exons_as))
save(overlap_int_not_empty_reads_exons_as, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_not_empty_reads_exons_as.rda"))
# INTRONS
overlap_int_not_empty_reads_introns <- summarizeOverlaps(sc68_introns, list_read_files, mode = "IntersectionNotEmpty", ignore.strand = FALSE)
colnames(overlap_int_not_empty_reads_introns) <- gsub('^/.*(A2|A|B).*','int_not_empty_\\1', colnames(overlap_int_not_empty_reads_introns))
save(overlap_int_not_empty_reads_introns, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_not_empty_reads_introns.rda"))
# INTRONS (ANTISENSE)
overlap_int_not_empty_reads_introns_as <- summarizeOverlaps(sc68_introns_as, list_read_files, mode = "IntersectionNotEmpty", ignore.strand = FALSE)
colnames(overlap_int_not_empty_reads_introns_as) <- gsub('^/.*(A2|A|B).*','int_not_empty_\\1', colnames(overlap_int_not_empty_reads_introns_as))
save(overlap_int_not_empty_reads_introns_as, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_not_empty_reads_introns_as.rda"))
###

### D. OVERLAP ANALYSIS CUTS/SUTS
# GenomicRanges::summarizeOverlaps

## D1. MODE: UNION
# CUTS
overlap_union_reads_xu_CUTs <- summarizeOverlaps(xu_CUTs, list_read_files, mode = "Union", ignore.strand = FALSE)
colnames(overlap_union_reads_xu_CUTs) <- gsub('^/.*(A2|A|B).*','union_\\1', colnames(overlap_union_reads_xu_CUTs))
save(overlap_union_reads_xu_CUTs, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_union_reads_xu_CUTs.rda"))
# SUTS
overlap_union_reads_xu_SUTs <- summarizeOverlaps(xu_SUTs, list_read_files, mode = "Union", ignore.strand = FALSE)
colnames(overlap_union_reads_xu_SUTs) <- gsub('^/.*(A2|A|B).*','union_\\1', colnames(overlap_union_reads_xu_SUTs))
save(overlap_union_reads_xu_SUTs, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_union_reads_xu_SUTs.rda"))
# ALL TRANSCRIPTS
overlap_union_reads_xu_all <- summarizeOverlaps(xu_all_transcripts, list_read_files, mode = "Union", ignore.strand = FALSE)
colnames(overlap_union_reads_xu_all) <- gsub('^/.*(A2|A|B).*','union_\\1', colnames(overlap_union_reads_xu_all))
save(overlap_union_reads_xu_all, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_union_reads_xu_all.rda"))
##

## D2. MODE: INTERSECTION STRICT
# CUTS
overlap_int_strict_reads_xu_CUTs <- summarizeOverlaps(xu_CUTs, list_read_files, mode = "IntersectionStrict", ignore.strand = FALSE)
colnames(overlap_int_strict_reads_xu_CUTs) <- gsub('^/.*(A2|A|B).*','int_strict_\\1', colnames(overlap_int_strict_reads_xu_CUTs))
save(overlap_int_strict_reads_xu_CUTs, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_strict_reads_xu_CUTs.rda"))
# SUTS
overlap_int_strict_reads_xu_SUTs <- summarizeOverlaps(xu_SUTs, list_read_files, mode = "IntersectionStrict", ignore.strand = FALSE)
colnames(overlap_int_strict_reads_xu_SUTs) <- gsub('^/.*(A2|A|B).*','int_strict_\\1', colnames(overlap_int_strict_reads_xu_SUTs))
save(overlap_int_strict_reads_xu_SUTs, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_strict_reads_xu_SUTs.rda"))
# ALL TRANSCRIPTS
overlap_int_strict_reads_xu_all <- summarizeOverlaps(xu_all_transcripts, list_read_files, mode = "IntersectionStrict", ignore.strand = FALSE)
colnames(overlap_int_strict_reads_xu_all) <- gsub('^/.*(A2|A|B).*','int_strict_\\1', colnames(overlap_int_strict_reads_xu_all))
save(overlap_int_strict_reads_xu_all, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_strict_reads_xu_all.rda"))
##

## D3. MODE: INTERSECTION NOT EMPTY
# CUTS
overlap_int_not_empty_reads_xu_CUTs <- summarizeOverlaps(xu_CUTs, list_read_files, mode = "IntersectionNotEmpty", ignore.strand = FALSE)
colnames(overlap_int_not_empty_reads_xu_CUTs) <- gsub('^/.*(A2|A|B).*','int_not_empty_\\1', colnames(overlap_int_not_empty_reads_xu_CUTs))
save(overlap_int_not_empty_reads_xu_CUTs, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_not_empty_reads_xu_CUTs.rda"))
# SUTS
overlap_int_not_empty_reads_xu_SUTs <- summarizeOverlaps(xu_SUTs, list_read_files, mode = "IntersectionNotEmpty", ignore.strand = FALSE)
colnames(overlap_int_not_empty_reads_xu_SUTs) <- gsub('^/.*(A2|A|B).*','int_not_empty_\\1', colnames(overlap_int_not_empty_reads_xu_SUTs))
save(overlap_int_not_empty_reads_xu_SUTs, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_not_empty_reads_xu_SUTs.rda"))
# ALL TRANSCRIPTS
overlap_int_not_empty_reads_xu_all <- summarizeOverlaps(xu_all_transcripts, list_read_files, mode = "IntersectionNotEmpty", ignore.strand = FALSE)
colnames(overlap_int_not_empty_reads_xu_all) <- gsub('^/.*(A2|A|B).*','int_not_empty_\\1', colnames(overlap_int_not_empty_reads_xu_all))
save(overlap_int_not_empty_reads_xu_all, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_not_empty_reads_xu_all.rda"))
###


### E. GROUP OBJECTS & CLEAN-UP

## E1. LIST OF GENOMIC TYPE OVERLAPS
# Remove all objects in environment
rm(list = ls())
# Load objects
lapply(list.files(pattern="^overlap.*genomic_type.rda$"), load, .GlobalEnv)
# Generate list from global environment
list_overlap_genomic_type_objects <- as.list.environment(.GlobalEnv)
# Save list
save(list_overlap_genomic_type_objects, file=file.path("~/riboRNA/R_files/", "list_overlap_genomic_type_objects.rda"))
##

## E2. LIST OF GENE TYPE OVERLAPS
# Remove all objects in environment
rm(list = ls())
# Load objects
lapply(list.files(pattern="^overlap.*gene_type.rda$"), load, .GlobalEnv)
# Generate list from global environment
list_overlap_gene_type_objects <- as.list.environment(.GlobalEnv)
# Save list
save(list_overlap_gene_type_objects, file=file.path("~/riboRNA/R_files/", "list_overlap_gene_type_objects.rda"))
##

## E3. LISTS OF REMAINING (NOT MUTUALLY EXCLUSIVE) SC68 FEATURES OVERLAPS
# UNION
# Remove all objects in environment
rm(list = ls())
# Load objects
lapply(list.files(pattern="^overlap_union.*(genes|exons|introns)(_as.|.)rda$"), load, .GlobalEnv)
# Generate list from global environment
list_overlap_union_other_objects <- as.list.environment(.GlobalEnv)
# Save list
save(list_overlap_union_other_objects, file=file.path("~/riboRNA/R_files/", "list_overlap_union_other_objects.rda"))
# INT STRICT
# Remove all objects in environment
rm(list = ls())
# Load objects
lapply(list.files(pattern="^overlap_int_strict.*(genes|exons|introns)(_as.|.)rda$"), load, .GlobalEnv)
# Generate list from global environment
list_overlap_int_strict_other_objects <- as.list.environment(.GlobalEnv)
# Save list
save(list_overlap_int_strict_other_objects, file=file.path("~/riboRNA/R_files/", "list_overlap_int_strict_other_objects.rda"))
# INT NOT EMPTY
# Remove all objects in environment
rm(list = ls())
# Load objects
lapply(list.files(pattern="^overlap_int_not_empty.*(genes|exons|introns)(_as.|.)rda$"), load, .GlobalEnv)
# Generate list from global environment
list_overlap_int_not_empty_other_objects <- as.list.environment(.GlobalEnv)
# Save list
save(list_overlap_int_not_empty_other_objects, file=file.path("~/riboRNA/R_files/", "list_overlap_int_not_empty_other_objects.rda"))
##

## E4. LISTS OF CUT/SUT OVERLAPS
# UNION
# Remove all objects in environment
rm(list = ls())
# Load objects
lapply(list.files(pattern="^overlap_union.*xu.*rda$"), load, .GlobalEnv)
# Generate list from global environment
list_overlap_union_xu_objects <- as.list.environment(.GlobalEnv)
# Save list
save(list_overlap_union_xu_objects, file=file.path("~/riboRNA/R_files/", "list_overlap_union_xu_objects.rda"))
# INT STRICT
# Remove all objects in environment
rm(list = ls())
# Load objects
lapply(list.files(pattern="^overlap_int_strict.*xu.*rda$"), load, .GlobalEnv)
# Generate list from global environment
list_overlap_int_strict_xu_objects <- as.list.environment(.GlobalEnv)
# Save list
save(list_overlap_int_strict_xu_objects, file=file.path("~/riboRNA/R_files/", "list_overlap_int_strict_xu_objects.rda"))
# INT NOT EMPTY
# Remove all objects in environment
rm(list = ls())
# Load objects
lapply(list.files(pattern="^overlap_int_not_empty.*xu.*rda$"), load, .GlobalEnv)
# Generate list from global environment
list_overlap_int_not_empty_xu_objects <- as.list.environment(.GlobalEnv)
# Save list
save(list_overlap_int_not_empty_xu_objects, file=file.path("~/riboRNA/R_files/", "list_overlap_int_not_empty_xu_objects.rda"))
##

## E5. LIST OF ALL OVERLAP OBJECTS
# Remove all objects in environment
rm(list = ls())
# Load objects
lapply(list.files(pattern="^overlap.*rda$"), load, .GlobalEnv)
# Generate list from global environment
list_overlap_all_objects <- as.list.environment(.GlobalEnv)
# Save list
save(list_overlap_all_objects, file=file.path("~/riboRNA/R_files/", "list_overlap_all_objects.rda"))
##

## E6. IMAGE FILE
# Load objects
lapply(list.files(pattern="^list_overlap.*rda$"), load, .GlobalEnv)
# Save image
save.image(file = "image_overlap_objects.Rdata")
##

## E7. CLEAN-UP
# Remove all objects in environment
rm(list = ls())
# Remove unused packages
detach(package:Rsamtools)
detach(package:GenomicRanges)
detach(package:Biostrings)
detach(package:IRanges)
detach(package:BiocGenerics)
###
