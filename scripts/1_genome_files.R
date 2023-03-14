#######
### A. PRE-REQUISITES
### B. GENOME IMPORT
### C. GENE LEVEL ANNOTATION
### D. GENOMIC LEVEL ANNOTATION
### E. GROUP OBJECTS & CLEAN-UP
#######

### A. PRE-REQUISITES
# Remove all objects
rm(list = ls())
# Load libraries
library(rtracklayer)
library("BSgenome.Scerevisiae.UCSC.sacCer3")
###

### B. GENOME IMPORT
# Download yeast genome GTF file from Ensembl
download.file("ftp://ftp.ensembl.org/pub/release-68/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.EF4.68.gtf.gz", destfile="~/riboRNA/genome_files/Saccharomyces_cerevisiae.EF4.68.gtf.gz")
# Use rtracklayer::import method to import yeast .gtf genome file to IRanges object and 
sc68_gtf <- import(gzfile("~/riboRNA/genome_files/Saccharomyces_cerevisiae.EF4.68.gtf.gz"), "gtf")
# ALTERNATIVELY: Unzip file first, then import decompressed file
# > system("gunzip ~/riboRNA/genome_files/Saccharomyces_cerevisiae.EF4.68.gtf.gz")
# > sc68_gtf <- import("~/riboRNA/genome_files/Saccharomyces_cerevisiae.EF4.68.gtf", "gtf")
# Detach rtracklayer package
detach(package:rtracklayer)
# Use methods::as method to cast IRanges object to GRanges object and sort
sc68_features <- sort(as(sc68_gtf, "GRanges"))
# Use base::gsub method to substitute chromosome labels to match those in the read files
seqlevels(sc68_features) <- gsub('Mito','M',seqlevels(sc68_features))
seqlevels(sc68_features) <- gsub('^(.*)$','chr\\1', seqlevels(sc68_features))
# Set chromosome lengths (from data present in BSgenome.Scerevisiae.UCSC.sacCer3 package)
seqlengths(sc68_features) <- seqlengths(Scerevisiae)
# Save substituted yeast genome features GRanges object 
save(sc68_features, file=file.path("~/riboRNA/R_files/", "sc68_features.rda"))
# Detach unused packages
detach(package:BSgenome.Scerevisiae.UCSC.sacCer3)
detach(package:BSgenome)
detach(package:Biostrings)
###

### C. GENE LEVEL ANNOTATION
## C1. EXTRACT EXONS
sc68_exons <- sc68_features[values(sc68_features)[["type"]] == "exon"]
save(sc68_exons, file=file.path("~/riboRNA/R_files/", "sc68_exons.rda"))
# Remove unused objects
rm(sc68_gtf, sc68_features)
##
## C2. INFER GENE REGIONS FROM EXONS
# Split "sc68_exons" to GenomicRangesList object via "gene_id" identifier
# Each list entry corresponds to a single gene and contains a GRangesList of its exons 
sc68_exons_list <- split(sc68_exons, mcols(sc68_exons)$gene_id)
# Use GenomicRanges::endoapply method to merge exons and introns for each gene (~7 min)
# Returns
 "sc68_genes_list" GRangesList object of genes (i.e. exons and introns merged)
# Algorithm:
# For each list entry
# 1. Assign first element (i.e. exon) to GRanges object "merged"
# 2. If multiple elements (i.e. exons) are present (if statement), compute min/max start/end values of all elements/exons and define as start/end of "merged"
# 3. Return "merged"
sc68_genes_list <- endoapply(sc68_exons_list, function(u) {
  merged <- u[1]
  if(length(u) > 1) {
    start(merged) <- min(start(u))
    end(merged) <- max(end(u))
  }
  return(merged)
})
# Unlist GRangesList to retrieve GRanges object of genes
sc68_genes <- unlist(sc68_genes_list)
# Save GRanges object "sc68_genes"
save(sc68_genes, file=file.path("~/riboRNA/R_files/", "sc68_genes.rda"))
##
## C3. INFER INTRONS FROM GENES AND EXONS
# Use GenomicRanges::psetdiff method to compare gene and exon lists (one list entry per gene, see above)
# Returns "sc68_introns_list" GRangesList object of introns (i.e. differences between genes and exons)
sc68_introns_list <- psetdiff(sc68_genes_list,sc68_exons_list)
# Obtain metadata from exons list (~17 min)
sc68_introns_list <- mendoapply(function(u,v) {
  if (length(v) > 0) {
    v$source <- u[1]$source
    v$type <- "intron"
    is.na(mcols(v)) <- c("score","phase")
    intron_no <- u$exon_number[1:(length(u$exon_number) - 1)]
    v$exon_number <- paste0(intron_no, "/", intron_no + 1)
    v$gene_biotype <- u[1]$gene_biotype
    v$gene_id <- u[1]$gene_id
    v$gene_name <- u[1]$gene_name
    is.na(mcols(v)) <- "protein_id"
    v$transcript_id <- u[1]$transcript_id
    v$transcript_name <- u[1]$transcript_name
  }
  return(v)
}, sc68_exons_list, sc68_introns_list)
# Unlist GRangesList to retrieve GRanges object of introns
sc68_introns <- unlist(sc68_introns_list)
# Save GRanges object "sc68_introns"
save(sc68_introns, file=file.path("~/riboRNA/R_files/", "sc68_introns.rda"))
# Remove GRangesList objects
rm(sc68_exons_list, sc68_genes_list, sc68_introns_list)
##
## C4. DEFINE ANTISENSE REGIONS OF GENES, EXONS & INTRONS
# Use base::ifelse method to swap strands
# GENES
sc68_genes_as <- sc68_genes
strand(sc68_genes_as)<- ifelse(strand(sc68_genes_as) == '+', '-', '+')
save(sc68_genes_as, file=file.path("~/riboRNA/R_files/", "sc68_genes_as.rda"))
# EXONS
sc68_exons_as <- sc68_exons
strand(sc68_exons_as)<- ifelse(strand(sc68_exons_as) == '+', '-', '+')
save(sc68_exons_as, file=file.path("~/riboRNA/R_files/", "sc68_exons_as.rda"))
# INTRONS
sc68_introns_as <- sc68_introns
strand(sc68_introns_as)<- ifelse(strand(sc68_introns_as) == '+', '-', '+')
save(sc68_introns_as, file=file.path("~/riboRNA/R_files/", "sc68_introns_as.rda"))
##
## C5. SUBSETS OF EXON TYPES
# Extract and save subsets of sequence gene biotypes from yeast genome exons GRanges object
# PROTEIN-CODING
sc68_type_protein_coding <- sc68_exons[mcols(sc68_exons)[["gene_biotype"]] == "protein_coding"]
save(sc68_type_protein_coding, file=file.path("~/riboRNA/R_files/", "sc68_type_protein_coding.rda"))
# PSEUDOGENES
sc68_type_pseudo <- sc68_exons[mcols(sc68_exons)[["gene_biotype"]] == "pseudogene"]
save(sc68_type_pseudo, file=file.path("~/riboRNA/R_files/", "sc68_type_pseudo.rda"))
# NON-CODING RNA
sc68_type_ncRNA <- sc68_exons[mcols(sc68_exons)[["gene_biotype"]] == "ncRNA"]
save(sc68_type_ncRNA, file=file.path("~/riboRNA/R_files/", "sc68_type_ncRNA.rda"))
# RIBSOMAL RNA
sc68_type_rRNA <- sc68_exons[mcols(sc68_exons)[["gene_biotype"]] == "rRNA"]
save(sc68_type_rRNA, file=file.path("~/riboRNA/R_files/", "sc68_type_rRNA.rda"))
# SMALL NUCLEOLAR RNA
sc68_type_snoRNA <- sc68_exons[mcols(sc68_exons)[["gene_biotype"]] == "snoRNA"]
save(sc68_type_snoRNA, file=file.path("~/riboRNA/R_files/", "sc68_type_snoRNA.rda"))
# SMULL NUCLEAR RNA
sc68_type_snRNA <- sc68_exons[mcols(sc68_exons)[["gene_biotype"]] == "snRNA"]
save(sc68_type_snRNA, file=file.path("~/riboRNA/R_files/", "sc68_type_snRNA.rda"))
# TRANSFER RNA
sc68_type_tRNA <- sc68_exons[mcols(sc68_exons)[["gene_biotype"]] == "tRNA"]
save(sc68_type_tRNA, file=file.path("~/riboRNA/R_files/", "sc68_type_tRNA.rda"))
###

### D. GENOMIC LEVEL
## D1. DEFINE AMBIGUOUS REGIONS
# Use IRanges::subsetByOverlaps method to define ambiguous regions based on overlaps between different annotated regions
sc68_ambig_exons_introns <- subsetByOverlaps(sc68_exons, sc68_introns)
sc68_ambig_exons_exons_as <- subsetByOverlaps(sc68_exons, sc68_exons_as)
sc68_ambig_exons_introns_as <- subsetByOverlaps(sc68_exons, sc68_introns_as)
sc68_ambig_introns_exons_as <- subsetByOverlaps(sc68_introns, sc68_exons_as)
sc68_ambig_introns_introns_as <- subsetByOverlaps(sc68_introns, sc68_introns_as)
sc68_ambig_exons_as_introns_as <- subsetByOverlaps(sc68_exons_as, sc68_introns_as)
# Save ambiguous regions (improve metadata annotation?!)
save(sc68_ambig_exons_introns, file=file.path("~/riboRNA/R_files/", "sc68_ambig_exons_introns.rda"))
save(sc68_ambig_exons_exons_as, file=file.path("~/riboRNA/R_files/", "sc68_ambig_exons_exons_as.rda"))
save(sc68_ambig_exons_introns_as, file=file.path("~/riboRNA/R_files/", "sc68_ambig_exons_introns_as.rda"))
save(sc68_ambig_introns_exons_as, file=file.path("~/riboRNA/R_files/", "sc68_ambig_introns_exons_as.rda"))
save(sc68_ambig_introns_introns_as, file=file.path("~/riboRNA/R_files/", "sc68_ambig_introns_introns_as.rda"))
save(sc68_ambig_exons_as_introns_as, file=file.path("~/riboRNA/R_files/", "sc68_ambig_exons_as_introns_as.rda"))
# Combine different ambigiuous regions
sc68_ambig <- c(sc68_ambig_exons_introns, sc68_ambig_exons_exons_as, sc68_ambig_exons_introns_as, sc68_ambig_introns_exons_as, sc68_ambig_introns_introns_as, sc68_ambig_exons_as_introns_as, ignore.mcols=TRUE)
# Remove overlapping regions with the GenomicRanges::disjoin method; save ambiguous regions object
sc68_ambig <- disjoin(sc68_ambig)
save(sc68_ambig, file=file.path("~/riboRNA/R_files/", "sc68_ambig.rda"))
##
## D2. DEFINE UNAMBIGUOUS REGIONS
# Define exclusive regions for exons, introns and their antisense strands
sc68_excl_exonic <- setdiff(sc68_exons, sc68_ambig)
sc68_excl_intronic <- setdiff(sc68_introns, sc68_ambig)
sc68_excl_opp_exonic <- setdiff(sc68_exons_as, sc68_ambig)
sc68_excl_opp_intronic <- setdiff(sc68_introns_as, sc68_ambig)
# Use GenomicRanges::gaps method to define intergenic regions; discard star (*) strand
sc68_excl_intergenic <- gaps(c(sc68_excl_exonic, sc68_excl_intronic, sc68_excl_opp_exonic, sc68_excl_opp_intronic, sc68_ambig))
sc68_excl_intergenic <- sc68_excl_intergenic[strand(sc68_excl_intergenic) != "*"]
# Save regions
save(sc68_excl_exonic, file=file.path("~/riboRNA/R_files/", "sc68_excl_exonic.rda"))
save(sc68_excl_intronic, file=file.path("~/riboRNA/R_files/", "sc68_excl_intronic.rda"))
save(sc68_excl_opp_exonic, file=file.path("~/riboRNA/R_files/", "sc68_excl_opp_exonic.rda"))
save(sc68_excl_opp_intronic, file=file.path("~/riboRNA/R_files/", "sc68_excl_opp_intronic.rda"))
save(sc68_excl_intergenic, file=file.path("~/riboRNA/R_files/", "sc68_excl_intergenic.rda"))
###

### E. GROUP OBJECTS & CLEAN-UP
# List of all objects
list_sc68_all_objects <- as.list.environment(.GlobalEnv)
save(list_sc68_all_objects, file=file.path("~/riboRNA/R_files/", "list_sc68_all_objects.rda"))
# GRangesList of all mutually exclusive genomic regions objects
list_sc68_genomic_type <- GRangesList()
list_sc68_genomic_type$exonic <- sc68_excl_exonic
list_sc68_genomic_type$intronic <- sc68_excl_intronic
list_sc68_genomic_type$exonic_as <- sc68_excl_opp_exonic
list_sc68_genomic_type$intronic_as <- sc68_excl_opp_intronic
list_sc68_genomic_type$intergenic <- sc68_excl_intergenic
list_sc68_genomic_type$ambiguous <- sc68_ambig
save(list_sc68_genomic_type, file=file.path("~/riboRNA/R_files/", "list_sc68_genomic_type.rda"))
# GRangesList of all gene types objects (mutually exclusive)
list_sc68_gene_type <- GRangesList() 
list_sc68_gene_type$protein_coding <- sc68_type_protein_coding
list_sc68_gene_type$pseudo <- sc68_type_pseudo
list_sc68_gene_type$ncRNA <- sc68_type_ncRNA
list_sc68_gene_type$rRNA <- sc68_type_rRNA
list_sc68_gene_type$snoRNA <- sc68_type_snoRNA
list_sc68_gene_type$snRNA <- sc68_type_snRNA
list_sc68_gene_type$tRNA <- sc68_type_tRNA
save(list_sc68_gene_type, file=file.path("~/riboRNA/R_files/", "list_sc68_gene_type.rda"))
# GRangesList of other types (genes, exons, introns and opposite strands; NOT mutually exclusive!)
list_sc68_other <- GRangesList()
list_sc68_other$genes <- sc68_genes
list_sc68_other$genes_as <- sc68_genes_as
list_sc68_other$exons <- sc68_exons
list_sc68_other$exons_as <- sc68_exons_as
list_sc68_other$introns <- sc68_introns
list_sc68_other$introns_as <- sc68_introns_as
save(list_sc68_other, file=file.path("~/riboRNA/R_files/", "list_sc68_other.rda"))
# Image file of all objects
save.image(file = "image_sc68_objects.Rdata")
# Remove all objects
rm(list = ls())
# Remove unused packages
detach(package:GenomicRanges)
detach(package:IRanges)
detach(package:BiocGenerics)
###
