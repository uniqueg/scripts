#######
### A. PRE-REQUISITES
### B. ANALYSIS GENOME ANNOTATION
### C. ANALYSIS READS
### D. COMBINE & CLEAN-UP
### E. GRAPHICAL OUTPUT
#######


### A. PRE-REQUISITES
# Remove all objects
rm(list = ls())
# Load libraries
library(GenomicRanges)
###


### B. ANALYSIS GENOME ANNOTATION
# List the fractions of different genome features

## B1. GENOMIC TYPES
# Load object list
load("list_sc68_genomic_type.rda")
# Use sapply to extract total interval width for each feature
genomic_type_annotation <- sapply(list_sc68_genomic_type, function(u) sum(width(u)))
# Coerce to matrix
genomic_type_annotation <- as.matrix(genomic_type_annotation)
# Assign column name
colnames(genomic_type_annotation) <- "total_nt"
##

## B2. GENE TYPES
# Load object list
load("list_sc68_gene_type.rda")
# Use sapply to extract total interval widths
gene_type_annotation <- sapply(list_sc68_gene_type, function(u) sum(width(u)))
# Coerce to matrix
gene_type_annotation <- as.matrix(gene_type_annotation)
# Assign column name
colnames(gene_type_annotation) <- "total_nt"
##

## B3. REMAINING (NOT MUTUALLY EXCLUSIVE)  SC68 FEATURES
# Load object list
load("list_sc68_other.rda")
# Use sapply to extract total interval widths
other_annotation <- sapply(list_sc68_other, function(u) sum(width(u)))
# Coerce to matrix
other_annotation <- as.matrix(other_annotation)
# Assign column name
colnames(other_annotation) <- "total_nt"
##

## B4. CUTS/SUTS
# Load object list
load("list_xu_objects.rda")
# Use sapply to extract total interval widths
xu_annotation <- sapply(list_xu_objects, function(u) sum(width(u)))
# Coerce to matrix
xu_annotation <- as.matrix(xu_annotation)
# Assign column name
colnames(xu_annotation) <- "total_nt"
##

## B5. COMBINE MATRICES, CALCULATE FREQUENCIES & CLEAN-UP
# Combine matrizes
feature_annotation <- rbind(genomic_type_annotation, gene_type_annotation, other_annotation, xu_annotation)
# Exctract total genome size
ds_genome_size <- 2 * sum(seqlengths(list_sc68_other$exons))
# Compute frequencies and add to matrix
feature_annotation <- cbind(feature_annotation, feature_annotation/ds_genome_size*100)
# Add column name
colnames(feature_annotation)[2] <- "total_nt_%"
# Remove unused objects
rm(list= ls()[!(ls() %in% 'feature_annotation')])
###


### C. ANALYSIS READS
# List the read counts within each region type

## C1. GENOMIC TYPES
# Load object list
load("list_overlap_genomic_type_objects.rda")
# Use lapply to obtain the total read counts for each genomic type in each sample
genomic_type_reads <- lapply(list_overlap_genomic_type_objects, function(u) assays(u)$counts)
# Coerce to matrix
genomic_type_reads <- do.call(cbind, genomic_type_reads)
# Sort rows alphabetically
genomic_type_reads <- genomic_type_reads[order(rownames(genomic_type_reads)), ]
##

## C2. GENE TYPES
# Load object list
load("list_overlap_gene_type_objects.rda")
# Use lapply to obtain the total read counts for each genomic type in each sample
gene_type_reads <- lapply(list_overlap_gene_type_objects, function(u) assays(u)$counts)
# Coerce to matrix
gene_type_reads <- do.call(cbind, gene_type_reads)
# Sort rows alphabetically
gene_type_reads <- gene_type_reads[order(rownames(gene_type_reads)), ]
##

## C3. REMAINING (NOT MUTUALLY EXCLUSIVE) SC68 FEATURES

# UNION
# Load object list
load("list_overlap_union_other_objects.rda")
# Use lapply to !!!!
other_union_reads <- lapply(seq_along(list_overlap_union_other_objects), function(u, n, i) {
  sums <- matrix(ncol=3)
  sums[,1] <- sum(assays(u[[i]][,1])$counts)
  sums[,2] <- sum(assays(u[[i]][,2])$counts)
  sums[,3] <- sum(assays(u[[i]][,3])$counts)
  colnames(sums) <- colnames(u[[i]])
  rownames(sums) <- n[[i]]
  return(sums)
}, u = list_overlap_union_other_objects, n = gsub('overlap_union_reads_', '', names(list_overlap_union_other_objects)))
# Combine list elements to matrix
other_union_reads <- do.call(rbind, other_union_reads)
# Sort rows alphabetically
other_union_reads <- other_union_reads[order(rownames(other_union_reads)), ]

# INT STRICT
# Load object list
load("list_overlap_int_strict_other_objects.rda")
# Use lapply to !!!!
other_int_strict_reads <- lapply(seq_along(list_overlap_int_strict_other_objects), function(u, n, i) {
  sums <- matrix(ncol=3)
  sums[,1] <- sum(assays(u[[i]][,1])$counts)
  sums[,2] <- sum(assays(u[[i]][,2])$counts)
  sums[,3] <- sum(assays(u[[i]][,3])$counts)
  colnames(sums) <- colnames(u[[i]])
  rownames(sums) <- n[[i]]
  return(sums)
}, u = list_overlap_int_strict_other_objects, n = gsub('overlap_int_strict_reads_', '', names(list_overlap_int_strict_other_objects)))
# Combine list elements to matrix
other_int_strict_reads <- do.call(rbind, other_int_strict_reads)
# Sort rows alphabetically
other_int_strict_reads <- other_int_strict_reads[order(rownames(other_int_strict_reads)), ]

# INT NOT EMPTY
# Load object list
load("list_overlap_int_not_empty_other_objects.rda")
# Use lapply to !!!!
other_int_not_empty_reads <- lapply(seq_along(list_overlap_int_not_empty_other_objects), function(u, n, i) {
  sums <- matrix(ncol=3)
  sums[,1] <- sum(assays(u[[i]][,1])$counts)
  sums[,2] <- sum(assays(u[[i]][,2])$counts)
  sums[,3] <- sum(assays(u[[i]][,3])$counts)
  colnames(sums) <- colnames(u[[i]])
  rownames(sums) <- n[[i]]
  return(sums)
}, u = list_overlap_int_not_empty_other_objects, n = gsub('overlap_int_not_empty_reads_', '', names(list_overlap_int_not_empty_other_objects)))
# Combine list elements to matrix
other_int_not_empty_reads <- do.call(rbind, other_int_not_empty_reads)
# Sort rows alphabetically
other_int_not_empty_reads <- other_int_not_empty_reads[order(rownames(other_int_not_empty_reads)), ]

# Combine matrices
other_reads <- cbind(other_union_reads, other_int_strict_reads, other_int_not_empty_reads)
##

## C4. CUTS/SUTS

# UNION
# Load object list
load("list_overlap_union_xu_objects.rda")
# Use lapply to !!!!
xu_union_reads <- lapply(seq_along(list_overlap_union_xu_objects), function(u, n, i) {
  sums <- matrix(ncol=3)
  sums[,1] <- sum(assays(u[[i]][,1])$counts)
  sums[,2] <- sum(assays(u[[i]][,2])$counts)
  sums[,3] <- sum(assays(u[[i]][,3])$counts)
  colnames(sums) <- colnames(u[[i]])
  rownames(sums) <- n[[i]]
  return(sums)
}, u = list_overlap_union_xu_objects, n = gsub('overlap_union_reads_', '', names(list_overlap_union_xu_objects)))
# Combine list elements to matrix
xu_union_reads <- do.call(rbind, xu_union_reads)
# Sort rows alphabetically
xu_union_reads <- xu_union_reads[order(rownames(xu_union_reads)), ]

# INT STRICT
# Load object list
load("list_overlap_int_strict_xu_objects.rda")
# Use lapply to !!!!
xu_int_strict_reads <- lapply(seq_along(list_overlap_int_strict_xu_objects), function(u, n, i) {
  sums <- matrix(ncol=3)
  sums[,1] <- sum(assays(u[[i]][,1])$counts)
  sums[,2] <- sum(assays(u[[i]][,2])$counts)
  sums[,3] <- sum(assays(u[[i]][,3])$counts)
  colnames(sums) <- colnames(u[[i]])
  rownames(sums) <- n[[i]]
  return(sums)
}, u = list_overlap_int_strict_xu_objects, n = gsub('overlap_int_strict_reads_', '', names(list_overlap_int_strict_xu_objects)))
# Combine list elements to matrix
xu_int_strict_reads <- do.call(rbind, xu_int_strict_reads)
# Sort rows alphabetically
xu_int_strict_reads <- xu_int_strict_reads[order(rownames(xu_int_strict_reads)), ]

# INT NOT EMPTY
# Load object list
load("list_overlap_int_not_empty_xu_objects.rda")
# Use lapply to !!!!
xu_int_not_empty_reads <- lapply(seq_along(list_overlap_int_not_empty_xu_objects), function(u, n, i) {
  sums <- matrix(ncol=3)
  sums[,1] <- sum(assays(u[[i]][,1])$counts)
  sums[,2] <- sum(assays(u[[i]][,2])$counts)
  sums[,3] <- sum(assays(u[[i]][,3])$counts)
  colnames(sums) <- colnames(u[[i]])
  rownames(sums) <- n[[i]]
  return(sums)
}, u = list_overlap_int_not_empty_xu_objects, n = gsub('overlap_int_not_empty_reads_', '', names(list_overlap_int_not_empty_xu_objects)))
# Combine list elements to matrix
xu_int_not_empty_reads <- do.call(rbind, xu_int_not_empty_reads)
# Sort rows alphabetically
xu_int_not_empty_reads <- xu_int_not_empty_reads[order(rownames(xu_int_not_empty_reads)), ]

# Combine matrices
xu_reads <- cbind(xu_union_reads, xu_int_strict_reads, xu_int_not_empty_reads)
##

## C5. COMBINE MATRICES & CLEAN-UP
# Combine matrices
feature_reads <- rbind(genomic_type_reads, gene_type_reads, other_reads, xu_reads)
# Load read objects
load("reads_A.rda")
load("reads_A2.rda")
load("reads_B.rda")
# Extract total read numbers
read_no_A <- length(reads_A)
read_no_A2 <- length(reads_A2)
read_no_B <- length(reads_B)
# Sort matrix columns alphabetically
feature_reads <- feature_reads[, order(colnames(feature_reads))]
# Compute matrix of frequencies
feature_reads_pct <- cbind(feature_reads[,c(1,4,7)]/read_no_A*100, feature_reads[,c(2,5,8)]/read_no_A2*100, feature_reads[,c(3,6,9)]/read_no_B*100)
# Set column names
colnames(feature_reads_pct) <- gsub('(.*)', '\\1_%', colnames(feature_reads_pct))
# Combine matrices
feature_reads <- cbind(feature_reads, feature_reads_pct)
# Sort matrix columns alphabetically
feature_reads <- feature_reads[, order(colnames(feature_reads))]
###


### D. COMBINE & CLEAN-UP
# Reorder read overlaps matrix according to rownames in annotation matrix by subsetting
feature_reads <- feature_reads[rownames(feature_annotation), , drop=FALSE]
# Combine matrices
mt_genome_analysis <- cbind(feature_annotation, feature_reads)
# Save resulting matrix
save(mt_genome_analysis, file=file.path("/home/akanitz/riboRNA/R_files/", "mt_genome_analysis.rda"))
# Remove unused objects
rm(list= ls()[!(ls() %in% 'mt_genome_analysis')])
###


### E. GRAPHICAL OUTPUT

## E1. PIE CHARTS GENOMIC TYPE
# Genome
pdf("pie_genomic_type_genome.pdf", width = 6, height = 6)
pie(mt_genome_analysis[1:6, 1], paste(gsub('_as', ' antisense', rownames(mt_genome_analysis)[1:6]), mt_genome_analysis[1:6, 1], sep = "\n"), main="S. cerevisiae genome (sc68)", clockwise = TRUE)
dev.off()
# Sample A, Mode: Union
pdf("pie_genomic_type_A_union.pdf", width = 6, height = 6)
pie(mt_genome_analysis[1:6, 15], paste(gsub('_as', ' antisense', rownames(mt_genome_analysis)[1:6]), mt_genome_analysis[1:6, 15], sep = "\n"), main="Translatome wt", sub = "Overlap mode: Union", clockwise = TRUE)
dev.off()
# Sample A2, Mode: Union
pdf("pie_genomic_type_A2_union.pdf", width = 6, height = 6)
pie(mt_genome_analysis[1:6, 17], paste(gsub('_as', ' antisense', rownames(mt_genome_analysis)[1:6]), mt_genome_analysis[1:6, 17], sep = "\n"), main="Transcriptome wt", sub = "Overlap mode: Union", clockwise = TRUE)
dev.off()
# Sample B, Mode: Union
pdf("pie_genomic_type_B_union.pdf", width = 6, height = 6)
pie(mt_genome_analysis[1:6, 19], paste(gsub('_as', ' antisense', rownames(mt_genome_analysis)[1:6]), mt_genome_analysis[1:6, 19], sep = "\n"), main="Translatome DTT", sub = "Overlap mode: Union", clockwise = TRUE)
dev.off()
# Sample A, Mode: Intersection Strict
pdf("pie_genomic_type_A_int_strict.pdf", width = 6, height = 6)
pie(mt_genome_analysis[1:6, 9], paste(gsub('_as', ' antisense', rownames(mt_genome_analysis)[1:6]), mt_genome_analysis[1:6, 9], sep = "\n"), main="Translatome wt", sub = "Overlap mode: Intersection Strict", clockwise = TRUE)
dev.off()
# Sample A2, Mode: Intersection Strict
pdf("pie_genomic_type_A2_int_strict.pdf", width = 6, height = 6)
pie(mt_genome_analysis[1:6, 11], paste(gsub('_as', ' antisense', rownames(mt_genome_analysis)[1:6]), mt_genome_analysis[1:6, 11], sep = "\n"), main="Transcriptome wt", sub = "Overlap mode: Intersection Strict", clockwise = TRUE)
dev.off()
# Sample B, Mode: Intersection Strict
pdf("pie_genomic_type_B_int_strict.pdf", width = 6, height = 6)
pie(mt_genome_analysis[1:6, 13], paste(gsub('_as', ' antisense', rownames(mt_genome_analysis)[1:6]), mt_genome_analysis[1:6, 13], sep = "\n"), main="Translatome DTT", sub = "Overlap mode: Intersection Strict", clockwise = TRUE)
dev.off()
# Sample A, Mode: Intersection Not Empty
pdf("pie_genomic_type_A_int_not_empty.pdf", width = 6, height = 6)
pie(mt_genome_analysis[1:6, 3], paste(gsub('_as', ' antisense', rownames(mt_genome_analysis)[1:6]), mt_genome_analysis[1:6, 3], sep = "\n"), main="Translatome wt", sub = "Overlap mode: Intersection Not Empty", clockwise = TRUE)
dev.off()
# Sample A2, Mode: Intersection Not Empty
pdf("pie_genomic_type_A2_int_not_empty.pdf", width = 6, height = 6)
pie(mt_genome_analysis[1:6, 5], paste(gsub('_as', ' antisense', rownames(mt_genome_analysis)[1:6]), mt_genome_analysis[1:6, 5], sep = "\n"), main="Transcriptome wt", sub = "Overlap mode: Intersection Not Empty", clockwise = TRUE)
dev.off()
# Sample B, Mode: Intersection Not Empty
pdf("pie_genomic_type_B_int_not_empty.pdf", width = 6, height = 6)
pie(mt_genome_analysis[1:6, 7], paste(gsub('_as', ' antisense', rownames(mt_genome_analysis)[1:6]), mt_genome_analysis[1:6, 7], sep = "\n"), main="Translatome DTT", sub = "Overlap mode: Intersection Not Empty", clockwise = TRUE)
dev.off()
##

## E2. PIE CHARTS GENE TYPE
# Genome
pdf("pie_gene_type_genome.pdf", width = 6, height = 6)
pie(mt_genome_analysis[7:13,1], paste(gsub('_as', ' antisense', rownames(mt_genome_analysis)[7:13]), mt_genome_analysis[7:13,1], sep = "\n"), main="S. cerevisiae genome (sc68)", clockwise = TRUE)
dev.off()
# Sample A, Mode: Union
pdf("pie_gene_type_A_union.pdf", width = 6, height = 6)
pie(mt_genome_analysis[7:13,15], paste(gsub('_as', ' antisense', rownames(mt_genome_analysis)[7:13]), mt_genome_analysis[7:13,15], sep = "\n"), main="Translatome wt", sub = "Overlap mode: Union", clockwise = TRUE)
dev.off()
# Sample A2, Mode: Union
pdf("pie_gene_type_A2_union.pdf", width = 6, height = 6)
pie(mt_genome_analysis[7:13,17], paste(gsub('_as', ' antisense', rownames(mt_genome_analysis)[7:13]), mt_genome_analysis[7:13,17], sep = "\n"), main="Transcriptome wt", sub = "Overlap mode: Union", clockwise = TRUE)
dev.off()
# Sample B, Mode: Union
pdf("pie_gene_type_B_union.pdf", width = 6, height = 6)
pie(mt_genome_analysis[7:13,19], paste(gsub('_as', ' antisense', rownames(mt_genome_analysis)[7:13]), mt_genome_analysis[7:13,19], sep = "\n"), main="Translatome DTT", sub = "Overlap mode: Union", clockwise = TRUE)
dev.off()
# Sample A, Mode: Intersection Strict
pdf("pie_gene_type_A_int_strict.pdf", width = 6, height = 6)
pie(mt_genome_analysis[7:13,9], paste(gsub('_as', ' antisense', rownames(mt_genome_analysis)[7:13]), mt_genome_analysis[7:13,9], sep = "\n"), main="Translatome wt", sub = "Overlap mode: Intersection Strict", clockwise = TRUE)
dev.off()
# Sample A2, Mode: Intersection Strict
pdf("pie_gene_type_A2_int_strict.pdf", width = 6, height = 6)
pie(mt_genome_analysis[7:13,11], paste(gsub('_as', ' antisense', rownames(mt_genome_analysis)[7:13]), mt_genome_analysis[7:13,11], sep = "\n"), main="Transcriptome wt", sub = "Overlap mode: Intersection Strict", clockwise = TRUE)
dev.off()
# Sample B, Mode: Intersection Strict
pdf("pie_gene_type_B_int_strict.pdf", width = 6, height = 6)
pie(mt_genome_analysis[7:13,13], paste(gsub('_as', ' antisense', rownames(mt_genome_analysis)[7:13]), mt_genome_analysis[7:13,13], sep = "\n"), main="Translatome DTT", sub = "Overlap mode: Intersection Strict", clockwise = TRUE)
dev.off()
# Sample A, Mode: Intersection Not Empty
pdf("pie_gene_type_A_int_not_empty.pdf", width = 6, height = 6)
pie(mt_genome_analysis[7:13,3], paste(gsub('_as', ' antisense', rownames(mt_genome_analysis)[7:13]), mt_genome_analysis[7:13,3], sep = "\n"), main="Translatome wt", sub = "Overlap mode: Intersection Not Empty", clockwise = TRUE)
dev.off()
# Sample A2, Mode: Intersection Not Empty
pdf("pie_gene_type_A2_int_not_empty.pdf", width = 6, height = 6)
pie(mt_genome_analysis[7:13,5], paste(gsub('_as', ' antisense', rownames(mt_genome_analysis)[7:13]), mt_genome_analysis[7:13,5], sep = "\n"), main="Transcriptome wt", sub = "Overlap mode: Intersection Not Empty", clockwise = TRUE)
dev.off()
# Sample B, Mode: Intersection Not Empty
pdf("pie_gene_type_B_int_not_empty.pdf", width = 6, height = 6)
pie(mt_genome_analysis[7:13,7], paste(gsub('_as', ' antisense', rownames(mt_genome_analysis)[7:13]), mt_genome_analysis[7:13,7], sep = "\n"), main="Translatome DTT", sub = "Overlap mode: Intersection Not Empty", clockwise = TRUE)
dev.off()
##

## E3. STACKED BAR CHART GENOMIC TYPE
dt <- mt_genome_analysis[1:6, c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20)]
# Infer unmapped reads, label and add to matrix
sums <- t(as.matrix(colSums(dt)))
unmapped <- 100 - sums
rownames(unmapped) <- "unmapped/overlapping"
dt_plot <- rbind(dt, unmapped)
# Labels and legend
lbls <- c("Genome", "Transl wt\nIntNotEmpty", "Transcr wt\nIntNotEmpty", "Transl DTT\nIntNotEmpty", "Transl wt\nIntStrict", "Transcr wt\nIntStrict", "Transl DTT\nIntStrict", "Transl wt\nUnion", "Transcr wt\nUnion", "Transl DTT\nUnion")
lgnd <- gsub('_as$', ' antisense', rownames(dt_plot))
# Plot
pdf("stacked_bar_genomic_type.pdf", width = 15, height = 6)
barplot(dt_plot, xlim = c(4, 100), ylim = c(0, 100), width = 6, space = 0.4, main = "S. cerevisiae genome (sc68)", names.arg = lbls, col = c(1:6, 8), cex.names=0.75, legend = lgnd)
dev.off()
##

## E4. STACKED BAR CHARTS GENE TYPE

# PROTEIN CODING, RIBOSOMAL RNA, UNANNOTATED & REST
# Subset data
dt <- mt_genome_analysis[7:13, c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20)]
# Calculate "rest" and unannotated
other <- t(as.matrix(colSums(dt[c(2, 3, 5:7), ])))
dt_plot <- rbind(dt[c(1, 4), ], other)
unannotated <- 100 - t(as.matrix(colSums(dt_plot)))
dt_plot <- rbind(dt_plot, unannotated)
rownames(dt_plot) <- c("protein coding", "rRNA", "other", "unannotated / no gene")
# Legend (labels as above; E3.)
lgnd <- rownames(dt_plot)
# Plot and print .pdf
pdf("stacked_bar_gene_type_all.pdf", width = 15, height = 6)
barplot(dt_plot, xlim = c(4, 100), ylim = c(0, 100), width = 6, space = 0.4, main = "S. cerevisiae genome (sc68)", names.arg = lbls, col = c(1:3, 8), cex.names = 0.75, legend = lgnd)
dev.off()

# "REST": ALL BUT PROTEIN CODING, RIBOSOMAL RNA & UNANNOTATED
# Subset data
dt <- mt_genome_analysis[c(8, 9, 11:13), c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20)]
# Normalize to 100
dt_plot <- apply(dt, 2, function(u) u / sum(u) * 100)
# Legend (labels as above; E3.)
lgnd <- rownames(dt_plot)
# Plot and print .pdf
pdf("stacked_bar_gene_type_rest.pdf", width = 15, height = 6)
barplot(dt_plot, xlim = c(4, 100), ylim = c(0, 100), width = 6, space = 0.4, main = "S. cerevisiae genome (sc68)", names.arg = lbls, col = c(1:6), cex.names = 0.75, legend = lgnd)
dev.off()
##

## E5. STACKED BAR CHART CUTS/SUTS
# Subset data
dt <- mt_genome_analysis[c(20, 21), c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20)]
# Normalize to 100
dt_plot <- dt
# Legend (labels as above; E3.)
lgnd <- rownames(dt_plot)
# Plot and print .pdf
pdf("stacked_bar_cuts_suts.pdf", width = 15, height = 6)
barplot(dt_plot, xlim = c(4, 100), ylim = c(0, ceiling(max(colSums(dt_plot)))), width = 6, space = 0.4, main = "S. cerevisiae genome (sc68)", names.arg = lbls, col = c(1:6), cex.names = 0.75, legend = lgnd)
dev.off()
##

## E6. RELATIVE CHANGES/ENRICHMENTS

# TRANSCRIPTOME WT VS GENOME
# Subset data
genome_distr <- mt_genome_analysis[c(1:9,11:22), 2]
transcr_wt <- mt_genome_analysis[c(1:9,11:22), c(6, 12, 18)]
# Calculate ratios
dt_plot <- t(transcr_wt / genome_distr)
# Labels and legend
# lgnd <- gsub('(u|i)(nion|nt)(|_)(|n|s)(|ot|trict)(|_)(|e)(|mpty)_A_%', '\\U\\1\\L\\2\\U\\4\\L\\5\\U\\7\\L\\8', rownames(dt_plot), perl=TRUE)
lgnd <- c("Int not empty", "Int strict", "Union")
lbls <- colnames(dt_plot)
# Plot to pdf
pdf("bar_changes_transcr_vs_genome.pdf", width = 21, height = 6)
barplot(dt_plot, ylim = c(0, ceiling(max(dt_plot))), main = "Relative changes (Transcriptome wt / Genome)", cex.names = 0.75, axis.lty = 1, legend = lgnd, beside = TRUE, names.arg = lbls)
dev.off()

# TRANSLATOME WT VS GENOME
# Subset data
transl_wt <- mt_genome_analysis[c(1:9,11:22), c(4, 10, 16)]
# Calculate ratios
dt_plot <- t(transl_wt / genome_distr)
# Labels and legend as above
# Plot to pdf
pdf("bar_changes_transl_vs_genome.pdf", width = 21, height = 6)
barplot(dt_plot, ylim = c(0, ceiling(max(dt_plot))), main = "Relative changes (Translatome wt / Genome)", cex.names = 0.75, axis.lty = 1, legend = lgnd, beside = TRUE, names.arg = lbls)
dev.off()

# TRANSLATOME DTT VS GENOME
# Subset data
transl_dtt <- mt_genome_analysis[c(1:9,11:22), c(8, 14, 20)]
# Calculate ratios
dt_plot <- t(transl_dtt / genome_distr)
# Labels and legend as above
# Plot to pdf
pdf("bar_changes_transl_dtt_vs_genome.pdf", width = 21, height = 6)
barplot(dt_plot, ylim = c(0, ceiling(max(dt_plot))), main = "Relative changes (Translatome DTT / Genome)", cex.names = 0.75, axis.lty = 1, legend = lgnd, beside = TRUE, names.arg = lbls)
dev.off()

# TRANSCRIPTOME WT VS TRANSLATOME WT
# Calculate ratios
dt_plot <- t(transl_wt / transcr_wt)
# Labels and legend as above
# Plot to pdf
pdf("bar_changes_transl_wt_vs_transcr_wt.pdf", width = 21, height = 6)
barplot(dt_plot, ylim = c(0, ceiling(max(dt_plot))), main = "Relative changes (Translatome wt / Transcriptome wt)", cex.names = 0.75, axis.lty = 1, legend = lgnd, beside = TRUE, names.arg = lbls)
dev.off()

# TRANSLATOME WT VS DTT
# Subset data
transl_dtt <- mt_genome_analysis[c(1:9,11:22), c(8, 14, 20)]
# Calculate ratios
dt_plot <- t(transl_dtt / transl_wt)
# Labels and legend as above
# Plot to pdf
pdf("bar_changes_transl_dtt_vs_wt.pdf", width = 21, height = 6)
barplot(dt_plot, ylim = c(0, ceiling(max(dt_plot))), main = "Relative changes (Translatome DTT / Translatome wt)", cex.names = 0.75, axis.lty = 1, legend = lgnd, beside = TRUE, names.arg = lbls)
dev.off()
###


### F. STATISTICS >>> DOES NOT WORK!?!

## GENOMIC TYPE

# TRANSCRIPTOME WT VS GENOME
# Subset data
chi_genome <- mt_genome_analysis[1:6, 2]
chi_transcr_wt_int_not_empty <- mt_genome_analysis[1:6, 6]
chi_transcr_wt_int_strict <- mt_genome_analysis[1:6, 12]
chi_transcr_wt_union <- mt_genome_analysis[1:6, 18]
# Calculate chi square
chisq.test(chi_genome, chi_transcr_wt_int_not_empty)
chisq.test(chi_genome, chi_transcr_wt_int_strict)
chisq.test(chi_genome, chi_transcr_wt_union)

# TRANSLATOME WT VS GENOME
# Subset data
chi_transl_wt_int_not_empty <- mt_genome_analysis[1:6, 4]
chi_transl_wt_int_strict <- mt_genome_analysis[1:6, 10]
chi_transl_wt_union <- mt_genome_analysis[1:6, 16]
# Calculate chi square
chisq.test(chi_genome, chi_transl_wt_int_not_empty)
chisq.test(chi_genome, chi_transl_wt_int_strict)
chisq.test(chi_genome, chi_transl_wt_union)

# TRANSCRIPTOME WT VS TRANSLATOME WT
# Calculate chi square
chisq.test(chi_transl_wt_int_not_empty, chi_transcr_wt_int_not_empty)
chisq.test(chi_transl_wt_int_strict, chi_transcr_wt_int_strict)
chisq.test(chi_transl_wt_union, chi_transcr_wt_union)

# TRANSLATOME WT VS DTT
# Subset data
chi_transl_dtt_int_not_empty <- mt_genome_analysis[1:6, 8]
chi_transl_dtt_int_strict <- mt_genome_analysis[1:6, 14]
chi_transl_dtt_union <- mt_genome_analysis[1:6, 20]
# Calculate chi square
chisq.test(chi_transl_dtt_int_not_empty, chi_transl_wt_int_not_empty)
chisq.test(chi_transl_dtt_int_strict, chi_transl_wt_int_strict)
chisq.test(chi_transl_dtt_union, chi_transl_wt_union)
##

## GENE TYPE

# TRANSCRIPTOME WT VS GENOME
# Subset data
chi_genome <- mt_genome_analysis[7:13, 2]
chi_transcr_wt_int_not_empty <- mt_genome_analysis[7:13, 6]
chi_transcr_wt_int_strict <- mt_genome_analysis[7:13, 12]
chi_transcr_wt_union <- mt_genome_analysis[7:13, 18]
# Calculate chi square
chisq.test(chi_genome, chi_transcr_wt_int_not_empty)
chisq.test(chi_genome, chi_transcr_wt_int_strict)
chisq.test(chi_genome, chi_transcr_wt_union)

# TRANSLATOME WT VS GENOME
# Subset data
chi_transl_wt_int_not_empty <- mt_genome_analysis[7:13, 4]
chi_transl_wt_int_strict <- mt_genome_analysis[7:13, 10]
chi_transl_wt_union <- mt_genome_analysis[7:13, 16]
# Calculate chi square
chisq.test(chi_genome, chi_transl_wt_int_not_empty)
chisq.test(chi_genome, chi_transl_wt_int_strict)
chisq.test(chi_genome, chi_transl_wt_union)

# TRANSCRIPTOME WT VS TRANSLATOME WT
# Calculate chi square
chisq.test(chi_transl_wt_int_not_empty, chi_transcr_wt_int_not_empty)
chisq.test(chi_transl_wt_int_strict, chi_transcr_wt_int_strict)
chisq.test(chi_transl_wt_union, chi_transcr_wt_union)

# TRANSLATOME WT VS DTT
# Subset data
chi_transl_dtt_int_not_empty <- mt_genome_analysis[7:13, 8]
chi_transl_dtt_int_strict <- mt_genome_analysis[7:13, 14]
chi_transl_dtt_union <- mt_genome_analysis[7:13, 20]
# Calculate chi square
chisq.test(chi_transl_dtt_int_not_empty, chi_transl_wt_int_not_empty)
chisq.test(chi_transl_dtt_int_strict, chi_transl_wt_int_strict)
chisq.test(chi_transl_dtt_union, chi_transl_wt_union)
###


### CLEAN-UP
# Remove all objects in environment
rm(list = ls())
# Remove unused packages
detach(package:GenomicRanges)
detach(package:IRanges)
detach(package:BiocGenerics)
###
