### A. PRE-REQUISITES
# Load library
suppressMessages(library(rtracklayer))
###

### B. GENOME IMPORT
# Use rtracklayer::import method to import human .gtf genome file to GRanges object 
hs71_features <- import("gencode.v16.annotation.gtf.gz", genome="hg19", asRangedData=FALSE) 
# Save HS71 features GRanges object
save(hs71_features, file="hs71_features.Rdata")
###

### C. CREATE GENE LIST
# Subset EXONS (discards several other categories, e.g. CDS, start_codon etc.)
hs71_exons <- hs71_features[values(hs71_features)[["type"]] == "exon"]
# Save exon GRanges object
save(hs71_exons, file="hs71_exons.Rdata")
# Split exons GRanges into GRangesList by 'gene_id' 
hs71_gene_list <- split(hs71_exons, hs71_exons$gene_id)
# For each gene (i.e. GRangesList element), merge overlapping exons into pseudo exons
hs71_gene_list <- endoapply(hs71_gene_list, function(gr) reduce(gr))
# Save GRangesList object
save(hs71_gene_list, file="hs71_gene_list_w_pseudo_exons.Rdata")
###

### D. EXPORT AS BED FILE
## BED file can only be created from GRangesList, but doing so directly blocks access to the individual pseudo exon ranges
## Detour via GRanges object, so that each pseudo exon keeps its ranges
# GRangesList to GRanges object, each pseudo exon in single row
gr <- unlist(hs71_gene_list)
# Assign additional metadata column 'pseudo_exon' (1: total number of pseudo exons) to GRanges
gr <- GRanges(seqnames=seqnames(gr), ranges=ranges(gr), strand=strand(gr), pseudo_exon=1:length(gr), seqlengths=seqlengths(gr))
# Split GRanges according to pseudo exon number, assign gene_id to each pseudo exon (multiple entries for each gene)
pseudo <- split(gr, gr$pseudo_exon)
names(pseudo) <- names(gr)
export.bed(pseudo, con="/import/bc2/home/zavolan/cherrmann/ALTSPLICE/GSE12946/info/hs71_pseudo_exons.bed")

