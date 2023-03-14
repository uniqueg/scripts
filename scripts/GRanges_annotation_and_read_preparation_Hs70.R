# Load required libraries
library(rtracklayer)				# also loads IRanges & GenomicRanges
library("BSgenome.Hsapiens.UCSC.hg19")		# only necessary if chromosome length information is required

# Set working directory

### Annotations

# Download annotation file, import as GRanges object, coerce to GRanges, sort and save
download.file("ftp://ftp.ensembl.org/pub/release-70/gtf/homo_sapiens/Homo_sapiens.GRCh37.70.gtf.gz", destfile="./Homo_sapiens.GRCh37.70.gtf.gz")
Hs70_gtf <- import(gzfile("./Homo_sapiens.GRCh37.70.gtf.gz"), "gtf")
Hs70_features <- sort(as(Hs70_gtf, "GRanges"))
# OPTIONAL: Add chromosome length info:		seqlengths(Hs70_features) <- seqlengths(Hsapiens)
# OPTIONAL: Subset features, e.g. exons:	Hs70_exons <- Hs70_features[mcols(Hs70_features)$type == "exon"]
# OPTIONAL: Group features, e.g. genes:		Hs70_genes <- split(Hs70_exons, mcols(Hs70_exons)$gene_id)
save(Hs70_features, file=("./Hs70_features.R")

### Reads (example: BED file reads.bed)

# Load reads table file, construct GRanges object and save
reads <- read.table("./reads.bed")
reads_gr <- GRanges(seqnames=Rle(reads[,1]), ranges=IRanges(start=reads[,2], end=reads[,3]), strand=reads[,6], seqinfo=Hsapiens@seqinfo)	# optionally: add metadata such as scores; see ?GRanges for details
save(reads_gr, file=("./reads.R")
