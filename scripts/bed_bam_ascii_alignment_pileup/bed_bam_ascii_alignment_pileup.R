#!/usr/bin/env Rscript

#==================#
#   HEADER START   #
#==================#
### Created: Sep 28, 2019
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
### Contact: alexander.kanitz@alumni.ethz.ch
#==================#
#    HEADER END    #
#==================#


#========================#
#   SESSION INFO START   #
#========================#
# R version 3.6.0 (2019-04-26)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)
#
# Matrix products: default
# BLAS/LAPACK: /scicore/soft/apps/OpenBLAS/0.3.1-GCC-7.3.0-2.30/lib/libopenblas_sandybridgep-r0.3.1.so
#
# locale:
#  [1] LC_CTYPE=en_US             LC_NUMERIC=C
#  [3] LC_TIME=de_CH.UTF-8        LC_COLLATE=en_US
#  [5] LC_MONETARY=de_CH.UTF-8    LC_MESSAGES=en_US
#  [7] LC_PAPER=de_CH.UTF-8       LC_NAME=C
#  [9] LC_ADDRESS=C               LC_TELEPHONE=C
# [11] LC_MEASUREMENT=de_CH.UTF-8 LC_IDENTIFICATION=C
#
# attached base packages:
# [1] stats4    parallel  stats     graphics  grDevices utils     datasets [8] methods   base
#
# other attached packages:
#  [1] rtracklayer_1.44.0          GenomicFeatures_1.36.0
#  [3] AnnotationDbi_1.46.0        GenomicAlignments_1.20.0
#  [5] Rsamtools_2.0.0             Biostrings_2.52.0
#  [7] XVector_0.24.0              SummarizedExperiment_1.14.0
#  [9] DelayedArray_0.10.0         BiocParallel_1.18.0
# [11] matrixStats_0.54.0          Biobase_2.44.0
# [13] GenomicRanges_1.36.0        GenomeInfoDb_1.20.0
# [15] IRanges_2.18.0              S4Vectors_0.22.0
# [17] BiocGenerics_0.30.0
# 
# loaded via a namespace (and not attached):
#  [1] Rcpp_1.0.1             compiler_3.6.0         prettyunits_1.0.2
#  [4] bitops_1.0-6           tools_3.6.0            zlibbioc_1.30.0
#  [7] progress_1.2.2         biomaRt_2.40.0         digest_0.6.18
# [10] bit_1.1-14             RSQLite_2.1.1          memoise_1.1.0
# [13] lattice_0.20-38        pkgconfig_2.0.2        rlang_0.3.4
# [16] Matrix_1.2-17          DBI_1.0.0              GenomeInfoDbData_1.2.1
# [19] httr_1.4.0             stringr_1.4.0          hms_0.4.2
# [22] bit64_0.9-7            grid_3.6.0             R6_2.4.0
# [25] XML_3.98-1.19          magrittr_1.5           blob_1.1.1
# [28] assertthat_0.2.1       stringi_1.4.3          RCurl_1.95-4.12       [31] crayon_1.3.4
#========================#
#    SESSION INFO END    #
#========================#


#==========================#
#   PRE-REQUISITES START   #
#==========================#
#---> LOAD OPTION PARSER <---#
if ( suppressWarnings(suppressPackageStartupMessages(require("optparse"))) == FALSE ) { stop("Package 'optparse' required!\nExecution aborted.") }

#---> GET SCRIPT NAME <---#
script <- sub("--file=", "", basename(commandArgs(trailingOnly=FALSE)[4]))

#---> DESCRIPTION <---#
description <- "Generates an ASCII-style pileup of read alignments in one or more BAM files against one or more regions specified in a BED file.\n"
author <- "Author: Alexander Kanitz, Biozentrum, University of Basel; alexander.kanitz@alumni.ethz.ch"
version <- "Version: 1.0.0 (28-SEP-2019)"
requirements <- c("optparse", "rtracklayer", "GenomicFeatures", "GenomicAlignments")
requirements_txt <- paste("Requires:", paste(requirements, collapse=", "), sep=" ")
msg <- paste(description, author, version, requirements_txt, sep="\n")

#---> COMMAND-LINE ARGUMENTS <---#
## List of allowed/recognized arguments
option_list <- list(
        make_option(
            "--bed",
            action="store",
            type="character",
            default=NULL,
            help="Required: Path to BED file of query regions. Consider the `--maximum-region-width` parameter.",
            metavar="file"
        ),
        make_option(
            "--bam",
            action="store",
            type="character",
            default=NULL,
            help="Required: Path to BAM file containing the alignments to be visualized against the query region(s). It is possible to use glob wildcards to select multiple BAM files at once. However, the feature is experimental and may lead to inconsistencies if a read contains an indel.",
            metavar="file"
        ),
        make_option(
            "--reference",
            action="store",
            type="character",
            default=NULL,
            help="Reference genome sequence in FASTA format. The file *MUST* be compressed with BGZIP. If supplied, the reference sequence for the query region(s) will be added to the output. Note that on the first run with a specific reference genome file, an FAI index is generated which will take some time.",
            metavar="file"
        ),
        make_option(
            "--annotations",
            action="store",
            type="character",
            default=NULL,
            help="Annotation file in GFF/GTF format used to annotate sequences. If supplied, features overlapping the query region(s) will be visualized in the output. Ensure that the argument to option `annotation-name-field` corresponds to a field in the annotations, otherwise the script will fail.",
            metavar="file"
        ),
        make_option(
            "--output-directory",
            action="store",
            type="character",
            default=getwd(),
            help="Output directory. One output file will be created for each region in `--bed` and the filenames will be generated from the basenames of the supplied BAM file(s) and the name field (4th column) of the BED file. [default \"%default\"]",
            metavar="dir"
        ),
        make_option(
            "--maximum-region-width",
            action="store",
            type="integer",
            default=200,
            help="Maximum input region width. Use with care as wide regions will use excessive resources. [default %default]",
            metavar="int"
        ),
        make_option(
            "--do-not-collapse-alignments",
            action="store_true",
            default=FALSE,
            help="Show alignments of reads with identical sequences individually."
        ),
        make_option(
            "--minimum-count",
            action="store",
            type="integer",
            default=1,
            help="Alignments of reads with less copies than the specified number will not be printed. Option is not considered if `do-not-collapse-alignments` is set. [default %default]",
            metavar="int"
        ),
        make_option(
            "--annotation-name-field",
            action="store",
            type="character",
            default="Name",
            help="Annotation field used to populate the `name` column in the output. [default \"%default\"]",
            metavar="str"
        ),
        make_option(
            "--padding-character",
            action="store",
            default=".",
            help="Character used for padding alignments. [default %default]",
            metavar="char"
        ),
        make_option(
            "--indel-character",
            action="store",
            default="-",
            help="Character to denote insertions and deletions in alignments. [default %default]",
            metavar="char"
        ),
        make_option(
            c("-h", "--help"),
            action="store_true",
            default=FALSE,
            help="Show this information and die."
        ),
        make_option(
            c("-v", "--verbose"),
            action="store_true",
            default=FALSE,
            help="Print log messages to STDOUT."
        )
)

## Parse command-line arguments
opt_parser <- OptionParser(usage="Usage: %prog [OPTIONS] --bed <PATH> --bam <PATH>\n", option_list = option_list, add_help_option=FALSE, description=msg)
opt <- parse_args(opt_parser)

## Die if any required arguments are missing...
if ( is.null(opt[["bed"]]) || is.null(opt[["bam"]]) ) {
    write("[ERROR] Required argument(s) missing!\n\n", stderr())
    stop(print_help(opt_parser))
}
#==========================#
#    PRE-REQUISITES END    #
#==========================#


#================#
#   MAIN START   #
#================#
#---> START MESSAGE <---#
if ( opt$verbose ) cat("Starting '", script, "'...\n", sep="")

#---> LOAD REQUIRED LIBRARIES <---#
if ( opt$verbose ) cat("Loading libraries...\n", sep="")
for (req in requirements) {
    if ( suppressWarnings(suppressPackageStartupMessages(require(req, character.only=TRUE))) == FALSE ) { cat("Package '", req, "' required!\nExecution aborted.") }
}

#---> RE-ASSIGN CLI ARGUMENTS <---#
if ( opt$verbose ) cat("Loading libraries...\n", sep="")
# Re-assign CLI arguments
fl.query <- opt[["bed"]]
fl.bam <- opt[["bam"]]
fl.ref <- opt[["reference"]]
fl.anno <- opt[["annotations"]]
dir.out <- opt[["output-directory"]]
width.max <- opt[["maximum-region-width"]]
collapse <- ! opt[["do-not-collapse-alignments"]]
count.min <- opt[["minimum-count"]]
char.pad <- opt[["padding-character"]]
char.indel <- opt[["indel-character"]]
field.name.anno <- opt[["annotation-name-field"]]

#---> IMPORT FILES <---#
# Print status message
if ( opt$verbose ) cat("Importing input files...\n", sep="")
# Load files
bed <- import(con=fl.query)
if (! is.null(fl.ref)) {ref <- FaFile(fl.ref)}
if (! is.null(fl.anno)) {anno <- import(con=fl.anno)}

#---> FIND BAM FILES <---#
# Print status message
if ( opt$verbose ) cat("Searching for BAM files...\n", sep="")
# Get list of BAM files
fl.bam <- dir(dirname(fl.bam), pattern=glob2rx(basename(fl.bam)), full.names=TRUE)
if (! length(fl.bam)) {stop("No BAM files found in input directory")}
fl.prefix <- paste(basename(tools::file_path_sans_ext(fl.bam)), collapse=".")

#--->   <---#
# Print status message
if ( opt$verbose ) cat("Iterating over regions in BED file...\n")
# Iterate over input regions
for(index in seq_along(bed)) {

    #---> PREPARE STACKING OF ALIGNMENTS  <---#
    # Assign current region
    region <- bed[index]
    # Print status message
    if ( opt$verbose ) cat("Processing region '", mcols(region)[["name"]], "'...\n", sep="")
    # Exit with error if region is too wide
    if (width(region) > width.max) {stop("Supplied region too large. Consider increasing the `width.max` parameter, but note that for very large regions, the memory footprint may be excessive.")}
    # Initialize DNAStringSet container object
    seq.out <- DNAStringSet()

    #---> ADD READ ALIGNMENTS  <---#
    # Print status message
    if ( opt$verbose ) cat("Iterating over BAM files...\n")
    # Iterate over BAM files
    for (bam in fl.bam) {
        # Print status message
        if ( opt$verbose ) cat("Adding read alignments for file '", bam,"'...\n", sep="")
        # Get alignments for current BAM file
        seq.stack <- stackStringsFromBam(
            bam,
            param=region,
            D.letter=char.indel,
            N.letter=char.indel,
            Lpadding.letter=char.pad,
            Rpadding.letter=char.pad
        )
        # Add to container
        seq.out <- append(seq.out, seq.stack)
    }
    # Get reverse complement if on minus strand
    if (as.character(strand(region))[[1]] == "-") {
        seq.out <- complement(seq.out)
    }
    # Convert to character
    seq.out <- as.character(seq.out)
    # Convert to dataframe for writing output
    df <- data.frame(seq=seq.out, count=rep(NA, length(seq.out)), stringsAsFactors=FALSE)

    #---> COLLAPSE READ ALIGNMENTS <---#
    if (collapse) {
        # Print status message
        if ( opt$verbose ) cat("Collapsing identical reads/alignments...\n")
        # Get unique alignments
        df <- data.frame(table(df[["seq"]]), stringsAsFactors=FALSE)
        if (! ncol(df) == 2) df <- data.frame(matrix(ncol = 2, nrow = 0))
        colnames(df) <- c("seq", "count")
        df[["seq"]] <- as.character(df[["seq"]])
        # Filter out any alignments that do not make the specified minimum count cutoff
        df[df[["count"]] >= count.min, ]
    }

    #---> SORTING ALIGNMENTS  <---#
    # Print status message
    if ( opt$verbose ) cat("Sorting alignments...\n")
    # Sort by position of first nucleotide, count and position of last nucleotide
    if (nrow(df)) {
        last_char <- nchar(df[["seq"]][[1]])
        pos.nuc.first <- regexpr(paste0("[^", char.pad, "\\.]"), df[["seq"]])
        pos.nuc.last <- last_char - regexpr(paste0("[^", char.pad, "\\.]"), unlist(lapply(df[["seq"]], reverse))) + 1
        df <- df[order(
            pos.nuc.first, df[["count"]], pos.nuc.last,
            decreasing=c(FALSE, TRUE, FALSE)
        ), ]
    }

    #---> ADD REFERENCE SEQUENCE  <---#
    if (! is.null(fl.ref)) {
        # Print status message
        if ( opt$verbose ) cat("Adding reference sequence...\n")
        # Generate name for reference sequence
        name.ref <- paste(
            seqnames(region),
            paste(start(region), end(region), sep="-"),
            strand(region), sep=":"
        )
        # Get sequence
        row.ref <- getSeq(ref, region)  # will create .fai index if not present
        # Get reverse if on minus strand
        if (as.character(strand(region))[[1]] == "-") {
            row.ref <- reverse(row.ref)
        }
        # Compile row and add to dataframe
        df.ref <- data.frame(seq=as.character(row.ref), count=name.ref, stringsAsFactors=FALSE)
        df <- rbind(df.ref, df)
    }

    #---> ADD ANNOTATIONS FOR REGION <---#
    # Add annotations
    if (! is.null(fl.anno)) {
        # Print status message
        if ( opt$verbose ) cat("Adding overlapping features...\n")
        # Find overlapping features
        features <- anno[to(findOverlaps(region, anno))]
        # Set order of addition from most upstream to most downstream
        if (as.character(strand(region)) == "+") {
            order.features <- order(start(features))
        } else {
            order.features <- rev(order(end(features)))
        }
        # Order features
        features <- features[order.features]
        # Iterate over features
        seqs <- sapply(seq_along(features), function(index) {
            feat <- features[index]
            # Set markup character depending on strand
            char.anno <- ifelse(strand(region) == "+", ">", "<")
            # Calculate lengths of feature (in region) & left/right padding
            diff.start <- start(feat) - start(region)
            n.padl <- max(0, diff.start)
            n.anno <- min(min(0, diff.start) + width(feat), width(region) - n.padl)
            n.padr <- width(region) - n.padl - n.anno
            # Generate string encompassing feature
            str.padl <- paste(rep(char.pad, n.padl), collapse="")
            str.anno <- paste(rep(char.anno, n.anno), collapse="")
            str.padr <- paste(rep(char.pad, n.padr), collapse="")
            str.final <- paste0(str.padl, str.anno, str.padr)
        })
        df.anno <- data.frame(seq=seqs, count=mcols(features)[[field.name.anno]], stringsAsFactors=FALSE)
        df <- rbind(df.anno, df)
    }

    #---> WRITE OUTPUT  <---#
    fl.out <- file.path(dir.out, paste(fl.prefix, mcols(region)[["name"]], "pileup", "tab", sep="."))
    # Print status message
    if ( opt$verbose ) cat("Writing output to file '", fl.out, "'...\n", sep="")
    # Write tab-separated output
    write.table(df, fl.out, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
}

#---> END MESSAGE <---#
if ( opt$verbose ) cat("Done.\n")
#================#
#    MAIN END    #
#================#
