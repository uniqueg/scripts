#!/usr/bin/env Rscript

# (c) 2017 Alexander Kanitz & Maria Katsantoni, Biozentrum, University of Basel
# (@) alexander.kanitz@unibas.ch


# TODO: implement an options '--regions' that is mutually exclusive with '--chromosome-sizes' and accepts a BED file of regions for which sliding windows shall be generated


#######################
###  PARSE OPTIONS  ###
#######################

# Load 'optparse' package
if ( suppressWarnings(suppressPackageStartupMessages(require("optparse"))) == FALSE ) { stop("[ERROR] Package 'optparse' required! Aborted.") }

# Get script name
script <- sub("--file=", "", basename(commandArgs(trailingOnly=FALSE)[4]))

# Build description message
description <- "Generates sliding windows of a genome and returns a BED file listing the number of BAM alignments overlapping each window.\n"
author <- "Author: Alexander Kanitz & Maria Katsantoni, Biozentrum, University of Basel"
version <- "Version: 1.1.0 (05-MAY-2017)"
requirements <- "Requires: optparse, Bioconductor, GenomicRanges, GenomicAlignments, Rsamtools, rtracklayer"
msg <- paste(description, author, version, requirements, sep="\n")

# Define list of arguments
option_list <- list(
    make_option(
        "--bam",
        action="store",
        type="character",
        default=NULL,
        help="BAM file. Required.",
        metavar="bam"
    ),
    make_option(
        "--bam-background",
        action="store",
        type="character",
        default=NULL,
        help="BAM file for calculating (optional) background counts. A separate output BED file with identical window names is generated.",
        metavar="bam"
    ),
    make_option(
        "--bai",
        action="store",
        type="character",
        default=NULL,
        help="BAM index (BAI) file. If not provided, the index file path is generated from the BAM file path by adding '.bai' as a (second) file extension.",
        metavar="bai"
    ),
    make_option(
        "--bai-background",
        action="store",
        type="character",
        default=NULL,
        help="BAM index (BAI) file for '--background-bam'. If not provided, the index file path is generated from the background BAM file path by adding '.bai' as a (second) file extension.",
        metavar="bai"
    ),
    make_option(
        "--ratios",
        action="store_true",
        default=FALSE,
        help="Indicate if an additional file with foreground over background count ratios shall be generated. Ignored if '--bam-background' is not supplied.",
    ),
    make_option(
        "--output-directory",
        action="store",
        type="character",
        default=getwd(),
        help="Directory where output files are written. If not provided, the current working directory is used.",
        metavar="dir"
    ),
    make_option(
        "--output-prefix",
        action="store",
        type="character",
        default="overlaps",
        help="Output filename prefix. Default: '%default'.",
        metavar="dir"
    ),
    make_option(
        "--paired",
        action="store",
        type="integer",
        default=NULL,
        help="For paired-end libraries, supply either 0, 1 or 2 to perform the counting on both mates, mates 1 or mates 2, respectively. By default, a single-end library is assumed.",
        metavar="int"
    ),
    make_option(
        "--discard-singletons",
        action="store_true",
        default=NA,
        help="For paired-end libraries, indicates whether unpaired reads shall be ignored. Option is ignored if '--paired' is unset.",
    ),
    make_option(
        "--proper-pairs-only",
        action="store_true",
        default=NA,
        help="For paired-end libraries, indicates whether only properly paired reads shall be considered. Option implies '--discard-singletons'. Option is ignored if '--paired' is unset.",
    ),
    make_option(
        "--chromosome-sizes",
        action="store",
        type="character",
        default=NULL,
        help="Headerless tab-separated file listing chromosome names and sizes in columns 1 and 2, respectively. If not provided, chromosome size information is obtained from the BAM file.",
        metavar="tsv"
    ),
    make_option(
        "--window-size",
        action="store",
        type="integer",
        default=1000,
        help="Width of each sliding window [default: %default].",
        metavar="int"
    ),
    make_option(
        "--step-size",
        action="store",
        type="integer",
        default=1000,
        help="Step size for the generation of sliding windows [default: %default]. If different from '--window-size', overlapping or gapped sliding windows will be generated.",
        metavar="int"
    ),
    make_option(
        "--window-size-background",
        action="store",
        type="integer",
        default=NULL,
        help="Width of each background sliding window [default: %default]. By default, assumes the same value as '--window-size' if '--bam-background' is supplied. Ignored if '--bam-background' is not supplied.",
        metavar="int"
    ),
    make_option(
        "--remove-incomplete-windows",
        action="store_true",
        default=FALSE,
        help="Do not consider incomplete windows at the ends of chromosomes. If '--bam-background' is supplied, incomplete windows will also be removed from the background and, moreover, it is ensured that only complete window:background window pairs remain."
    ),
    make_option(
        "--strand",
        action="store",
        type="character",
        default="+-",
        help="Determines the strand of the generated sliding windows. One of '+-' (counts for both plus and minus strand reported in same file; default), '+/-' (counts for both plus and minus strand reported in separate files; default), '*' (sum of counts for plus and minus strand), '+' and '-' (counts for plus and minus strands only, respectively).",
        metavar="string"
    ),
    make_option(
        "--overlap-type",
        action="store",
        type="character",
        default="any",
        help="One of 'any' (default), '5p', '3p', 'left', 'right', or 'within'. By default, all overlaps are accepted. Set this option to count overlaps towards a given window only if a read's 5'-end ('5p'), 3'-end ('3p'), left side ('left'), or right side ('right') overlap with it, or if it is fully contained in it ('within').",
        metavar="string"
    ),
    make_option(
        c("-h", "--help"),
        action="store_true",
        default=FALSE,
        help="Show this information and die."
    ),
    make_option(
        c("-u", "--usage"),
        action="store_true",
        default=FALSE,
        dest="help",
        help="Show this information and die."
    ),
    make_option(
        c("-v", "--verbose"),
        action="store_true",
        default=FALSE,
        help="Print log messages to STDOUT."
    )
)

# Parse command-line arguments
opt_parser <- OptionParser(usage=paste("Usage:", script, "[OPTIONS] --bam <BAM>\n", sep=" "), option_list = option_list, add_help_option=FALSE, description=msg)
opt <- parse_args(opt_parser)

# Re-assign variables
chr.sizes.file <- opt[["chromosome-sizes"]]
bam.file <- opt[["bam"]]
bam.bg.file <- opt[["bam-background"]]
bai.file <- opt[["bai"]]
bai.bg.file <- opt[["bai-background"]]
ratios <- opt[["ratios"]]
out.dir <- opt[["output-directory"]]
out.prefix <- opt[["output-prefix"]]
paired <- opt[["paired"]]
is.paired <- opt[["discard-singletons"]]
is.proper.pair <- opt[["proper-pairs-only"]]
win.size <- opt[["window-size"]]
step.size <- opt[["step-size"]]
win.size.bg <- opt[["window-size-background"]]
remove.incomplete <- opt[["remove-incomplete-windows"]]
str <- opt[["strand"]]
type <- opt[["overlap-type"]]
verb <- opt[["verbose"]]

# Validate required arguments
if ( is.null(bam.file) ) {
    print_help(opt_parser)
    stop("[ERROR] Required argument missing! Aborted.")
}

# Validate allowed values
str.allowed <- c("+-", "+/-", "*", "+", "-")
type.allowed <- c("any", "5p", "3p", "left", "right", "within")
paired.allowed <- c(0, 1, 2)
if ( ! str %in% str.allowed ) {
    print_help(opt_parser)
    stop("[ERROR] Illegal argument for option '--strand'! Aborted.")
}
if ( ! type %in% type.allowed ) {
    print_help(opt_parser)
    stop("[ERROR] Illegal argument for option '--overlap-type'! Aborted.")
}
if ( ! paired %in% paired.allowed && ! is.null(paired) ) {
    print_help(opt_parser)
    stop("[ERROR] Illegal argument for option '--paired'! Aborted.")
}

# Set dependent options
if ( is.null(paired) ) is.paired <- is.proper.pair <- NA
if ( ! is.na(is.proper.pair) && is.proper.pair ) is.paired <- TRUE
if ( ! is.null(bam.bg.file) && is.null(win.size.bg) ) win.size.bg <- win.size
out.prefix <- file.path(out.dir, out.prefix)


#################
###  IMPORTS  ###
#################

# Import required packages
if ( verb ) cat("Loading required packages...\n", sep="'")
if ( suppressWarnings(suppressPackageStartupMessages(require("GenomicRanges"))) == FALSE ) { stop("[ERROR] Package 'GenomicRanges' required! Aborted.") }
if ( suppressWarnings(suppressPackageStartupMessages(require("GenomicAlignments"))) == FALSE ) { stop("[ERROR] Package 'GenomicAlignments' required! Aborted.") }
if ( suppressWarnings(suppressPackageStartupMessages(require("Rsamtools"))) == FALSE ) { stop("[ERROR] Package 'Rsamtools' required! Aborted.") }
if ( suppressWarnings(suppressPackageStartupMessages(require("rtracklayer"))) == FALSE ) { stop("[ERROR] Package 'rtracklayer' required! Aborted.") }


##############
###  MAIN  ###
##############

# Import BAM file as GAlignments object
if ( verb ) cat("Importing alignments...\n", sep="'")
if ( is.null(paired) ) {
    param <- NULL
} else if ( paired == 0 ) {
    param <- ScanBamParam(
        flag=scanBamFlag(
            isPaired = is.paired,
            isProperPair = is.proper.pair,
            isFirstMateRead = NA,
            isSecondMateRead = NA,
        )
    )
} else if ( paired == 1 ) {
    param <- ScanBamParam(
        flag=scanBamFlag(
            isPaired = is.paired,
            isProperPair = is.proper.pair,
            isFirstMateRead = TRUE,
            isSecondMateRead = FALSE,
        )
    )
} else if ( paired == 2 ) {
    param <- ScanBamParam(
        flag=scanBamFlag(
            isPaired = is.paired,
            isProperPair = is.proper.pair,
            isFirstMateRead = FALSE,
            isSecondMateRead = TRUE,
        )
    )
}
if ( is.null(bai.file) ) bai.file <- bam.file
gr <- as(readGAlignments(file=bam.file, index=bai.file, param=param), "GRanges")
if ( type == "5p" ) {
    start(gr) <- ifelse(strand(gr) == "+", start(gr), end(gr))
    end(gr) <- ifelse(strand(gr) == "+", start(gr), end(gr))
}
if ( type == "3p" ) {
    start(gr) <- ifelse(strand(gr) == "+", end(gr), start(gr))
    end(gr) <- ifelse(strand(gr) == "+", end(gr), start(gr))
}
if ( type == "left" ) {
    start(gr) <- end(gr) <- pmin(start(gr), end(gr))
}
if ( type == "right" ) {
    start(gr) <- end(gr) <- pmax(start(gr), end(gr))
}
if ( ! is.null(bam.bg.file) ) {
    if ( ! is.null(bam.bg.file) && is.null(bai.bg.file) ) bai.bg.file <- bam.bg.file
    gr.bg <- as(readGAlignments(file=bam.bg.file, index=bai.bg.file, param=param), "GRanges")
    if ( type == "5p" ) {
        start(gr.bg) <- ifelse(strand(gr.bg) == "+", start(gr.bg), end(gr.bg))
        end(gr.bg) <- ifelse(strand(gr.bg) == "+", start(gr.bg), end(gr.bg))
    }
    if ( type == "3p" ) {
        start(gr.bg) <- ifelse(strand(gr.bg) == "+", end(gr.bg), start(gr.bg))
        end(gr.bg) <- ifelse(strand(gr.bg) == "+", end(gr.bg), start(gr.bg))
    }
    if ( type == "left" ) {
        start(gr.bg) <- end(gr.bg) <- pmin(start(gr.bg), end(gr.bg))
    }
    if ( type == "right" ) {
        start(gr.bg) <- end(gr.bg) <- pmax(start(gr.bg), end(gr.bg))
    }
}

# Obtain chromosome sizes as Seqinfo object
if ( verb ) cat("Retrieving chromosome sizes...\n", sep="'")
if ( is.null(chr.sizes.file) ) {
    seq.info <- seqinfo(gr)
} else {
    chr.sizes <- read.delim(chr.sizes.file, stringsAsFactors=FALSE, header=FALSE)
    seq.info <- Seqinfo(chr.sizes[, 1], seqlengths=chr.sizes[, 2])
}

# Generate sliding windows
if ( verb ) cat("Generating sliding windows...\n", sep="'")
wins <- tileGenome(seqlengths=seq.info, tilewidth=step.size, cut.last.tile.in.chrom=TRUE)
end(wins) <- pmin(start(wins) + win.size - 1, seqlengths(wins)[as.character(seqnames(wins))])
if ( str != "*" ) {
    wins.plus <- wins.minus <- wins
    strand(wins.plus) <- "+"
    strand(wins.minus) <- "-"  
    if ( str == "+-" || str == "+/-" ) wins <- c(wins.plus, wins.minus)
    if ( str == "+" ) wins <- wins.plus
    if ( str == "-" ) wins <- wins.minus
}
wins$name <- paste(as.character(seqnames(wins)), paste(as.character(start(wins) - 1), as.character(end(wins)), sep="-"), as.character(strand(wins)), sep=":")
if ( ! is.null(bam.bg.file) ) {
    wins.bg <- wins
    shift <- (win.size.bg - win.size) / 2
    start(wins.bg) <- pmax(start(wins.bg) - shift, 1)
    end(wins.bg) <- pmin(end(wins.bg) + shift, seqlengths(wins.bg)[as.character(seqnames(wins.bg))])
}

# Remove incomplete windows
if ( remove.incomplete ) {
    if ( verb ) cat("Removing incomplete windows...\n", sep="'")
    wins <- wins[width(wins) == win.size]
    if ( ! is.null(bam.bg.file) ) {
        wins.bg <- wins.bg[width(wins.bg) == win.size.bg]
        wins <- wins[wins$name %in% wins.bg$name]
        wins.bg <- wins.bg[wins.bg$name %in% wins$name]
    }
}

# Count overlaps
if ( verb ) cat("Counting overlaps...\n", sep="'")
if ( type == "within" ) {
    hits <- findOverlaps(gr, wins, type="within")
    cts <- merge(1:length(wins), as.data.frame(table(subjectHits(hits))), by=1, all=TRUE)[, 2]
    cts[is.na(cts)] <- 0
    wins$score <- cts
    if ( ! is.null(bam.bg.file) ) {
        hits.bg <- findOverlaps(gr.bg, wins.bg, type="within")
        cts.bg <- merge(1:length(wins.bg), as.data.frame(table(subjectHits(hits.bg))), by=1, all=TRUE)[, 2]
        cts.bg[is.na(cts.bg)] <- 0
        wins.bg$score <- cts.bg
    }
} else {
    wins$score <- countOverlaps(wins, gr, type="any")
    if ( ! is.null(bam.bg.file) ) wins.bg$score <- countOverlaps(wins.bg, gr.bg, type="any")
}
if ( ! is.null(bam.bg.file) && ratios ) {
    wins.ratios <- wins
    wins.ratios$score <- wins$score / wins.bg$score
    wins.ratios$score[is.nan(wins.ratios$score) | is.na(wins.ratios$score) | is.infinite(wins.ratios$score)] <- 0
}
if ( str == "+/-" ) {
    wins.plus <- wins[strand(wins) == "+"]
    wins.minus <- wins[strand(wins) == "-"]
    if ( ! is.null(bam.bg.file) ) {
        wins.bg.plus <- wins.bg[strand(wins.bg) == "+"]
        wins.bg.minus <- wins.bg[strand(wins.bg) == "-"]
        if ( ratios ) {
            wins.ratios.plus <- wins.ratios[strand(wins.ratios) == "+"]
            wins.ratios.minus <- wins.ratios[strand(wins.ratios) == "-"]
        }
    }
}

# Write output files in BED format
if ( verb ) cat("Writing output...\n", sep="'")
if ( str == "+/-" ) {
    export(wins.plus, paste(out.prefix, "fg", "plus", "bed", sep="."), format="bed")
    export(wins.minus, paste(out.prefix, "fg", "minus", "bed", sep=".") , format="bed")
    if ( ! is.null(bam.bg.file) ) {
        export(wins.bg.plus, paste(out.prefix, "bg", "plus", "bed", sep=".") , format="bed")
        export(wins.bg.minus, paste(out.prefix, "bg", "minus", "bed", sep=".") , format="bed")
        if ( ratios ) {
            export(wins.ratios.plus, paste(out.prefix, "fg_over_bg", "plus", "bed", sep=".") , format="bed")
            export(wins.ratios.minus, paste(out.prefix, "fg_over_bg", "minus", "bed", sep=".") , format="bed")
        }
    }
} else {
    export(wins, paste(out.prefix, "fg", "bed", sep=".") , format="bed")
    if ( ! is.null(bam.bg.file) ) {
        export(wins.bg, paste(out.prefix, "bg", "bed", sep=".") , format="bed")
        if ( ratios ) {
            export(wins.ratios, paste(out.prefix, "fg_over_bg", "bed", sep=".") , format="bed")
        }
    }
}

# Write log
cat("Done.\n")
