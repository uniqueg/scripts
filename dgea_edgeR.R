#!/usr/bin/env Rscript

# (c) 2016 Alexander Kanitz, Biozentrum, University of Basel
# (@) alexander.kanitz@unibas.ch


#################
###  IMPORTS  ###
#################

# Import required packages
if ( suppressWarnings(suppressPackageStartupMessages(require("optparse"))) == FALSE ) { stop("[ERROR] Package 'optparse' required! Aborted.") }
if ( suppressWarnings(suppressPackageStartupMessages(require("edgeR"))) == FALSE ) { stop("[ERROR] Package 'edgeR' required! Aborted.") }


#######################
###  PARSE OPTIONS  ###
#######################

# Get script name
script <- sub("--file=", "", basename(commandArgs(trailingOnly=FALSE)[4]))

# Build description message
description <- "Run differential (gene) expression analyses with edgeR.\n"
author <- "Author: Alexander Kanitz, Biozentrum, University of Basel"
version <- "Version: 1.0.0 (13-JAN-2017)"
requirements <- "Requires: edgeR, optparse"
msg <- paste(description, author, version, requirements, sep="\n")

# Define list of arguments
option_list <- list(
                make_option(
                    "--count-table",
                    action="store",
                    type="character",
                    default=NULL,
                    help="Count table containing feature/gene IDs in the first column, sequencing read counts in the remaining columns and a header line with sample IDs. The header line must contain one field less than the remaining lines (i.e. no field for the first/feature ID column). Required.",
                    metavar="tsv"
                ),
                make_option(
                    "--output-directory",
                    action="store",
                    type="character",
                    default=".",
                    help="Directory where output files shall be written. Default: \".\".",
                    metavar="directory"
                ),
                make_option(
                    "--run-id",
                    action="store",
                    type="character",
                    default="experiment",
                    help="String used as analysis identifier prefix for output files. Default: 'experiment'.",
                    metavar="file"
                ),
                make_option(
                    "--exclude",
                    action="store",
                    type="character",
                    default=NULL,
                    help="Text file containing identifiers of samples to be excluded from the analysis. Expected format: One identifier per line. Identifiers have to match column headers of '--count-table'. Default: NULL.",
                    metavar="file"
                ),
                make_option(
                    "--id-field-separator",
                    action="store",
                    type="character",
                    default=".",
                    help="Character for splitting sample identifiers (i.e. column names of '--count-table') for the extraction of group identifiers together with '--group-id-fields'. If multiple sample identifier components are used to construct group identifiers, the individual components are joined together by the same character. Ignored if '--annotation' is supplied. Default: '.'.",
                    metavar="char"
                ),
                make_option(
                    "--group-id-fields",
                    action="store",
                    type="character",
                    default="1",
                    help="String for identifying the components of sample identifiers that are to be used for the construction of group identifiers (together with '--id-field-separator'). The format is similar to that of field selection option of the Bash 'cut' command, i.e. a string of integers, separted by either commas or dashes. Example: '1,7,3-5' would use sample identifier fields 1, 7, 3, 4, and 5 (in that order) to build the group identifier. Ignored if '--annotation' is supplied. Default: 1 (i.e. all characters of the sample identifier until the first occurence of '--id-field-separator').",
                    metavar="string"
                ),
                make_option(
                    "--annotation",
                    action="store",
                    type="character",
                    default=NULL,
                    help="Annotation table with one row for each sample, a sample identifier column whose values match the column names of '--count-table' (set with '--sample-id-column') and a column with sample group information for each sample (set with '--group-column'). If not provided, sample groups are derived from the column names of '--count-table'. If provided, samples without annotation are not considered. Specify '--has-header-annotation' if table includes a header line. Default: NULL.",
                    metavar="tsv"
                ),
                make_option(
                    "--has-header-annotation",
                    action="store_true",
                    default=FALSE,
                    help="Indicates whether the annotation table includes a header line. Default: FALSE."
                ),
                make_option(
                    "--sample-id-column",
                    action="store",
                    type="integer",
                    default=1,
                    help="Annotation table column/field number containing sample identifiers. Ignored if '--annotation' is not specified. Default: 1.",
                    metavar="int"
                ),
                make_option(
                    "--group-column",
                    action="store",
                    type="integer",
                    default=2,
                    help="Annotation table column/field number containing category/group information for each sample. Samples belonging to the same group are analyzed together, by comparing them to samples of another group. Default: 2.",
                    metavar="int"
                ),
                make_option(
                    "--comparisons",
                    action="store",
                    type="character",
                    default=NULL,
                    help="Table with one row for each comparison to be made, a column each for the reference and query groups (set with '--reference-column' and '--query-column', respectively) and optionally a (short) descriptive name for the contrast/comparison to be used in plotting, filenames etc (set with '--comparison-name-column'). Group names listed in the reference and query columns have to match the ones extracted from the column/sample names of '--count-table' or, if provided, those specified in '--annotation'. If provided, only the comparisons listed in this table are analyzed. If not provided, all sample groups are compared with each other. Specify '--has-header-comparisons' if table includes a header line. Default: NULL.",
                    metavar="tsv"
                ),
                make_option(
                    "--has-header-comparisons",
                    action="store_true",
                    default=FALSE,
                    help="Indicates whether the comparison table includes a header line. Default: FALSE."
                ),
                make_option(
                    "--reference-column",
                    action="store",
                    type="integer",
                    default=1,
                    help="Comparison table column/field number containing sample group identifiers matching those extracted from the column/sample names of '--count-table' or, if provided, specified in '--annotation'. Groups specified in this field/column serve as references in the comparisons, with 'log fold changes = expression of query / expression of reference'. Ignored if '--comparisons' is not specified. Default: 1.",
                    metavar="int"
                ),
                make_option(
                    "--query-column",
                    action="store",
                    type="integer",
                    default=2,
                    help="Comparison table column/field number containing sample group identifiers matching those extracted from the column/sample names of '--count-table' or, if provided, specified in '--annotation'. Groups specified in this field/column serve as queries in the comparisons, i.e. 'log fold changes = expression of query / expression of reference'. Ignored if '--comparisons' is not specified. Default: 2.",
                    metavar="int"
                ),
                make_option(
                    "--comparison-name-column",
                    action="store",
                    type="integer",
                    default=NULL,
                    help="Optional comparison table column/field containing short descriptive names for each comparison/contrast to be considered. If provided, these short names are used for plotting and (after removal of whitespace etc.) in output filenames. If not provided, descriptive names are generated from the query and reference group identifiers. Default: NULL.",
                    metavar="int"
                ),
                make_option(
                    "--minimum-sample-number-per-comparison",
                    action="store",
                    type="integer",
                    default=3,
                    help="Minimum number of samples required for a comparison to be considered. Be aware that results obtained from analyses with no or very few replicates may not be trustworthy. Default: 3.",
                    metavar="int"
                ),
                make_option(
                    "--minimum-sample-number-per-endpoint",
                    action="store",
                    type="integer",
                    default=0.1,
                    help="Minimum number of samples per group/comparison endpoint required for a comparison to be considered. Be aware that results obtained from analyses with no or very few replicates may not be trustworthy. Default: 1.",
                    metavar="int"
                ),
                make_option(
                    "--common-dispersion",
                    action="store",
                    type="numeric",
                    default=0.1,
                    help="Common dispersion value to be used when no replicates are available for both reference and query (i.e. total sample number n = 2). Default: 0.1.",
                    metavar="FLOAT"
                    ),
                make_option(
                    "--minimum-cpm-global",
                    action="store",
                    type="numeric",
                    default=0,
                    help="For all comparisons, consider only features/genes with the specified minimum read count per million reads (CPM) across at least '--samples-with-minimum-cpm-global'. Set to 0 to disable filtering. Default: 0.",
                    metavar="float"
                ),
                make_option(
                    "--samples-with-minimum-cpm-global",
                    action="store",
                    type="numeric",
                    default=0,
                    help="Number (integer >1) or fraction (float <=1) of samples that have to have at least '--minimum-cpm-global' for a feature/gene to be considered for any comparison. Set to 0 to disable filtering. Default: 0.",
                    metavar="float|int"
                ),
                make_option(
                    "--minimum-cpm-per-comparison",
                    action="store",
                    type="numeric",
                    default=0,
                    help="For a given comparison, consider only features/genes with the specified minimum read count per million reads (CPM) across at least '--samples-with-minimum-cpm-per-comparison'. Set to 0 to disable filtering. Default: 0.",
                    metavar="float"
                ),
                make_option(
                    "--samples-with-minimum-cpm-per-comparison",
                    action="store",
                    type="numeric",
                    default=0,
                    help="Number (integer >1) or fraction (float <=1) of samples that have to have at least '--minimum-cpm-per-comparison' for a feature/gene to be considered for a given comparison. Set to 0 to disable filtering. Default: 0.",
                    metavar="float|int"
                ),
                make_option(
                    "--p-adjust-method",
                    action="store",
                    type="character",
                    default="BH",
                    help="Method for multiple testing correction of P values. One of 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', or 'none'. Default: 'BH'.",
                    metavar="string"
                ),
                make_option(
                    "--adjusted-p-value-threshold",
                    action="store",
                    type="numeric",
                    default=0.05,
                    help="Adjusted P value threshold below which features/genes are considered as differentially expressed (given that the '--log-fold-change-threshold' is also met). Default: 0.05.",
                    metavar="float"
                ),
                make_option(
                    "--log-fold-change-threshold",
                    action="store",
                    type="numeric",
                    default=0,
                    help="Log fold change threshold above which features/genes are considered as differentially expressed (given that the '--adjusted-p-value-threshold' is also met). Default: 0.",
                    metavar="float"
                ),
                make_option(
                    "--mds-plot-gene-number",
                    action="store",
                    type="numeric",
                    default=500,
                    help="Number (integer >1) or fraction (float <=1) of genes/features to be considered for generating the multidimensional scaling plots for each comparison. Default: 500.",
                    metavar="float|int"
                ),
                make_option(
                    "--mds-plot-symbol-expansion",
                    action="store",
                    type="numeric",
                    default=1.8,
                    help="Expansion factor for MDS plot symbols. Default: 1.8.",
                    metavar="float"
                ),
                make_option(
                    "--plot-width",
                    action="store",
                    type="integer",
                    default=17,
                    help="Plot width in inches. Default: 17.",
                    metavar="int"
                ),
                make_option(
                    "--plot-width-legend",
                    action="store",
                    type="integer",
                    default=3,
                    help="Width of plot legend in inches. Default: 3.",
                    metavar="int"
                ),
                make_option(
                    "--plot-height",
                    action="store",
                    type="numeric",
                    default=14.5,
                    help="Plot height in inches. Default: 14.5.",
                    metavar="float"
                ),
                make_option(
                    "--plot-height-title",
                    action="store",
                    type="numeric",
                    default=0.5,
                    help="Height of plot title in inches. Default: 0.5.",
                    metavar="float"
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
opt_parser <- OptionParser(usage=paste("Usage:", script, "[OPTIONS] --count-table <TSV>\n", sep=" "), option_list = option_list, add_help_option=FALSE, description=msg)
opt <- parse_args(opt_parser)

# Re-assign variables
cnts <- opt[["count-table"]]
out.dir.root <- opt[["output-directory"]]
run.id <- opt[["run-id"]]
excl <- opt[["exclude"]]
id.sep <- opt[["id-field-separator"]]
id.flds <- opt[["group-id-fields"]]
anno <- opt[["annotation"]]
anno.header <- opt[["has-header-annotation"]]
id.col <- opt[["sample-id-column"]]
grp.col <- opt[["group-column"]]
comp <- opt[["comparisons"]]
comp.header <- opt[["has-header-comparisons"]]
ref.col <- opt[["reference-column"]]
query.col <- opt[["query-column"]]
name.col <- opt[["comparison-name-column"]]
cutoff.min_sample_no.total <- opt[["minimum-sample-number-per-comparison"]]
cutoff.min_sample_no.endpoint <- opt[["minimum-sample-number-per-endpoint"]]
common.disp <- opt[["common-dispersion"]]
cutoff.cpm.global <- opt[["minimum-cpm-global"]]
cutoff.sample_no.global <- opt[["samples-with-minimum-cpm-global"]]
cutoff.cpm.comp <- opt[["minimum-cpm-per-comparison"]]
cutoff.sample_no.comp <- opt[["samples-with-minimum-cpm-per-comparison"]]
method.p.adj <- opt[["p-adjust-method"]]
cutoff.p.adj <- opt[["adjusted-p-value-threshold"]]
cutoff.abs.log.fc <- opt[["log-fold-change-threshold"]]
mds.gene.no <- opt[["mds-plot-gene-number"]]
mds.sym.cex <- opt[["mds-plot-symbol-expansion"]]
plot.width <- opt[["plot-width"]]
plot.width.legend <- opt[["plot-width-legend"]]
plot.height <- opt[["plot-height"]]
plot.height.title <- opt[["plot-height-title"]]
verb <- opt[["verbose"]]

# Validate required arguments
if ( is.null(cnts) ) {
        print_help(opt_parser)
        stop("[ERROR] Required argument missing! Aborted.")
}


###################
###  FUNCTIONS  ###
###################

# Get group identifiers from sample identifiers
get.groups.from.sample.ids <- function(x, split=".", fields=1) {
    return(sapply(lapply(strsplit(x, split=split, fixed=TRUE), "[", fields), paste, collapse=split))
}

# Returns a vector of integers from a Bash 'cut'-like field selection string containing commas and dashes
get.field.selection.vector <- function(x) {

    # Remove whitespace
    x <- gsub("\\s", "", x)

    # Return NULL if string is not of required format
    if ( ! grepl("^\\d+([,-]\\d+)*$", x) ) return(NULL)

    # Split string by special characters
    x <- strsplit(unlist(strsplit(x, ",")), "-")

    # Return NULL if multiple dashes are used in one expression
    if ( any(sapply(x, length) > 2) ) return(NULL)

    # Iterate over terms
    x <- lapply(x, function(term) {

        # Convert to integer
        term <- as.integer(term)

        # Return integer or sequence of integers
        if ( length(term) == 1 ) return(term) else return(seq(term[1], term[2]))

    })

    # Obtain vector of sorted unique integers
    x <- unique(unlist(x))

    # Return vector
    return(x)

}

# Get comparison identifiers from group names
get.comparison.identifier.from.group.names <- function(query, ref, join="over", sep=".") {
    return(paste(query, join, ref, sep=sep))
}

# Plot MDS, BCV and smear plot
plot.dgea <- function(
    dge = NULL,
    dge.ex = NULL,
    title = NULL,
    out.dir = ".",
    prefix = "plots",
    suffix = "svg",
    sep = ".",
    mds.gene.no = 500,
    width = 17,
    width.legend = 3,
    height = 14.5,
    height.title = 0.5,
    cex.sym = 2,
    legend.pos = "top"
) {

    # Initialize results container
    r <- list()

    # Set parameters for multidimensional scaling plots
    col <- c("black", "red")[as.integer(dge$samples$group)]
    pch <- c(0, 5)[as.integer(dge$samples$group)]
    if ( mds.gene.no <= 1 ) mds.gene.no <- ceiling(mds.gene.no * nrow(dge))

    # Build output filename
    out.file <- file.path(out.dir, paste(prefix, suffix, sep=sep))

    # Open graphics device
    svg(out.file, width=width, height=height)

    # Get default margins
    mar <- mar.bak <- par("mar")

    # Set layout grid
    width.plot <- (width - width.legend) / 2
    height.plot <- (height - height.title) / 2
    layout(
        matrix(
            c(
                rep(1, width),
                rep(2, width.plot),
                rep(3, width.legend),
                rep(4, width.plot),
                rep(5, width.plot),
                rep(6, width.legend),
                rep(7, width.plot)
            ),
            nrow=3,
            byrow=TRUE
        ),
        heights=c(
            height.title,
            rep(height.plot, 2)
        )
    )

    # Plot title
    par(mar=c(0, 0, 0, 0), oma=c(0, 0, 0, 0))
    plot.new()
    if ( ! is.null(title) ) text(0.5, 0.5, title, cex=1.5, font=2)

    # Generate multidimensional scaling plot, gene selection: common
    mar[4] <- 0.1
    par(mar=mar)
    r$mds.common <- plotMDS(
        dge,
        main=paste("MDS (features:", mds.gene.no, "; gene selection: common)", sep=""),
        top=mds.gene.no,
        gene.selection="common",
        col=col,
        pch=pch
    )

    # Plot legend for multidimensional scaling plots
    mar[2] <- 0.1
    par(mar=mar)
    plot(1, type="n", axes=FALSE, xlab="", ylab="")
    legend(legend.pos, legend=dge.ex$comparison, pch=unique(pch), col=unique(col), bty="n")

    # Generate multidimensional scaling plot, gene selection: pairwise
    mar[2] <- mar.bak[2]
    mar[4] <- mar.bak[4]
    par(mar=mar)
    r$mds.pairwise <- plotMDS(
        dge,
        main=paste("MDS (features:", mds.gene.no, "; gene selection: pairwise)", sep=""),
        top=mds.gene.no,
        gene.selection="pairwise",
        col=col,
        pch=pch
    )

    # Generate BCV plot
    mar[4] <- 0.1
    par(mar=mar)
    r$bcv <- plotBCV(dge, main="BCV")

    # Plot placeholder
    mar[2] <- 0.1
    par(mar=mar)
    plot.new() 

    # Generate smear plot
    mar[2] <- mar.bak[2]
    mar[4] <- mar.bak[4]
    par(mar=mar)
    r$smear <- plotSmear(dge.ex, main="Smear plot", de.tags=rownames(dge)[as.logical(dge.ex$table$diffExpr)])

    # Close graphics device
    invisible(dev.off())

    # Return results container
    return(r)

}


##############
###  MAIN  ###
##############

# Initialize global variables
warn <- list()

# Write log
if ( verb ) cat("Creating output directory...\n", sep="'")

# Create output directory
dir.create(out.dir.root, recursive=TRUE, showWarnings=FALSE)

# Write log
if ( verb ) cat("Importing data...\n", sep="'")

# Import sample data
cnts <- read.delim(cnts, stringsAsFactors=FALSE, row.names=1)

# Import sample exclude list
if ( ! is.null(excl) ) excl <- readLines(excl)

# Import annotation table
if ( ! is.null(anno) ) anno <- read.delim(anno, header=anno.header, stringsAsFactors=FALSE)

# Import comparison table
if ( ! is.null(comp) ) comp <- read.delim(comp, header=comp.header, stringsAsFactors=FALSE)

# Write log
if ( verb ) cat("Processing data...\n", sep="")

# Filter data for samples to be excluded
if ( ! is.null(excl) ) cnts <- cnts[, ! colnames(cnts) %in% excl]

# Get group identifiers
if ( ! is.null(anno) ) {
    if ( ncol(anno) < max(id.col, grp.col) ) {
        print_help(opt_parser)
        stop("[ERROR] Annotation table misses required columns! Aborted.")
    }
    cnts <- cnts[, colnames(cnts) %in% anno[, id.col]]
    anno <- anno[anno[, id.col] %in% colnames(cnts), ]
    anno <- anno[match(colnames(cnts), anno[, id.col]), ]
    grps <- anno[, grp.col]
} else {
    grps <- get.groups.from.sample.ids(colnames(cnts), split=id.sep, fields=id.flds)
}

# Get comparisons
if ( ! is.null(comp) ) {
    if ( ncol(comp) < max(ref.col, query.col, name.col) ) {
        print_help(opt_parser)
        stop("[ERROR] Annotation table misses required columns! Aborted.")
    }
    nms <- if ( ! is.null(name.col) ) comp[, name.col] else get.comparison.identifier.from.group.names(query=comp[, query.col], ref=comp[, ref.col])
    filter.group_missing <- as.logical(comp[, ref.col] %in% grps & comp[, query.col] %in% grps)
    warn$group_missing <- mapply(function(name, keep) {
        if ( ! keep ) cat("[WARNING] One or both groups missing for comparison ", name, ". Skipped.\n", sep="'")
    }, nms, filter.group_missing)
    comp <- comp[filter.group_missing, ]
    nms <- nms[filter.group_missing]
    comp <- t(comp[, c(ref.col, query.col)])
} else {
    comp <- combn(unique(grps), 2)
    comp <- cbind(comp, apply(comp, 2, rev))
    nms <- get.comparison.identifier.from.group.names(query=comp[2, ], ref=comp[1, ])
}

# Discard comparisons with too few samples
n.ref <- table(grps)[match(comp[1, ], names(table(grps)))]
n.query <- table(grps)[match(comp[2, ], names(table(grps)))]
n.comp <- rowSums(data.frame(n.ref, n.query))
filter.too_few_samples <- as.logical(n.comp >= cutoff.min_sample_no.total & n.ref >= cutoff.min_sample_no.endpoint & n.query >= cutoff.min_sample_no.endpoint)
warn$too_few_samples <- mapply(function(name, keep) {
    if ( ! keep ) cat("[WARNING] Not enough samples available for comparison ", name, ". Skipped.\n", sep="'")
}, nms, filter.too_few_samples)
comp <- as.data.frame(comp[, filter.too_few_samples], stringsAsFactors=FALSE)
nms <- nms[filter.too_few_samples]

# Generate DGEList object
dge <- DGEList(counts=cnts, group=grps)

# Consider only features that are expressed (CPM > x) in at least n samples
if ( cutoff.cpm.global > 0 | cutoff.sample_no.global > 0 ) {
    if ( cutoff.sample_no.global < 1 ) cutoff.sample_no.global <- floor(ncol(dge) * cutoff.sample_no.global)
    dge <- dge[rowSums(cpm(dge) > cutoff.cpm.global) >= cutoff.sample_no.global, , keep.lib.sizes=FALSE]
}

# Write log
if ( verb ) cat("Starting differential expression analyses...\n", sep="")

# Iterate over contrasts
results <- mapply(function(contrast, name) {

    # Assign fields to individual variables
    ref <- contrast[1]
    query <- contrast[2]

    # Generate safe and "pretty" versions of name
    name.safe <- gsub("_+", "_", gsub(":", "_", gsub("_*\\((.*)\\)", ".\\1", gsub(",", "", gsub(" ", "_", name)))))
    name.pretty <- gsub(".", " ", gsub("_+", " ", name), fixed=TRUE)

    # Write log
    cat("Processing comparison ", name.pretty, "...\n", sep="'")

    # Initialize results container
    r <- list()

    # Subset data for current comparison
    r$dge <- dge[, dge$samples$group %in% c(ref, query), keep.lib.sizes=TRUE]

    # Filter expressed features in current comparison
    if ( cutoff.cpm.comp > 0 | cutoff.sample_no.comp > 0 ) {
        if ( cutoff.sample_no.comp < 1 ) cutoff.sample_no.comp <- floor(ncol(r$dge) * cutoff.sample_no.comp)
        r$dge <- r$dge[rowSums(cpm(r$dge) > cutoff.cpm.comp) >= cutoff.sample_no.comp, , keep.lib.sizes=FALSE]
    }

    # Ensure that enough data is available
    if ( nrow(r$dge) > 1 ) {

        # Calculate normalization factors for RNA composition bias
        r$dge <- calcNormFactors(r$dge)

        # Estimate (or set) dispersion
        if ( length(r$dge$samples$group) > 2 ) r$dge <- estimateDisp(r$dge) else r$dge$common.dispersion <- common.disp

        # Calculate fold changes and P values
        r$dge.ex <- exactTest(r$dge, pair=c(ref, query))

        # Add false discovery rate
        r$dge.ex$table$FDR <- p.adjust(p=r$dge.ex$table$PValue, method=method.p.adj)

        # Add differential expression flag
        r$dge.ex$table$diffExpr <- as.numeric(decideTestsDGE(r$dge.ex, adjust.method=method.p.adj, p.value=cutoff.p.adj, lfc=cutoff.abs.log.fc))

        # Create output directory for current contrast
        out.dir <- file.path(out.dir.root, name.safe)
        dir.create(out.dir, recursive=TRUE, showWarnings=FALSE)

        # Get output filename prefix
        out.prefix <- file.path(out.dir, paste(run.id, name.safe, sep="."))

        # Write out fold change table
        out.file <- paste(out.prefix, "fold_changes", "tsv", sep=".")
        write.table(r$dge.ex$table, out.file, row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")

        # Write out lists of all, differentially expressed, "upregulated" and "downregulated" features/genes
        out.file.all <- paste(out.prefix, "all", "ids", "tsv", sep=".")
        write.table(rownames(r$dge.ex$table), out.file.all, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
        out.file.de <- paste(out.prefix, "differentially_expressed", "ids", "tsv", sep=".")
        write.table(rownames(r$dge.ex$table[abs(r$dge.ex$table$diffExpr) == 1, ]), out.file.de, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
        out.file.up <- paste(out.prefix, "up", "ids", "tsv", sep=".")
        write.table(rownames(r$dge.ex$table[r$dge.ex$table$diffExpr == 1, ]), out.file.up, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
        out.file.down <- paste(out.prefix, "down", "ids", "tsv", sep=".")
        write.table(rownames(r$dge.ex$table[r$dge.ex$table$diffExpr == -1, ]), out.file.down, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

        # Generate MDS, BCV and smear plots
        r$plots <- plot.dgea(
            dge          = r$dge,
            dge.ex       = r$dge.ex,
            title        = name.pretty,
            out.dir      = out.dir,
            prefix       = paste(run.id, name.safe, "plots", sep="."),
            mds.gene.no  = mds.gene.no,
            cex.sym      = mds.sym.cex,
            width        = plot.width,
            width.legend = plot.width.legend,
            height       = plot.height,
            height.title = plot.height.title
        )

    } else {

        # Issue warning
        cat("[WARNING] No features in data table. Skipped.\n", sep="'")

    }

    # Return results
    return(r)

}, comp, nms)

# Build output filename for R image
out.file.rimage <- file.path(out.dir.root, paste(run.id, "RData", sep="."))

# Write log
if ( verb ) cat("Writing R session image to file ", out.file.rimage ,"...\n", sep="'")

# Save image
save.image(out.file.rimage)

# Write log
cat("Done.\n")
