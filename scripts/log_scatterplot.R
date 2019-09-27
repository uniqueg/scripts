#!/usr/bin/Rscript

#==================#
#   HEADER START   #
#==================#
### Created: Apr 01, 2014
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
#==================#
#    HEADER END    #
#==================#


#==========================#
#   PRE-REQUISITES START   #
#==========================#
#---> LOAD LIBRARIES <---#
if ( suppressWarnings(suppressPackageStartupMessages(require("optparse"))) == FALSE ) { stop("Package 'optparse' required!\nExecution aborted.") }

#---> GET SCRIPT NAME <---#
script <- sub("--file=", "", basename(commandArgs(trailingOnly=FALSE)[4]))

#---> COMMAND-LINE ARGUMENTS <---#
## List of allowed/recognized arguments
option_list <- list(
	make_option(c("-x", "--x-values"), action="store", type="character", default="", help="REQUIRED: File containing x values (format: ID (TAB) value).", metavar="file"),
	make_option(c("-y", "--y-values"), action="store", type="character", default="", help="REQUIRED: File containing y values (format see --x-values).", metavar="file"),
	make_option(c("-p", "--prefix"), action="store", type="character", default="", help="REQUIRED: Output filename prefix/basename (suffix '.svg' is added).", metavar="file"),
	make_option(c("-f", "--filter"), action="store", type="character", default="", help="File indicating which IDs to plot (format: one ID per line).", metavar="file"),
	make_option("--header", action="store_true", default=FALSE, help="Indicate if '--x-values' and '--y-values' files have a header line."), 
	make_option("--input-log", action="store_true", default=FALSE, help="Specify if input values are in log space."),
	make_option("--modify-na", action="store", type="numeric", default="0", help="Zero and unavailable '--x/y-values' are set to the specified (via --(x/y-)min) or data-derived minimum value minus this number (default: 0). This is a convenience option to set both the modifier for x and y values together.", metavar="num"),
    make_option("--x-modify-na", action="store", type="numeric", default=NA, help="Zero and unavailable '--x-values' are set to the specified (via --(x-)min) or data-derived minimum value minus this number (default: 0).", metavar="num"),
    make_option("--y-modify-na", action="store", type="numeric", default=NA, help="Zero and unavailable '--y-values' are set to the specified (via --(y-)min) or data-derived minimum value minus this number (default: 0).", metavar="num"),
	make_option("--min", action="store", type="numeric", default=NA, help="Set minimum log2 value to plot (default: derive from data, separately for the x and y values). Overruled by --x-min and --y-min. This is a convenience option to set the minima for the x and y axis together when expecting symmetrical plots.", metavar="num"),
	make_option("--max", action="store", type="numeric", default=NA, help="Set maximum log2 value to plot (default: derive from data, separately for the x and y values). Overruled by --x-max and --y-max. This is a convenience option to set the maxima for the x and y axis together when expecting symmetrical plots.", metavar="num"),		
    make_option("--x-min", action="store", type="numeric", default=NA, help="Set minimum log2 x value to plot (default: derive from data).", metavar="num"),
    make_option("--x-max", action="store", type="numeric", default=NA, help="Set maximum log2 x value to plot (default: derive from data).", metavar="num"),
    make_option("--y-min", action="store", type="numeric", default=NA, help="Set minimum log2 y value to plot (default: derive from data).", metavar="num"),
    make_option("--y-max", action="store", type="numeric", default=NA, help="Set maximum log2 y value to plot (default: derive from data).", metavar="num"),
	make_option("--title", action="store", type="character", default="", help="Plot title.", metavar="file"),
	make_option("--x-label", action="store", type="character", default="x", help="Label for x-axis.", metavar="file"),
	make_option("--y-label", action="store", type="character", default="y", help="Label for y-axis.", metavar="file"),
	make_option("--corr", action="store", type="character", default="", help="Calculate x,y correlations for 'union', 'intersection' or 'both' and print to the plot (default: no correlation printed).", metavar="type"),
    make_option("--corr-method", action="store", type="character", default="pearson", help="Method for calculating correlations. One of 'pearson', 'spearman' or 'kendall' (default: 'pearson').", metavar="type"),
	make_option("--counts", action="store_true", default=FALSE, help="Add feature counts (x values, y values, union, intersection) to the plot."),
	make_option("--diagonal", action="store_true", default=FALSE, help="Add a diagonal with slope 1 and y-intercept 0 to the plot."),
	make_option(c("-h", "--help"), action="store_true", default=FALSE, help="Show this information and die!"),
	make_option(c("-u", "--usage"), action="store_true", default=FALSE, dest="help", help="Show this information and die!"),
	make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print log messages [default: TRUE]."),
	make_option(c("-q", "--quiet"), action="store_false", dest="verbose", help="Shut up!")
)

## Description
description <- "Generates a log2 density scatterplot for the given x- and y-values.\nNote: No input file validation performed.\n"
author <- "Author: Alexander Kanitz & Andreas R. Gruber, Biozentrum, University of Basel"
version <- "Version: 1.2 (21-MAY-2014)"
requirements <- "Requires: optparse, geneplotter, RColorBrewer"
msg <- paste(description, author, version, requirements, sep="\n")

## Parse command-line arguments
opt_parser <- OptionParser(usage="Usage: %prog (OPTIONS) --x-values [FILE] --y-values [FILE] --prefix [FILE]\n", option_list = option_list, add_help_option=FALSE, description=msg)
opt <- parse_args(opt_parser)

## Die if any required arguments are missing...
if 	( opt$`x-values` == "" || opt$`y-values` == "" || opt$prefix == "" ) { 
	write("[ERROR] Required argument(s) missing!\n\n", stderr())	
	stop(print_help(opt_parser))
}

## Die if illegal --corr value
if (! opt$corr %in% c("union", "intersection", "both", "")) {
    write("[ERROR] Illegal argument to the '--corr' option! Use one of 'union', 'intersection', 'both' or '' (empty string; default).\n\n", stderr())    
    stop(print_help(opt_parser))
}

## Die if illegal --corr-method value
if (! opt$`corr-method` %in% c("pearson", "spearman", "kendall")) {
    write("[ERROR] Illegal argument to the '--corr-method' option! Use one of 'pearson' (default), 'spearman' or 'kendall'.\n\n", stderr())                             
    stop(print_help(opt_parser))
}

## Set --x-modify-na and --y-modify-na, if applicable
if ( is.na(opt$`x-modify-na`) ) opt$`x-modify-na` <- opt$`modify-na`
if ( is.na(opt$`y-modify-na`) ) opt$`y-modify-na` <- opt$`modify-na` 

## Set --x-min, --x-max, --y-min and --y-max based on --min and --max, if applicable
if ( is.na(opt$`x-min`) && ! is.na(opt$min) ) opt$`x-min` <- opt$min
if ( is.na(opt$`y-min`) && ! is.na(opt$min) ) opt$`y-min` <- opt$min
if ( is.na(opt$`x-max`) && ! is.na(opt$max) ) opt$`x-max` <- opt$max
if ( is.na(opt$`y-max`) && ! is.na(opt$max) ) opt$`y-max` <- opt$max

## Build output filename
outfile <- paste(opt$prefix, "svg", sep=".")

#==========================#
#    PRE-REQUISITES END    #
#==========================#

#================#
#   MAIN START   #
#================#

#---> START MESSAGE <---#
if ( opt$verbose ) cat("Starting '", script, "'...\n", sep="")

#---> LOAD LIBRARIES <---#
# Print status message
if ( opt$verbose ) cat("Loading libraries...\n", sep="")
# Load libraries
if ( suppressWarnings(suppressPackageStartupMessages(require("geneplotter"))) == FALSE ) { stop("Package 'geneplotter' required!\nExecution aborted.") }
if ( suppressWarnings(suppressPackageStartupMessages(require("RColorBrewer"))) == FALSE ) { stop("Package 'RColorBrewer' required!\nExecution aborted.") }

#---> LOAD X VALUES FILE <---#
# Print status message
if ( opt$verbose ) cat("Reading file '", basename(opt$`x-values`), "'...\n", sep="")
# Load data
df_x <- as.data.frame(read.table(opt$`x-values`, sep="\t", header=opt$header, row.names=1, col.names=c("id", "value", rep("NULL", max(count.fields(opt$`x-values`, sep = "\t")) - 2 )), colClasses=c("character", "numeric", rep("NULL", max(count.fields(opt$`x-values`, sep = "\t")) - 2 ))))

#---> LOAD Y VALUES FILE <---#
# Print status message
if ( opt$verbose ) cat("Reading file '", basename(opt$`y-values`), "'...\n", sep="")
# Load data
df_y <- as.data.frame(read.table(opt$`y-values`, sep="\t", header=opt$header, row.names=1, col.names=c("id", "value", rep("NULL", max(count.fields(opt$`y-values`, sep = "\t")) - 2 )), colClasses=c("character", "numeric", rep("NULL", max(count.fields(opt$`y-values`, sep = "\t")) - 2 ))))

#---> FILTER OUT NA & ZERO VALUES <---#
# Print status message
if ( opt$verbose ) cat("Filtering out zero values...\n")
# Filter data
df_x <- df_x[! is.na(df_x$value), , drop=FALSE]
df_y <- df_y[! is.na(df_y$value), , drop=FALSE]
# Check if values are in logspace
if ( ! opt$`input-log` ) {
	df_x <- df_x[df_x$value != 0, , drop=FALSE]
	df_y <- df_y[df_y$value != 0, , drop=FALSE]
}

#---> APPLY FILTER (IF APPLICABLE) <---#
if ( opt$filter != "" ) {
	# Print status message
	if ( opt$verbose ) cat("Reading filter file '", basename(opt$filter), "'...\n", sep="")
	# Load data
	filter <- read.table(opt$filter, header=FALSE, stringsAsFactors=FALSE, col.names="id", colClasses="character")
	# Print status message
	if ( opt$verbose ) cat("Filtering IDs provided in '", basename(opt$filter), "'...\n", sep="")
	# Filter data
	df_x <- df_x[row.names(df_x) %in% filter$id, , drop=FALSE]
	df_y <- df_y[row.names(df_y) %in% filter$id, , drop=FALSE]
}

#---> CONVERT TO LOGSPACE (IF APPLICABLE) <---#
# Check if values are in logspace
if ( ! opt$`input-log` ) {
	# Print status message
	if ( opt$verbose ) cat("Converting values to logspace (log2)...\n")
	# Convert to logspace
	df_x <- log2(df_x)
	df_y <- log2(df_y)
}

#---> SET MIN AND MAX VALUES <---#	
# Check if --x-min is available
if ( is.na(opt$`x-min`) ) {
	# Print status message
	if ( opt$verbose ) cat("Deriving x-axis minimum plot value from data...\n")
	# Derive minimum plot value
	opt$`x-min` <- min(df_x)
}
# Check if --y-min is available
if ( is.na(opt$`y-min`) ) {
    # Print status message
    if ( opt$verbose ) cat("Deriving y-axis minimum plot value from data...\n")
    # Derive minimum plot value
    opt$`y-min` <- min(df_y)
}
# Check if --max was specified
if ( is.na(opt$`x-max`) ) {
	# Print status message
	if ( opt$verbose ) cat("Deriving x-axis maximum plot value from data...\n")
	# Derive minimum plot value
	opt$`x-max` <- max(df_x)
}
# Check if --max was specified
if ( is.na(opt$`y-max`) ) {
    # Print status message
    if ( opt$verbose ) cat("Deriving y-axis maximum plot value from data...\n")
    # Derive minimum plot value
    opt$`y-max` <- max(df_x)
}
## Die if minimum value is greater than maximum value
if ( opt$`x-min` > opt$`x-max` ) stop("Minimum x value greater than maximum x value!\nExecution aborted.")
if ( opt$`y-min` > opt$`y-max` ) stop("Minimum y value greater than maximum y value!\nExecution aborted.")
# Calculate distance between minimum (floor) and maximum (ceiling) plot values
abs_x <- ceiling(opt$`x-max`) - floor(opt$`x-min`)
abs_y <- ceiling(opt$`y-max`) - floor(opt$`y-min`)

#---> SUBSET VALUES IN MIN / MAX RANGE <---#
# Print status message
if ( opt$verbose ) cat("Removing 'out of range' data...\n")
# Subset data in allowed range
df_x <- df_x[df_x >= opt$`x-min`, , drop=FALSE]
df_y <- df_y[df_y >= opt$`y-min`, , drop=FALSE]
df_x <- df_x[df_x <= opt$`x-max`, , drop=FALSE]
df_y <- df_y[df_y <= opt$`y-max`, , drop=FALSE]

#---> GENERATE XY UNION AND INTERSECTION <---#
# Print status message
if ( opt$verbose ) cat("Performing set operations (union and intersection) on xy values...\n")
# Set operations
df_xy_union <- merge(df_x, df_y, by=0, all=TRUE, sort=FALSE)
df_xy_intersection <- merge(df_x, df_y, by=0, sort=FALSE)

#---> SORT DATAFRAMES BY IDS <---#
# Print status message
if ( opt$verbose ) cat("Sorting union and intersection data by IDs...\n")
# Convert IDs to factors
df_xy_union$Row.names <- as.factor(df_xy_union$Row.names)
df_xy_intersection$Row.names <- as.factor(df_xy_intersection$Row.names)
# Sort data frames by IDs
df_xy_union <- df_xy_union[order(df_xy_union$Row.names), ]
df_xy_intersection <- df_xy_intersection[order(df_xy_intersection$Row.names), ]

#---> XY UNION: SET NA TO MODIFIED MIN VALUES <---#
# Print status message
if ( opt$verbose ) cat("Setting NA values in the union of x and y IDs to the adjusted (see --(x/y-)modify-na) x/y minima: ", opt$`x-min`, " and ", opt$`y-min` , "\n", sep="")
# Replace NAs
df_xy_union$value.x[is.na(df_xy_union$value.x)] <- opt$`x-min` - opt$`x-modify-na`
df_xy_union$value.y[is.na(df_xy_union$value.y)] <- opt$`y-min` - opt$`y-modify-na`

#---> COMPUTE COUNTS <---#
# Print status message
if ( opt$verbose ) cat("Counting features...\n")
# Compute counts
n_x <- nrow(df_x)
n_y <- nrow(df_y)
n_union <- nrow(df_xy_union)
n_intersection <- nrow(df_xy_intersection)

#---> SET PLOT COLORS <---#
# Print status message
if ( opt$verbose ) cat("Setting plot colors...\n")
# Set plot colors and order
colors  <- densCols(df_xy_intersection$value.x, df_xy_intersection$value.y, colramp=colorRampPalette(c("grey80","black")))
order <- rev(order(colors))

#---> PLOT <---#
# Print status message
if ( opt$verbose ) cat("Plotting to file '", outfile, "'...\n", sep="")
# Plot
svg(outfile)
par(mar = c(5.1, 5.1, 4.1, 2.1))
plot(0, type="n", xlim=c(floor(opt$`x-min`),ceiling(opt$`x-max`)), ylim=c(floor(opt$`y-min`),ceiling(opt$`y-max`)), axes=FALSE, xlab=NA, ylab=NA, main=opt$title)
box()
axis(side=1, labels=seq(floor(opt$`x-min`), ceiling(opt$`x-max`), abs_x/10), at=seq(floor(opt$`x-min`), ceiling(opt$`x-max`), abs_x/10))
axis(side=2, labels=seq(floor(opt$`y-min`), ceiling(opt$`y-max`), abs_y/10), at=seq(floor(opt$`y-min`), ceiling(opt$`y-max`), abs_y/10))
mtext(side=1, opt$`x-label`, line=3.5)
mtext(side=2, opt$`y-label`, line=2.5)
points(df_xy_intersection$value.x[order], df_xy_intersection$value.y[order], col=colors[order], pch=16)
if (opt$corr == "both" || opt $corr == "intersection") {
	text(floor(opt$`x-min`) + 0.0 * abs_x, floor(opt$`y-min`) + 0.98 * abs_y, adj=c(0,NA), paste("Ri =", signif(cor(df_xy_intersection$value.x, df_xy_intersection$value.y, method=opt$`corr-method`), digits=4)), cex=1.2)
}
if (opt$corr == "both" || opt $corr == "union") {
	text(floor(opt$`x-min`) + 0.0 * abs_x, floor(opt$`y-min`) + 0.93 * abs_y, adj=c(0,NA), paste("Ru =", signif(cor(df_xy_union$value.x, df_xy_union$value.y, method=opt$`corr-method`), digits=4)), cex=1.2)
}
if (opt$counts) {
	text(floor(opt$`x-min`) + 0.7 * abs_x, floor(opt$`y-min`) + 0.13 * abs_y, adj=c(0,NA), paste("Features x:", n_x))
	text(floor(opt$`x-min`) + 0.7 * abs_x, floor(opt$`y-min`) + 0.09 * abs_y, adj=c(0,NA), paste("Features y:", n_y))
	text(floor(opt$`x-min`) + 0.7 * abs_x, floor(opt$`y-min`) + 0.05 * abs_y, adj=c(0,NA), paste("Intersection:", n_intersection))
	text(floor(opt$`x-min`) + 0.7 * abs_x, floor(opt$`y-min`) + 0.01 * abs_y, adj=c(0,NA), paste("Union:", n_union))
}
if (opt$diagonal) {
	abline(0,1)
}
graphics.off()

#---> END MESSAGE <---#
if ( opt$verbose ) cat("Done.\n")
