#!/usr/bin/env Rscript

##############################
###  GLOBAL PARAMETERS //  ###
##############################
region_types <- c("three_prime_utr", "five_prime_utr", "CDS", "none")
def_region_types <- c("three_prime_utr", "five_prime_utr", "CDS")
reg_colors_def <- c("#33ccff", "#666666", "#ff8000")
##############################
###  // GLOBAL PARAMETERS  ###
##############################


######################
###  FUNCTIONS //  ###
######################
formatTypes <- function(chr) {
    chr <- gsub("five_prime_utr", "5'UTR", chr)
    chr <- gsub("three_prime_utr", "3'UTR", chr)
    chr <- gsub("none", "Undefined", chr)
    return(chr)
}
#-----------------------#
pieChart <- function(x, col, legend, main=NULL, cex.main=1.6, legend.pos="topright", cex.leg=1.4) {
    pie(x, clockwise=TRUE, col=col, labels=NA, main=main, cex.main=cex.main)
    legend(x=legend.pos, legend=legend, fill=col, bty="n", cex=cex.leg)
    return(NULL)
}
#-----------------------#
chiSquare <- function(x, p) {
#
    chsq <- chisq.test(x=x, p=p)
    meth <- chsq$method
    sep <- paste(rep("-", nchar(meth)), collapse="")
    stat <- paste("Statistic:", chsq$statistic)
    df <- paste("Degrees of freedom:", chsq$parameter)
    p <- paste("P-value:", chsq$p.value)
    obs <- paste("Observed:", paste(chsq$observed, collapse=", "))
    exp <- paste("Expected:", paste(chsq$expected, collapse=", "))
    res <- paste("Pearson residuals:", paste(chsq$residuals, collapse=", "))
    stdres <- paste("Standardized residuals:", paste(chsq$stdres, collapse=", "))
    chsq$print <- paste(meth, sep, stat, df, p, obs, exp, res, stdres, sep="\n")
    return(chsq)
}
#-----------------------#
categoryVectorListToBED <- function(vectorList, nms, categories, start, prefix) {
#
    grl <- GRangesList(mapply(function(name, start) {
        site <- vectorList[[name]]
        site_coll <- rle(site)
        type <- site_coll$values %in% categories
        site_coll <- setNames(site_coll$lengths[type], site_coll$values[type])
        start <- start + c(0, cumsum(site_coll[-length(site_coll)]))
        ranges <- IRanges(start=start, width=site_coll)
        seqnames <- Rle(rep(name, length(ranges)))
        gr <- GRanges(seqnames=seqnames, ranges=ranges, name=names(site_coll))
    }, nms, start))
    gr <- unlist(grl)
    names(gr) <- NULL
    outFile <- paste(prefix, "regions.bed", sep=".")
    write(paste("Writing region types to BED file ", outFile, "...", sep="'"), stdout())
    export(gr, outFile, format="bed")
    return(gr)
}
#-----------------------#
plotFormats <- function(FUN, formats, prefix) {
    if ("pdf" %in% formats) {
        outFile <- paste(prefix, "pdf", sep=".")
        write(paste("Plotting pie chart to file ", outFile, "...", sep="'"))
        pdf(outFile)
        dump <- FUN
        dev.off()
    }
#    if ("png" %in% formats) {
#        outFile <- paste(prefix, "png", sep=".")
#        write(paste("Plotting pie chart to file ", outFile, "...", sep="'"))
#        png(outFile)
#        dump <- FUN
#        dev.off()
#    }
#    if ("svg" %in% formats) {
#        outFile <- paste(prefix, "svg", sep=".")
#        write(paste("Plotting pie chart to file ", outFile, "...", sep="'"))
#        svg(outFile)
#        dump <- FUN
#        dev.off()
#    }
}
#-----------------------#
processAnnotations <- function(gtf, outDir) {

    # Load packages
    write("Loading package 'rtracklayer'...", stdout())
    library("rtracklayer")

    # Make output directory & prepare output filenames
    write("Generating output directory...", stdout())
    dir.create(outDir, recursive=TRUE, showWarnings = FALSE)
    gtfBase <- unlist(strsplit(basename(gtf), ".gtf"))[1]
    out_prefix <- file.path(outDir, gtfBase)

    # Import GTF annotations
    write(paste("Importing GTF annotation data from file ", gtf, "...", sep="'"), stdout())
    gtf <- import(gtf, format="gtf", asRangedData=FALSE)
    outFile <- paste(out_prefix, "GRanges.R", sep=".")
    write(paste("Saving GTF object in R file ", outFile, "...", sep="'"), stdout())
    save(gtf, file=outFile)

    # Subset exons and split by transcript identifier
    write("Subsetting exons...", stdout())
    exons_gr <- gtf[gtf$type == "exon"]
    mcols(exons_gr) <- list(transcript_id=factor(exons_gr$transcript_id))
    exons_grl <- split(exons_gr, exons_gr$transcript_id)

    # Subset & process region annotations (5' UTR, CDS, 3' UTR)
    write("Subsetting/processing region annotations...", stdout())
    regions_gr <- gtf[gtf$type %in% c("five_prime_utr", "CDS", "stop_codon", "three_prime_utr")]
    regions_gr$type[regions_gr$type == "stop_codon"] <- "CDS"
    mcols(regions_gr) <- list(transcript_id=factor(regions_gr$transcript_id), type=factor(regions_gr$type))
    regions_grl <- split(regions_gr, regions_gr$transcript_id)

    # Get exons with and without region annotation
    write("Identify exons with/without region annotations...", stdout())
    exons_w_regions_grl <- exons_grl[intersect(names(exons_grl), names(regions_grl))]
    exons_wo_regions_grl <- exons_grl[setdiff(names(exons_grl), names(exons_w_regions_grl))]
    exons_w_partial_regions_grl <- psetdiff(exons_w_regions_grl, regions_grl)

    # Update region annotations with type "none"
    write("Add annotation type 'none' to unavailable/undefined region annotations...", stdout())
    no_regions_gr <- c(unlist(exons_wo_regions_grl, use.names=FALSE), unlist(exons_w_partial_regions_grl, use.names=FALSE))
    mcols(no_regions_gr)$type <- factor(rep("none", length(no_regions_gr)))
    regions_gr <- sort(c(regions_gr, no_regions_gr))
    regions_grl <- split(regions_gr, regions_gr$transcript_id)

    # Get list of vectors of region types
    # List contains one vector for each transcript
    # Each vector is composed of region types for each position
    write("Obtaining region type information per nucleotide (this may take long)...", stdout())
    reg_vec_all_ls <- lapply(regions_grl, function(trx) {
        strand <- unique(strand(trx))
        if ( length(strand) != 1 | ! strand %in% c("+", "-") ) {
            write("[WARNING] Strand information unclear.", stderr())
            return(NULL)
        }
        vec <- as.character(unlist(mapply(rep, x=trx$type, each=width(trx))))
        if ( strand == "+" ) return(vec) else return(rev(vec))
    })
    outFile <- paste(out_prefix, "regionByNucleotide.R", sep=".")
    write(paste("Saving nucleotide-level composition information in file ", outFile, "...", sep="'"), stdout())
    save(reg_vec_all_ls, file=outFile)

    # Summarize region type nucleotide composition
    write("Counting nucleotides per region type...", stdout())
    reg_cts_all <- table(unlist(reg_vec_all_ls))
    reg_cts_all <- setNames(as.numeric(reg_cts_all), names(reg_cts_all))
    outFile <- paste(out_prefix, "regionCounts.R", sep=".")
    write(paste("Saving counts in file ", outFile, "...", sep="'"), stdout())
    save(reg_cts_all, file=outFile)

    # Generate BED file of regions
    names_found <- names(reg_vec_all_ls[! sapply(reg_vec_all_ls, is.null)])
    start <- rep(1, length(names_found))
    gr <- categoryVectorListToBED(reg_vec_sites_ls, names_found, region_types, start, out_prefix)

    # Generate pie chart
    reg_cts_all_def <- reg_cts_all[def_region_types]
    plotFormats(pieChart(x=reg_cts_all_def[def_region_types], col=reg_colors_def, legend=formatTypes(def_region_types), main="control"), formats=c("pdf", "png", "svg"), prefix=paste(out_prefix, "regionCounts.pie", sep="."))

    # Return list of objects
    obj_ls <- list(gtf=gtf, regions_gr=regions_gr, reg_vec_all_ls=reg_vec_all_ls, reg_cts_all=reg_cts_all, gr=gr, reg_cts_all_def=reg_cts_all_def)
    return(obj_ls)

}
#-----------------------#
processSample <- function (csv, regionPerNt, reg_cts_all, outDir) {

    # Make output directory & prepare output filenames
    write("Generating output directory...", stdout())
    dir.create(outDir, recursive=TRUE, showWarnings = FALSE)
    csvBase <- unlist(strsplit(basename(csv), ".csv"))[1]
    out_prefix <- file.path(outDir, csvBase)

    # Loading annotation data
    write(paste("Obtaining annotation R objects...", sep="'"), stdout())
    if ( mode(regionPerNt) == "character" ) {
        load(regionPerNt)
    } else {
        reg_vec_all_ls <- regionPerNt
    }
    if ( mode(regionCounts) == "character" ) {
        load(regionCounts)
    } else {
        reg_cts_all <- regionCounts
    }

    # Importing CSV file of sites
    write(paste("Importing sites from file ", csv, "...", sep="'"), stdout())
    sites <- read.delim(csv, stringsAsFactors=FALSE)

    # Get list of vectors of region types
    # List contains one vector for each transcript
    # Each vector is composed of region types for each position
    reg_vec_sites_ls <- apply(sites, 1, function(site) {
        reg_vec_all_ls[[site["seqnames"]]][site["start"]:site["end"]]
    })
    names(reg_vec_sites_ls) <- sites$seqnames
    outFile <- paste(out_prefix, "regionByNucleotide.R", sep=".")
    write(paste("Saving nucleotide-level composition information in file ", outFile, "...", sep="'"), stdout())
    save(reg_vec_sites_ls, file=outFile)

    # Summarize region type nucleotide composition
    write("Counting nucleotides per region type...", stdout())
    reg_cts_sites <- table(unlist(reg_vec_sites_ls))
    reg_cts_sites <- setNames(as.numeric(reg_cts_sites), names(reg_cts_sites))
    outFile <- paste(out_prefix, "regionCounts.R", sep=".")
    write(paste("Saving counts in file ", outFile, "...", sep="'"), stdout())
    save(reg_cts_sites, file=outFile)

    # Generate BED file of regions
    names_found <- names(reg_vec_sites_ls[! sapply(reg_vec_sites_ls, is.null)])
    start <- sites$start[sites$seqnames %in% names_found]
    gr <- categoryVectorListToBED(reg_vec_sites_ls, names_found, region_types, start, out_prefix)

    # Generate pie chart
    reg_cts_sites_def <- reg_cts_sites[def_region_types]
    plotFormats(pieChart(x=reg_cts_sites_def, col=reg_colors_def, legend=formatTypes(def_region_types), main="sample"), formats=c("pdf", "png", "svg"), prefix=paste(out_prefix, "regionCounts.pie", sep="."))

    # Run Pearson's Chi-squared test
    write("Running Pearson's Chi-squared test...", stdout())
    reg_cts_all_def <- reg_cts_all[def_region_types]
    chsq <- chiSquare(x=reg_cts_sites_def, p=reg_cts_all_def/sum(reg_cts_all_def))
    outFile <- paste(out_prefix, "chiSquare.txt", sep=".")
    write(paste("Writing Chi-square summary to file", outFile, "...", sep="'"), stdout())
    write(chsq$print, file=outFile)

    # Return list of objects
    obj_ls <- list(reg_vec_all_ls=reg_vec_all_ls, reg_cts_all=reg_cts_all, sites=sites, reg_vec_sites_ls=reg_vec_sites_ls, reg_cts_sites=reg_cts_sites, reg_cts_sites_def=reg_cts_sites_def, reg_cts_all_def=reg_cts_all_def, gr=gr, chsq=chsq)
    return(obj_ls)

}
######################
###  // FUNCTIONS  ###
######################


#################
###  MAIN //  ###
#################

# Initiate objects
obj_ls <- NULL

# Process annotations
if ( ! is.null(gtf) ) {
    annot_obj_ls <- processAnnotations(gtf, outDir)
    if ( is.null(regionPerNt) ) regionPerNt <- annot_obj_ls$reg_vec_all_ls
    if ( is.null(regionCounts) ) regionCounts <- annot_obj_ls$reg_cts_all
}

# Process sample
if ( ! any(is.null(c(csv, regionPerNt, regionCounts))) ) {
    sample_obj_ls <- processSample(csv, regionPerNt, regionCounts, outDir)
}

# Save session
outFile <- file.path(outDir, "session.R")
write(paste("Saving R session in file ", outFile, "...", sep="'"), stdout())
save.image(file=outFile)

#################
###  // MAIN  ###
#################
