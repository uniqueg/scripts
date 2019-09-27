#!/usr/bin/env R

# Alexander Kanitz, alexander.kanitz@alumni.ethz.ch
# Biozentrum, University of Basel


################
## PARAMETERS ##
################

# Working directory
workingDir <- "/scicore/home/zavolan/kanitz/PROJECTS/SpliceFactorsReprogramming"

# Flags
plot <- TRUE
log <- TRUE

# Sample/group names
names <- list()
names$reference <- "NotInfected"
names$query <- "Esrp2"
names$comparison <- paste(names$query, names$reference, sep="_vs_")

# Output directory
outDir <- file.path("analyzedData/kallistoQuant_2.analysis", names$comparison)
dir.create(outDir, showWarnings = FALSE)

# Input filenames
files <- list()
files$abundances <- "analyzedData/kallistoQuant_2/kallistoQuant.tpm"
files$t2g <- "publicResources/GRCm38_75/Mus_musculus.GRCm38.75.no_patches.lookup_table.trxID_geneID.tab"
files$id2name <- "publicResources/GRCm38_75/Mus_musculus.GRCm38.75.no_patches.lookup_table.geneID_geneName.tab"

# Columns to process/use as reference
cols <- list()
cols$reference <- c(5, 15, 20, 25, 20, 35, 10)
cols$query <- c(2, 12, 17, 22, 27, 32, 7)
cols$commonRef <- 1

# Axis labels
axisLabels <- list()
axisLabels$x <- "Day"
axisLabels$yAbund <- "Abundance (TPM)"
axisLabels$yFldCh <- paste("Log2 fold change (", names$query, " - ", names$reference, ")", sep="")
axisLabels$yFract <- "Fraction of total gene expression"
axisLabels$yFrcDf <- "Difference in fraction of total gene expression"
axisLabels$yFrcDf <- paste("Difference in fraction of total gene expression (", names$query, " - ", names$reference, ")", sep="")
axisLabels$yIntAbund <- "Abundance relative to T-1 (TPM)"
axisLabels$yIntFldCh <- paste("Log2 fold change relative to T-1 (", names$query, " - ", names$reference, ")", sep="")
axisLabels$yIntFract <- "Fraction of total gene expression relative to T-1"
axisLabels$yIntFrcDf <- paste("Difference in fraction of total gene expression relative to T-1 (", names$query, " - ", names$reference, ")", sep="")
axisLabels$yNullAbund <- "Abundance relative to T0 (TPM)"
axisLabels$yNullFract <- "Fraction of total gene expression relative to T0"
axisLabels$xMarks <- c(0,1,2,3,4,5,15)
axisLabels$xMarksAlt <- c(1,2,3,4,5,15)

# Output filename prefixes
prefixes <- list()
prefixes$abund <- file.path(outDir, paste("abund", names$comparison, sep="."))
prefixes$fldCh <- file.path(outDir, paste("fldCh", names$comparison, sep="."))
prefixes$fract <- file.path(outDir, paste("fract", names$comparison, sep="."))
prefixes$frcDf <- file.path(outDir, paste("frcDf", names$comparison, sep="."))
prefixes$intAbund <- file.path(outDir, paste("intervalAbund", names$comparison, sep="."))
prefixes$intFldCh <- file.path(outDir, paste("intervalFldCh", names$comparison, sep="."))
prefixes$intFract <- file.path(outDir, paste("intervalFract", names$comparison, sep="."))
prefixes$intFrcDf <- file.path(outDir, paste("intervalFrcDf", names$comparison, sep="."))
prefixes$nullAbund <- file.path(outDir, paste("nullAbund", names$comparison, sep="."))
prefixes$nullFract <- file.path(outDir, paste("nullFract", names$comparison, sep="."))
prefixes$tables <- file.path(outDir, paste("tables", names$comparison, sep="."))
prefixes$session <- file.path(outDir, paste("session.R", names$comparison, sep="."))


###############
## FUNCTIONS ##
###############
loadData <- function(files, log=TRUE) {

    # Print log message
    if (log) print("Loading data...")

    # Initialize results object
    results <- list()

    # Load transcript -> gene ID lookup table
    t2g <- read.delim(files$t2g, col.names=c("trx_id", "gene_id"), stringsAsFactors=FALSE)
    results$t2g <- t2g[order(t2g$trx_id),]

    # Load gene ID -> gene name lookup table
    results$id2name <- read.delim(files$id2name, col.names=c("gene_id", "gene_name"), stringsAsFactors=FALSE)

    # Load & process transcript abundances
    results$df <- read.delim(files$abundances, stringsAsFactors=FALSE)

    # Return results object
    return(results)

}
# ---------- #
processData <- function(obj, log=TRUE) {

    # Print log message
    if (log) print("Pre-prrocessing data...")

    # Initialize results object
    results <- list()

    # Process transcript abundances
    obj$df <- merge(obj$t2g, obj$df, by.x=1, by.y=0)
    obj$df <- merge(obj$id2name, obj$df)
    results$trx <- list(count=length(obj$df$trx_id), trx_id=obj$df$trx_id, gene_id=obj$df$gene_id, gene_name=obj$df$gene_name, abundances=list(values=obj$df[,4: ncol(obj$df)]))

    # Aggregate & process gene abundances
    obj$df <- aggregate(results$trx$abundances$values, list(results$trx$gene_id), sum)
    obj$df <- merge(obj$id2name, obj$df, by=1)
    results$gen <- list(count=length(obj$df[,1]), gene_id=obj$df[,1], gene_name=obj$df[,2], abundances=list(values=obj$df[, 3:ncol(obj$df)]))

    # Get transcript fractions
    results$trx$fractions <- list(values=results$trx$abundances$values / merge(results$trx$gene_id, cbind(results$gen$gene_id, results$gen$abundances$values), by=1)[,-1])

    # Return results object
    return(results)

}
# ---------- #
selectSamples <- function(data, cols, names, log=TRUE) {

    # Print log message
    if (log) print("Selecting samples...")

    # Transcript abundances
    data$trx$abundances$reference <- data$trx$abundances$values[, cols$reference]
    data$trx$abundances$query <- data$trx$abundances$values[, cols$query]

    # Gene abundances
    data$gen$abundances$reference <- data$gen$abundances$values[, cols$reference]
    data$gen$abundances$query <- data$gen$abundances$values[, cols$query]

    # Transcript fractions
    data$trx$fractions$reference <- data$trx$fractions$values[, cols$reference]
    data$trx$fractions$query <- data$trx$fractions$values[, cols$query]

    # Sample information
    data$reference <- list(name=names$reference, length=ncol(data$trx$abundances$reference))
    data$query <- list(name=names$query, length=ncol(data$trx$abundances$query))
    data$comparison <- names$comparison

    # Return results object
    return(data)

}
# ---------- #
getCommonRefData <-function(data, col=1, log=TRUE) {

    # Print log message
    if (log) print("Calculating data using common reference...")

    # Get common reference from reference sample
    data$trx$abundances$null <- data$trx$abundances$reference[,col]
    data$gen$abundances$null <- data$gen$abundances$reference[,col]
    data$trx$fractions$null <- data$trx$fractions$reference[,col]

    # Transcript abundances
    data$trx$abundances$nullReference <- log2(data$trx$abundances$reference[,-col] / data$trx$abundances$null)
    data$trx$abundances$nullQuery <- log2(data$trx$abundances$query[,-col] / data$trx$abundances$null)

    # Gene abundances
    data$gen$abundances$nullReference <- log2(data$gen$abundances$reference[,-col] / data$gen$abundances$null)
    data$gen$abundances$nullQuery <- log2(data$gen$abundances$query[,-col] / data$gen$abundances$null)

    # Transcript fractions
    data$trx$fractions$nullReference <- data$trx$fractions$reference[,-col] - data$trx$fractions$null
    data$trx$fractions$nullQuery <- data$trx$fractions$query[,-col] - data$trx$fractions$null

    # Return results object
    return(data)

}
# ---------- #
getPreviousRefData <- function(data, log=TRUE) {

    # Print log message
    if (log) print("Calculating data using previous data point as reference...")

    # Transcript abundances
    data$trx$abundances$intervalReference <- cbind(data$trx$abundances$nullReference[, 1, drop=FALSE], data$trx$abundances$nullReference[, -1] - data$trx$abundances$nullReference[, -ncol(data$trx$abundances$nullReference)])
    data$trx$abundances$intervalQuery <- cbind(data$trx$abundances$nullQuery[, 1, drop=FALSE], data$trx$abundances$nullQuery[, -1] - data$trx$abundances$nullQuery[, -ncol(data$trx$abundances$nullQuery)])

    # Gene abundances
    data$gen$abundances$intervalReference <- cbind(data$gen$abundances$nullReference[, 1, drop=FALSE], data$gen$abundances$nullReference[, -1] - data$gen$abundances$nullReference[, -ncol(data$gen$abundances$nullReference)])
    data$gen$abundances$intervalQuery <- cbind(data$gen$abundances$nullQuery[, 1, drop=FALSE], data$gen$abundances$nullQuery[, -1] - data$gen$abundances$nullQuery[, -ncol(data$gen$abundances$nullQuery)])

    # Transcript fractions
    data$trx$fractions$intervalReference <- cbind(data$trx$fractions$nullReference[, 1, drop=FALSE], data$trx$fractions$nullReference[, 2:ncol(data$trx$fractions$nullReference)] - data$trx$fractions$nullReference[, 1:(ncol(data$trx$fractions$nullReference) - 1)])
    data$trx$fractions$intervalQuery <- cbind(data$trx$fractions$nullQuery[, 1, drop=FALSE], data$trx$fractions$nullQuery[, 2:ncol(data$trx$fractions$nullQuery)] - data$trx$fractions$nullQuery[, 1:(ncol(data$trx$fractions$nullQuery) - 1)])

    # Return results object
    return(data)

}
# ---------- #
getGroupComparisonData <- function(data, log=TRUE) {

    # Print log message
    if (log) print("Calculating data for intergroup comparisons...")

    # Get abundance fold changes: transcripts
    data$trx$abundances$foldChanges <- log2(data$trx$abundances$query / data$trx$abundances$reference)
    data$trx$abundances$intervalFoldChanges <- data$trx$abundances$intervalQuery - data$trx$abundances$intervalReference

# Get abundance fold changes: genes
    data$gen$abundances$foldChanges <- log2(data$gen$abundances$query / data$gen$abundances$reference)
    data$gen$abundances$intervalFoldChanges <- data$gen$abundances$intervalQuery - data$gen$abundances$intervalReference

    # Get fraction differences
    data$trx$fractions$differences <- data$trx$fractions$query - data$trx$fractions$reference
    data$trx$fractions$intervalDifferences <- data$trx$fractions$intervalQuery - data$trx$fractions$intervalReference

    # Return results object
    return(data)

}
# ---------- #
writeOutputTables <- function(obj, prefix, log=TRUE) {

    # Print log message
    if (log) print("Writing output tables...")

    # Transcript abundances
    for (i in 1:length(obj$trx$abundances)) {
        outFile <- paste(prefix, "trx", "abundances", names(obj$trx$abundances)[[i]], "csv", sep=".")
        df <- cbind(trx_id=obj$trx$trx_id, gene_id=obj$trx$gene_id, gene_name=obj$trx$gene_name, obj$trx$abundances[[i]])
        write.table(df, outFile, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
    }

    # Gene abundances
    for (i in 1:length(obj$gen$abundances)) {
        outFile <- paste(prefix, "gen", "abundances", names(obj$gen$abundances)[[i]], "csv", sep=".")
        df <- cbind(gene_id=obj$gen$gene_id, gene_name=obj$gen$gene_name, obj$gen$abundances[[i]])
        write.table(df, outFile, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
    }

    # Transcript fractions
    for (i in 1:length(obj$trx$fractions)) {
        outFile <- paste(prefix, "trx", "fractions", names(obj$trx$fractions)[[i]], "csv", sep=".")
        df <- cbind(trx_id=obj$trx$trx_id, gene_id=obj$trx$gene_id, gene_name=obj$trx$gene_name, obj$trx$fractions[[i]])
        write.table(df, outFile, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
    }

}
# ---------- #
plotAbund <- function(obj, ref="none", prefix, xlab, ylab, xmarks, log=TRUE) {

    # Check reference type and set parameters accordingly
    if ( ref == "none" ) {

        logMsg <- "Plotting abundances for gene"
        refGAll <- obj$gen$abundances$reference
        queryGAll <- obj$gen$abundances$query
        refTAll <- obj$trx$abundances$reference
        queryTAll <- obj$trx$abundances$query
        yHigh <- 1
        yLow <- 0

    } else if ( ref == "common" ) {

        logMsg <- "Plotting abundances (common reference) for gene"
        refGAll <- obj$gen$abundances$nullReference
        queryGAll <- obj$gen$abundances$nullQuery
        refTAll <- obj$trx$abundances$nullReference
        queryTAll <- obj$trx$abundances$nullQuery
        yHigh <- 1
        yLow <- -1

    } else if ( ref == "previous" ) {

        logMsg <- "Plotting abundances (previous data as reference) for gene"
        refGAll <- obj$gen$abundances$intervalReference
        queryGAll <- obj$gen$abundances$intervalQuery
        refTAll <- obj$trx$abundances$intervalReference
        queryTAll <- obj$trx$abundances$intervalQuery
        yHigh <- 1
        yLow <- -1

    } else {

        print("[WARNING] Unknown reference type. Nothing is plotted.")
        return(NULL)

    }

    # Get generic x values
    x <- 1:length(xmarks)

    # Iterate over genes
    for (i in 1:obj$gen$count) {

        # Get gene ID and name
        id <- obj$gen$gene_id[[i]]
        name <- obj$gen$gene_name[[i]]

        # Print log message
        if (log) print(paste(logMsg, " ", name, " (", id, ")...", sep=""))

        # Get plotting data: gene
        refG <- as.matrix(refGAll[i, , drop=FALSE])
        queryG <- as.matrix(queryGAll[i, , drop=FALSE])
        refG[is.infinite(refG)] <- NA
        queryG[is.infinite(queryG)] <- NA

        # Get plotting data: transcripts
        idxT <- grep(id, obj$trx$gene_id)
        refT <- as.matrix(refTAll[idxT, , drop=FALSE])
        queryT <- as.matrix(queryTAll[idxT, , drop=FALSE])
        refT[is.infinite(refT)] <- NA
        queryT[is.infinite(queryT)] <- NA

        # Get number and IDs of transcripts
        numT <- length(idxT)
        idsT <- obj$trx$trx_id[idxT]

        # Get y-axis limits
        maxY <- max(refG, refT, queryG, queryT, yHigh, na.rm=TRUE)
        minY <- min(refG, refT, queryG, queryT, yLow, na.rm=TRUE)

        # Build output filename
        file <- paste(prefix, paste(name, id, sep="_"), "pdf", sep=".")

        # Open graphics device
        pdf(file)

        # Plot canvas, frame/box & axes
        plot(x, rep(0, length(x)), type="n", xlim=c(min(x), max(x)*1.3), ylim=c(minY, maxY), main=paste(id, name, sep=" | "), xlab=xlab, ylab=ylab, xaxt="n")
        axis(1, at=x, labels=xmarks)

        # Plot transcript data
        for (j in 1:numT) {
            lines(x, refT[j, ], type="b", col=j+1, lty=2)
            lines(x, queryT[j, ], type="b", col=j+1, lty=1)
        }

        # Plot gene data
        lines(x, refG, type="b", col=1, lty=2)
        lines(x, queryG, type="b", col=1, lty=1)

        # Plot legend
        legend("bottomright", legend=c(id, idsT), col=c(1,1:numT+1), lty=c(1,rep(1, numT)), bty="n", cex=0.55)
        legend("topright", legend=c(obj$query$name, obj$reference$name), col=c(1, 1), lty=c(1,2), box.col=0, bg=0, cex=0.75)

        # Close graphics device
        dev.off()

    }
}
# ---------- #
plotFldCh <- function(obj, ref="none", prefix, xlab, ylab, xmarks, log=TRUE) {

    # Check reference type and set parameters accordingly
    if ( ref == "none" ) {

        logMsg <- "Plotting fold changes for gene"
        fcGAll <- obj$gen$abundances$foldChanges
        fcTAll <- obj$trx$abundances$foldChanges

    } else if ( ref == "previous" ) {

        logMsg <- "Plotting fold changes (previous data as reference) for gene"
        fcGAll <- obj$gen$abundances$intervalFoldChanges
        fcTAll <- obj$trx$abundances$intervalFoldChanges

    } else {

        print("[WARNING] Unknown reference type. Nothing is plotted.")
        return(NULL)

    }

    # Get generic x values
    x <- 1:length(xmarks)

    # Iterate over genes
    for (i in 1:obj$gen$count) {

        # Get gene ID and name
        id <- obj$gen$gene_id[[i]]
        name <- obj$gen$gene_name[[i]]

        # Print log message
        if (log) print(paste(logMsg, " ", name, " (", id, ")...", sep=""))

        # Get plotting data: gene
        fcG <- as.matrix(fcGAll[i, , drop=FALSE])
        fcG[is.infinite(fcG)] <- NA

        # Get plotting data: transcripts
        idxT <- grep(id, obj$trx$gene_id)
        fcT <- as.matrix(fcTAll[idxT, , drop=FALSE])
        fcT[is.infinite(fcT)] <- NA

        # Get number and IDs of transcripts
        numT <- length(idxT)
        idsT <- obj$trx$trx_id[idxT]

        # Get y-axis limits
        minY <- min(fcG, fcT, -1, na.rm=TRUE)
        maxY <- max(fcG, fcT, 1, na.rm=TRUE)

        # Build output filename
        file <- paste(prefix, paste(name, id, sep="_"), "pdf", sep=".")

        # Open graphics device
        pdf(file)

        # Plot canvas, frame/box & axes
        plot(x, rep(0, length(x)), type="n", xlim=c(min(x), max(x)*1.3), ylim=c(minY, maxY), main=paste(id, name, sep=" | "), xlab=xlab, ylab=ylab, xaxt="n")
        axis(1, at=x, labels=xmarks)

        # Plot transcript data
        for (j in 1:numT) {
            lines(x, fcT[j, ], type="b", col=j+1, lty=1)
        }

        # Plot gene data
        lines(x, fcG, type="b", col=1, lty=1)

        # Plot legend
        legend("topright", legend=c(id, idsT), col=c(1,1:numT+1), lty=c(1,rep(1, numT)), bty="n", cex=0.55)

        # Close graphics device
        dev.off()

    }
}
# ---------- #
plotFract <- function(obj, ref="none", prefix, xlab, ylab, xmarks, log=TRUE) {

    # Check reference type and set parameters accordingly
    if ( ref == "none" ) {

        logMsg <- "Plotting transcript fractions for gene"
        refTAll <- obj$trx$fractions$reference
        queryTAll <- obj$trx$fractions$query

    } else if ( ref == "common" ) {

        logMsg <- "Plotting transcript fractions (common reference) for gene"
        refTAll <- obj$trx$fractions$nullReference
        queryTAll <- obj$trx$fractions$nullQuery

    } else if ( ref == "previous" ) {

        logMsg <- "Plotting transcript fractions (previous data as reference) for gene"
        refTAll <- obj$trx$fractions$intervalReference
        queryTAll <- obj$trx$fractions$intervalQuery

    } else {

        print("[WARNING] Unknown reference type. Nothing is plotted.")
        return(NULL)

    }

    # Get generic x values
    x <- 1:length(xmarks)

    # Iterate over genes
    for (i in 1:obj$gen$count) {

        # Get gene ID and name
        id <- obj$gen$gene_id[[i]]
        name <- obj$gen$gene_name[[i]]

        # Print log message
        if (log) print(paste(logMsg, " ", name, " (", id, ")...", sep=""))

        # Get plotting data: transcripts
        idxT <- grep(id, obj$trx$gene_id)
        refT <- as.matrix(refTAll[idxT, , drop=FALSE])
        queryT <- as.matrix(queryTAll[idxT, , drop=FALSE])
        refT[is.infinite(refT)] <- NA
        queryT[is.infinite(queryT)] <- NA

        # Get number and IDs of transcripts
        numT <- length(idxT)
        idsT <- obj$trx$trx_id[idxT]

        # Build output filename
        file <- paste(prefix, paste(name, id, sep="_"), "pdf", sep=".")

        # Open graphics device
        pdf(file)

        # Plot canvas, frame/box & axes
        plot(x, rep(0, length(x)), type="n", xlim=c(min(x), max(x)*1.3), ylim=c(0, 1), main=paste(id, name, sep=" | "), xlab=xlab, ylab=ylab, xaxt="n")
        axis(1, at=x, labels=xmarks)

        # Plot transcript data
        for (j in 1:numT) {
            lines(x, refT[j, ], type="b", col=j+1, lty=2)
            lines(x, queryT[j, ], type="b", col=j+1, lty=1)
        }

        # Plot legend
        legend("bottomright", legend=c(idsT), col=c(1:numT+1), lty=c(rep(1, numT)), bty="n", cex=0.55)
        legend("topright", legend=c(obj$query$name, obj$reference$name), col=c(1, 1), lty=c(1,2), box.col=0, bg=0, cex=0.75)

        # Close graphics device
        dev.off()

    }
}
# ---------- #
plotFrcDf <- function(obj, ref="none", prefix, xlab, ylab, xmarks, log=TRUE) {

    # Check reference type and set parameters accordingly
    if ( ref == "none" ) {

        logMsg <- "Plotting transcript fraction differences for gene"
        fdTAll <- obj$trx$fractions$differences

    } else if ( ref == "previous" ) {

        logMsg <- "Plotting transcript fraction differences (previous data as reference) for gene"
        fdTAll <- obj$trx$fractions$intervalDifferences

    } else {

        print("[WARNING] Unknown reference type. Nothing is plotted.")
        return(NULL)

    }

    # Get generic x values
    x <- 1:length(xmarks)

    # Iterate over genes
    for (i in 1:obj$gen$count) {

        # Get gene ID and name
        id <- obj$gen$gene_id[[i]]
        name <- obj$gen$gene_name[[i]]

        # Print log message
        if (log) print(paste(logMsg, " ", name, " (", id, ")...", sep=""))

        # Get plotting data: transcripts
        idxT <- grep(id, obj$trx$gene_id)
        fdT <- as.matrix(fdTAll[idxT, , drop=FALSE])
        fdT[is.infinite(fdT)] <- NA

        # Get number and IDs of transcripts
        numT <- length(idxT)
        idsT <- obj$trx$trx_id[idxT]

        # Build output filename
        file <- paste(prefix, paste(name, id, sep="_"), "pdf", sep=".")

        # Open graphics device
        pdf(file)

        # Plot canvas, frame/box & axes
        plot(x, rep(0, length(x)), type="n", xlim=c(min(x), max(x)*1.3), ylim=c(-1, 1), main=paste(id, name, sep=" | "), xlab=xlab, ylab=ylab, xaxt="n")
        axis(1, at=x, labels=xmarks)

        # Plot transcript data
        for (j in 1:numT) {
            lines(x, fdT[j, ], type="b", col=j+1, lty=1)
        }

        # Plot legend
        legend("topright", legend=c(idsT), col=c(1:numT+1), lty=c(rep(1, numT)), bty="n", cex=0.55)

        # Close graphics device
        dev.off()

    }
}
# ---------- #
plotAll <- function(obj, prefixes, labels, log=TRUE) {

    # Print log message
    if (log) print("Generating plots...")

    # Attach variables in containers to environment for ease of use
    with(c(labels, prefixes), {

        # No reference
        plotAbund(obj=obj, ref="none", prefix=abund, xlab=x, ylab=yAbund, xmarks=xMarks, log=log)
        plotFldCh(obj=obj, ref="none", prefix=fldCh, xlab=x, ylab=yFldCh, xmarks=xMarks, log=log)
        plotFract(obj=obj, ref="none", prefix=fract, xlab=x, ylab=yFract, xmarks=xMarks, log=log)
        plotFrcDf(obj=obj, ref="none", prefix=frcDf, xlab=x, ylab=yFrcDf, xmarks=xMarks, log=log)

        # Comparisons against previous
        plotAbund(obj=obj, ref="previous", prefix=intAbund, xlab=x, ylab=yIntAbund, xmarks=xMarksAlt, log=log)
        plotFldCh(obj=obj, ref="previous", prefix=intFldCh, xlab=x, ylab=yIntFldCh, xmarks=xMarksAlt, log=log)
        plotFract(obj=obj, ref="previous", prefix=intFract, xlab=x, ylab=yIntFract, xmarks=xMarksAlt, log=log)
        plotFrcDf(obj=obj, ref="previous", prefix=intFrcDf, xlab=x, ylab=yIntFrcDf, xmarks=xMarksAlt, log=log)

        # Comparisons against common reference
        plotAbund(obj=obj, ref="common", prefix=nullAbund, xlab=x, ylab=yNullAbund, xmarks=xMarksAlt, log=log)
        plotFract(obj=obj, ref="common", prefix=nullFract, xlab=x, ylab=yNullFract, xmarks=xMarksAlt, log=log)

    })

}


##########
## MAIN ##
##########

# Set working directory
setwd(workingDir)

# Load data
raw <- loadData(files=files, log=log)

# Process data
ls <- processData(obj=raw, log=log)

# Select samples
ls <- selectSamples(data=ls, cols=cols, names=names, log=log)

# Get data for comparisons against common reference
ls <- getCommonRefData(data=ls, col=cols$commonRef, log=log)

# Get data for comparisons against previous column
ls <- getPreviousRefData(data=ls, log=log)

# Calculate abundance fold changes and fraction differences
ls <- getGroupComparisonData(data=ls, log=log)

# Write output tables
writeOutputTables(obj=ls, prefix=prefixes$tables, log=log)

# Save session
save.image(prefixes$session)

# Generate plots
if (plot) plotAll(obj=ls, prefixes=prefixes, labels=axisLabels, log=log)
