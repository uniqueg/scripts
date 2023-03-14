#!/usr/bin/env Rscript

#===========#
#  HEAD //  #
#===========#

# (c) 2016, Alexander Kanitz, Biozentrum, University of Basel

#===========#
#  // HEAD  #
#===========#


#===========#
#  TODO //  #
#===========#

# Dependencies
# Check for collision of option names and argument keys
# Re-organize 'obj': gene/transcript/biotypes
# Option --include-cDNA-coordinates?
# Option --output-format:
#       BED
#       BED12
#       both
# Option --region-types:
#       GROUPS
#       ======
#       all: gene-level, transcript-level, translation-level, intergenic, annotations
#       gene-level: genes, transcripts, isoforms-gene-level, protein-coding-gene-level, biotypes-gene-level
#       transcript-level: transcripts, isoforms-transcript-level, protein-coding-transcript-level, biotypes-transcript-level
#       isoforms-gene-level: exonic, intronic, untranscribed_gene (reduced)
#       isoforms-transcript-level: exons, introns, untranscribed_transcript
#       protein-coding-gene-level: 5'-UTR_gene, CDS_gene, 3'-UTR_gene
#       protein-coding-transcript-level: 5'-UTR_transcript, CDS_transcript, 3'-UTR_transcript
#       biotypes: biotypes-gene-level, biotypes-transcript-level
#       INDIVIDUAL
#       ==========
#               intergenic
#               genes
#               transcripts*
#               exons*
#               introns*
#               untranscribed*
#       5'-UTR*
#               CDS*
#       3'-UTR*
#       transcription_start_sites*
#       cleavage_sites*
#       translation_start_sites*
#       translation_stop_sites*
#       biotypes-transcript-level
#       exons_reduced*
#       introns_reduced*
#       untranscribed_reduced*
#       5'-UTR_reduced*
#       CDS_reduced*
#       3'-UTR_reduced*
#       transcription_start_sites_reduced*
#       cleavage_sites_reduced*
#       translation_start_sites*
#       translation_stop_sites*
#       biotypes-gene-level
#*: BED12 output possible

#===========#
#  // TODO  #
#===========#


#=====================#
#  PRE-REQUISITES //  #
#=====================#


#=====================#
#  // PRE-REQUISITES  #
#=====================#


#================#
#  FUNCTIONS //  #
#================#

loadPackages <- function(packages, install = FALSE, repo = "http://cran.us.r-project.org", status = FALSE) {

# [Description]
#     Loads one or more packages and either dies when packages are unavailable or tries to install
#     them.
# [Parameters]
#     packages       : A vector of package names
#     install        : Boolean indicating whether it shall be attempted to install missing libraries
#                      default: FALSE, i.e. execution is aborted when missing packages are
#                      encountered)
#     repo           : Repository for downloading packages, when install = TRUE (default:
#                      http://cran.us.r-project.org)
#     status         : Boolean indicating whether status messages shall be printed (default: FALSE)
# [Return value]
#     N/A
# [Dependencies]
#     printStatus(), printError()

    #---> STATUS MESSAGE <---#
    printStatus("Loading dependencies...", status)

    #---> BODY // <---#

        #---> Iterate over packages <---#
        for ( package in packages ) {

            #---> Status message <---#
            printStatus(paste("Loading package '", package, "'...", sep=""), status)

            #---> Load package and attach to environment <---#
            success <- suppressMessages(suppressWarnings(
                           require(
                               package,
                               quietly = TRUE,
                               warn.conflicts = FALSE,
                               character.only = TRUE
                           )
                       ))

            #---> If package was not successfully attached... <---#
            if ( ! success ) {

                #---> Either try to install package... <---#
                if ( install ) {

                    #---> Die if package is unavailable <---#
                    pkgNames <- available.packages(contrib.url("http://cran.us.r-project.org", type="source"))[ , 1]
                        if ( ! package %in% pkgNames ) {
                            printError(
                                c("Package '", package, "' is not available in repository '",
                                  repo, "'! Verify the package name and repository.")
                            )
                            quit(save = "no", status = 1, runLast = FALSE)
                        }

                    #---> Install package <---#
                    install.packages(package, repos=repo, dependencies = TRUE, quiet = TRUE)
                    require(package, quietly = TRUE, warn.conflicts = FALSE, character.only = TRUE)

                #---> ...or die! <---#
                } else {
                    printError(
                        c("Package '", package,
                          "' could not be loaded/attached! Verify whether it is installed.")
                    )
                    quit(save = "no", status = 1, runLast = FALSE)
                }
            }
        }

    #---> // BODY <---#

    #---> STATUS MESSAGE <---#
    printStatus("Dependencies loaded...", status)

}
#-----------------------#
getLogDateTime <- function(template = '%Y/%m/%d %H:%M:%S') {

# [Description]
#     Returns current date/time stamp in a format suitable for log entries.
# [Parameters]
#     template       : POSIX time template string (default: "%Y/%m/%d %H:%M:%S")
# [Return value]
#     Date/time string
# [Dependencies]
#     N/A

    #---> BODY // <---#

        #---> Build output filename <---#
        dateTime <- format.Date(Sys.time(), template)

    #---> // BODY <---#

    #---> RETURN VALUE <---#
    return(dateTime)

}
#-----------------------#
printError <- function(message, toScreen = TRUE) {

# [Description]
#     Print error message to STDERR.
# [Parameters]
#     message        : Error message string
#     print          : Boolean whether message shall be printed, e.g. value of '--verbose' option
#                      (default: TRUE)
# [Return value]
#     N/A
# [Dependencies]
#     N/A

    #---> BODY // <---#

        #---> Build message <---#
        message <- c("[ERROR] ", message, "\n",
                     "[ERROR] Execution aborted!")
        message <- paste(message, collapse="")

        #---> Print message <---#
        if (toScreen) { write(message, stderr()) }

    #---> // BODY <---#

}
#-----------------------#
printWarning <- function(message, toScreen = TRUE) {

# [Description]
#     Print warning to STDERR.
# [Parameters]
#     message        : Warning message string
#     toScreen       : Boolean whether message shall be printed, e.g. value of '--verbose' option
#                      (default: TRUE)
# [Return value]
#     N/A
# [Dependencies]
#     N/A

    #---> BODY // <---#

        #---> Build message <---#
        message <- c("[WARNING] ", message)

        #---> Print message <---#
        if (toScreen) { write(message, stderr()) }

    #---> // BODY <---#

}
#-----------------------#
printStatus <- function(message, toScreen = TRUE) {

# [Description]
#     Print status message to STDERR.
# [Parameters]
#     message        : Status message string
#     toScreen       : Boolean whether message shall be printed, e.g. value of '--verbose' option
#                      (default: TRUE)
# [Return value]
#     N/A
# [Dependencies]
#     getLogDateTime()

    #---> BODY // <---#

        #---> Build message <---#
        message <- paste(c("[", getLogDateTime(), "] ", message), collapse="")

        #---> Print message <---#
        if (toScreen) { write(message, stderr()) }

    #---> // BODY <---#

}
#-----------------------#
getScriptName <- function() {

# [Description]
#     Extract name of executing script.
# [Parameters]
#     N/A
# [Return value]
#     Name of executing script
# [Dependencies]
#     N/A

    #---> BODY // <---#

        #---> Extract script name <---#
        scriptName <- sub("--file=", "", basename(commandArgs(trailingOnly=FALSE)[4]))

    #---> // BODY <---#

    #---> RETURN VALUE <---#
    return(scriptName)

}
#-----------------------#
formatOptions <- function(parsedOptions, displayWidth = 100, minFlagWidth = 22, flagPrefix = 4, flagSuffix = 2) {

# [Description]
#     Formats list of option for '--usage' output based on an optparse::OptionParser object.
# [Parameters]
#     parsedOptions  : An optparse::OptionParser object
#     displayWidth   : Absolute width after which option descriptions should be wrapped (default:
#                      100)
#     minFlagWidth   : Minimum width of spaced reserved for option flags (default: 22). The width
#                      actually used for plotting is the maximum of this value and the number of
#                      characters of the longest flag (plus 'metavalue').
#     flagPrefix     : Number of empty columns before option flags (default: 4)
#     flagSuffix     : Number of empty columns after option flags (default: 2)
# [Return value]
#     Formatted options string
# [Dependencies]
#     N/A

    #---> BODY // <---#

        #---> Apply over list of options <---#
        formattedOptions <- lapply(parsedOptions@`options`, function(opt) {

            #---> Format flags & 'metavar' <---#
            if ( is.na(opt@`short_flag`) ) {
                flags <- opt@`long_flag`
            } else {
                flags <- paste(opt@`long_flag`, " | ", flags <- opt@`short_flag`, sep="")
            }
            if ( length(opt@`metavar`) > 0 && nchar(opt@`metavar`) > 0 ) {
                flags <- paste(flags, " <", opt@`metavar`, ">", sep = "")
            }

            #---> Format option description <---#
            lastChar <- substr(opt@`help`, nchar(opt@`help`), nchar(opt@`help`))
            if ( lastChar %in% c(".", "?", "!", ";", ",") ) {
                helpPrefix <- substr(opt@`help`, 1, nchar(opt@`help`) - 1)
            } else {
                helpPrefix <- opt@`help`
                lastChar <- "."
            }

            #---> Format default value <---#
            if (
               is.null(opt@`default`)                 ||
               is.na(opt@`default`)                   ||
               opt@`default` == ""                    ||
               opt@`default` == FALSE                 ||
               identical(character(0), opt@`default`) ||
               identical(numeric(0), opt@`default`)   ||
               identical(integer(0), opt@`default`)   ||
               identical(logical(0), opt@`default`)
               ) {
                default <- ""
            } else if ( is.character(opt@`default`) ) {
                default <- paste(" Default: '", opt@`default`, "'.", sep="")
            } else {
                default <- paste(" Default: ", opt@`default`, ".", sep="")
            }
            help <- paste(helpPrefix, lastChar, default, sep="")

            #---> Build string <---#
            optString <- c(flags, help)

            #---> Return string <---#
            return(optString)

        })

        #---> Determine width for option flags <---#
        flagWidth <- max(nchar(sapply(formattedOptions, "[", 1)), minFlagWidth)

        #---> Build pre- and suffix for flags <---#
        flagPrefix <- paste(rep(" ", flagPrefix), collapse="")
        flagSuffix <- paste(rep(" ", flagSuffix), collapse="")

        #---> Define template string for flag formatting <---#
        templateString <- paste("%s%-", flagWidth, "s%s", sep="")

        #---> Determine width for option description <---#
        descrWidth <- displayWidth - nchar(flagPrefix) - flagWidth - nchar(flagSuffix)

        #---> Set prefix for additional description lines <---#
        descrPrefix <- paste("\n", paste(rep(" ", displayWidth - descrWidth), collapse=""), sep="")

        #---> Format all options and collapse into string <---#
        optString <- paste(sapply(formattedOptions, function(opt) {
            flags <- sprintf(templateString, flagPrefix, opt[[1]], flagSuffix)
            descr <- paste(strwrap(opt[[2]], width=descrWidth), collapse=descrPrefix)
            string <- paste(flags, descr, sep="")
        }), collapse="\n")

    #---> // BODY <---#

    #---> RETURN VALUE <---#
    return(optString)

}
#-----------------------#
userLicense <- function() {

# [Description]
#     Build and return license string.
# [Parameters]
#     N/A
# [Return value]
#     License string
# [Dependencies]
#     N/A

    #---> BODY // <---#

        #---> Build string <---#
        license <- paste(c(
"{{LICENSE}}"
        ), collapse="")

    #---> // BODY <---#

    #---> RETURN VALUE <---#
    return(license)

}
#-----------------------#
userVersion <- function(scriptName) {

# [Description]
#     Build and return version string.
# [Parameters]
#     scriptName     : Name of the script for which a usage string shalle be generated
# [Return value]
#     Version string
# [Dependencies]
#     N/A

    #---> BODY // <---#

        #---> Build string <---#
        version <- paste(c(
scriptName, ", v1.0 (01-MAR-2016)
(c) 2016 Alexander Kanitz, University of Basel"
        ), collapse="")

    #---> // BODY <---#

    #---> RETURN VALUE <---#
    return(version)

}
#-----------------------#
usage <- function(scriptName, parsedOptions) {

# [Description]
#     Build and return help/usage string.
# [Parameters]
#     scriptName     : Name of the script for which a usage string shalle be generated
#     parsedOptions  : An object of class OptionParser, returned by optparse::OptionParser()
# [Return value]
#     Usage string
# [Dependencies]
#     userVersion(), formatOptions()

    #---> BODY // <---#

        #---> Get version string <---#
        versionString <- paste(userVersion(scriptName), collapse="")
        versionString <- unlist(strsplit(versionString, "\n"))
        versionString <- paste('    ', versionString, sep="")
        versionString <- paste(versionString, collapse="\n")

        #---> Build string <---#
        usage <- paste(c(
"[VERSION INFORMATION]
", versionString, "

[CONTACT INFORMATION]
    Alexander Kanitz, alexander.kanitz@alumni.ethz.ch
    Biozentrum, University of Basel

[USAGE]
    ", scriptName, " [OPTIONS] <GTF>

[DESCRIPTION]
    Extracts regions in BED/BED12 format from a GTF gene annotation file. See notes for the
    different region types that can be extracted.

[DEPENDENCIES]
    optparse, GenomicRanges, rtracklayer

[OPTIONS]
", formatOptions(parsedOptions), "

[USAGE NOTES]
    (1) {{}}."
        ), collapse="")

    #---> // BODY <---#

    #---> RETURN VALUE <---#
    return(usage)

}
#-----------------------#
validateArgs <- function(arguments, options, parsedOptions, keys=NULL, exactly=NULL, atLeast=exactly, atMost=exactly) {

# [Description]
#     Checks for unrecognized options, verifies count of positional command-line arguments, adds 
#     names/keys to arguments and adds arguments to options; missing names will be set to "posN" 
#     where N is the position of the argument; surplus names will be ignored
# [Parameters]
#     arguments      : Vector of positional arguments, such as is produced by optparse::parse_args()
#     options        : List of option values, such as is produced by optparse::parse_args()
#     parsedOptions  : An object of class OptionParser, returned by optparse::OptionParser()
#     keys           : Character vector of argument names; must be of equal length as 'args'
#     exactly        : Exact count of required positional arguments
#     atLeast        : Minimum required count of positional arguments
#     atMost         : Maximum allowed count of positional arguments
# [Return value]
#     N/A
# [Dependencies]
#     printError()

    #---> BODY // <---#

        # Check for unrecognized options
        idx <- grep('^--', arguments)
        if ( length(idx) ) {
            printError(c("Unrecognized option: ", arguments[idx]))
            write(usage(options$`scriptName`, parsedOptions), stderr())
            quit(save = "no", status=1, runLast=FALSE)
        }

        # Validate count of positional arguments
        if ( ! is.null(atLeast) && length(arguments) < atLeast ) {
            printError(c("At least ", atLeast, " arguments required."))
            write(usage(options$`scriptName`, parsedOptions), stderr())
            quit(save = "no", status=1, runLast=FALSE)
        }
        if ( ! is.null(atMost) && length(arguments) > atMost ) {
            printError(c("At most ", atMost, " arguments allowed."))
            write(usage(options$`scriptName`, parsedOptions), stderr())
            quit(save = "no", status=1, runLast=FALSE)
        }

        # Convert argument vector to list
        arguments <- as.list(arguments)

        # Add names to arguments list
        if ( length(arguments) - length(keys) > 0 ) {
            keys <- c(keys, paste("pos", (length(keys) + 1) : length(arguments), sep=""))
        } else if ( length(arguments) - length(keys) < 0 ) {
            keys <- keys[0:length(arguments)]
        }
        names(arguments) <- keys

        # Add arguments to options
        options <- c(options, arguments)

    #---> // BODY <---#

    #---> RETURN VALUE <---#
    return(options)

}
#-----------------------#
validateOptions <- function(options, parsedOptions) {

# [Description]
#     Validates command-line options.
# [Parameters]
#     options        : List of option values, such as is produced by optparse::parse_args()
#     parsedOptions  : An object of class OptionParser, returned by optparse::OptionParser()
# [Return value]
#     N/A
# [Dependencies]
#     printError()

    #---> BODY // <---#

        # Print usage information and die if required argument is missing
        if ( is.null(options$gtf) ) {
            printError("Required argument missing.")
            write(usage(options$`scriptName`, parsedOptions), stderr())
            quit(save = "no", status=1, runLast=FALSE)
        }

        # Print usage information and die if file is not found
        if ( ! file.exists(options$gtf) ) {
            printError(c("File '", options$gtf, "' was not found."))
            write(usage(options$`scriptName`, parsedOptions), stderr())
            quit(save = "no", status=1, runLast=FALSE)
        }
        if ( ! is.null(options$`chromosome-sizes`) && ! file.exists(options$`chromosome-sizes`) ) {
            printError(c("File '", options$`chromosome-sizes`, "' was not found."))
            write(usage(options$`scriptName`, parsedOptions), stderr())
            quit(save = "no", status=1, runLast=FALSE)
        }

    #---> // BODY <---#

    #---> RETURN VALUE <---#
    return(options)

}
#-----------------------#
parseOptions <- function(argNames=NULL, argCount=NULL, minArgs=argCount, maxArgs=argCount) {

# [Description]
#     Returns validated command-line options.
# [Parameters]
#     N/A
# [Return value]
#     List of CLI options, as well as 'scriptName'
# [Dependencies]
#     optparse, loadPackages(), usage(), userVersion(), userLicense(), validateOptions()

    #---> LOAD FUNCTION DEPENDENCIES <---#
    loadPackages("optparse", FALSE)

    #---> BODY // <---#

        #---> Generate list of options <---#
        option_list <- list(
            make_option("--chromosome-sizes", action="store", type="character", default=NULL, help="Tab-separated chromsome sizes table of the format <CHR_NAME>\t<SIZE>, without header. Required for extracting intergenic regions.", metavar="FILE"),
            make_option("--verbose", action="store_true", default=FALSE, help="Print log messages to STDOUT."),
            make_option("--version", action="store_true", default=FALSE, help="Show version information and exit."),
            make_option("--license", action="store_true", default=FALSE, help="Show license information and exit."),
            make_option("--help", action="store_true", default=FALSE, help="Show this screen and exit."),
            make_option("--usage", dest="help", action="store_true", default=FALSE, help="Show this screen and exit.")
        )

        #---> Parse options <---#
        parsedOptions <- OptionParser(option_list = option_list, add_help_option=FALSE)
        opts_args <- parse_args(parsedOptions, print_help_and_exit=FALSE, positional_arguments=TRUE)
        options <- opts_args$options

        #---> Get script name <---#
        options$`scriptName` <- getScriptName()

        #---> Show usage / version / license <---#
        if ( options$`help` ) {
            write(usage(options$`scriptName`, parsedOptions), stderr())
            quit(save = "no", status=0, runLast=FALSE)
        }
        if ( options$`version` ) {
            write(userVersion(options$`scriptName`), stderr())
            quit(save = "no", status=0, runLast=FALSE)
        }
        if ( options$`license` ) {
            write(userLicense(), stderr())
            quit(save = "no", status=0, runLast=FALSE)
        }

        #---> Add arguments to options & validation <---#
        options <- validateArgs(opts_args$args, options, parsedOptions, keys=argNames, exactly=argCount, atLeast=minArgs, atMost=maxArgs)
        options <- validateOptions(options, parsedOptions)

        #---> Add arguments to options <---#
        options <- c(options, args)

    #---> // BODY <---#

    #---> RETURN VALUE <---#
    return(options)

}
#-----------------------#
getGenicRegions <- function(obj, verbose) {

# [Description]
#     Given GTF-derived GRanges object 'obj$all_gr', extracts genic regions,
#     i.e. all ranges for which the 'type' metadata column equals 'gene'
# [Parameters]
#     obj            : List of range objects (must include element 'all_gr')
#     verbose        : Print status messages
# [Return value]
#     Updated objects list
# [Dependencies]
#     printStatus(), GenomicRanges

    #---> STATUS MESSAGE <---#
    printStatus("Extracting genic regions...", verbose)

    #---> BODY // <---#

        #---> Extract genic ranges <---#
        obj$genic_gr <- obj$all_gr[obj$all_gr$type == "gene"]
        obj$genic_gr$name <- obj$genic_gr$gene_id

    #---> // BODY <---#

    #---> RETURN VALUE <---#
    return(obj)

}
#-----------------------#
getIntergenicRegions <- function(obj, verbose) {

# [Description]
#     Given GTF-derived GRanges object 'obj$genic_gr', extracts intergenic
#     regions, i.e. the complement of genic regions
# [Parameters]
#     obj            : List of range objects (must include element genic_gr')
#     verbose        : Print status messages
# [Return value]
#     Updated objects list
# [Dependencies]
#     printStatus(), GenomicRanges

    #---> STATUS MESSAGE <---#
    printStatus("Extracting intergenic regions...", verbose)

    #---> BODY // <---#

        #---> Extract intergenic ranges <---#
        obj$intergenic_gr <- gaps(obj$genic_gr)
        obj$intergenic_gr$name <- paste("INTERGENIC", sprintf("%011d", 1:length(obj$intergenic_gr)), sep="")

    #---> // BODY <---#

    #---> RETURN VALUE <---#
    return(obj)

}
#-----------------------#
getTranscribedRegions <- function(obj, verbose) {

# [Description]
#     Given GTF-derived GRanges object 'obj$all_gr', extracts transcribed
#     regions, i.e. the combination of all exons and introns of each transcript
# [Parameters]
#     obj            : List of range objects (must include element 'all_gr')
#     verbose        : Print status messages
# [Return value]
#     Updated objects list
# [Dependencies]
#     printStatus(), GenomicRanges

    #---> STATUS MESSAGE <---#
    printStatus("Extracting transcribed regions...", verbose)

    #---> BODY // <---#

        #---> Extract untranscribed regions <---#
        obj$transcribed_gr <- obj$all_gr[obj$all_gr$type == "transcript"]
        obj$transcribed_gr$name <- paste(
            obj$transcribed_gr$gene_id,
            obj$transcribed_gr$transcript_id,
            sep="|"
        )

    #---> // BODY <---#

    #---> RETURN VALUE <---#
    return(obj)

}
#-----------------------#
getUntranscribedRegions <- function(obj, verbose) {

# [Description]
#     Given GTF-derived GRanges objects 'obj$genic_gr' and
#     'obj$transcribed_gr', extracts untranscribed regions, i.e. regions of a
#     transcript that are part of a gene but are not transcribed in an
#     individual transcript
# [Parameters]
#     obj            : List of range objects (must include element 'all_gr')
#     verbose        : Print status messages
# [Return value]
#     Updated objects list
# [Dependencies]
#     printStatus(), GenomicRanges

    #---> STATUS MESSAGE <---#
    printStatus("Extracting untranscribed regions...", verbose)

    #---> BODY // <---#

        #---> Extract untranscribed regions <---#
        tmp_gr <- obj$transcribed_gr
        id_idx <- match(obj$transcribed_gr$gene_id, obj$genic_gr$gene_id)
        ranges(tmp_gr) <- ranges(obj$genic_gr[id_idx])
        obj$untranscribed_gr <- unlist(
            psetdiff(
                split(tmp_gr, tmp_gr$transcript_id),
                split(obj$transcribed_gr, obj$transcribed_gr$transcript_id)
            )
        )
        obj$untranscribed_gr$transcript_id <- names(obj$untranscribed_gr)
        names(obj$untranscribed_gr) <- NULL
        obj$untranscribed_gr$gene_id <- obj$transcribed_gr$gene_id[
            match(
                obj$untranscribed_gr$transcript_id,
                obj$transcribed_gr$transcript_id
            )
        ]
        tmp_uniq_name <- paste("UNTRANSCRIBED", sprintf("%011d", 1:length(obj$untranscribed_gr)), sep="")
        obj$untranscribed_gr$name <- paste(
            obj$untranscribed_gr$gene_id,
            obj$untranscribed_gr$transcript_id,
            tmp_uniq_name,
            sep="|"
        )

    #---> // BODY <---#

    #---> RETURN VALUE <---#
    return(obj)

}
#-----------------------#
getExons <- function(obj, verbose) {

# [Description]
#     Given GTF-derived GRanges object 'obj$all_gr', extracts the exons of all
#     transcripts
# [Parameters]
#     obj            : List of range objects (must include element 'all_gr')
#     verbose        : Print status messages
# [Return value]
#     Updated objects list
# [Dependencies]
#     printStatus(), GenomicRanges

    #---> STATUS MESSAGE <---#
    printStatus("Extracting exons...", verbose)

    #---> BODY // <---#

        #---> Extract exons <---#
        obj$exons_gr <- obj$all_gr[obj$all_gr$type == "exon"]
        obj$exons_gr$name <- paste(
            obj$exons_gr$gene_id,
            obj$exons_gr$transcript_id,
            obj$exons_gr$exon_id,
            sep="|"
        )

    #---> // BODY <---#

    #---> RETURN VALUE <---#
    return(obj)

}
#-----------------------#
getExonicRegions <- function(obj, verbose) {

# [Description]
#     Given GTF-derived GRanges object 'obj$exons_gr', extracts the exonic
#     regions of all genes, i.e. the merge of all exons of each gene
# [Parameters]
#     obj            : List of range objects (must include element 'exons_gr')
#     verbose        : Print status messages
# [Return value]
#     Updated objects list
# [Dependencies]
#     printStatus(), GenomicRanges

    #---> STATUS MESSAGE <---#
    printStatus("Extracting exonic regions...", verbose)

    #---> BODY // <---#

        #---> Extract exons <---#
        obj$exonic_gr <- unlist(reduce(split(obj$exons_gr, obj$exons_gr$gene_id)))
        obj$exonic_gr$gene_id <- names(obj$exonic_gr)
        names(obj$exonic_gr) <- NULL

    #---> // BODY <---#

    #---> RETURN VALUE <---#
    return(obj)

}
#-----------------------#
getIntrons <- function(obj, verbose) {

# [Description]
#     Given several GTF-derived GRanges objects in 'obj', extracts the introns
#     of all transcripts
# [Parameters]
#     obj            : List of range objects (must include elements
#                      'transcribed_gr', 'untranscribed_gr', 'genic_gr' and
#                      'exons_gr')
#     verbose        : Print status messages
# [Return value]
#     Updated objects list
# [Dependencies]
#     printStatus(), GenomicRanges

    #---> STATUS MESSAGE <---#
    printStatus("Extracting introns...", verbose)

    #---> BODY // <---#

        #---> Extract exons <---#
        tmp_big_gr <- obj$transcribed_gr
        id_idx <- match(obj$transcribed_gr$gene_id, obj$genic_gr$gene_id)
        ranges(tmp_big_gr) <- ranges(obj$genic_gr[id_idx])
        tmp_small_gr <- c(
            GRanges(
                seqnames(obj$exons_gr),
                IRanges(ranges(obj$exons_gr), names=obj$exons_gr$transcript_id),
                strand(obj$exons_gr)
            ),
            GRanges(
                seqnames(obj$untranscribed_gr),
                IRanges(ranges(obj$untranscribed_gr), names=obj$untranscribed_gr$transcript_id),
                strand(obj$untranscribed_gr)
            )
        )
        obj$introns_gr <- unlist(
            psetdiff(
                split(tmp_big_gr, tmp_big_gr$transcript_id),
                reduce(split(tmp_small_gr, names(tmp_small_gr)))
            )
        )
        obj$introns_gr$transcript_id <- names(obj$introns_gr)
        names(obj$introns_gr) <- NULL
        obj$introns_gr$gene_id <- obj$transcribed_gr$gene_id[
            match(
                obj$introns_gr$transcript_id,
                obj$transcribed_gr$transcript_id
            )
        ]
        tmp_uniq_name <- paste("INTRON", sprintf("%011d", 1:length(obj$introns_gr)), sep="")
        obj$introns_gr$name <- paste(
            obj$introns_gr$gene_id,
            obj$introns_gr$transcript_id,
            tmp_uniq_name,
            sep="|"
        )

    #---> // BODY <---#

    #---> RETURN VALUE <---#
    return(obj)

}
#-----------------------#
getIntronicRegions <- function(obj, verbose) {

# [Description]
#     Given GTF-derived GRanges object 'obj$introns_gr', extracts the exonic
#     regions of all genes, i.e. the merge of all exons of each gene
# [Parameters]
#     obj            : List of range objects (must include element 'introns_gr')
#     verbose        : Print status messages
# [Return value]
#     Updated objects list
# [Dependencies]
#     printStatus(), GenomicRanges

    #---> STATUS MESSAGE <---#
    printStatus("Extracting intronic regions...", verbose)

    #---> BODY // <---#

        #---> Extract introns <---#
        obj$intronic_gr <- unlist(reduce(split(obj$introns_gr, obj$introns_gr$gene_id)))
        obj$intronic_gr$gene_id <- names(obj$intronic_gr)
        names(obj$intronic_gr) <- NULL

    #---> // BODY <---#

    #---> RETURN VALUE <---#
    return(obj)

}
#-----------------------#
getCodingRegions <- function(obj, verbose) {

# [Description]
#     Given GTF-derived GRanges object 'obj$all_gr', extracts the coding
#     regions/CDS of all transcripts
# [Parameters]
#     obj            : List of range objects (must include element 'all_gr')
#     verbose        : Print status messages
# [Return value]
#     Updated objects list
# [Dependencies]
#     printStatus(), GenomicRanges

    #---> STATUS MESSAGE <---#
    printStatus("Extracting coding regions...", verbose)

    #---> BODY // <---#

        #---> Extract exons <---#
        obj$cds_gr <- obj$all_gr[obj$all_gr$type == "CDS"]
        tmp_uniq_name <- paste("CDS", sprintf("%011d", 1:length(obj$cds_gr)), sep="")
        obj$cds_gr$name <- paste(
            obj$cds_gr$gene_id,
            obj$cds_gr$transcript_id,
            tmp_uniq_name,
            sep="|"
        )

    #---> // BODY <---#

    #---> RETURN VALUE <---#
    return(obj)

}
getFiveUTRs <- function(obj, verbose) {

# [Description]
#     Given GTF-derived GRanges object 'obj$all_gr', extracts the coding
#     regions/CDS of all transcripts
# [Parameters]
#     obj            : List of range objects (must include element 'all_gr')
#     verbose        : Print status messages
# [Return value]
#     Updated objects list
# [Dependencies]
#     printStatus(), GenomicRanges

    #---> STATUS MESSAGE <---#
    printStatus("Extracting coding regions...", verbose)

    #---> BODY // <---#

        #---> Extract exons <---#
        obj$cds_gr <- obj$all_gr[obj$all_gr$type == "CDS"]
        tmp_uniq_name <- paste("CDS", sprintf("%011d", 1:length(obj$cds_gr)), sep="")
        obj$cds_gr$name <- paste(
            obj$cds_gr$gene_id,
            obj$cds_gr$transcript_id,
            tmp_uniq_name,
            sep="|"
        )

    #---> // BODY <---#

    #---> RETURN VALUE <---#
    return(obj)

}



#================#
#  // FUNCTIONS  #
#================#



#================#
#  // FUNCTIONS  #
#================#


#===========#
#  MAIN //  #
#===========#

#---> PROCESS GLOBAL OPTIONS <---#
opt <- parseOptions(argNames="gtf")

#---> STATUS MESSAGE <---#
printStatus(c("Starting '", opt$`scriptName`, "'..."), opt$`verbose`)

#---> LOAD GLOBAL DEPENDENCIES <---#
loadPackages(c("GenomicRanges", "rtracklayer"))

#---> INITIALIZE GLOBAL VARIABLES <---#
obj <- list()

# TIME
t1 <- proc.time()
# /TIME

#---> BODY // <---#


    #---> Import GTF file <---#
    printStatus(c("Importing GTF file '", opt$`gtf`, "'..."), opt$`verbose`)
    obj$all_gr <- import(con=opt$gtf, format="gtf")

    #---> If available, import chromosome sizes table and update GRanges object  <---#
    if ( ! is.null(opt$`chromosome-sizes`) ) {
        printStatus(c("Reading chromosome size table '", opt$`chromosome-sizes`, "'..."), opt$`verbose`)
        obj$seqlen <- read.delim(file=opt$`chromosome-sizes`, header=FALSE, stringsAsFactors=FALSE, col.names=c("chr", "size"))
        seqlengths(obj$all_gr) <- obj$seqlen[,2][match(names(seqlengths(obj$all_gr)), obj$seqlen[,1])]
    }

    #---> Get genic regions <---#
    obj <- getGenicRegions(obj, opt$`verbose`)

    #---> Get intergenic regions <---#
    if ( ! is.null(opt$`chromosome-sizes`) ) {
        obj <- getIntergenicRegions(obj, opt$`verbose`)
    }

    #---> Get transcribed regions <---#
    obj <- getTranscribedRegions(obj, opt$`verbose`)

    #---> Get untranscribed regions <---#
    obj <- getUntranscribedRegions(obj, opt$`verbose`)

    #---> Get exons <---#
    obj <- getExons(obj, opt$`verbose`)

    #---> Get exonic regions <---#
    obj <- getExonicRegions(obj, opt$`verbose`)

    #---> Get introns <---#
    obj <- getIntrons(obj, opt$`verbose`)

    #---> Get intronic regions <---#
    obj <- getIntronicRegions(obj, opt$`verbose`)

    #---> Get coding regions <---#
    obj <- getCodingRegions(obj, opt$`verbose`)

    #---> Get 5'-UTRs <---#
    obj <- getFiveUTRs(obj, opt$`verbose`)

    #---> Get 3'-UTRs <---#
    obj <- getThreeUTRs(obj, opt$`verbose`)

    #---> {{}} <---#



    #---> {{}} <---#



#---> // BODY <---#

#---> SAVE SESSION <---#
printStatus("Saving session...", opt$`verbose`)
save.image()

#---> STATUS MESSAGE <---#
printStatus("Done.", opt$`verbose`)

#---> PROGRAM EXIT <---#
quit(save = "no", status=0, runLast=FALSE)

#===========#
#  // MAIN  #
#===========#
