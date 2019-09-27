#!/usr/bin/env Rscript


#=====================#
#  ISSUES & IDEAS //  #
#=====================#
# TODO: Output: add summary, smear plot (with feature highlighting?)
# TODO: Put query column derivation in function (together with tests: are columns left?)
# TODO: Add support for secondary (i.e. gene) feature summing
# TODO: Add functionality to run more than one contrast (read parameters from file)
#=====================#
#  // ISSUES & IDEAS  #
#=====================#


#========================#
#  GENERIC FUNCTIONS //  #
#========================#
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
            success <- suppressWarnings(
                           require(
                               package,
                               quietly = TRUE,
                               warn.conflicts = FALSE,
                               character.only = TRUE
                           )
                       )

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
    printStatus("Dependencies loaded.", status)

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
        if (toScreen) write(message, stderr())

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
        message <- paste(message, collapse="")

        #---> Print message <---#
        if (toScreen) write(message, stderr())

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
        message <- c("[", getLogDateTime(), "] ", message)
        message <- paste(message, collapse="")

        #---> Print message <---#
        if (toScreen) write(message, stderr())

    #---> // BODY <---#

}
#========================#
#  // GENERIC FUNCTIONS  #
#========================#


#===================================#
#  CLI OPTION-RELATED FUNCTIONS //  #
#===================================#
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

            #---> Format default value <---#
            if (
               is.null(opt@`default`)                 ||
               is.na(opt@`default`)                   ||
               opt@`default` == ""                    ||
               identical(character(0), opt@`default`) ||
               identical(numeric(0), opt@`default`)   ||
               identical(integer(0), opt@`default`)   ||
               identical(logical(0), opt@`default`)
               ) {
                default <- " (no default)"
            } else if ( is.character(opt@`default`) ) {
                default <- paste(" (default: '", opt@`default`, "')", sep="")
            } else {
                default <- paste(" (default: ", opt@`default`, ")", sep="")
            }

            #---> Format option description <---#
            lastChar <- substr(opt@`help`, nchar(opt@`help`), nchar(opt@`help`))
            if ( lastChar %in% c(".", "?", "!", ";", ",") ) {
                helpPrefix <- substr(opt@`help`, 1, nchar(opt@`help`) - 1)
            } else {
                helpPrefix <- opt@`help`
                lastChar <- "."
            }
            help <- paste(helpPrefix, default, lastChar, sep="")

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
'The MIT License (MIT)

Copyright (c) <<COPYRIGHT_YEAR>> <<AUTHOR_NAME>>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.'
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
scriptName, ", v1.0 (Sep 29th, 2015)
(c) 2015 Alexander Kanitz"
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
    Alexander Kanitz <alexander.kanitz@alumni.ethz.ch>
    Biozentrum, University of Basel

[USAGE]
    ", scriptName, " [OPTIONS] --counts <FILE> --reference <STRING>

[DESCRIPTION]
    Perform differential expression analyses with the R/Bioconductor package edgeR.

[DEPENDENCIES]
    optparse, edgeR

[OPTIONS]
", formatOptions(parsedOptions), "

[USAGE NOTES]
    (1) Strings supplied to '--reference' (required) and '--query' (optional) can consist of
        positive comma-separated integers and may contain dashes to indicate ranges. Duplicate
        values are explicitly allowed and will select a column multiple times. Spaces may be
        inserted at any place to improve readability, but mind the quoting rules of your shell."
        ), collapse="")

    #---> // BODY <---#

    #---> RETURN VALUE <---#
    return(usage)

}
#-----------------------#
validateColumnSelection <- function(string) {

# [Description]
#     Returns a list of integers from a string containing spaces, commas and dashes.
# [Parameters]
#     string         : List of options and values, such as is produced by optparse::parse_args()
# [Return value]
#     Integer vector
# [Dependencies]
#     N/A

    #---> INITIALIZE FUNCTION VARIABLES <---#
    allowed <- c(",", "-", 0:9)
    results <- list()
    results$results <- NULL
    results$error <- list()
    results$warnings <- list()

    #---> BODY // <---#

        #---> Remove whitespace <---#
        string <- gsub("\\s", "", string)

        #---> Return error if empty string <---#
        if ( string == "" ) {
            err <- list(type="EmptyString", message="String is empty.")
            results$error[[length(results$error) + 1]] <- err
        }

        #---> Return error if forbidden characters are found <---#
        charVec <- unlist(strsplit(string, ""))
        illegal <- ! charVec %in% allowed
        if ( any(illegal) ) {
            err <- list(type="IllegalChar", message=paste("String contains illegal characters:", paste(unique(charVec[illegal]), collapse=", ")))
            results$error[[length(results$error) + 1]] <- err
        }

        #---> Return error if string does not begin or end with a digit <---#
        if ( grepl("^\\D", string) ) {
            err <- list(type="IllegalFormat", message=paste("First character ('", charVec[1], "') is not a digit (whitespace ignored).", sep=""))
            results$error[[length(results$error) + 1]] <- err
        }
        if ( grepl("\\D$", string) ) {
            err <- list(type="IllegalFormat", message=paste("Last character ('", charVec[length(charVec)], "') is not a digit (whitespace ignored).", sep=""))
            results$error[[length(results$error) + 1]] <- err
        }

        #---> Split by commas <---#
        terms <- unlist(strsplit(string, ","))

        #---> Return error if empty strings are found <---#
        if ( any(terms == "") ) {
            terms <- terms[terms != ""]
            err <- list(type="IllegalFormat", message="No valid number or range between commas.")
            results$error[[length(results$error) + 1]] <- err
        }

        #---> Split by dashes <---#
        intVec <- suppressWarnings(as.integer(unlist(lapply(terms, function(term) {
            if ( grepl("-", term) ) {
                if ( grepl("^\\d+-\\d+$", term) ) {
                    startEnd <- unlist(strsplit(term, "-"))
                    return(startEnd[1]:startEnd[2])
                } else {
                    return(NA)
                }
            } else {
                return(term)
            }
        }))))

        #---> Return error if dash used illegaly <---#
        if ( any(is.na(intVec)) ) {
            err <- list(type="IllegalFormat", message="Illegal dash syntax. Use 'start-end'.")
            results$error[[length(results$error) + 1]] <- err
        }

        #---> Return warning if duplicated numbers <---#
        if ( any(duplicated(na.omit(intVec))) ) {
            warn <- list(type="Duplicate", message=paste("Columns selected more than once:", paste(unique(intVec[duplicated(intVec)]), collapse=", ")))
            results$warnings[[length(results$warnings) + 1]] <- warn
        }

        #---> Add index vector to results list
        results$results <- intVec

    #---> // BODY <---#

    #---> RETURN VALUE <---#
    return(results)

}
#-----------------------#
validateOptions <- function(options, parsedOptions) {

# [Description]
#     Validates command-line options.
# [Parameters]
#     options        : List of options and values, such as is produced by optparse::parse_args()
#     parsedOptions  : An object of class OptionParser, returned by optparse::OptionParser()
# [Return value]
#     N/A
# [Dependencies]
#     printError(), printWarning()

    #---> BODY // <---#

        #---> Ensure that required options/arguments are specified <---#
        if ( is.null(options$`counts`) ) {
            printError(c("Required option missing. Please specify option '--count' and specify a valid arguments."))
            stop(usage(options$`scriptName`, parsedOptions))
        }

        if ( is.null(options$`reference`) ) {
            printError(c("Required option missing. Please specify option '--reference' and specify a valid arguments."))
            stop(usage(options$`scriptName`, parsedOptions))
        }

        #---> Set defaults if not specified <---#
        if ( is.null(options$`output-directory`) ) {
            options$`output-directory` <- getwd()
        }

        #---> Ensure that input files exist <---#
        if ( ! file.exists(options$`counts`) ) {
            printError(c("The file '", options$`counts`, "', specified as an argument to option '--counts' could not be found."))
            stop(usage(options$`scriptName`, parsedOptions))
        }

        if ( ! is.null(options$`abundances`) && ! file.exists(options$`abundances`) ) {
            printError(c("The file '", options$`abundances`, "', specified as an argument to option '--abundances' could not be found."))
            stop(usage(options$`scriptName`, parsedOptions))
        }

        if ( ! is.null(options$`lookup-table`) && ! file.exists(options$`lookup-table`) ) {
            printError(c("The file '", options$`lookup-table`, "', specified as an argument to option '--lookup-table' could not be found."))
            stop(usage(options$`scriptName`, parsedOptions))
        }

        if ( ! is.null(options$`aliases-primary`) && ! file.exists(options$`aliases-primary`) ) {
            printError(c("The file '", options$`aliases-primary`, "', specified as an argument to option '--aliases-primary' could not be found."))
            stop(usage(options$`scriptName`, parsedOptions))
        }

        if ( ! is.null(options$`aliases-secondary`) && ! file.exists(options$`aliases-secondary`) ) {
            printError(c("The file '", options$`aliases-secondary`, "', specified as an argument to option '--aliases-secondary' could not be found."))
            stop(usage(options$`scriptName`, parsedOptions))
        }

        #---> Try to generate output directory if it does not yet exist <---#
        dir.create(options$`output-directory`, recursive=TRUE, showWarnings = FALSE)

        #---> Validate column selector format and return vector of selected columns  <---#
        options$`reference` <- validateColumnSelection(options$`reference`)
        if ( length(options$`reference`[["error"]]) > 0 ) {
            invisible(lapply(options$`reference`[["error"]], function(error) {
                printError(error[["message"]])
            }))
            stop("\n", usage(options$`scriptName`, parsedOptions))
        } else { options$`reference` <- options$`reference`[["results"]] }

        if ( ! is.null(options$`query`) ) {
            options$`query` <- validateColumnSelection(options$`query`)
            if ( length(options$`query`[["error"]] ) > 0 ) {
                invisible(lapply(options$`query`[["error"]], function(error) {
                    printError(error[["message"]])
                }))
                stop(usage(options$`scriptName`, parsedOptions))
           } else { options$`query` <- options$`query`[["results"]] }
        }

        #---> Ensure user-specified values are legal <---#
        if ( ! options$`adjustP` %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr") ) {
            printError(c("Illegal argument specified for option '--adjustP'. Use one of the following: 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', or 'fdr'."))
            stop(usage(options$`scriptName`, parsedOptions))
        }

        if ( ! options$`common-dispersion` > 0 ) {
            printError(c("Illegal argument specified for option '--mds-top-n'. Only positive numbers allowed."))
            stop(usage(options$`scriptName`, parsedOptions))
        }

        if ( ! options$`mds-top-n` > 0 ) {
            printError(c("Illegal argument specified for option '--mds-top-n'. Only positive integers allowed."))
            stop(usage(options$`scriptName`, parsedOptions))
        }

    #---> // BODY <---#

    #---> RETURN VALUE <---#
    return(options)

}
#-----------------------#
parseOptions <- function() {

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
            make_option("--counts", action="store", type="character", default=NULL, help="Tab-separated data table containing integer count data for each feature (rows) and condition (column). Required.", metavar="FILE"),
            make_option("--abundances", action="store", type="character", default=NULL, help="Tab-separated data table containing abundance data for each feature (rows) and condition (column).", metavar="FILE"),
            make_option("--output-directory", action="store", type="character", default=NULL, help="Absolute or relative filename of output directory. It will be attempted to create the directory. If not specified, the current working directory will be used.", metavar="PATH"),
            make_option("--name", action="store", type="character", default="edgeR_comparison", help="Descriptive name for analysis.", metavar="STRING"),
            make_option("--reference", action="store", type="character", default=NULL, help="Fields of counts data table that shall be used as the reference sample group. See notes for input format. Required.", metavar="FIELDS"),
            make_option("--query", action="store", type="character", default=NULL, help="Fields of counts data table that shall be used as the query sample group. See notes for input format. By default, all columns but the ones specified to '--reference' are used.", metavar="FIELDS"),
            make_option("--lookup-table", action="store", type="character", default=NULL, help="Tab-separated lookup table indicating a secondary feature for each primary feature identifier (those present in the data tables supplied to '--counts' and '--abundances'). A header is expected. Data for primary features are summed per secondary feature and, wherever possible, analyses are going to be performed for secondary features too. A typical example is a transcript ID -> gene ID lookup table.", metavar="FILE"),
            make_option("--aliases-primary", action="store", type="character", default=NULL, help="Tab-separated lookup table indicating main IDs for the primary features and mappings to one or more corresponding aliases (one per column). A header is expected.", metavar="FILE"),
            make_option("--aliases-secondary", action="store", type="character", default=NULL, help="Tab-separated lookup table indicating main IDs for the secondary features and mappings to one or more corresponding aliases (one per column). A header is expected. Ignored if '--lookup-table' is not specified.", metavar="FILE"),
            make_option("--common-dispersion", action="store", type="numeric", default=0.1, help="Common dispersion value to be used when no replicates are available for reference and query (total N = 2).", metavar="FLOAT"),
            make_option("--adjustP", action="store", type="character", default="BH", help="Method for adjusting P-values for multiple testing. One of 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', or 'fdr'.", metavar="STRING"),
            make_option("--mds-top-n", action="store", type="integer", default=500, help="Number of features to be used for the multi-dimensional scaling plot.", metavar="INT"),
            make_option("--mds-pairwise", action="store_true", default=FALSE, help="When preparing multi-dimensional scaling plots, use the top N features (see '--mds-top-n') for each pairwise comparison. By default, the common top N features are used."),
            make_option("--verbose", action="store_true", default=FALSE, help="Print log messages to STDOUT."),
            make_option("--version", action="store_true", default=FALSE, help="Show version information and exit."),
            make_option("--license", action="store_true", default=FALSE, help="Show license information and exit."),
            make_option("--usage", action="store_true", default=FALSE, help="Show this screen and exit.")
        )

        #---> Parse options <---#
        parsedOptions <- OptionParser(option_list = option_list)
        options <- parse_args(parsedOptions)

        #---> Get script name <---#
        options$`scriptName` <- getScriptName()

        #---> Show usage / version / license <---#
        if ( options$`usage` ) {
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

        #---> Validate options <---#
        options <- validateOptions(options, parsedOptions)

    #---> // BODY <---#

    #---> RETURN VALUE <---#
    return(options)

}
#===================================#
#  // CLI OPTION-RELATED FUNCTIONS  #
#===================================#


#======================#
#  OTHER FUNCTIONS //  #
#======================#
loadData <- function(counts, abundances=NULL, lookup=NULL, aliasesPrimary=NULL, aliasesSecondary=NULL, verbose=FALSE) {

# [Description]
#     Load required data tables.
# [Parameters]
#     counts           : Tab-separated data table containing integer count data for each feature
#                        (rows) and condition (column)
#     abundances       : Tab-separated data table containing abundance data for each feature (rows)
#                        and condition (column)
#     lookup           : Tab-separated lookup table indicating a secondary feature for each primary
#                        feature identifier (those present in the data tables supplied to 'counts'
#                        and 'abundances')
#     aliasesPrimary   : Tab-separated lookup table indicating main IDs for the primary features
#                        and mappings to one or more corresponding aliases (one per column)
#     aliasesSecondary : Primary : Tab-separated lookup table indicating main IDs for the secondary
#                        features and mappings to one or more corresponding aliases (one per column)
#     verbose          : Print status messages (default: FALSE)
# [Return value]
#     List of data tables.
# [Dependencies]
#     printStatus()

    #---> STATUS MESSAGE <---#
    printStatus(c("Loading data..."), verbose)

    #---> INITIALIZE FUNCTION VARIABLES <---#
    data <- list()

    #---> BODY // <---#

        #---> Load counts data table <---#
        data$`counts` <- read.delim(counts)

        #---> Load abundances data table <---#
        if ( ! is.null(abundances) ) data$`abundances` <- read.delim(abundances)

        #---> Load lookup table <---#
        if ( ! is.null(lookup) ) data$`lookup` <- read.delim(lookup)

        #---> Load aliases for primary features <---#
        if ( ! is.null(aliasesPrimary) ) data$`primary` <- read.delim(aliasesPrimary)

        #---> Load aliases for primary features <---#
        if ( ! is.null(aliasesSecondary) ) data$`secondary` <- read.delim(aliasesSecondary)

    #---> // BODY <---#

    #---> STATUS MESSAGE <---#
    printStatus(c("Data loaded."), verbose)

    #---> RETURN VALUE <---#
    return(data)

}
#-----------------------#
aggregateCountsAbundances <- function(data, verbose=FALSE) {

# [Description]
#     Aggregate Generate MDS plot (`limma` package) in PDF format.
# [Parameters]
#     data           : A list of data objects containing a count table 'counts', a lookup table
#                      'lookup' that maps primary to secondary IDs, and optionally a data table
#                      'abundances'
#     verbose        : Print status messages (default: FALSE)
# [Return value]
#     Updated data object list.
# [Dependencies]
#     printStatus()

    #---> STATUS MESSAGE <---#
    printStatus(c("Aggregating counts and abundances..."), verbose)

    #---> // BODY <---#

        #---> Group primary IDs by secondary ID <---#
        featGroups <- aggregate(as.character(data$`lookup`[,1])~as.character(data$`lookup`[,2]), data$`lookup`, c)
        featGroups <- setNames(featGroups[,2], featGroups[,1])

        #---> Aggregate counts <---#
        data$`countsGroup` <- t(sapply(names(featGroups), function(groupName) {
            colSums(data$`counts`[match(featGroups[[groupName]], rownames(data$`counts`)),])
        }))

        #---> Aggregate abundances <---#
        if ( ! is.null(data$`abundances`) ) {
            data$`abundancesGroup` <- t(sapply(names(featGroups), function(groupName) {
                colSums(data$`abundances`[match(featGroups[[groupName]], rownames(data$`abundances`)),])
            }))
        }

        #---> Add reverse lookup table to data object list <---#
        data$`lookupRev` <- featGroups

    #---> // BODY <---#

    #---> STATUS MESSAGE <---#
    printStatus(c("Counts and abundances aggregated."), verbose)

    #---> RETURN VALUE <---#
    return(data)

}
#-----------------------#
edgeRCompareTwoGroups <- function(counts, reference, query=setdiff(1:ncol(counts), reference), adjustP="BH", name=match.call()[[1]], abundances=NULL, commonDisp=0.1, verbose=FALSE) {

# [Description]
#     Identify differentially expressed features between two sample groups with `edgeR`
# [Parameters]
#     counts           : Matrix or data frame containing count data for each feature (rows) and
#                        condition (column)
#     reference        : Vector of column names or indices for the selection of reference data
#                        columns from data supplied to 'counts'
#     query            : Vector of column names or indices for the selection of query data columns
#                        from data supplied to 'counts' (default: all columns but the ones specified
#                        to parameter 'reference')
#     adjustP          : Method for adjusting P values for multiple testing; refer to
#                        `stats::p.adjust` for allowed values and further information (default:
#                        "BH", i.e. 'Benjamini-Hochberg').
#     name             : Identifier for analysis / sample contrast (default: generic name derived
#                        from the name of this function)
#     abundances       : Matrix or data frame containing abundance data foe each feature (rows) and
#                        condition (column); if supplied, assumes the exact same layout as object
#                        supplied to 'counts' (default: NULL)
#     commonDisp       : Dispersion value that is set when no replicates are available for either
#                        sample (default: 0.1).
#     verbose          : Print status messages (default: FALSE)
# [Return value]
#     A list of tables, summaries and information about the analysis. Inspect resulting object with 
#     str(object).
# [Dependencies]
#     edgeR, printStatus(), loadPackages()

    #---> STATUS MESSAGE <---#
    printStatus(c("Determining differentially expressed features for contrast '", name, "'...'"), verbose)

    #---> LOAD FUNCTION DEPENDENCIES <---#
    loadPackages("edgeR")

    #---> INITIALIZE FUNCTION VARIABLES <---#
    results <- list()

    #---> BODY // <---#

        #---> Add contrast name to results object <---#
        results$`contrast_name` <- name

        #---> Add column names for reference and query to results object <---#
        results$`reference` <- colnames(counts)[reference]
        results$`query` <- colnames(counts)[query]

        #---> Select reference and query data and set grouping factor <---#
        df <- counts[ , c(reference, query) ]
        group <- factor(c(rep(1, length(query)), rep(2, length(reference))))

        #---> Add data identifier to count table <---#
        colnames(df) <- paste(colnames(df), "count", sep=".")

        #---> Add grouping factor, counts and, if available, abundances to results object <---#
        results$`group` <- group
        results$`counts` <- df
        if ( ! is.null(abundances) ) {
            colnames(abundances) <- paste(colnames(abundances), "abundance", sep=".")
            results$`abundances` <- abundances[, c(reference, query) ]
        }

        #---> Initialize DGEList data object <---#
        dge_ls <- DGEList(df, group=group)

        #---> Calculate normalization and dispersion factors <---#
        dge_ls <- calcNormFactors(dge_ls)
        if (length(results$group) == 2) {
            dge_ls$common.dispersion <- commonDisp
        } else {
            dge_ls <- estimateCommonDisp(dge_ls)
            dge_ls <- estimateTagwiseDisp(dge_ls)
        }

        #---> Add DGEList object to results objects <---#
        results$`DGEList` <- dge_ls

        #---> Compute differences in the means of counts between two groups <---#
        dge_ex <- exactTest(dge_ls)

        #---> Correct P-values for multiple testing and determine DE features <---#
        dge_ex$table$`fdr` <- p.adjust(dge_ex$table$PValue, method=adjustP, n=nrow(df))
        dge_ex$table$`isDE` <- decideTestsDGE(dge_ex)

        #---> Add DGEExact object to results objects <---#
        results$`DGEExact` <- dge_ex

        #---> Add feature IDs & sample overview to results object <---#
        results$`featIDs` <- rownames(dge_ex$table)
        results$`sample_overview` <- cbind(features=rownames(dge_ls$samples), dge_ls$samples)
        rownames(results$`sample_overview`) <- NULL

        #---> Add numbers of differentially expressed features to results object <---#
        results$`n_de` <- sum(! results$`isDE` == 0)
        results$`n_enr_depl` <- summary(results$`isDE`)

    #---> // BODY <---#

    #---> STATUS MESSAGE <---#
    printStatus("Determination of differentially expressed features completed.", verbose)

    #---> RETURN VALUE <---#
    return(results)

}
#-----------------------#
plot.MDS <- function(counts, outDir, topN=500, pairwise=TRUE, names=NULL, group=NULL, title=NULL, xlab="Dimension 1", ylab="Dimension 2", verbose=FALSE) {

# [Description]
#     Generate MDS plot (`limma` package) in PDF format.
# [Parameters]
#     counts         : Matrix of counts (data.frame or DGEList object)
#     outDir         : Absolute or relative filename of output directory
#     topN           : Top N features are considered (default: 500)
#     pairwise       : Choose top genes separately for each comparison (default: TRUE)
#     names          : Sample names (default: 1 to N)
#     group          : Factor for sample groups (default: no grouping)
#     title          : Plot title and filename prefix
#     xlab           : Label for x axis (default: "Dimension 1")
#     ylab           : Label for y axis (default: "Dimension 2")
#     verbose        : Print status messages (default: FALSE)
# [Return value]
#     Return value of plot function
# [Dependencies]
#     edgeR, printStatus(), loadPackages()

    #---> STATUS MESSAGE <---#
    printStatus(c("MDS plot is being generated..."), verbose)

    #---> LOAD FUNCTION DEPENDENCIES <---#
    loadPackages("edgeR")

    #---> BODY // <---#

        # Set gene selection method
        if ( pairwise ) method <- "pairwise" else method <- "common"

        # Set default group if not specified
        if ( is.null(group) ) group <- factor(rep(1, ncol(counts)))

        # Convert count table to `DGEList` if specified as `data.frame`
        if (! class(counts) == "DGEList") counts <- DGEList(counts, group)

        # Set output filename
        outFile <- file.path(outDir, paste(title, "MDS", "pdf", sep="."))

        # Open plotting device
        pdf(file=outFile)

        # Generate plot
        return <- plotMDS(counts, top=topN, gene.selection=method, labels=names, col=as.integer(group), main=title, xlab=xlab, ylab=ylab)

        # Close plotting device
        dev.off()

    #---> // BODY <---#

    #---> STATUS MESSAGE <---#
    printStatus(paste("MDS plot written to file '", outFile, "'.", sep=""), verbose)

    #---> RETURN VALUE <---#
    return(return)

}
#-----------------------#
plot.Smear <- function(object, deTags, outDir, title=NULL, xlab="Average logCPM", ylab="logFC", verbose=FALSE) {

# [Description]
#     Generate MDS plot (`limma` package) in PDF format.
# [Parameters]
#     object         : limma/edgeR DGEList or DGEExact object
#     deTags         : Names of differentially expressed features
#     outDir         : Absolute or relative filename of output directory
#     group          : Factor for sample groups (default: no grouping)
#     title          : Plot title and filename prefix
#     xlab           : Label for x axis (default: "Dimension 1")
#     ylab           : Label for y axis (default: "Dimension 2")
#     verbose        : Print status messages (default: FALSE)
# [Return value]
#     Return value of plot function
# [Dependencies]
#     edgeR, printStatus(), loadPackages()

    #---> STATUS MESSAGE <---#
    printStatus(c("Smear plot is being generated..."), verbose)

    #---> LOAD FUNCTION DEPENDENCIES <---#
    loadPackages("edgeR")

    #---> BODY // <---#

        # Set output filename
        outFile <- file.path(outDir, paste(title, "smear", "pdf", sep="."))

        # Open plotting device
        pdf(file=outFile)

        # Generate plot
        return <- plotSmear(object, de.tags=deTags, main=title, xlab=xlab, ylab=ylab)

        # Close plotting device
        dev.off()

    #---> // BODY <---#

    #---> STATUS MESSAGE <---#
    printStatus(paste("Smear plot written to file '", outFile, "'.", sep=""), verbose)

    #---> RETURN VALUE <---#
    return(return)

}
#-----------------------#
edgeRWriteResults <- function(results, outDir, topMDS=500, pairwiseMDS=FALSE, aliasesPrimary=NULL, lookup=NULL, lookupRev=NULL, aliasesSecondary=NULL, verbose=FALSE) {

# [Description]
#     Print summary, result tables and plot graphs for edgeR analysis.
# [Parameters]
#     results          : List containing element 'primary' and, optionally, 'secondary', both of
#                        representing results lists returned by edgeRCompareTwoGroups()
#     outDir           : Absolute or relative filename of output directory
#     topMDS           : Number of top features to consider for MDS plot (default: 500)
#     pairwiseMDS      : Choose top genes separately for each comparison in MDS plot (default:
#                        FALSE)
#     aliasesPrimary   : Matrix or data frame containing main IDs for primary features and mappings
#                        to one or more corresponding aliases (one per column)
#     lookup           : Matrix or data frame containing a mapping from primary to secondary IDs
#                        (one column each)
#     lookupRev        : List of vectors representing secondary to primary ID mappings
#     aliasesSecondary : Matrix or data frame containing main IDs for primary features and mappings
#                        to one or more corresponding aliases (one per column)
#     verbose          : Print status messages (default: FALSE)
# [Return value]
#     N/A
# [Dependencies]
#     printStatus()

    #---> STATUS MESSAGE <---#
    printStatus(c("Results are being written..."), verbose)

    #---> BODY // <---#

        #---> Write sample summary tables <---#
        outFile <- file.path(outDir, paste(results$primary$contrast_name, "samples", "tab", sep="."))
        write.table(results$primary$`sample_overview`, outFile, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
        if ( ! is.null(results$secondary) ) {
            outFile <- file.path(outDir, paste(results$secondary$contrast_name, "samples", "tab", sep="."))
            write.table(results$secondary$`sample_overview`, outFile, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
        }

        #---> Get primary feature IDs <---#
        if ( ! is.null(aliasesPrimary) ) {
            primFeatIDs <- merge(results$primary$featIDs, aliasesPrimary, by=1)
            colnames(primFeatIDs) <- colnames(aliasesPrimary)
        } else {
            primFeatIDs <- data.frame(feature=results$primary$featIDs)
        }

        #---> Get secondary feature IDs & merge with primary IDs <---#
        if ( ! is.null(results$secondary) ) {
            if ( ! is.null(aliasesSecondary) ) {
                secFeatIDs <- merge(results$secondary$featIDs, aliasesSecondary, by=1)
                colnames(secFeatIDs) <- colnames(aliasesSecondary)
            } else {
                secFeatIDs <- data.frame(feature=results$secondary$featIDs)
            }
        }

        #---> Assemble feature table: Primary <---#

            #---> Compile results table for all primary features <---#
            tmp <- merge(primFeatIDs, results$primary$DGEExact$table, by.x=1, by.y=0)
            rownames(tmp) <- NULL

            #---> Add counts and abundances <---#
            tmp <- merge(tmp, results$primary$`counts`, by.x=1, by.y=0)
            if ( ! is.null(results$primary$abundances) ) tmp <- merge(tmp, results$primary$`abundances`, by.x=1, by.y=0)

            #---> Add to results object <---#
            results$primary$`all_feats` <- tmp

        #---> Assemble feature table: Secondary <---#

        if ( ! is.null(results$secondary) ) {

            #---> Compile results table for all secondary features <---#
            tmp <- merge(secFeatIDs, results$secondary$DGEExact$table, by.x=1, by.y=0)
            rownames(tmp) <- NULL

            #---> Add counts and abundances <---#
            tmp <- merge(tmp, results$secondary$`counts`, by.x=1, by.y=0)
            if ( ! is.null(results$secondary$abundances) ) tmp <- merge(tmp, results$secondary$`abundances`, by.x=1, by.y=0)

            #---> Add to results object <---#
            results$secondary$`all_feats` <- tmp

        }

        #---> Add secondary feature stats to primary feature table <---#

            #---> Add secondary feature IDs & DE flag to primary feature table <---#
            if ( ! is.null(results$secondary) ) {
                tmp <- merge(data$lookup, secFeatIDs, by.x=2, by.y=1)
                tmp <- merge(tmp, setNames(as.integer(results$secondary$all_feats$isDE), secFeatIDs[,1]), by.x=1, by.y=0)[,c(2,1,3,4)]
                colnames(tmp)[4] <- "isGroupDE"
                results$primary$`all_feats` <- merge(results$primary$`all_feats`, tmp, by=1)
            }

        #---> Add primary feature stats to secondary feature table <---#

            #---> Add number of primary features and number of DE primary features to secondary feature table <---#
            if ( ! is.null(results$secondary) ) {
                primFeats <- sapply(lookupRev, length)
                primFeatsDE <- sapply(names(lookupRev), function(groupName) {
                    sum(abs(results$primary$all_feats$isDE[match(lookupRev[[groupName]], results$primary$all_feats[,1])]))
                })
                tmp <- merge(primFeats, primFeatsDE, by=0)
                colnames(tmp) <- c(colnames(results$secondary$`all_feats`)[1], "inGroup", "inGroupDE")
                results$secondary$`all_feats` <- merge(results$secondary$`all_feats`, tmp, by=1)
            }

        #---> Write tables with all features <---#
        outFile <- file.path(outDir, paste(results$primary$contrast_name, "all_features", "tab", sep="."))
        write.table(results$primary$`all_feats`, outFile, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

        outFile <- file.path(outDir, paste(results$secondary$contrast_name, "all_features", "tab", sep="."))
        write.table(results$secondary$`all_feats`, outFile, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

        #---> Subset DE features, order by P-value, add to results object & write table <---#
        results$primary$`de_feats` <- results$primary$all_feats[! results$primary$all_feats$isDE == 0, ]
        results$primary$`de_feats` <- results$primary$de_feats[order(results$primary$de_feats$PValue), ]
        outFile <- file.path(outDir, paste(results$primary$contrast_name, "de_features", "tab", sep="."))
        write.table(results$primary$`de_feats`, outFile, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

        results$secondary$`de_feats` <- results$secondary$all_feats[! results$secondary$all_feats$isDE == 0, ]
        results$secondary$`de_feats` <- results$secondary$de_feats[order(results$secondary$de_feats$PValue), ]
        outFile <- file.path(outDir, paste(results$secondary$contrast_name, "de_features", "tab", sep="."))
        write.table(results$secondary$`de_feats`, outFile, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

        #---> Generate MDS plots <---#
        if ( length(results$primary$group) > 2 ) {

            results$mds <- with(results$primary, plot.MDS (
                counts=DGEList,
                outDir=outDir,
                topN=topMDS,
                pairwise=pairwiseMDS,
                names=sample_overview$features,
                group=group,
                title=contrast_name,
                verbose=verbose
            ))
            if ( ! is.null(results$secondary) ) {
                results$mdsGroup <- with(results$secondary, plot.MDS (
                    counts=DGEList,
                    outDir=outDir,
                    topN=topMDS,
                    pairwise=pairwiseMDS,
                    names=sample_overview$features,
                    group=group,
                    title=contrast_name,
                    verbose=verbose
                ))
            }
        }

        #---> Generate smear plots <---#
        results$smear <- with(results$primary, plot.Smear (
            object=DGEExact,
            deTags=results$primary$de_feats[,1],
            outDir=outDir,
            title=contrast_name,
            verbose=verbose
        ))
        if ( ! is.null(results$secondary) ) {
            results$smearGroup <- with(results$secondary, plot.Smear (
                object=DGEExact,
                deTags=results$secondary$de_feats[,1],
                outDir=outDir,
                title=contrast_name,
                verbose=verbose
            ))
        }

        #---> Write summaries <---#
# TODO

        #---> Save results <---#
        outFile <- file.path(outDir, paste(results$primary$contrast_name, "results", "R", sep="."))
        save(results, file=outFile)

    #---> // BODY <---#

    #---> STATUS MESSAGE <---#
    status <- paste("Output files have been written to directory '", outDir, "'.", sep="")
    printStatus(status, verbose)

}
#======================#
#  // OTHER FUNCTIONS  #
#======================#


#===========#
#  MAIN //  #
#===========#

#---> PROCESS GLOBAL OPTIONS <---#
opt <- parseOptions();

#---> STATUS MESSAGE <---#
printStatus(c("Starting '", opt$`scriptName`, "'..."), opt$`verbose`);

#---> INITIALIZE GLOBAL VARIABLES <---#
data <- list()
results <- list()

#---> BODY // <---#

    #---> Load data <---#
    data <- with(opt, loadData(
        counts=counts,
        abundances=opt$abundances,
        lookup=opt$`lookup-table`,
        aliasesPrimary=opt$`aliases-primary`,
        aliasesSecondary=opt$`aliases-secondary`,
        verbose=verbose
    ))

    #---> Aggregate data <---#
    if ( ! is.null (opt$`lookup-table`) ) data <- aggregateCountsAbundances(data=data, verbose=opt$`verbose`)

    #----> Set query columns if not specified by user <---#
    if ( is.null(opt$`query`) ) {
        opt$`query` <- setdiff(1:ncol(data$`counts`), opt$`reference`)
        printStatus(c("No columns specified for query group. Using columns:\n", toString(opt$`query`)), opt$`verbose`)
    }

    #---> Run edgeR comparison on primary features <---#
    results$primary <- with(opt, edgeRCompareTwoGroups(
        counts=data$`counts`,
        reference=reference,
        query=query,
        adjustP=adjustP,
        name=name,
        abundances=data$`abundances`,
        commonDisp=`common-dispersion`,
        verbose=verbose
    ))

    #---> Run edgeR comparison on secondary features <---#
    if ( ! is.null(opt$`lookup-table`) ) {
        results$secondary <- with(opt, edgeRCompareTwoGroups(
            counts=data$`countsGroup`,
            reference=reference,
            query=query,
            adjustP=adjustP,
            name=paste(name, "group", sep="-"),
            abundances=data$`abundancesGroup`,
            commonDisp=`common-dispersion`,
            verbose=verbose
        ))
    }

    #---> Write results <---#
    with(opt, edgeRWriteResults(
        results=results,
        outDir=`output-directory`,
        topMDS=`mds-top-n`,
        pairwiseMDS=`mds-pairwise`,
        aliasesPrimary=data$`primary`,
        lookup=data$`lookup`,
        lookupRev=data$`lookupRev`,
        aliasesSecondary=data$`secondary`,
        verbose=verbose
    ))

#---> // BODY <---#

#---> STATUS MESSAGE <---#
printStatus("Done.", opt$`verbose`);

#---> SAVE SESSION <---#
outFile <- file.path(opt$`output-directory`, paste(opt$`name`, "session", "R", sep="."))
save.image(file=outFile)

#---> PROGRAM EXIT <---#
quit(save = "no", status=0, runLast=FALSE)

#===========#
#  // MAIN  #
#===========#
