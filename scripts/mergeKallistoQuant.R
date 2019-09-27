#!/usr/bin/env Rscript


#=====================#
#  ISSUES & IDEAS //  #
#=====================#
# TODO: Alternative source for sample names
# TODO: Resolve sample name conflicts
# TODO: Validate options
# TODO: Make script generic (remove/change defaults, remove references to Kallisto)
#=====================#
#  // ISSUES & IDEAS  #
#=====================#


#========================#
#  GENERIC FUNCTIONS //  #
#========================#
getLogDateTime <- function(template = '%Y/%m/%d %H:%M:%S') {

# [Description]
#     Returns current date/time stamp in a format suitable for log entries.
# [Options]
#     template       : POSIX time template string (default: "%Y/%m/%d %H:%M:%S")
# [Return value]
#     Date/time string
# [Dependencies]

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
# [Options]
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
# [Options]
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
# [Options]
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

        #---> Print message <---#
        if (toScreen) { write(message, stderr()) }

    #---> // BODY <---#

}
#-----------------------#
loadPackages <- function(packages, install = FALSE, repo = "http://cran.us.r-project.org", status = FALSE) {

# [Description]
#     Loads one or more packages and either dies when packages are unavailable
#     or tries to install them.
# [Options]
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
            if ( ! success )

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

    #---> // BODY <---#

    #---> STATUS MESSAGE <---#
    printStatus("Dependencies loaded...", status)

}
#-----------------------#
dataFrameFromNamedVectors <- function(ls, nms=NULL, all=TRUE) {

# [Description]
#     Create merged dataframe from list of named vectors
# [Options]
#     ls             : List of named vectors
#     nms            : Vector providing column names for the resulting table (length of vector and 'ls' has to be identical -> currently not validated!) (default: use serial number 1 to n
#     all            : Whether NAs shall be introduced when not all vectors have identical names (see '?merge')
# [Return value]
#     Merged dataframe
# [Dependencies]
#     N/A

    #---> BODY // <---#

        #---> Generate column names if not specified <---#
        if (is.null(nms)) nms <- seq(1, length(ls))

        #---> Coerce first vector to dataframe and set name <---#
        df <- as.data.frame(ls[[1]])
        colnames(df) <- nms[[1]]

        #---> Iterate over remaining vectors <---#
        for (i in 2:length(ls)) {

            #---> Merge current vector with existing dataframe using the row names (dataframe) and names (vector) <---#
            df <- merge(df, ls[[i]], by=0, all=all)

            #---> Reset dataframe row names, remove column & set name <---#
            row.names(df) <- df$Row.names
            df <- df[,-1]
            colnames(df)[i] <- nms[i]

        }

    #---> // BODY <---#

    #---> RETURN VALUE <---#
    return(df)

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
# [Options]
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
# [Options]
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
# [Options]
#     N/A
# [Return value]
#     License string
# [Dependencies]
#     N/A

    #---> BODY // <---#

        #---> Build string <---#
        license <- paste(c(
'
The MIT License (MIT)

Copyright (c) 2015 Alexander Kanitz

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
SOFTWARE.
'
        ), collapse="")

    #---> // BODY <---#

    #---> RETURN VALUE <---#
    return(license)

}
#-----------------------#
userVersion <- function(scriptName) {

# [Description]
#     Build and return version string.
# [Options]
#     scriptName     : Name of the script for which a usage string shalle be generated
# [Return value]
#     Version string
# [Dependencies]
#     N/A

    #---> BODY // <---#

        #---> Build string <---#
        version <- paste(c(
scriptName, ", v1.0 (Sep 2, 2015)
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
# [Options]
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
    ", scriptName, " [OPTIONS]

[DESCRIPTION]
    Recursively finds Kallisto-based abundance estimates and merges them into a single table.

[DEPENDENCIES]
    optparse

[OPTIONS]
", formatOptions(parsedOptions), "

[USAGE NOTES]
    (1) When supplied with the '--plaintext' option, the current Kallisto version (v0.42.2)
        produces an output directory with a file 'abundance.tsv' containing the inferred
        abundances. By default, the names of the Kallisto output directories are used as sample
        names when producing the merged table of abundance estimates. Use '--sampleNamePrefix'
        and '--sampleNameSuffix' to only extract parts of these directory names for sample name
        generation. Note that conflicts in sample names are not resolved (currently untested)."
        ), collapse="")

    #---> // BODY <---#

    #---> RETURN VALUE <---#
    return(usage)

}
#-----------------------#
validateOptions <- function(options, parsedOptions) {

# [Description]
#     Validates command-line options.
# [Options]
#     options        : List of options and values, such as is produced by optparse::parse_args()
#     parsedOptions  : An object of class OptionParser, returned by optparse::OptionParser()
# [Return value]
#     N/A
# [Dependencies]
#     printError(), printWarning()


    #---> BODY // <---#

        # Currently no options are validated

    #---> // BODY <---#

}
#-----------------------#
parseOptions <- function() {

# [Description]
#     Returns validated command-line options.
# [Options]
#     N/A
# [Return value]
#     List of CLI options, as well as 'scriptName'
# [Dependencies]
#     optparse, loadPackages(), usage(), userVersion(), userLicense(), validateOptions()

    #---> LOAD SUBROUTINE DEPENDENCIES <---#
    loadPackages("optparse", FALSE)

    #---> BODY // <---#

        #---> Generate list of options <---#
        option_list <- list(
            make_option("--inputDir", action="store", type="character", default=".", help="Root directory containing one or more Kallisto output directories.", metavar="DIR"),
            make_option("--outFile", action="store", type="character", default="./kallistoQuant.tab", help="Output filename.", metavar="FILE"),
            make_option("--filenameGlob", action="store", type="character", default="abundance.tsv", help="Name of abundance estimate files in Kallisto output directories.", metavar="GLOB"),
            make_option("--colSelect", action="store", type="character", default="est_counts", help="Name of column whose values shall be extracted (default: 'est_counts').", metavar="STRING"),
            make_option("--round", action="store_true", default=FALSE, help="Whether values shall be rounded to the next integer."),
            make_option("--sampleNamePrefix", action="store", type="character", default="", help="Sample names are extracted from Kallisto output directories. The specified prefix is ignored. See also '--sampleNameSuffix'.", metavar="STRING"),
            make_option("--sampleNameSuffix", action="store", type="character", default="", help="The specified prefix is ignored during sample name extraction. See '--sampleNamePrefix'.", metavar="STRING"),
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

        #---> Verify options <---#
        validateOptions(options, parsedOptions)

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
main <- function(opt) {

# [Description]
#     Generic 'main' function
# [Options]
#     opt            : Options list
# [Return value]
#     N/A
# [Dependencies]
#     N/A

    #---> BODY // <---#

        #---> Recursively find files of specified name <---#
        filePaths <- sort(dir(opt$`inputDir`, pattern=glob2rx(opt$`filenameGlob`), recursive=TRUE, full.names=TRUE))

        #---> Extract sample names by removing specified pre- and suffix from 'dirnames' <---#
        sampleNames <- gsub(opt$`sampleNameSuffix`, "", gsub(opt$`sampleNamePrefix`, "", basename(dirname(filePaths))))

        #---> Read dataframes <---#
        df_ls <- lapply(filePaths, read.table, header=TRUE, row.names=1, sep="\t", stringsAsFactors=FALSE)

        #---> Select column of interest and return as named vector <---#
        vec_ls <- lapply(df_ls, function(df) {
            setNames(df[, opt$`colSelect`], row.names(df))
        })

        #---> Round to nearest integer if desired <---#
        if ( opt$`round` ) { vec_ls <- lapply(vec_ls, round) }

        #---> Merge named vecor list into dataframe <---#
        df <- dataFrameFromNamedVectors(vec_ls, sampleNames)

        #---> Write dataframe <---#
        write.table(df, opt$`outFile`, col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE)

    #---> // BODY <---#

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

#---> BODY // <---#

    #---> Call 'main' function <---#
    main(opt)

#---> // BODY <---#

#---> STATUS MESSAGE <---#
printStatus("Done.", opt$`verbose`);

#---> PROGRAM EXIT <---#
quit(save = "no", status=0, runLast=FALSE)

#===========#
#  // MAIN  #
#===========#
