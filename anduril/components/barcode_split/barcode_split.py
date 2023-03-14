#!/usr/bin/env python

## COMMENTS
# 

import component_skeleton.main
import sys, os
import subprocess

### Script to be executed from the anduril component_skeleton
def execute(cf):

### A. PREREQUISITES
# Get in and output arguments as specified in component.xml
    infolder = cf.get_input("infolder")
    infile = cf.get_input("infile")
    outfolder = cf.get_output("outfolder")
    barcode_split = cf.get_parameter("barcode_split_dir", "string")
    barcodes = cf.get_parameter("barcodes", "string")
    mismatches = cf.get_parameter("mismatches", "string")
    eol = cf.get_parameter("eol", "string")


# Check for presence of infolder(Anduril output from previous component) or infile (new user input)
# If infolder: replace the variable infile with the file in the infolder
    if infolder:
        if infile:
            print "\n[ERROR] infile and infolder given:\nPlease give only one input\n"
            return -1
    
        d = os.listdir(infolder)
        infile = os.path.join(infolder, d[0])
    else:
        if not infile:
            print "\n[ERROR] no infile or infolder given\n"
            return -1

# Create outfile names
    #base_name = str.rsplit(os.path.basename(infile),".")[0]
    #outfile = os.path.join(outfolder, base_name + "_no_bar.fastq")

# Make output directories
    os.mkdir(outfolder)
    unmatched = os.path.join(outfolder, "unmatched")
    os.mkdir(unmatched)
    print "\nReads without barcodes stored in folder\n %s.\n These will not be considered for further analysis\n." % unmatched
    
### B. CALL TO FUNCTIONS/SCRIPTS
# Make strings for optional part of system call, depending on specified parameters
# Allowed mismatches
    if isinstance(mismatches, str):
        if mismatches == "":
            mis = ""
            print "Use default setting for mismatches %s\n" % mis
        else:
            mis = "--mismatches " + mismatches
            print "Mismatches set to %s\n" % mis
    else:
        print "[ERROR] Option mismatches not specified correctly, must be string. \nEmpty string calls default specified in barcode_split.pl\n"
        return -1
        
# Barcode at 3' end (eol = "3_END") or 5' end of sequence
    if eol == "5_END":
        e = ""
        print "Barcodes at 5'end, option --eol DISabled: %s\n" % e
    else:
        if eol == "3_END":
            e = "--eol"
            print "Barcodes at 3'end, option --eol ENabled: %s\n" % e
        else:
            print "[ERROR] Option eol not specified correctly, set eol to \"5_END\" OR \"3_END\".\n"
            return -1

# Call actual component script with necessary input parameters
    proc = subprocess.Popen("perl %s --infile %s --bc \"%s\" --outdir %s --outdir_um %s %s --trim %s" % (barcode_split, infile, barcodes, outfolder, unmatched, e, mis),
                                                        stdout=subprocess.PIPE,
                                                        stderr=subprocess.PIPE,
                                                        shell=True
                                                        )
# Wait until script has been executed
    stdout_value, stderr_value = proc.communicate()
    print stdout_value
    print stderr_value

# Return error
    if proc.poll() > 0:
        print '\tstderr:', repr(stderr_value.rstrip())
        return -1

# Return success
    return 0

### C. Execute specified script from Anduril
component_skeleton.main.main(execute)

