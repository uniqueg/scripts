#!/usr/bin/env python

import component_skeleton.main
import sys, os
import subprocess

### A. Script to be executed from the anduril component_skeleton
def execute(cf):

# Get in and output arguments as specified in component.xml
    infolder = cf.get_input("infolder")
    infile = cf.get_input("infile")
    outfolder = cf.get_output("outfolder")
    gtf = cf.get_parameter("gtf", "string")
    index = cf.get_parameter("index", "string")
    genome = cf.get_parameter("genome", "string")
    mode = cf.get_parameter("mode", "string")
    ignore_strand = cf.get_parameter("ignore_strand", "string")

# Check for presence of infolder(Anduril output from previous component) or infile (new user input)
# If infolder: replace the variable infile with the file in the infolder
    if infolder:
        if infile:
            print "infile and infolder given:\nPlease give only one input\n"
            return -1
    
        d = os.listdir(infolder)
        infile = os.path.join(infolder, d[0])
    else:
        if not infile:
            print "no infile of infolder given\n"
            return -1
# Set outfile name
    base_name = str.rsplit(os.path.basename(infile),".")[0]
    outfile = os.path.join(outfolder, base_name + "_counttable" + ".tab")

# Make output directory
    os.mkdir(outfolder)

# Call summarize_overlaps with --ignore_strand=false
    if ignore_strand == "FALSE":
        proc = subprocess.Popen("Rscript summarize_overlaps_GTF_BAM.R --gtf %s --bam %s --out %s --index %s --genome %s --mode %s" % (gtf, infile, outfile, index, genome, mode),
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            shell=True
                            )
# Call actual component script with necessary input parameters
    else:
        if ignore_strand == "TRUE":
            proc = subprocess.Popen("Rscript summarize_overlaps_GTF_BAM.R --gtf %s --bam %s --out %s --index %s --genome %s --mode %s --ignore-strand" % (gtf, infile, outfile, index, genome, mode),
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            shell=True
                            )
        else:
            print "option ignore_strand ambiguous\n"
            return -1

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

### B. Execute specified script from Anduril
component_skeleton.main.main(execute)
