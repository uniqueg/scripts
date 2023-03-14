#!/usr/bin/env python

import component_skeleton.main
import sys, os
import subprocess

### A. Script to be executed from the anduril component_skeleton
def execute(cf):

# Get in and output arguments as specified in component.xml
    infolder = cf.get_input("infolder")
    dupli_dir = cf.get_input("duplicates")
    infile = cf.get_input("infile")
    outfolder = cf.get_output("outfolder")
    picard_dir = cf.get_parameter("picard_dir", "string")
    jar_dir = cf.get_parameter("jar_dir", "string")
    filter = cf.get_parameter("filter", "string")
    validation_stringency = cf.get_parameter("validation_stringency", "string")

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
            print "no infile or infolder given\n"
            return -1
        
# fetch the duplicates file from the duplicates input folder
    dup = os.listdir(dupli_dir)
    duplicates = os.path.join(dupli_dir, dup[0])
        
# Set outfile name
    base_name = str.rsplit(os.path.basename(infile),".")[0]
    outfile = os.path.join(outfolder, base_name + "_trx_uniq_gen.sam")

# Make output directory
    os.mkdir(outfolder)

# Make sorted SAM file
    proc = subprocess.Popen("%s -Xms2g -jar %s/FilterSamReads.jar FILTER=%s READ_LIST_FILE=%s WRITE_READS_FILES=false VALIDATION_STRINGENCY=%s INPUT=%s OUTPUT=%s" % (jar_dir, picard_dir, filter, duplicates, validation_stringency, infile, outfile),
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

### B. Execute specified script from Anduril
component_skeleton.main.main(execute)

