#!/usr/bin/env python

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
    fastq_dir = cf.get_parameter("fastq_dump_dir", "string")
    minReadLen = cf.get_parameter("minReadLen", "int")

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

# Create outfile names
    base_name = str.rsplit(os.path.basename(infile),".")[0]
    outfile = os.path.join(outfolder, base_name + ".fastq")

# Make output directories
    os.mkdir(outfolder)
    
### B. CALL TO FUNCTIONS/SCRIPTS
# Call actual component script with necessary input parameters
    proc = subprocess.Popen("%s --minReadLen %s --stdout %s > %s" % (fastq_dir, minReadLen, infile, outfile),
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

