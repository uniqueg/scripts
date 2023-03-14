#!/usr/bin/env python

## COMMENTS
# adapter options not yet included in call to segemehl
# polyA and hardclip not yet accessible by user

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
    unmatched_folder = cf.get_output("unmatched_folder")
    segemehl = cf.get_parameter("segemehl", "string")
    index = cf.get_parameter("index", "string")
    genome = cf.get_parameter("genome", "string")
    threads = cf.get_parameter("threads", "int")
    differences = cf.get_parameter("differences", "int")
    accuracy = cf.get_parameter("accuracy", "int")
    clipacc = cf.get_parameter("clipacc", "int")
    prime3 = cf.get_parameter("prime3", "string")
    prime5 = cf.get_parameter("prime5", "string")


# Check for presence of infolder(Anduril output from previous component) or infile (new user input)
# If infolder: replace the variable infile with the file in the infolder
    if infolder:
        if infile:
            print "\n[ERROR] Infile and infolder given:\nPlease give only one input\n"
            return -1
    
        d = os.listdir(infolder)
        infile = os.path.join(infolder, d[0])
    else:
        if not infile:
            print "\n[ERROR] No infile or infolder given\n"
            return -1

# Create outfile names
    base_name = str.rsplit(os.path.basename(infile),".")[0]
    outfile = os.path.join(outfolder, base_name + "_segemehl.sam")
    unmatched_file = os.path.join(unmatched_folder, base_name + "_segemehl_unmatched")

# Make output directories
    os.mkdir(outfolder)
    os.mkdir(unmatched_folder)
    
### B. CALL TO FUNCTIONS/SCRIPTS
# Make strings for calling optional parameters
    if prime3 != "N/A":
        prime = "--prime3 " + prime3
    else:
        if prime5 != "N/A":
            prime = "--prime5 " + prime5
        else:
            prime = ""
    
    print "\nAdapter removal set to %s\n" % prime

# Call actual component script with necessary input parameters
    proc = subprocess.Popen("%s --index %s --database %s --query %s %s --differences %i --accuracy %i --threads %i --polyA --hardclip --clipacc %i --silent --outfile %s --nomatchfilename %s" % (segemehl, index, genome, infile, prime, differences, accuracy, threads, clipacc, outfile, unmatched_file),
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

