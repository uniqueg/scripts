#!/usr/bin/env python

import component_skeleton.main
import sys, os
import subprocess

### A. Script to be executed from the anduril component_skeleton
def execute(cf):

# Get in and output arguments as specified in component.xml
    infolder = cf.get_input("infolder")
    infile = cf.get_input("infile")
    unique_dir = cf.get_output("unique_mappers")
    dupl_dir = cf.get_output("duplicates")

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
    unique_file = os.path.join(unique_dir, base_name + "_unique" + ".sam")
    dupl_file = os.path.join(dupl_dir, base_name + "_dupl" + ".txt")

# Make output directory
    unique_mappers = unique_dir
    os.mkdir(unique_mappers)
    duplicates = dupl_dir
    os.mkdir(duplicates)

# Call actual component script with necessary input parameters
    proc = subprocess.Popen("perl remove_multi_mappers.pl --in %s --uniq %s --dupl %s" % (infile, unique_file, dupl_file),
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
