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
    samtools_dir = cf.get_parameter("samtools_dir", "string")

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
# Set outfile name
    base_name = str.rsplit(os.path.basename(infile),".")[0]
    outfile = os.path.join(outfolder, base_name)

# Make output directory
    os.mkdir(outfolder)

# Make sorted BAM file
    proc = subprocess.Popen("%s view -bS %s | %s sort - %s" % (samtools_dir, infile, samtools_dir, outfile),
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

# Make name for outfile.bam
    #os.mkdir(os.path.join(os.path.split(outfolder)[0],"bai"))
    #bai_folder = os.path.join(os.path.split(outfolder)[0],"bai")
    bamfile = os.path.join(outfile + ".bam")
    #bai_file = os.path.join(bai_folder, base_name + ".bam.bai")

# Make bam.bai index file
    proc2 = subprocess.Popen("%s index %s" % (samtools_dir, bamfile),
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            shell=True
                            )

# Wait until script has been executed
    stdout_value, stderr_value = proc2.communicate()
    print stdout_value
    print stderr_value

# Return error
    if proc2.poll() > 0:
        print '\tstderr:', repr(stderr_value.rstrip())
        return -1

# Return success
    return 0

### B. Execute specified script from Anduril
component_skeleton.main.main(execute)

