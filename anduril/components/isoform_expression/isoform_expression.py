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
    prog = cf.get_parameter("prog", "string")
    gene_trx_tab = cf.get_parameter("gene_trx_tab", "string")
    likelihood_fct = cf.get_parameter("likelihood_fct", "string")
    var = cf.get_parameter("var", "int")

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
    outfile = os.path.join(outfolder, base_name + "_iso" + ".tab")

# Make output directory
    os.mkdir(outfolder)


# Call actual component script with necessary input parameters
    proc = subprocess.Popen("export LD_LIBRARY_PATH='/import/bc2/soft/app/gcc/4.6.3/Linux/lib64':$LD_LIBRARY_PATH; %s -g %s -e %s -s %i %s > %s" % (prog, gene_trx_tab, likelihood_fct, var, infile, outfile),
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
