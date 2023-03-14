#!/bin/env/user python

import jinja2
import yaml
import os, sys
import subprocess

### A. FUNCTIONS

def main():
    
    ### 1. GET VARIABLES AND LOAD PARAMETERS

    # Check for correct usage
    if not len(sys.argv) == 2:
        print 'Usage: ./setup_pipeline.py params.yaml\n'
        sys.exit(0)

    # Get params file and load parameters with yaml
    parameter_file = sys.argv[1]
    p = open(parameter_file, 'r')
    params = yaml.load(p)
    p.close()
    
    # Load template dir (with network_template.and, hosts_template.and etc.) into jinja2
    krinihome = params["KRINI_HOME"]
    andfolder = os.path.join(krinihome, "anduril")
    temp_env = jinja2.Environment(loader=jinja2.FileSystemLoader(andfolder))

    # Get output, log and network folders
    outfolder = params["DATA"]
    logfolder = params["LOG"]
    netfolder = params["NETWORK"]
    
    ### 2. RENDER TEMPLATE FILES
 
    ## B.1 HOSTS.conf
    # Import template_hosts.conf
    hosts_temp = temp_env.get_template('template_hosts.conf')
    # Render template file with parameters from params.yaml file
    rendered_hosts_template = hosts_temp.render(params)
    # Write the completed hosts.conf file
    o = open(os.path.join(netfolder, 'hosts.conf'), 'w')
    for l in rendered_hosts_template:
        o.write(l)
    o.close()    
    
     # Variable for hosts.conf file for anduril call:
    hosts = os.path.join(netfolder, 'hosts.conf')
    
    
    ## B.2 NETWORK.and
    # Import template_network.and
    network_temp = temp_env.get_template('template_network.and')
    # Fill template file with parameters from params.yaml file
    rendered_network_template = network_temp.render(params)
    # Write the completed network.and file
    o = open(os.path.join(netfolder, 'network.and'), 'w')
    for l in rendered_network_template:
        o.write(l)
    o.close()
    
    # Variable for network.and file for anduril call:
    network = os.path.join(netfolder, 'network.and')
    
    
    ## B.3 PART_1.and
    # Import template_network.and
    part_1_temp = temp_env.get_template('template_processFASTQ_part_1.and')
    # Fill template file with parameters from params.yaml file
    rendered_part_1_template = part_1_temp.render(params)
    # Write the completed network.and file
    o = open(os.path.join(netfolder, 'processFASTQ_part_1.and'), 'w')
    for l in rendered_part_1_template:
        o.write(l)
    o.close()
    
    ## B.4 PART_2.and
    # Import template_network.and
    part_2_temp = temp_env.get_template('template_processFASTQ_part_2.and')
    # Fill template file with parameters from params.yaml file
    rendered_part_2_template = part_2_temp.render(params)
    # Write the completed network.and file
    o = open(os.path.join(netfolder, 'processFASTQ_part_2.and'), 'w')
    for l in rendered_part_2_template:
        o.write(l)
    o.close()
    
    
    '''
    ### 3. START ANDURIL PIPELINE
    
    # TO DO: remove shell=True, then command and options has to be a sequence, not a string!    
    proc = subprocess.call("anduril run %s -c components -d %s --hosts %s --log %s" % (network, outfolder, hosts, logfolder), shell=True) 
    
    
    proc = subprocess.Popen("anduril run yaml_1st_test.and -c components -d import/bc2/home/zavolan/cherrmann/PIPELINE/TEST_OUTPUT --hosts yaml_hosts_1.conf",
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

    '''

    ### 4. Print some feedback to STDOUT
    
    print 'The following parameters were found in "%s":' % parameter_file
    for key in params:
        print key + ": " + `params[key]`
    
    print

    print 'Output will be stored in "%s".' % outfolder
    print 'Logs will be stored in "%s".' % logfolder

    print

    print 'The following network files were created in %s":' % netfolder
    print os.path.join(netfolder, 'hosts.conf')
    print os.path.join(netfolder, 'network.and')
    print os.path.join(netfolder, 'processFASTQ_part_1.and')
    print os.path.join(netfolder, 'processFASTQ_part_2.and')

### B. MAIN

if __name__ == "__main__":
    main()
