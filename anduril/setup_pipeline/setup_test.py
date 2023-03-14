#!/bin/env/user python

import jinja2
import yaml
import os, sys
import subprocess

def main():
    
### A. PREREQUISITES    
    # Check for correct usage
    if not len(sys.argv) == 2:
        print 'Usage: ./setup_pipeline.py params.yaml\n'
        sys.exit(0)

    # Get params file and load parameters with yaml
    parameter_file = sys.argv[1]
    print parameter_file
    p = open(parameter_file, 'r')
    params = yaml.load(p)
    print params
    p.close()
    
    print params[OUTPUT]
    
    # load directory with all templates (network_template.and and hosts_template.and) into jinja2:
    temp_env = jinja2.Environment(loader=jinja2.FileSystemLoader('/import/bc2/home/zavolan/cherrmann/PIPELINE'))


### B. RENDER TEMPLATE FILES
 
    ## B.1 NETWORK.and
    # Import network_template.and
    network_temp = temp_env.get_template('network_template.and')
    # Fill template file with parameters from params.yaml file
    rendered_network_template = network_temp.render(params)
    # Write the completed network.and file
    o = open('yaml_1st_test.and', 'w')
    for l in rendered_network_template:
        o.write(l)
    o.close()
    
    ## B.2 HOSTS.conf
    # Import hosts_template.conf
    '''hosts_temp = temp_env.get_template('hosts_template.conf')
    # Render template file with parameters from params.yaml file
    rendered_hosts_template = hosts_temp.render(params)
    # Write the completed hosts.conf file
    o = open('yaml_hosts_1.conf', 'w')
    for l in rendered_hosts_template:
        o.write(l)
    o.close()    
'''
### C. START ANDURIL PIPELINE
    # TO DO: remove shell=True, then command and options has to be a sequence, not a string!    
    '''
    proc = subprocess.call("anduril run yaml_1st_test.and -c components -d import/bc2/home/zavolan/cherrmann/PIPELINE/TEST_OUTPUT --hosts yaml_hosts_1.conf", shell=True)
    
    
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
if __name__ == "__main__":
    main()
