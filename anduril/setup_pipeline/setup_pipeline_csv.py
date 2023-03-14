#!/bin/env/user python

import jinja2
import yaml
import os, sys

def extractFormat(infile):
    format = str.upper(str.rsplit(os.path.basename(infile),".")[1])

    return format

def createCSVfile(input_files_list):

    o = open('csvfile', 'w')
    o.write('filepath\tformat\n')

    for f in input_files_list:

        format = extractFormat(f)
        
        o.write('%s\t%s\n' %(f, format))


    o.close()

    return 'csvfile'


def main():

    if not len(sys.argv) == 2:
        print 'Usage: ./setup_pipeline.py params.yaml\n'
        sys.exit(0)

    parameter_file = sys.argv[1]
    print parameter_file
    p = open(parameter_file, 'r')
    params = yaml.load(p)
    print params
    p.close()
    input_dir = params['INPUT_FILES_DIR']
    input_dir_files = os.listdir(input_dir)

    infiles_abspath_list = []
    for f in input_dir_files:
        infiles_abspath_list.append(os.path.join(input_dir, f))

    #input_files_list = params['INPUT_FILES'].split()
    #input_files_list += infiles_abspath_list

    csvfilepath = createCSVfile(infiles_abspath_list)
    params['CSVFILE'] = csvfilepath


    # load directory with all templates (network.and template) into jinja2:
    temp_env = jinja2.Environment(loader=jinja2.FileSystemLoader('/import/bc2/home/zavolan/cherrmann/PIPELINE'))

    # render output html files
    # peaks page
    network_temp = temp_env.get_template('network_template.and')
    rendered_network_template = network_temp.render(params)
    o = open('yaml_1st_test.and', 'w')
    for l in rendered_network_template:
        o.write(l)
    o.close()
    
    hosts_temp = temp_env.get_template('hosts_template.conf')
    rendered_hosts_template = hosts_temp.render(params)
    o = open('yaml_hosts_1.conf', 'w')
    for l in rendered_hosts_template:
        o.write(l)
    o.close()    

if __name__ == "__main__":
    main()
