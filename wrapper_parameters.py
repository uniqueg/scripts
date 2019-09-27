#!/bin/env/user python

## Alexander Kanitz, Biozentrum, University of Basel
## 30-JUL-2014

### PREREQUISITES
import yaml
import os, sys

### FUNCTIONS
def main():
    
    ## Check for correct usage
    if not len(sys.argv) == 11:
        print 'Usage: ./wrapper_parameters.py <template_parameters.yaml> <input_file_directory> <characters_to_remove_from_the_end_for_sample_id> <base_execution_folder> <out_dir> <dir_to_setup_pipeline_batch.py_script> <bundle_dir> <template_folder> <filename_workflow_template_without_folder> <filename_host_conf_template_without_folder>\n'
        sys.exit(0)

    # Get files
    infiles = [ os.path.join(sys.argv[2], f) for f in os.listdir(sys.argv[2]) if os.path.isfile(os.path.join(sys.argv[2], f)) ]

    ## Get parameters
    parameter_file = sys.argv[1]
    p = open(sys.argv[1], 'r')
    params = yaml.load(p)
    p.close()

    ## Open pipeline setup wrapper file for writing and add shell bang
    setup_wrapper = open(os.path.join(os.path.dirname(os.path.realpath(__file__)), "wrapper_setup.sh"), 'w')
    setup_wrapper.write("#!/bin/sh\n")
    
    ## Open anduril run wrapper file for writing and add shell bang
    run_wrapper = open(os.path.join(os.path.dirname(os.path.realpath(__file__)), "wrapper_run.sh"), 'w')
    run_wrapper.write("#!/bin/sh\n")

    # Iterate over all files
    for infile in infiles:
        
        # Make a file-specific copy of params
        params_this = params.copy()
        
        # Get sample ID
        sample_this = os.path.basename(infile)[:-int(sys.argv[3])]
        
        # Update SAMPLE_ID
        params_this["SAMPLE_ID"] = sample_this
        
        # Update INPUT_FILE_PATH
        params_this["INPUT_FILE_PATH"] = infile
        
        # Update EXECUTION_DIR
        params_this["EXECUTION_DIR"] = os.path.join(sys.argv[4], sample_this)
        
        # Update TMP_FOLDER
        params_this["TMP_FOLDER"] = os.path.join(sys.argv[4], sample_this, "tmp")
        
        ## Write parameters to file
        outfile_this = os.path.join(sys.argv[5], ("parameters_" + sample_this + ".yaml"))
        fp = open(outfile_this, 'w')
        fp.write(yaml.dump(params_this))
        fp.close()
        
        # Write to setup wrapper script
        setup_wrapper.write("/import/bc2/home/zavolan/krini/bin/python " + sys.argv[6] + " " + outfile_this + " " + sys.argv[5] + " " + sys.argv[8] + " " + sys.argv[9] + " " + sys.argv[10] + "\n")

        # Write to anduril run wrapper script
        run_wrapper.write("/usr/bin/time /import/bc2/home/zavolan/krini/bin/anduril run " + os.path.join(sys.argv[5], ("workflow_" + sample_this + ".and")) + "--exec-mode remote --execution-dir " + params_this["EXECUTION_DIR"] + " --bundle " + sys.argv[7] + " --hosts " + os.path.join(sys.argv[5], ("hosts_" + sample_this + ".conf")) + " --log " + params_this["EXECUTION_DIR"] + "_log\n")
        
    # Close pipeline setup wrapper file
    setup_wrapper.close()
    
    # Close anduril run wrapper file
    run_wrapper.close()

### MAIN
if __name__ == "__main__":
    main()
