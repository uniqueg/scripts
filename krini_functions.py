#!/usr/bin/env python

## Alexander Kanitz, Biozentrum, University of Basel (alexander.kanitz@unibas.ch)
## 26-MAR-2015
## Version 1.1 (20-AUG-2015)

##########################
#<--- GLOBAL MODULES --->#
##########################
import sys
import os


####################################
#<--- GENERIC HELPER FUNCTIONS --->#
####################################

#<-- Print header -->#
def printHeader(header, prefix="== ", suffix=" ==", print_top_sep=False, print_bottom_sep=False, top_sep_symbol="=", bottom_sep_symbol="=", top_sep_length=15, bottom_sep_length=15, use_header_length=True):

    # Build header
    header = str(prefix) + str(header) + str(suffix)

    # Print top separator
    if print_top_sep:

        # Determine separator length
        if use_header_length:
            top_sep_length = len(header)

        # Build separator
        top_sep = str(top_sep_symbol) * top_sep_length

        # Print separator
        print top_sep

    # Print header
    print header

    # Print bottom separator
    if print_bottom_sep:

        # Determine separator length
        if use_header_length:
            bottom_sep_length = len(header)

        # Build separator
        bottom_sep = str(bottom_sep_symbol) * bottom_sep_length

        # Print separator
        print bottom_sep

#<-- Print title / value -->#
def printTitleValue(title, value):

    # Print header
    printHeader(title)

    # Print time
    print value

#<-- Print key / value pair -->#
def printKeyValue(key, value):

    # Print time
    print key, value

#<-- Print dictionary -->#
def printDict(dict, header, format="{key:17s}: {value}"):

    # Print header
    printHeader(header)

    # Print dictionary
    for key in sorted(dict):
        string = format.format(key=str(key), value=str(dict[key]))
        print string

#<-- Get runtime in seconds -->#
def getRuntimeSeconds(runtime):

    # Split h:mm:ss string
    runtime_list = runtime.split(':')

    # Calculate seconds
    seconds = int(runtime_list[0]) * 3600 + int(runtime_list[1]) * 60 + int(runtime_list[2])

    # Return seconds
    return(seconds)


#############################################
#<--- SET LINE-BUFFER FOR STDOUT/STDERR --->#
#############################################

#<-- Line buffer STDOUT/STDERR -->#
def lineBufferStdOutErr():

    # Set line buffer by reopening 'stdout'/'stderr' file descriptors in write mode and set 1 as buffer size
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 1)
    sys.stderr = os.fdopen(sys.stderr.fileno(), 'w', 1)


###################################
#<--- GET EXECUTION DIRECTORY --->#
###################################
def getExecPath(tempdir):

    # Get execution path
    execDir = os.path.dirname(tempdir)

    # Die if execution path is not accessible/writable
    if not os.access(execDir, os.W_OK):
        sys.stderr.write("[ERROR] Execution path '%s' is not writable!\n[ERROR] Execution aborted.\n" % execDir)
        sys.exit(1)    

    # Get execution directory from temporary directory
    return(execDir)


###############################################
#<--- VALIDATE INPUT, OUTPUT & PARAMETERS --->#
###############################################

#<-- Wrapper -->#
def validateParameters(component, execDir):

    # Validate input/output/parameters
    component = validateInputPaths(component)
    component = validateOutputPaths(component, execDir)
    component = validateRequiredParams(component)

    # Return 'component' object
    return(component)

#<-- Input -->#
def validateInputPaths(component):

    # Import regex module
    import re

    # Get dictionary of input ports
    inputDict = component.input.to_dict()

    # Compile regular expressions
    reInFile = re.compile(r'^INFILE_')
    reInDir = re.compile(r'^INDIR_')

    # Loop over dictionary of input ports
    for port, destination in inputDict.items():

        # For input ports for which a destination is indicated...
        if destination is not None:

            # Convert to strings
            port = str(port)
            destination = str(destination)

            # Die if destination does not exist
            if not os.path.exists(destination):
                sys.stderr.write('[ERROR] Input destination "{destination}" to port "{port}" not found!\n[ERROR] Execution aborted.\n'.format(destination=destination, port=port))
                sys.exit(1)

            # If input port requests FILE:
            if reInFile.match(port):

                # Die if destination is not a file
                if not os.path.isfile(destination):
                    sys.stderr.write('[ERROR] Input destination "{destination}" to port "{port}" is not a file!\n[ERROR] Execution aborted.\n'.format(destination=destination, port=port))
                    sys.exit(1)

                # Die if destination file is not readable
                if not os.access(destination, os.R_OK):
                    sys.stderr.write('[ERROR] Input destination file "{destination}" to port "{port}" is not readable!\n[ERROR] Execution aborted.\n'.format(destination=destination, port=port))
                    sys.exit(1)

            # If input port requests DIRECTORY:
            elif reInDir.match(port):

                # Die if destination is not a directory
                if not os.path.isdir(destination):
                    sys.stderr.write('[ERROR] Input destination "{destination}" to port "{port}" is not a directory!\n[ERROR] Execution aborted.\n'.format(destination=destination, port=port))
                    sys.exit(1)

                # Die if destination directory is not accessible
                if not os.access(destination, os.X_OK):
                    sys.stderr.write('[ERROR] Input destination directory "{destination}" to port "{port}" is not accessible!\n[ERROR] Execution aborted.\n'.format(destination=destination, port=port))
                    sys.exit(1)

                # Die if destination directory is empty
                if len(os.listdir(destination)) == 0:
                    sys.stderr.write('[ERROR] Input destination directory "{destination}" to port "{port}" is not empty!\n[ERROR] Execution aborted.\n'.format(destination=destination, port=port))
                    sys.exit(1)

            # Die if input port has illegal formatting:
            else:
                sys.stderr.write('[ERROR] Input port "{port}" has illegal format!\n[ERROR] Execution aborted.\n'.format(port=port))
                sys.exit(1)

    # Update and return component
    component.input.update(inputDict)
    return(component)

#<-- Output -->
def validateOutputPaths(component, execDir):

    # Import regex module
    import re

    # Get dictionary of output ports
    outputDict = component.output.to_dict()

    # Compile regular expressions
    reOutFile = re.compile(r'^OUTFILE_')
    reOutDir = re.compile(r'^OUTDIR_')
    reOutDirMake = re.compile(r'^OUTDIRMAKE_')

    # Loop over output list
    for port, destination in outputDict.items():

        # Convert to strings
        port = str(port)
        destination = str(destination)

        # Assert that the destination output file path equals the execution path plus the keyword
        if not destination == os.path.join(execDir, port):
            sys.stderr.write('[ERROR] Path to output port "{port}" not within execution directory "{execDir}"!\n[ERROR] Execution aborted.\n'.format(port=port, execDir=execDir))
            sys.exit(1)

        # If output port requests FILE:
        if reOutFile.match(port):

            # Die if destination file exists
            if os.path.exists(destination):
                sys.stderr.write('[ERROR] Output destination file "{destination}" from port "{port}" already exists!\n[ERROR] Execution aborted.\n'.format(destination=destination, port=port))
                sys.exit(1)

            # Die if destination file is not writable
            try:
                open(destination, 'w').close()
            except OSError, errorCode:
                sys.stderr.write('[ERROR] Output destination file "{destination}" from port "{port}" is not writable!\n[ERROR] Execution aborted.\n'.format(destination=destination, port=port))
                sys.exit(1)

        # If output port requests ABSENT DIRECTORY (created by the executable):
        elif reOutDir.match(port):

            # Die if destination folder exists
            if os.path.exists(destination):
                sys.stderr.write('[ERROR] Output destination directory "{destination}" from port "{port}" already exists!\n[ERROR] Execution aborted.\n'.format(destination=destination, port=port))
                sys.exit(1)

        # If output port requests PRESENT DIRECTORY:
        elif reOutDirMake.match(port):

            # Ensure that destination is interpreted as directory by executing program (append '/')
            destination = outputDict[port] = outputDict[port] + '_TEST/'
            destination = outputDict[port]

            # Try to create directory if it does not exist
            if not os.path.exists(destination):
                try:
                    os.mkdir(destination)
                except OSError, errorCode:
                    sys.stderr.write('[ERROR] Output destination directory "{destination}" from port "{port}" does not exist and cannot be created!\n[ERROR] Execution aborted.\n'.format(destination=destination, port=port))
                    sys.exit(1)

            # Die if destination is not a directory
            if not os.path.isdir(destination):
                sys.stderr.write('[ERROR] Output destination directory "{destination}" to port "{port}" is not a directory!\n[ERROR] Execution aborted.\n'.format(destination=destination, port=port))
                sys.exit(1)

            # Die if destination directory is not accessible
            if not os.access(destination, os.X_OK):
                sys.stderr.write('[ERROR] Output destination directory "{destination}" to port "{port}" is not accessible!\n[ERROR] Execution aborted.\n'.format(destination=destination, port=port))
                sys.exit(1)

        # Die if output port has illegal formatting:
        else:
            sys.stderr.write('[ERROR] Output port "{port}" has illegal format!\n[ERROR] Execution aborted.\n'.format(port=port))
            sys.exit(1)

    # Update and return component
    component.output.update(outputDict)
    return(component)

#<-- Parameters -->#
def validateRequiredParams(component):

    # Import modules
    import subprocess
    import re

    # Get dictionary of parameters
    paramDict = component.param.to_dict()

    # Check whether executable is installed and accessible (using shell function 'which')
    executableList = paramDict['_executable'].split(None, 1)
    process = subprocess.Popen(["which", str(executableList[0])], stdout=subprocess.PIPE)
    location = process.stdout.read().rstrip()
    if not os.path.exists(location):
        sys.stderr.write('[ERROR] Executable "{executable}" not found!\n[ERROR] Execution aborted.\n'.format(executable=executableList[0]))
        sys.exit(1)
    else:
        executableList[0] = location
        paramDict['_executable'] = ' '.join(executableList)

    # Check whether '_execMode' has allowed value
    if paramDict['_execMode'] not in ['remote', 'local', 'none']:
        sys.stderr.write('[ERROR] Illegal execution mode "{_execMode}"!\n[ERROR] Execution aborted.\n'.format(_execMode=paramDict['_execMode']))
        sys.exit(1)

    # Check if '_cores' is integer
    try:
        int(paramDict['_cores'])
    except ValueError:
        sys.stderr.write('[ERROR] Value of parameter "_cores" ({_cores}) is not an integer!\n[ERROR] Execution aborted.\n'.format(_cores=paramDict['_cores']))
        sys.exit(1)

    # Check if '_membycore' is valid memory string
    membycore_format = re.compile(r'^\d+[KMG]?$')
    if not membycore_format.match(str(paramDict['_membycore'])):
        sys.stderr.write('[ERROR] Value of parameter "_membycore" ({_membycore}) is of the wrong format (integer, optionally followed by "K", "M" or "G", expected)!\n[ERROR] Execution aborted.\n'.format(_membycore=paramDict['_membycore']))
        sys.exit(1)

    # Check if '_runtime' is valid time string
    runtime_format = re.compile(r'^\d+:\d\d\:\d\d$')
    if not runtime_format.match(paramDict['_runtime']):
        sys.stderr.write('[ERROR] Value of parameter "_runtime" ({{_runtime}}) is of the wrong format (hh:mm:ss expected)!\n[ERROR] Execution aborted.\n'.format(_runtime=paramDict['_runtime']))
        sys.exit(1)

    # Update and return component
    component.param.update(paramDict)
    return(component)


############################
#<--- VALIDATE COMMAND --->#
############################

#<-- Remove command expression markup and handle optional/multi/switch parameters & input/output ports
def renderCommand (component, command, tempdir, execdir, separator='$$$', spacerNonPositional='^^^', spacerPositional='###', supportedRedirectors=['>', '2>', '&>']):

    # Import regular expression module
    import re

    # Replace reserved metavalues/placeholders
    command = replacePlaceholders(component, command, tempdir, execdir)

    # Get dictionaries of output ports and parameters
    inputDict = component.input.to_dict()
    outputDict = component.output.to_dict()
    paramDict = component.param.to_dict()

    # Compile regular expressions
    regexExpr = re.compile(r'([^\s]+?)' + r'(' + re.escape(spacerNonPositional) + r'|' + re.escape(spacerPositional) + r')' + r'(.+)')
    regexArgMeta = re.compile(r'\{\{' + r'(.*)' + r'\}\}')
    regexInt = re.compile(r'^\d+$')
    regexMultiList = re.compile(r'\[\[' + r'(.*)' + r'\]\]')  
    regexEmptyOrSpace = re.compile(r'^[^\s]*$')  

    # Split command by separator
    commandList = command.split(separator)

    # Initialize 'stdout' and 'stderr'
    stdout = stderr = None

    # Initialize redirector flag (after the first redirector, only further redirectors are allowed)
    redirectorEncountered = False

    # Iterate over command list/index
    for index, expression in enumerate(commandList):

        # Initialize replacement value
        replace = ''

        # If expression matches expected format
        matchExpr = regexExpr.match(expression)
        if matchExpr:

            # Extract option, opt/arg spacer & argument
            option = str(matchExpr.group(1))
            spacer = str(matchExpr.group(2))
            argument = str(matchExpr.group(3))

            # Die if illegal spacer is encountered
            if spacer not in [spacerNonPositional, spacerPositional]:
                sys.stderr.write('[ERROR] The spacer "{spacer}" in command expression "{expression}" is not allowed!\n[ERROR] Execution aborted.\n'.format(spacer=spacer, expression=expression))
                sys.exit(1)

            # If option is redirector
            if option in supportedRedirectors:
                redirectorEncountered = True
                if '_' + argument not in paramDict.keys() or paramDict['_' + argument]:
                    validateRedirectionTarget(component, outputDict[argument])
                    stdout, stderr = setStdOutErr(stdout, stderr, option, outputDict[argument])
                replace = ''

            # Die if a redirector expression has been processed before
            elif redirectorEncountered:
                sys.stderr.write('[ERROR] Non-redirector expression "{expression}" encountered after redirector in command "{command}"!\n[ERROR] Execution aborted.\n'.format(expression=expression, command=command))
                sys.exit(1)

            # Check for inputs, other outputs and parameters
            else:

                # If argument is a metavalue of the form {{VALUE}}
                matchArgMeta = regexArgMeta.match(argument)
                if matchArgMeta:
                    argContent = str(matchArgMeta.group(1))

                    # Metavalue: {{FALSE}} / unset switch
                    if argContent in ['FALSE', '']:
                        replace = ''

                    # Metavalue: {{TRUE}} / set switch
                    elif argContent == 'TRUE':
                        replace = option

                    # Metavalue: {{<int>}} / repeat switch
                    elif regexInt.match(argContent):
                        replace = ' '.join([option] * int(argContent))

                    # Metavalue: {{[[item1//item2//.../itemN]]}} / repeat option with arguments
                    elif regexMultiList.match(argContent):
                        argsList = argContent[2:-2].split('//')
                        replaceList = []
                        for arg in argsList:
                            if regexEmptyOrSpace.match(arg):
                                replaceList.append(option + ' ' + arg)
                        if replaceList:
                            replace = ' '.join(replaceList)
                        else:
                            replace = ''

                    # Die if unknown metavalue
                    else:
                        sys.stderr.write('[ERROR] The metavalue "{argContent}" in expression "{expression}" was not recognized!\n[ERROR] Execution aborted.\n'.format(argContent=argContent, expression=expression))
                        sys.exit(1)

                # If argument is input port
                elif argument in inputDict.keys():

                    if inputDict[argument] is None:
                        replace = ''
                    elif spacer == spacerNonPositional:
                        replace = option + ' ' + inputDict[argument]
                    elif spacer == spacerPositional:
                        replace = inputDict[argument]

                # If argument is output port
                elif argument in outputDict.keys():

                    if '_' + argument in paramDict.keys() and not paramDict['_' + argument]:
                        replace = ''
                    elif spacer == spacerPositional:
                        replace = outputDict[argument]
                    elif spacer == spacerNonPositional:
                        replace = option + ' ' + outputDict[argument]
                    else:
                        sys.stderr.write('[ERROR] Unexpected output port expression "{expression}"!\n[ERROR] Execution aborted.\n'.format(expression=expression))
                        sys.exit(1)

                # If option is non-positional
                elif spacer == spacerNonPositional:
                    replace = option + ' ' + argument

                # If option is positional
                elif spacer == spacerPositional:
                    replace = argument

                # Catch exception
                else:
                    sys.stderr.write('[ERROR] Unexpected parameter expression "{expression}"!\n[ERROR] Execution aborted.\n'.format(expression=expression))
                    sys.exit(1)

            # Replace expression
            commandList[index] = replace

        # Ignore first expression (command executable)
        elif index == 0:
            pass

        # Die if expression is not of expected format
        else:
            sys.stderr.write('[ERROR] The expression "{expression}" is not of the expected format!\n[ERROR] Execution aborted.\n'.format(expression=expression))
            sys.exit(1)

    # Remove empty strings from processed command list
    commandList = [item for item in commandList if item]

    # Ensure executable is addressed by full absolute path
    commandList[0] = paramDict['_executable']

    # Join command list items by spaces
    command = ' '.join(commandList)

    # Return command string
    return(command, stdout, stderr)
    pass

#<-- Substitute placeholder variables with corresponding values -->#
def replacePlaceholders(component, command, tempdir, execdir):

    # Import regex module
    import re

    # Build regulard expression for finding placeholder variables (markup: {{VARIABLE}})
    placeholder = re.compile(r'\{\{' + r'([A-Z_]+)' + r'\}\}')

    # Find all placeholder variables
    for placeholder in placeholder.finditer(command):

        # Initialize replacement variable
        replace = False

        # Get placeholder variable without markup
        var = placeholder.group(1)

        # Get replacement values
        if   var == 'CORES':
            replace = str(component.param._cores)

        elif var == 'MEMBYCORE':
            replace = str(component.param._membycore)

        elif var == 'TEMPDIR':
            replace = str(tempdir)

        elif var == 'EXECDIR':
            replace = str(execdir)

        # Replace placeholder with value
        if replace:
            command = command.replace(placeholder.group(), replace)

    # Return modified command
    return(command)

#<-- Validate whether redirection target is an output port -->#
def validateRedirectionTarget(component, target):

    # Get output dictionary
    output = component.output.to_dict()

    # Die if target is not an output port
    if target not in output.values():
        sys.stderr.write('[ERROR] Redirection target "{target}" is not an output port!\n[ERROR] Execution aborted.\n'.format(target=target))
        sys.exit(1)

#<-- Parse supported output redirectors -->#
def setStdOutErr(stdout, stderr, type, target):

    # Initialize 'stdout' and 'stderr'
    stdout = stderr = None

    # Set target according to type
    if   type ==  ">" and stdout is None:
        stdout = target

    elif type == "2>" and stderr is None:
        stderr = target

    elif type == "&>" and stdout is None and stderr is None:
        stdout = stderr = target

    else:
        sys.stderr.write('[ERROR] The command contains ambiguous redirectors!\n[ERROR] Execution aborted.\n'.format())
        sys.exit(1)

    # Return 'command', 'stdout' & 'stderr'
    return(stdout, stderr)


######################################################################
#<--- PRINT INPUT/OUTPUT FILES/DIRECTORIES, PARAMETERS & COMMAND --->#
######################################################################

#<-- Wrapper -->#
def printParameters(component, command, execDir, tempdir):

    # Get dictionaries from Anduril's 'component' Struct
    meta = component.meta.to_dict()
    input = component.input.to_dict()
    output = component.output.to_dict()
    param = component.param.to_dict()

    # Print input/output files/directories and parameters
    printMeta(meta, execDir, tempdir)
    printDict(input, "Input ports", format="{key:31s}: {value}")
    printDict(output, "Output ports", format="{key:31s}: {value}")
    printDict(param, "Parameters", format="{key:31s}: {value}")

    # Print command
    printTitleValue("Command", command)

#<-- Print select metadata -->#
def printMeta(metaDict, execDir, tempdir):

    # Print header
    printHeader("Component/instance information")

    # Print values
    printKeyValue("{0:31s}:".format("Component name"), str(metaDict['componentName']))
    printKeyValue("{0:31s}:".format("Instance name"), str(metaDict['instanceName']))
    printKeyValue("{0:31s}:".format("Component path"), str(metaDict['componentPath']))
    printKeyValue("{0:31s}:".format("Execution directory"), str(execDir))
    printKeyValue("{0:31s}:".format("Temporary directory"), str(tempdir))


###########################
#<--- EXECUTE COMMAND --->#
###########################

#<-- Wrapper -->#
def executeCommand(component, execDir, command, stdout, stderr):

    # Import shlex module
    import shlex

    # Split command into arguments
    command = shlex.split(command)

    # Get STDOUT/STDERR filenames
    stdout_filename, stderr_filename = getStdOutErrFiles(stdout, stderr, component, execDir)

    # Execute command remotely, locally or not at all
    if   component.param._execMode == "remote":
        exit_status = executeRemotelyDRMAAtoSGE(component, command, stdout_filename, stderr_filename)
    elif component.param._execMode == "local":
        exit_status = executeLocally(command, stdout_filename, stderr_filename)
    elif component.param._execMode == "none":
        exit_status = -1
    else:
        sys.stderr.write('[ERROR] Illegal execution mode {mode}!\n[ERROR] Execution aborted.\n'.format(mode=component.param._execMode))
        sys.exit(1)

    # Print command output
    printStdOutErr(stdout, stdout_filename, stderr, stderr_filename)

    # Print exit status
    printTitleValue("Exit status", exit_status)

    # Return exit status
    return(exit_status)

#<-- Get STDOUT/STDERR filenames -->#
def getStdOutErrFiles(stdout, stderr, component, execDir):

    # Get basic file path from execution directory and instance name
    base_path = os.path.join(execDir, component.meta.instanceName)

    # Set 'stdout'
    if stdout is None:
        stdout = base_path + ".stdout"

    # Set 'stderr'
    if stderr is None:
        stderr = base_path + ".stderr"

    # Return 'stdout' and 'stderr' filehandles
    return(stdout, stderr)

#<-- Local execution -->#
def executeLocally(command, stdout, stderr):

    # Import modules
    import subprocess
    import datetime

    # Open 'stdout' and 'stderr' filehandles
    stdout_handle = open(stdout, 'w')
    stderr_handle = open(stderr, 'w')

    # Print progress section
    printHeader("Progress")

    # Get & print submit/start time
    time_submit = time_start = datetime.datetime.now()
    printKeyValue("{0:31s}:".format("Started"), "{:%Y-%b-%d, %H:%M:%S}".format(time_start))

    # Execute command
    exit_status = subprocess.call(command, stdout=stdout_handle, stderr=stderr_handle)

    # Get & print end time
    time_end = datetime.datetime.now()
    printKeyValue("{0:31s}:".format("Finished"), "{:%Y-%b-%d, %H:%M:%S}".format(time_end))

    # Print time stats
    printTimeStats(submit=time_submit, start=time_start, end=time_end)

    # Close 'stdout' and 'stderr' filehandles
    stdout_handle.close()
    stderr_handle.close()

    # Return exit status
    return(exit_status)

#<-- DRMAA/SGE execution -->#
def executeRemotelyDRMAAtoSGE(component, command, stdout, stderr):

    # Import modules
    import drmaa
    import datetime

    # Start/initialize DRMAA session
    session = drmaa.Session()
    session.initialize()

    # Create job template
    job_template = session.createJobTemplate()

    # Render job template
    job_template.jobName = component.meta.instanceName
    job_template.outputPath = ":" + stdout     # ':' required by DRMAA
    job_template.errorPath = ":" + stderr      # ':' required by DRMAA
    job_template.remoteCommand = command[0]
    job_template.args = command[1:]
    job_template.nativeSpecification = '-v LD_LIBRARY_PATH="{}"'.format(os.environ["LD_LIBRARY_PATH"])
    job_template.nativeSpecification += ' -pe smp {cpus} -l membycore={mem} -l runtime={rt} -b yes -shell no -p 0 -w e'.format(cpus=component.param._cores, mem=component.param._membycore, rt=component.param._runtime)

    # Submit job
    try:
        job_id = session.runJob(job_template)
    except drmaa.errors.DeniedByDrmException, errorCode:
        sys.stderr.write('[ERROR] Job was denied by DRM with the following error code:\n[ERROR] {errorCode}\n[ERROR] Execution aborted.\n'.format(errorCode=errorCode))
        sys.exit(1)

    # Print job ID
    printTitleValue("Job ID", job_id)

    # Print progress section header
    printHeader("Progress")

    # Get & print job submit time
    time_submit = datetime.datetime.now()
    printKeyValue("{0:31s}:".format("Submitted"), "{:%Y-%b-%d, %H:%M:%S}".format(time_submit))

    # Check job status until job runs; die in case of problems
    checkDRMAAJobStatus(session, job_id)

    # Get & print job start time
    time_start = datetime.datetime.now()
    printKeyValue("{0:31s}:".format("Started"), "{:%Y-%b-%d, %H:%M:%S}".format(time_start))

    # Wait for job to end
    job_info = waitForDRMAAJobToFinish(session, job_id, component)

    # Get & print job end time
    time_end = datetime.datetime.now()
    printKeyValue("{0:31s}:".format("Finished"), "{:%Y-%b-%d, %H:%M:%S}".format(time_end))

    # Delete job template
    session.deleteJobTemplate(job_template)

    # Exit DRMAA session
    session.exit()

    # Print time stats
    printTimeStats(submit=time_submit, start=time_start, end=time_end)

    # Print resource stats
    printDRMAAResourceStats(job_info.resourceUsage)

    # Return exit status
    return(job_info.exitStatus)

#<-- DRMAA: Check status of jobs that are not running -->#
def checkDRMAAJobStatus(session, job_id, timeout_active=0.5, tries_inactive=60, timeout_inactive=60):
### Continuously checks the status of a DRMAA submitted job based on the job ID and the DRMAA session.
### Function waits (timeout_active seconds) as long as job is in active queue and returns successfully (no return value) if job execution commences.
### Function dies with error if (a) job fails or (b) stays in an inactive state continuously (tries_inactive * timeout_inactive seconds). 

    # Import modules
    import drmaa, time

    # Set tries
    tries = tries_inactive

    # While there are tries left...
    while tries > 0:

        # Get job status
        status = session.jobStatus(job_id)

        # If job...
        # ...is in active queue, reset tries, sleep and try again
        if   status == 'queued_active':
            tries = tries_inactive
            time.sleep(timeout_active)
        # ...is running, break out of loop
        elif status == 'running':
            break
        # ...failed, die with error
        elif status == 'failed':
            sys.stderr.write("[ERROR] Job '%s' failed!\n[ERROR] Execution aborted.\n" % job_id)
            sys.exit(1)
        # Else sleep and try again for a specified amount of times, then die
        else:
            if tries == 1:
                sys.stderr.write("[ERROR] Job '%s' is in a prolonged inactive state.\n[ERROR] Execution aborted.\n" % job_id)
                session.control(job_id, drmaa.JobControlAction.TERMINATE)
                sys.exit(1)
            else:
		tries -= 1
                print "Job '{0}' not in active queue (job status: '{1}'). Trying again in {2} seconds ({3} tries left).".format(job_id, status, timeout_inactive, tries)
                time.sleep(timeout_inactive)

#<--- DRMAA: Wait for running jobs to finish --->#
def waitForDRMAAJobToFinish(session, job_id, component):

    # Import DRMAA module
    import drmaa

    # Wait for job to finish...
    try:
        job_info = session.wait(job_id, timeout=int(getRuntimeSeconds(component.param._runtime)))

    # ...but die if it times out
    except drmaa.errors.ExitTimeoutException:
        sys.stderr.write("[ERROR] Job '%s' exceeded specified runtime (%s seconds).\n[ERROR] Execution aborted.\n" % (job_id, component.param._runtime))
        sys.exit(1)

    # Die if job was aborted
    if job_info.wasAborted:
        sys.stderr.write("[ERROR] Job '%s' was aborted.\n[ERROR] Execution aborted.\n" % job_id)
        sys.exit(1)

    # Die if job was killed
    if job_info.hasSignal:
        sys.stderr.write("[ERROR] Job '%s' was killed (signal: %s).\n[ERROR] Execution aborted.\n" % (job_id, job_info.terminatedSignal))
        sys.exit(1)

    # Die if job ended prematurely for any other reason
    if not job_info.hasExited:
        sys.stderr.write("[ERROR] Job '%s' exited prematurely.\n[ERROR] Execution aborted.\n" % job_id)
        sys.exit(1)

    # Return job info
    return(job_info)

#<-- Print time statistics -->#
def printTimeStats(submit, start, end):

    # Calculate durations
    queue_time = start - submit
    run_time = end - start
    total_time = end - submit

    # Print header
    printHeader("Time statistics")

    # Print durations
    printKeyValue("{0:31s}:".format("Queue time"), "{0:.3f} s".format(float(queue_time.total_seconds()), "s"))
    printKeyValue("{0:31s}:".format("Runtime"), "{0:.3f} s".format(float(run_time.total_seconds()), "s"))
    printKeyValue("{0:31s}:".format("Total time"), "{0:.3f} s".format(float(total_time.total_seconds()), "s"))

#<-- Print resource statistics -->#
def printDRMAAResourceStats(dict):

   # Print header
   printHeader("Resource statistics")

   # Print resources of interest
   printKeyValue("{0:31s}:".format("Real time"), "{0:.3f} s".format(float(dict['ru_wallclock'])))
   printKeyValue("{0:31s}:".format("User time"), "{0:.3f} s".format(float(dict['ru_utime'])))
   printKeyValue("{0:31s}:".format("System time"), "{0:.3f} s".format(float(dict['ru_stime'])))
   printKeyValue("{0:31s}:".format("CPU time"), "{0:.3f} s".format(float(dict['cpu'])))
   printKeyValue("{0:31s}:".format("Integral memory"), "{0:.3f} Gbs".format(float(dict['mem'])))
   printKeyValue("{0:31s}:".format("Maximum virtual memory"), "{0:.0f} bytes".format(float(dict['maxvmem'])))
   printKeyValue("{0:31s}:".format("Maximum resident set size"), "{0:.3f} kb".format(float(dict['ru_maxrss'])))
   printKeyValue("{0:31s}:".format("Accumulated I/O usage"), "{0:.3f} Gb".format(float(dict['io'])))
   printKeyValue("{0:31s}:".format("I/O wait time"), "{0:.3f} s".format(float(dict['iow'])))
   printKeyValue("{0:31s}:".format("Block input operations"), "{0:.0f}".format(float(float(dict['ru_inblock']))))
   printKeyValue("{0:31s}:".format("Block output operations"), "{0:.0f}".format(float(dict['ru_oublock'])))
   printKeyValue("{0:31s}:".format("Soft page faults"), "{0:.0f}".format(float(dict['ru_minflt'])))
   printKeyValue("{0:31s}:".format("Hard page faults"), "{0:.0f}".format(float(dict['ru_majflt'])))
   printKeyValue("{0:31s}:".format("Voluntary context switches"), "{0:.0f}".format(float(dict['ru_nvcsw'])))
   printKeyValue("{0:31s}:".format("Involuntary context switches"), "{0:.0f}".format(float(dict['ru_nivcsw'])))

#<-- Print STDOUT/STDERR -->#
def printStdOutErr(stdout, stdout_filename, stderr, stderr_filename):

    # STDOUT: Print header
    printHeader("Command STDOUT")

    # If STDOUT was redirected to output file write only link to STDOUT file
    if stdout is not None:
        sys.stdout.write("< [STDOUT] was redirected to output file '%s'. >\n" % stdout_filename)

    # Else write entire content of STDOUT file
    else:
        if os.path.exists(stdout_filename):
            if os.path.getsize(stdout_filename) > 0:
                sys.stdout.write('< [STDOUT] as saved in file "{file}": >\n'.format(file=stdout_filename))
                stdout_handle = open(stdout_filename, 'r')
                sys.stdout.write(stdout_handle.read())
            else:
                sys.stdout.write('< [STDOUT] as saved in file "{file}" is empty. >\n'.format(file=stdout_filename))
        else:
            sys.stdout.write('< [STDOUT] file "{file}" was not produced. >\n'.format(file=stdout_filename))

    # STDERR: Print header
    printHeader("Command STDERR")

    # If STDERR was redirected to output file write only link to STDERR file
    if stderr is not None:
        sys.stdout.write("< [STDERR] was redirected to output file '%s'. >\n" % stderr_filename)

    # Else write entire content of STDERR file
    else:
        if os.path.exists(stderr_filename):
            if os.path.getsize(stderr_filename) > 0:
                sys.stdout.write('< [STDERR] as saved in file "{file}": >\n'.format(file=stderr_filename))
                stderr_handle = open(stderr_filename, 'r')
                sys.stdout.write(stderr_handle.read())
            else:
                sys.stdout.write('< [STDERR] as saved in file "{file}" is empty. >\n'.format(file=stderr_filename))
        else:
            sys.stdout.write('< [STDERR] file "{file}" was not produced. >\n'.format(file=stderr_filename))


#######################################
#<--- CREATE MISSING OUTPUT FILES --->#
#######################################
def createMissingOutputFiles(component):

    # Import regex module
    import re

    # Get directory of output ports
    outputDict = component.output.to_dict()

    # Compile regular expressions
    reOutFile = re.compile(r'^OUTFILE_')
    reOutDir = re.compile(r'^OUTDIR_')
    reOutDirMake = re.compile(r'^OUTDIRMAKE_')

    # Loop over output list
    for port, destination in outputDict.items():

        # Convert to strings
        port = str(port)
        destination = str(destination)

        # If destination file does not exist
        if not os.path.exists(destination):

            # If output port requests FILE:
            if reOutFile.match(port):

                # Try to create (open/close) file
                try:
                    open(destination, 'w').close()
                except OSError, errorCode:
                    sys.stderr.write('[ERROR] Trying to create output destination file "{destination}" from port "{port}" raised an OSError!\n[ERROR] {error}\n[ERROR] Execution aborted.\n'.format(destination=destination, port=port, errorCode=errorCode))
                    sys.exit(1)

            # If output port requests ABSENT DIRECTORY (created by the executable):
            elif reOutDir.match(port):

                # Try to create directory
                if not os.path.exists(destination):
                    try:
                        os.mkdir(destination)
                    except OSError, errorCode:
                        sys.stderr.write('[ERROR] Trying to create output destination directory "{destination}" from port "{port}" raised an OSError!\n[ERROR] {error}\n[ERROR] Execution aborted.\n'.format(destination=destination, port=port, errorCode=errorCode))
                        sys.exit(1)

            # Die if output port has illegal formatting:
            else:
                sys.stderr.write('[ERROR] Output port "{port}" has illegal format!\n[ERROR] Execution aborted.\n'.format(port=port))
                sys.exit(1)


################
#<--- MAIN --->#
################
def main(component, command, tempdir):

    # Line buffering STDOUT/STDERR
    lineBufferStdOutErr()

    # Get execution path
    execDir = getExecPath(tempdir)

    # Validate/modify input and output directories/files
    component = validateParameters(component, execDir)

    # Validate and return sanitized command and command STDOUT/STDERR ports
    command, stdout, stderr = renderCommand(component, command, tempdir, execDir)

    # Print job metadata, input/output files, parameters
    printParameters(component, command, execDir, tempdir)

    # Execute command
    exit_status = executeCommand(component, execDir, command, stdout, stderr)

    # Create missing output files
    createMissingOutputFiles(component)

    # Return exit status
    return(exit_status)
