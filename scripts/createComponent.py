#!/usr/bin/env python

##################
#<--- HEADER --->#
##################
## Alexander Kanitz, Biozentrum, University of Basel (alexander.kanitz@unibas.ch)
## 08-APR-2015
## Version 1.0 (08-APR-2015)

## TODO: Write documentation
## TODO: Input / output folders & prefixes (how can Anduril find the right output for the next step?!); test folder behavior (in worst case: parse component.xml, datatypes.xml, component.json...)
## TODO: Include short parameter synopsis for homepage (for all DOCSTRINGS)
## TODO: Include manual URL as tag or attribute
## TODO: Unicode support
## TODO: After execution, report which output files/folders are empty/non-empty
## TODO: Re-factor ugly parts of script
## TODO: Change section dictionary structure from [{tagName: <CLASS: I/O/P>, tagAttributes: {...}, tagValue: <DOCSTRING>, ...] to ???
             

#################################
#<--- IMPORT GLOBAL MODULES --->#
#################################
import sys, os


################
#<--- MAIN --->#
################
def main():

    # Import modules
    import xlrd
    import re

    # Assign CLI arguments
    xlsx = sys.argv[1]
    rootdir = sys.argv[2]

    # Open XLSX workbook file
    workbook = xlrd.open_workbook(xlsx)

    # Iterate over worksheets (one sheet = one component)
    for worksheet in workbook.sheet_names():

        # Ignore worksheets names 'Template' and 'Example'
        if worksheet not in ['TEMPLATE', 'SANDBOX', 'EXAMPLE'] and not re.match(r'^BUILD', worksheet):

            # Read worksheet and convert to dictionary of sections
            print '==\nProcessing worksheet "{}"...'.format(worksheet)
            lines = readWorksheet(workbook, worksheet)
            sectionDict = fileToDict(lines)

            # Generate boolean parameters for optional outputs
            sectionDict = addOptionalOutputBools(sectionDict)

            # Extract component name
            componentName = getComponentName(sectionDict)
            print '--\nComponent name: {}'.format(componentName)

            # Build command
            command = buildCommand(sectionDict)
            print '--\nThe following command template was built:\n{}'.format(command)

            # Create output directory
            outdir = os.path.join(rootdir, componentName)
            createOutDir(outdir)
            print '--\nOutput directory: "{}"'.format(outdir)

            # Enforce docstring format
            sectionDict = formatDocStrings(sectionDict)

            # Write JSON output file
            outJSON = os.path.join(outdir, 'component.json')
            writeJSON(sectionDict, outJSON)
            print '--\nComplete parameter file in JSON format written to "{}".'.format(outJSON)

            # Write XML output file
            outXML = os.path.join(outdir, 'component.xml')
            writeXML(sectionDict, outXML)
            print 'Anduril-compatible parameter file in XML format written to "{}".'.format(outXML)

            # Write Python wrapper
            outPy = os.path.join(outdir, 'component.py')
            writePyWrapper(command, outPy)
            print 'Generic wrapper script for Anduril-based execution of the component written to "{}".'.format(outXML)

            # Print status message
            print '--\nFinished processing worksheet "{}".'.format(worksheet)

    # Print status message and exit
    print '==\nDone.'
    sys.exit(0)
    pass


###################################
#<--- WORKSHEET TO DICTIONARY --->#
###################################
def readWorksheet(workbook, worksheet_name):
# Given a workbook object and a worksheet name, reads all rows of the worksheet, returning the content in a list of lists (cells per row)
# Requires the 'xlrd' module to handle Excel-style XLSX workbooks

    # Import XLRD module
    import xlrd

    # Open worksheet or die
    try:
        worksheet = workbook.sheet_by_name(worksheet_name)
    except xlrd.biffh.XLRDError, error_code:
        sys.stderr.write('[ERROR] {code}\n[ERROR] Execution aborted.\n'.format(code = error_code))
        sys.exit(1)

    # Initiate outer list (1 item = 1 row), iterate over rows and add them to outer list
    rows_list = []
    for row in range(0, worksheet.nrows):

        # Initiate inner list (1 item = 1 cell), iterate over cells and add them to inner list
        row_list = []
        for cell in worksheet.row(row):
            row_list.append(str(cell.value))

        rows_list.append(row_list)

    # Return list of lists
    return(rows_list)
    pass

def fileToDict(lines):
# Given a list of lists (outer list: rows, inner lists: fields), returns a dictionary of the form {section1: [entry1, ..., entryN], ..., section 3: [entry1, ..., entryN] }. List items in turn XML tags of the form {tagName: ..., tagAttributes: {...}, tavValue: ...}. 
# Sections, section headers and bodies as well as comment lines are recognized by specific markup values in the first fields of a given line.

    # Initialize dictionary of sections, current section name, (sub)section header, list of entries of current section
    sectionDict = {}
    sectionName = None
    sectionHeader = []
    sectionList = []

    # Iterate over lines
    for line in lines:

        # Complete current section & start new section
        if line[0] == '##':
            sectionDict = finishSection(sectionName, sectionList, sectionDict)
            sectionName = setSection(line)
            sectionHeader = []
            sectionList = []

        # Get (sub)section header
        elif line[0] == '#':
            sectionHeader = line[2:]

        # Ignore comment and empty lines
        elif line[0] == '!' or not any(line):
            pass

        # Process section body line and add to section list
        elif sectionName is not None:
            lineDict = processLine(sectionHeader, line, sectionName)
            sectionList.append(lineDict)

        # Die if section body before section declaration
        else:
            sys.stderr.write('[ERROR] Section body line "{}" encountered, but no section set!\nExecution aborted.\n'.format(line))
            sys.exit(1)

    # Finish last section
    sectionDict = finishSection(sectionName, sectionList, sectionDict)

    # Return dictionary of sections
    return(sectionDict)
    pass


def setSection(line):
# Extracts and returns section name from section initiation line (list of field values)
    try:
        section = line[1]
    except IndexError:
        sys.stderr.write('[ERROR] No section name provided in line "{}"!\nExecution aborted.\n'.format(line))
        sys.exit(1)        
    if not section or section is None:
        sys.stderr.write('[ERROR] No valid section name provided in line "{}"!\nExecution aborted.\n'.format(line))
        sys.exit(1)
    return(section)
    pass

def processLine(ids, line, section):
# Given a section body line (list of field values), generates and returns XML tag structure dictionary of the form {tagName: ..., tagAttributes: {...}, tavValue: ...}
# Adds a valid option name for each entry of sections 'inputs', 'outputs' & 'parameters'
# Reserved XML characters in XML tag values are escaped
    lineDict = {}
    ids = list(ids)
    values = line[2:]
    if 'optionName' in ids:
        keys, values = addValidOptionName(ids, values, section)
    if not len(ids) == len(values):
        sys.stderr.write('[ERROR] Number of IDs ({}) and values ({}) differ!\nExecution aborted.\n'.format(ids, values))
        sys.exit(1)
    lineDict['tagName'] = line[0]
    lineDict['tagAttributes'] = {key:value for key,value in zip(ids,values) if value}
    lineDict['tagValue'] = escapeReservedXML(line[1])
    return(lineDict)
    pass

def addValidOptionName(keys, values, section):
# Given two lists of keys of values, as well as a section name, the program-specific option names (usually of the form --option) are edited for use in Anduril
# Option-specific markup is added to 

    # Import regular expression module
    import re

    # Get index of 'optionName' bareword in list 'keys' and get corresponding value (die if it does not exist)
    index = keys.index('optionName')
    try:
        string = values[index]
    except IndexError:
        sys.stderr.write('[ERROR] No "optionName" provided in line "{}"!\nExecution aborted.\n'.format(line))
        sys.exit(1)

    # Generate Anduril-compatible version of option name
    string = string.lstrip('-')
    string = string.replace('-', '_')
    string = string.replace('@', 'at')
    string = string.replace('#', 'sharp')
    string = string.replace('?', 'quMark')
    string = string.replace('!', 'exclMark')
    string = string.translate(None, '"$%^&*()+={}[]\/<>,.`~|')
    string = string.translate(None, "'")
    if re.match(r'^\d', string):
        string = '_' + string
    if string == 'parameter':
        string = 'parameter_'

    # Add option-specific markup
    if 'fileClass' in keys:
        classIndex = keys.index('fileClass')
        try:
            className = values[classIndex]
        except IndexError:
            sys.stderr.write('[ERROR] No input/output file class (file, directory) provided in line "{}"!\nExecution aborted.\n'.format(line))
            sys.exit(1)
        if section == 'inputs':
            if className == 'directory':
                string = 'INDIR_' + string
            else:
                string = 'INFILE_' + string
        if section == 'outputs':
            if className == 'directory':
                string = 'OUTDIR_' + string
            elif className == 'directoryMake':
                string = 'OUTDIRMAKE_' + string
            else:
                string = 'OUTFILE_' + string

    # Add Anduril-compatible option name and bareword 'name' to lists values and keys, respectively
    keys.append('name')
    values.append(string)

    # Return modified lists
    return(keys, values)
    pass

def escapeReservedXML(string):
# Escapes the XML reserved characters (&,<,>,',") in a string and returns it
    string = string.replace('&', '&amp;')
    string = string.replace('<', '&lt;')
    string = string.replace('>', '&gt;')
    string = string.replace("'", '&apos;')
    string = string.replace('"', '&quot;')
    return(string)

def finishSection(sectionName, sectionList, sectionDict):
# Adds a list of section items to the section dictionary
# Returns an empty dictionary if no section name was previously defined (first instance)
    if sectionName is None:
        return({})
    sectionDict[sectionName] = sectionList
    return(sectionDict)
    pass


############################################################
#<--- ADD BOOLEAN PARAMETERS FOR OPTIONAL OUTPUT PORTS --->#
############################################################
def addOptionalOutputBools(dictionary):
# As Anduril does not support optional outputs, boolean switches are generated for every optional output port and added to the Anduril XML component descriptor file as parameters
    dictionary['optOutBoolParameters'] = []
    for element in dictionary['outputs']:
        if element['tagAttributes']['optional'] == 'true':
            tmpDict = {}
            tmpDict['name'] = '_' + element['tagAttributes']['name']
            tmpDict['type'] = 'boolean'
            tmpDict['default'] = 'false'
            elementDict = {}
            elementDict['tagName'] = 'parameter'
            elementDict['tagAttributes'] = tmpDict
            elementDict['tagValue'] = 'Switch for enabling the optional output port "{}". Set to true in workflow file if this output is desired.'.format(element['tagAttributes']['optionName'])
            dictionary['optOutBoolParameters'].append(elementDict)
    return(dictionary)
    pass


##############################
#<--- GET COMPONENT NAME --->#
##############################
def getComponentName(sectionDict):
# Given the section dictionary, returns the component name from the value field of section 'header'
    componentName = sectionDict['header'][0]['tagValue']
    return(componentName)
    pass


#########################
#<--- BUILD COMMAND --->#
#########################
def buildCommand(dictionary):

    # Die if parameters have identical names
    checkForDuplicates(dictionary)

    # Get option name dictionaries
    internalDict = getParamNameDict(dictionary['internalParameters'], all=True)
    iDict = getParamNameDict(dictionary['inputs'], all=False)
    oDict = getParamNameDict(dictionary['outputs'], all=False)
    pDict = getParamNameDict(dictionary['parameters'], all=False)
    posDict = getParamPositionDict(dictionary)

    # Add executables, parameters, inputs, outputs, positional parameters/inputs/outputs & redirected outputs
    prefix, suffix = addExecutable(internalDict)
    prefix, suffix = addParameters(pDict, prefix, suffix, spacer='^^^')
    prefix, suffix = addParameters(iDict, prefix, suffix, spacer='^^^')
    prefix, suffix = addParameters(oDict, prefix, suffix, spacer='^^^')
    prefix, suffix = addParameters(posDict, prefix, suffix, spacer='###')

    # Merge prefix (command with placeholders) and suffix (replacement dictionary) into valid python <STRING>.format(<DICT>) expression
    command = mergePrefixSuffix(prefix, suffix)

    # Return command
    return(command)
    pass

def checkForDuplicates(dictionary):
    paramsList = dictionary['inputs'] + dictionary['outputs'] + dictionary['parameters'] + dictionary['internalParameters']
    paramNameDict = getParamNameDict(paramsList, all=True)
    if not len(paramNameDict) == len(paramsList):
        sys.stderr.write('[ERROR] Parameters with identical names encountered!\nExecution aborted.\n'.format())
        sys.exit(1)

def getParamNameDict(paramsList, all=True):
    dictionary = {}
    for param in paramsList:
        if all or param['tagAttributes']['positional'] == "false":
            dictionary[param['tagAttributes']['name']] = param['tagAttributes']
    return(dictionary)

def getParamPositionDict(inDict):
    paramsList = inDict['inputs'] + inDict['outputs'] + inDict['parameters']
    outDict = {}
    for param in paramsList:
        key = param['tagAttributes']['positional']
        if key != "false":
            if key in outDict.keys():
                sys.stderr.write('[ERROR] Parameters with identical positions ("{pos}") encountered!\nExecution aborted.\n'.format(pos=key))
                sys.exit(1)
            else:
                outDict[key] = param['tagAttributes']
    return(outDict)

def addExecutable(dictionary):
    executable = dictionary['_executable']['name']
    prefix = '{' + executable + '}'
    suffix = executable + '=' + executable
    return(prefix, suffix)
    pass

def addParameters(dictionary, prefix, suffix, separator='$$$', spacer='^^^'):
    for key in sorted(dictionary.keys()):
        option = dictionary[key]['optionName']
        if 'redirect' in dictionary[key].keys() and dictionary[key]['redirect'] != "false":
            option = dictionary[key]['redirect']
        if 'default' in dictionary[key].keys():
            prefix += separator + option + spacer + '{' + dictionary[key]['name'] + '}'
            suffix += ', ' + dictionary[key]['name'] + '=' + dictionary[key]['name']
        else:
            prefix += separator + option + spacer + dictionary[key]['name']
    return(prefix, suffix)
    pass

def mergePrefixSuffix(prefix, suffix):
    command = "'" + prefix + "'" + '.format' + '(' + suffix + ')'
    return(command)
    pass


###################################
#<--- CREATE OUTPUT DIRECTORY --->#
###################################
def createOutDir(path):
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise
    pass


#############################
#<--- FORMAT DOCSTRINGS --->#
#############################
def formatDocStrings(dictionary):
    dictionary = formatDocSection(dictionary, 'header', ['doc'])
    dictionary = formatDocSection(dictionary, 'type-parameters')
    dictionary = formatDocSection(dictionary, 'inputs')
    dictionary = formatDocSection(dictionary, 'outputs')
    dictionary = formatDocSection(dictionary, 'parameters')
    dictionary = formatDocSection(dictionary, 'internalParameters')
    return(dictionary)
    pass

def formatDocSection(dictionary, key, tagFilter=[]):
    for element in dictionary[key]:
        if tagFilter and element['tagName'] not in tagFilter:
            continue
        doc = element['tagValue']
        doc = doc[0].upper() + doc[1:]
        if doc[-1] not in '.!?':
            doc += '.'
        element['tagValue'] = doc
    return(dictionary)
    pass


######################
#<--- WRITE JSON --->#
######################
def writeJSON(dictionary, outfile):
    import json
    with open(outfile, 'w') as file:
        json.dump(dictionary, file)
    pass


#####################
#<--- WRITE XML --->#
#####################
def writeXML(dictionary, outfile):
    filehandle = open(outfile, 'w')
    filehandle.write('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n')
    filehandle.write(XMLtag('component'))
    printTagOneLiners(filehandle, 'header', dictionary)
    printTagOneLiners(filehandle, 'credits', dictionary, tagFilter=['reference'])
    printTagOneLiners(filehandle, 'categories', dictionary)
    printLauncherSection(filehandle, dictionary)
    printTagOneLiners(filehandle, 'requires', dictionary, ['name', 'URL', 'optional', 'version', 'type'])
    filehandle.write(XMLtag('type-parameters', indent=1))
    printTagDocThreeLiners(filehandle, 'type-parameters', dictionary, ['name', 'extends'])
    filehandle.write(XMLtag('type-parameters', indent=1, close=True))
    filehandle.write(XMLtag('inputs', indent=1))
    printTagDocThreeLiners(filehandle, 'inputs', dictionary, ['name', 'type', 'optional', 'array'])
    filehandle.write(XMLtag('inputs', indent=1, close=True))
    filehandle.write(XMLtag('outputs', indent=1))
    printTagDocThreeLiners(filehandle, 'outputs', dictionary, ['name', 'type', 'array'])
    filehandle.write(XMLtag('outputs', indent=1, close=True))
    filehandle.write(XMLtag('parameters', indent=1))
    printTagDocThreeLiners(filehandle, 'internalParameters', dictionary, ['name', 'type', 'default'])
    printTagDocThreeLiners(filehandle, 'optOutBoolParameters', dictionary, ['name', 'type', 'default'])
    printTagDocThreeLiners(filehandle, 'parameters', dictionary, ['name', 'type', 'default'])
    filehandle.write(XMLtag('parameters', indent=1, close=True))
    filehandle.write(XMLtag('component', close=True))
    filehandle.close()
    pass

def printTagOneLiners(filehandle, key, dictionary, tagFilter=False, attrFilter=False, indent=1):
    for element in dictionary[key]:
        if tagFilter and element['tagName'] in tagFilter:
            continue 
        filehandle.write(XMLtag(element['tagName'], attributes=element['tagAttributes'], attributesFilter=attrFilter, indent=indent, newline=False))
        if not element['tagValue'] == 'n/a':
            filehandle.write(element['tagValue'])
        filehandle.write(XMLtag(element['tagName'], indent=0, close=True))
    pass

def printLauncherSection(filehandle, dictionary, indent=1):
    filehandle.write(XMLtag(dictionary['launcher'][0]['tagName'], attributes=dictionary['launcher'][0]['tagAttributes'], indent=indent))    
    for element in dictionary['launcherArguments']:
        filehandle.write(XMLtag(element['tagName'], attributes=element['tagAttributes'], indent=indent+1, singleton=True))
    filehandle.write(XMLtag(dictionary['launcher'][0]['tagName'], indent=indent, close=True))    
    pass

def printTagDocThreeLiners(filehandle, key, dictionary, attrFilter=False, indent=2):
    for element in dictionary[key]:
        filehandle.write(XMLtag(element['tagName'], attributes=element['tagAttributes'], attributesFilter=attrFilter, indent=indent))
        filehandle.write(XMLtag('doc', indent=indent+1, newline=False))
        docstring = element['tagValue'].capitalize()
        
        filehandle.write(element['tagValue'])
        filehandle.write(XMLtag('doc', indent=0, close=True))
        filehandle.write(XMLtag(element['tagName'], indent=indent, close=True))
    pass

def XMLtag(string, attributes=False, attributesFilter=False, filterExclude=False, indent=0, separator='    ', close=False, singleton=False, newline=True):
    attr = ''
    if singleton:
        suffix = ' />'
    else:
        suffix = '>'
    if close:
        prefix = '</'
    else:
        prefix = '<'
        if attributes:
            attributes = dict(attributes)
            if attributesFilter:
                if filterExclude:
                    for key in attributesFilter:
                        del attributes[key]
                    for key, value in sorted(attributes.items()):
                        attr += ' {0}="{1}"'.format(key, value)
                else:
                    for key in (set(attributes) - set(attributesFilter)):
                        del attributes[key]
                    attributesList = sorted(attributes.items(), key=lambda i:attributesFilter.index(i[0]))
                    for element in attributesList:
                        attr += ' {0}="{1}"'.format(element[0], element[1]) 
            else:
                for key, value in sorted(attributes.items()):
                    attr += ' {0}="{1}"'.format(key, value)
    if newline:
        newline = '\n'
    else:
        newline = ''
    return(str(separator) * int(indent) + prefix + str(string) + attr + suffix + newline)
    pass


################################
#<--- WRITE PYTHON WRAPPER --->#
################################
def writePyWrapper(command, outfile):
    with open(outfile, 'w') as file:
        file.write('#!/usr/bin/env python\n')
        file.write('\n')
        file.write('# Import modules\n')
        file.write('import sys\n')
        file.write('import anduril\n')
        file.write('import krini_functions\n')
        file.write('from anduril.args import *\n')
        file.write('\n')
        file.write('# Command template\n')
        file.write('command = ' + command + '\n')
        file.write('\n')
        file.write('# Execute command\n')
        file.write('exit_status = krini_functions.main(component, command, tempdir)\n')
        file.write('\n')
        file.write('# Return exit status\n')
        file.write('sys.exit(exit_status)\n')
    pass


##################
#<--- FOOTER --->#
##################
if __name__ == '__main__':
    sys.exit(main())
