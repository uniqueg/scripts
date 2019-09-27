#!/usr/bin/env python

#===========#
#  HEAD //  #
#===========#

# (c) 2016, Alexander Kanitz, Biozentrum, University of Basel

#===========#
#  // HEAD  #
#===========#


#===========#
#  TODO //  #
#===========#

# re-do score calculation (implement scale, information content, ...)
# update description
# add notes:
#  - expected matrix format
#  - matrix types and differences
#  - applied pseudocounts and matrix conversions
#  - scaling
#  - information content
# add license
# add author info
# handle formatting of usage, notes, license, author info, description, ...

#===========#
#  // TODO  #
#===========#


#=====================#
#  PRE-REQUISITES //  #
#=====================#

from __future__ import print_function

import sys, os
import time
import re
import math
import string

#=====================#
#  // PRE-REQUISITES  #
#=====================#


#================#
#  FUNCTIONS //  #
#================#

# Parse CLI arguments/options
def parseOpts():
    import optparse
    description = "Scans sequences in a FASTA file for the presence of binding motifs as represented by a position matrix. Calculates probabilities for each position of each sequence and reports them, together with the corresponding coordinates, if log probabilities are above a specified threshold."
    notes = formatUsageNotes(["Position weight matrix must contain a line 'Pos \t A \t C \t G \t T', followed by lines listing the nucleotide position and nucleotide probabilities (0 to 1).", "Some other note."])
    version = '%prog version 1.0'
    parser = optparse.OptionParser(description=description, version=version, epilog=notes)
    parser.add_option("--matrix-type", type="string", dest="mtype", default="ppm", help="Type of matrix. Allowed values are 'pfm', 'ppm' (default) and 'pwm'. See notes for details.", metavar="STRING")
    parser.add_option("--ppm-pseudocount", type="float", dest="pseudo", default=0.001, help="The specified value is used to replace 0 values in position probability matrices (default: 0.001). See notes for details.", metavar="FLOAT")
    parser.add_option("--scale", type="string", dest="scale", default="minmax", help="Report scaled scores in addition to probability products. Allowed values are 'max', 'minmax' (default) and 'none'. See notes for details.", metavar="STRING")
    parser.add_option("--threshold", type="float", dest="threshold", default=0.75, help="Report only those sites whose *scaled* score equals or exceeds the specified value (default: 0.75). Ignored if scaling is turned off.", metavar="FLOAT")
    parser.add_option("--information-content", action="store_true", dest="entropy", default=False, help="Consider positional information content when calculating scores. See notes for details.")
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False, help="Print log messages.")
    parser.add_option("--debug", action="store_true", dest="debug", default=False, help="Print debugging information.")
    #parser.add_option("--usage", action="store_true", dest="help", default=False, help="Print this help screen.")
    #parser.add_option("--version", action="store_true", dest="license", default=False, help="Print version.")
    parser.add_option("--license", action="store_true", dest="verbose", default=False, help="Print license information.")
    license = ""
    (opts, args) = parser.parse_args()
    # TODO: Call 'help', 'usage', 'license' and 'version' functions
    return(opts, args)

# Format notes string
def formatUsageNotes(stringList, prefix='[NOTES]', spacer='  ', bullet='- '):
    stringList.insert(0, prefix)
    string = "{spacer}{bullet}".format(spacer=spacer, bullet=bullet).join(stringList)
    return(string)

# Get formatted time string
def getLogTime():
    return(time.strftime("%d-%b-%Y %H:%M:%S"))

# Print arguments
def printArgs(args):
    for arg in args:
        print("  - {}".format(arg), file=sys.stderr)

# Print options
def printOpts(opts):
    for key, value in opts.__dict__.items():
        print("  - {key}: {value}".format(key=key, value=value), file=sys.stderr)

# Ensure that the correct number of arguments is given
def verifyArgsCount(args, exactly=None, atLeast=None, atMost=None):
    if exactly:
        if atLeast or atMost:
            sys.exit("[{time}] [ERROR] Specify either 'exactly' or one or both of 'atLeast' and 'atMost'. Aborted.".format(time=getLogTime()))
        atLeast = atMost = exactly
    if len(args) < atLeast:
        sys.exit("{time}] [ERROR] At least {minVal} arguments required. Aborted.".format(time=getLogTime(), minVal=atLeast))
    if atMost and len(args) > atMost:
        sys.exit("{time}] [ERROR] At most {maxVal} arguments allowed. Aborted.".format(time=getLogTime(), maxVal=atMost))

# Add arguments to options dictionary
def addArgsToOpts(opts, args, keys=None):
    if keys:
        if not len(keys) is len(args):
            sys.exit("{time}] [ERROR] Number of keywords ({keys}) does not equal number of arguments ({args}). Aborted.".format(time=getLogTime(), keys=len(keys), args=len(args)))
    else:
        keys = []
        for i in range(1, len(args)):
            keys.append("pos_{}".format(i))
    opts.__dict__.update(dict(zip(keys, args)))
    return(opts)

# Verify and process options
def processOpts(opts):
    verifyFiles([opts.fasta, opts.mt])
    verifyStringSet(opts.mtype, ['pfm', 'ppm', 'pwm'], context='option --matrix-type')
    verifyStringSet(opts.scale, ['none', 'minmax', 'max'], context='option --scale')
    verifyFloatRange(opts.threshold, atLeast=0.0, atMost=1.0, context='option --threshold')
    verifyFloatRange(opts.pseudo, atLeast=1e-10, atMost=1e-1, context='option --ppm-pseudocount')
    return(opts)

# Verify that files are available and read permissions are set
def verifyFiles(fileList):
    for f in fileList:
        if not os.path.exists(f):
            sys.exit("{time}] [ERROR] Path '{path}' is not available. Aborted.".format(time=getLogTime(), path=f))
        if not os.path.isfile(f):
            sys.exit("{time}] [ERROR] '{path}' is not a file. Aborted.".format(time=getLogTime(), path=f))
        try:
            with open(f, 'r') as fh:
                pass
        except IOError:
            sys.exit("{time}] [ERROR] '{path}' cannot be opened. Aborted.".format(time=getLogTime(), path=f))

# Verify if an expression is or can be coerced to type 'string' and it is included in a set of 
# allowed strings
def verifyStringSet(value, stringSet, context='N/A'):
    if not isinstance(value, str):
        try:
            value = str(value)
        except ValueError:
            sys.exit("{time}] [ERROR] Value '{value}' in context '{context}' is not of and cannot be coerced to type 'string'. Aborted".format(time=getLogTime(), value=value, context=context))
    if not value in stringSet:
        sys.exit("[{time}] [ERROR] The specified value '({value})' in context '{context}' is not allowed. See usage. Aborted.".format(time=getLogTime(), value=value, minVal=atLeast, context=context))

# Verify if an expression is or can be coerced to type 'float' and its value is in allowed range
def verifyFloatRange(value, atLeast=None, atMost=None, context='N/A'):
    if not isinstance(value, float):
        try:
            value = float(value)
        except ValueError:
            sys.exit("{time}] [ERROR] Value '{value}' in context '{context}' is not of and cannot be coerced to type 'float'. Aborted".format(time=getLogTime(), value=value, context=context))
    if value < atLeast:
        sys.exit("[{time}] [ERROR] The specified value '({value})' in context '{context}' is smaller than the allowed minimum value ({minVal}). Aborted.".format(time=getLogTime(), value=value, minVal=atLeast, context=context))
    if value > atMost:
        sys.exit("[{time}] [ERROR] The specified value '({value})' in context '{context}' is bigger than the allowed maximum value ({minVal}). Aborted.".format(time=getLogTime(), value=value, minVal=atMost, context=context))

# Read matrix
def readMatrix(matrixFile):
    header = True
    reHeader = re.compile(r"^Pos\tA\tC\tG\t[TU]$")
    mt = []
    with open(matrixFile, 'r') as mtHandle:
        for line in mtHandle:
            if reHeader.match(line):
                header = False
                continue
            line = line.rstrip()
            if line is "":
                continue
            if not header:
                line = line.rstrip()
                mt.append(map(float, line.split("\t")[1:]))
    if not mt:
        sys.exit("[{time}] [ERROR] No (valid) entries in position weight matrix. Aborted.".format(time=getLogTime()))
    return(mt)

# Convert PFM to PPM in place
def ppmFromPfm(mt):
    flat = [nuc for pos in mt for nuc in pos]
    for num in flat:
        if num < 0:
            sys.exit("[{time}] [ERROR] Matrix contains negative value '{num}'. Check matrix and matrix type. Aborted.".format(time=getLogTime(), num=num))
    if 0 in flat:
        for position in mt:
            position[:] = [ num + 1 for num in position ]
    for position in mt:
        position[:] = [ num / sum(position) for num in position ]

# Convert PPM to PWM in place
def pwmFromPpm(mt, pseudo=0.001):
    flat = [nuc for pos in mt for nuc in pos]
    for num in flat:
        if not 0 <= num <= 1:
            sys.exit("[{time}] [ERROR] Matrix contains illegal probability value '{num}'. Check matrix and matrix type. Aborted.".format(time=getLogTime(), num=num))
    if 0 in flat:
        for position in mt:
            position[:] = [ num if num > 0 else pseudo for num in position]
    for position in mt:
        position[:] = [ math.log(num * 4, 2) for num in position ]

# Get PWM information content vector
def getPwmInfoContent(mt):
    infCont = []
    for position in mt:
        posInfCont = 0
        for nuc in position:
            posInfCont += ( math.pow(2, nuc) / 4 ) * nuc
        posInfCont *= -1
        infCont.append(posInfCont)
    return(infCont)

# Get minimum possible probability
def getMinLikelihood(mt, entropy=False):
    minL = 0
    if not entropy:
        entropy = len(mt) * [float(-1)]
    for mtPos, entPos in zip(mt, entropy):
        minL += -entPos * min(mtPos)
    return(minL)

# Get maximum possible probability
def getMaxLikelihood(mt, entropy=False):
    maxL = 0
    if not entropy:
        entropy = len(mt) * [float(-1)]
    for mtPos, entPos in zip(mt, entropy):
        maxL += -entPos * max(mtPos)
    return(maxL)

# Scan FASTA file
def scanFasta(fasta, pwm, scale='none', threshold=0.0, minL=None, maxL=None, entropy=False, verbose=True, chunksize=1000):
    pwmLength = len(pwm)
    allResults = []
    with open(fasta, 'r') as fastaHandle:
        for count, (name, seq) in enumerate(readFasta(fastaHandle)):
            seqNum = processSeq(seq)
            seqResults = scanSeq(name=name, seqNum=seqNum, seq=seq, pwm=pwm, pwmLength=pwmLength, scale=scale, threshold=threshold, minL=minL, maxL=maxL, entropy=entropy)
            allResults += seqResults
            if not (count + 1) % chunksize:
                for r in allResults:
                    print("\t".join(r), file=sys.stdout)
                allResults = []
                if verbose:
                    print("[{time}] Processed {count} sequences...".format(time=getLogTime(), count=count + 1), file=sys.stderr)
        for r in allResults:
            print("\t".join(r), file=sys.stdout)
    fastaHandle.close()

# Read FASTA file line by line
def readFasta(handle):
    name, seq = None, []
    for line in handle:
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                yield (name, ''.join(seq).upper())
            name, seq = line[1:], []
        else:
            seq.append(line)
    if name:
        yield (name, ''.join(seq).upper())

# Process sequence
def processSeq(seq):
    if re.match(r'[^ACGTN]', seq):
        print("{time}] [WARNING] Sequence '{seq}' contains unrecognized characters. Converted to 'N'.".format(time=getLogTime(), seq=seq), file=sys.stderr)
        seq = re.sub(r'[^ACGTN]', 'N', seq)
    seq = seq.translate(string.maketrans("ACGTN", "01234"))
    return(seq)

# Scan sequence
def scanSeq(name, seqNum, seq, pwm, pwmLength, scale='none', threshold=0.0, minL=None, maxL=None, entropy=False):
    seqNum = map(int, list(seqNum))
    if not entropy:
        entropy = len(pwm) * [float(-1)]
    resultsAll = []
    for winStart in range(len(seqNum) - pwmLength + 1):
        win = seqNum[winStart:winStart + pwmLength]
        winSeq = seq[winStart:winStart + pwmLength]
        resultsWin = [name, str(winStart), str(winStart + pwmLength), winSeq]
        logSum = 0
        for pos, nuc in enumerate(win):
            if nuc is 4:
                continue
            logSum += -entropy[pos] * pwm[pos][nuc]
        resultsWin.append(str(logSum))
        if not 'scale' is 'none':
            logSumScaled = scaleLogLikelihood(scale=scale, logL=logSum, minL=minL, maxL=maxL)
            resultsWin.append(str(logSumScaled))
            if logSumScaled > threshold:
                resultsAll.append(resultsWin)
        else:
            resultsAll.append(resultsWin)
    return(resultsAll)

# Scale sum of log likelihoods
def scaleLogLikelihood(scale, logL, minL, maxL):
    if str(scale) == 'minmax':
        scaledL = (logL - minL) / (maxL - minL)
    elif str(scale) == 'max':
        scaledL = logL / maxL
    else:
        sys.exit("[{time}] [ERROR] Illegal scaling type '{scale}'. See usage. Aborted.".format(time=getLogTime(), scale=str(scale)))
    return(scaledL)

#================#
#  // FUNCTIONS  #
#================#


#===========#
#  MAIN //  #
#===========#

def main():

    opts, args = parseOpts()
    if opts.verbose:
        print("[{time}] Processing arguments and options...".format(time=getLogTime()), file=sys.stderr)
    if opts.debug:
        print("[{time}] [DEBUG] Arguments and options before processing:".format(time=getLogTime()), file=sys.stderr)
        print("Positional arguments", file=sys.stderr)
        printArgs(args)
        print("Options and values", file=sys.stderr)
        printOpts(opts)
    verifyArgsCount(args, exactly=2)
    opts = addArgsToOpts(opts, args, keys=['fasta', 'mt'])
    opts = processOpts(opts)
    if opts.debug:
        print("[{time}] [DEBUG] Options after processing:".format(time=getLogTime()), file=sys.stderr)
        printOpts(opts)

    if opts.verbose:
        print("[{time}] Reading and processing matrix...".format(time=getLogTime()), file=sys.stderr)
    opts.mt = readMatrix(opts.mt)
    if opts.debug:
        print("[{time}] [DEBUG] Type of matrix: {type}".format(time=getLogTime(), type=opts.mtype), file=sys.stderr)
        print("[{time}] [DEBUG] Length of matrix: {length}".format(time=getLogTime(), length=len(opts.mt)), file=sys.stderr)
        print("[{time}] [DEBUG] Matrix:".format(time=getLogTime()), file=sys.stderr)
        print(opts.mt, file=sys.stderr)
    if opts.mtype is 'pfm':
        ppmFromPfm(mt=opts.mt)
        opts.mtype = 'ppm'
        if opts.debug:
            print("[{time}] [DEBUG] Type of matrix: {type}".format(time=getLogTime(), type=opts.mtype), file=sys.stderr)
            print("[{time}] [DEBUG] Length of matrix: {length}".format(time=getLogTime(), length=len(opts.mt)), file=sys.stderr)
            print("[{time}] [DEBUG] Matrix:".format(time=getLogTime()), file=sys.stderr)
            print(opts.mt, file=sys.stderr)
    if opts.mtype is 'ppm':
        pwmFromPpm(mt=opts.mt, pseudo=opts.pseudo)
        opts.mtype = 'pwm'
        if opts.debug:
            print("[{time}] [DEBUG] Type of matrix: {type}".format(time=getLogTime(), type=opts.mtype), file=sys.stderr)
            print("[{time}] [DEBUG] Length of matrix: {length}".format(time=getLogTime(), length=len(opts.mt)), file=sys.stderr)
            print("[{time}] [DEBUG] Matrix:".format(time=getLogTime()), file=sys.stderr)
            print(opts.mt, file=sys.stderr)
    if opts.entropy:
        opts.entropy = getPwmInfoContent(opts.mt)
        if opts.debug:
            print("[{time}] [DEBUG] Matrix information content:".format(time=getLogTime()), file=sys.stderr)
            print(opts.entropy, file=sys.stderr)
    if str(opts.scale) == 'minmax':
        opts.minL = getMinLikelihood(opts.mt, opts.entropy)
        if opts.debug:
            print("[{time}] [DEBUG] Minimum likelihood: {prob}".format(time=getLogTime(), prob=opts.minL), file=sys.stderr)
    else:
        opts.minL = None
    if not str(opts.scale) == 'none':
        opts.maxL = getMaxLikelihood(opts.mt, opts.entropy)
        if opts.debug:
            print("[{time}] [DEBUG] Maximum likelihood: {prob}".format(time=getLogTime(), prob=opts.maxL), file=sys.stderr)
    else:
        opts.maxL = None

    if opts.verbose:
        print("[{time}] Scanning sequences...".format(time=getLogTime()), file=sys.stderr)
    scanFasta(fasta=opts.fasta, pwm=opts.mt, threshold=opts.threshold, scale=opts.scale, minL=opts.minL, maxL=opts.maxL, entropy=opts.entropy, verbose=opts.verbose, chunksize=1000)
    if opts.verbose:
        print("[{time}] Done.".format(time=getLogTime()), file=sys.stderr)
    sys.exit(0)

if __name__ == "__main__":
    main()

#===========#
#  // MAIN  #
#===========#
