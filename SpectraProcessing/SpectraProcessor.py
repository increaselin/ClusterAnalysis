#! /usr/bin/env python
# 
# Program: SpectraProcessor.py
#
# Author: Mike Lum
#
# Usage: ./SpectraProcessor.py [-vh]
#
# Description: Task waits for machine idle time (non-working hours), scans the 
#               input directory for new spectra files, runs the eqw measuring program,
#               storing the log in my Dropbox folder, and then moving the measured file into
#               a "processing completed" directory.
#
# Revision History:
#    Date        Vers.    Author        Description
#    2016-03-01  1.0f2    Lum        Relocated some code, had to change the include list.
#    2015-07-14  1.0f2    Lum        Change from command line call to a module-based system
#    2015-06-20  1.0f1    Lum        Produce custom line lists for each spectrum
#    2014        0.0a0    Lum        First checked in
#
# To Do:
#    
#

from collections import deque
import os
import sys
import datetime as dt
import time
import glob
import random
import shutil
#import subprocess as sub

from SpectraProcessing import SpectraData as SD
from constants import constants as k
from constants import directories as d
from constants import DBFormats as DB
from Databases import LineLookup as ll
from EqWidthMeasure import eqwm
from MOOGInterface import MOOGInterface as MI

thisScriptName = 'SpectraProcessor.py'
interpreterName = 'python'

# Constants
kWeekdays=['Monday','Tuesday','Wednesday','Thursday','Friday','Saturday','Sunday']
#kEQWMCommand = 'python -m EqWidthMeasure/eqwm'
#kEQWMFlags = '-a -v'


# Globals
gVerboseMode = False
gListFiles = []

# Functions
def printHelpText():
    print('Program: SpectraProcessor.py\n')
    print('Task waits for machine idle time (non-working hours), scans the ')
    print('input directory for new spectra files, runs the eqw measuring program,')
    print('storing the log in my Dropbox folder, and then moving the measured file into')
    print('a "processing completed" directory.\n')
    # Flags and options here
    print('Usage: ./SpectraProcessor.py [-hv]\n')
    print('Options:')
    print('-h: Print this help text.')
    print('-v: Use verbose progress and error messages')

# Utility functions
def dow():
    # dow = Day of the Week
    return dt.datetime.today().weekday()

def IsTimeToSleep():
    hour = dt.datetime.now().hour
    if hour > 9 and hour < 16 and dow() < 5:
    #    return True
    # Always run when I say to, now
        return False

    else:
        return False

def errorOutput(errorStr):
    if gVerboseMode:
        timestamp = dt.datetime.now().strftime("%y/%m/%d %H:%M:%S")
        print(timestamp + ": Error in processing: {0}.".format(errorStr))
    return
    
def GoToSleep():
    # Wakeup time is 8:00PM, tonight
    now = dt.datetime.now()
    wakeupTime = dt.datetime(year=now.year, month=now.month, day=now.day, hour=20)
    sleepTime = (wakeupTime - now).seconds
    time.sleep(sleepTime)
#    print('Sleep time would be {0:d} seconds. Sleeping for 1 second, instead.'.format(sleepTime))
#    time.sleep(1)
#    print('Waking.')
    return
    
def GetFileToProcess():
    fileList = sorted(glob.glob(d.ToBeProcessed+'/*.fits'))
    if len(fileList) == 0:
        return ''
    
    now=time.time()
    oldestFile = fileList[0], now-os.path.getctime(fileList[0])
    if len(fileList) > 1:
        for thisFile in fileList[1:]:
            age = now-os.path.getctime(thisFile)
            if age > oldestFile[1]:
                oldestFile = thisFile, age
    return oldestFile[0]

def GetListFiles():

    global gListFiles
    
    lisFileList = glob.glob(k.ListDirectory+'/*.lis')
    
    gListFiles = []
    
    lineList = []
    
    for lisFile in lisFileList:
        lineList.extend(MI.ReadMOOGLineFile(lisFile))
    
    # Pyspeckit has a serious memory leak, due to circular references in object definitions.
    # It will run out of memory after processing about 80 lines on a 2GB machine. So, we need
    # to break out line list into <80 line segments, and call the processor on each list. 
    minIndex = 0
    maxIndex = minIndex+k.MaxLineListFileLength
    listFileCount = 1
    
    tempLisFileLead = k.ListDirectory+k.TempListDirectory+'/'+k.TempListFileHead
    
    while True:
        maxIndex = minIndex+k.MaxLineListFileLength
        if maxIndex > len(lineList):
            maxIndex = len(lineList)
        newFileName = "{0}{1:03d}{2}".format(tempLisFileLead, listFileCount, k.TempListFileExt)
        with open(newFileName,'w') as listFile:
            for line in lineList[minIndex:maxIndex]:
                listFile.write("{0:8.3f} {1:10.1f} {2:11.4f} {3:11.5f} {4:4.1f} {5:21.1f}          {6}\n".format(line[0], line[1], line[2], line[3], line[4], line[5], line[6]))
        
        gListFiles.append(newFileName)
        
        if maxIndex == len(lineList):
            break
        else:
            minIndex = maxIndex
            listFileCount += 1
        
        
def makeListFiles(spectraFile, listDirectory, remeasure = False):
# Function takes the passed spectra, and builds a list of lines
# which should lie within the spectra's orders. Due to a problem
# with the automated measuring process, the list is split into
# multiple files, with no more than 200 in each list.
# Output lists are stored in the passed directory, with the full
# path name(s) returned.
# The remeasure parameter tells whether to re-measure a line for
# a given spectra, if it has been measured _FOR THE PASSED FILE_
# previously.
    lineListFileList = []
# Determine wavelength range coverage of this spectra
    objectName, numApertures, apertureSize, apBoundaries = SD.ReadSpectraInfo(spectraFile)
    
# Remove the 'bad' regions
    newRegions = []
    for mask in k.BadSpectralRegions:
        for actual in apBoundaries:
            if mask[1] < actual[0] or mask[0] > actual[1]:
            # No overlap
                pass
            elif mask[0] < actual[0] and mask[1] > actual[1]:
            # Complete overlap
                apBoundaries.remove(actual)
            elif mask[1] < actual[1]:
            # Overlap on the low side of the actual:
                if mask[0] > actual[0]:
                # Mask is internal to the region - need to split
                    apBoundaries.append((actual[0], mask[0]))
                apBoundaries.append((mask[1], actual[1]))
                apBoundaries.remove(actual)
            else: # mask[0] > actual[0]
            # Overlap on the high side of the actual:
                apBoundaries.append((actual[0], mask[0]))
                apBoundaries.remove(actual)
                
# Do the lookup (return format is in 'MOOG' text format)
    lineData, refList = ll.getLines([], apBoundaries, dataFormat='')
    
    # We do not want to re-measure lines (or maybe we do...?)
    if remeasure:
        measuredLines = []
    else:
        measuredLines = ll.getMeasuredLines(spectraFile)
    unmeasuredLines = [line for line in lineData if (line[1], line[0]) not in measuredLines]
    
    # Eliminate duplicate lines from overlapping windows and 
    # sort the returned lines by wavelength
    sortedLines = sorted(set(unmeasuredLines), key=lambda line: line[0])
   
# Write the list out in 75 line segments. This is done to prevent the memory leak
# issue in pyspeckit, when a large number of lines are measured.
    indexList = range(0, len(sortedLines), k.MaxLineListFileLength)
    lineListFileList = []
    
    for index, lineIndex in enumerate(indexList):
        # Note: For mp runs, we'll want to add a time or proc.# stamp to the log file name
        # File enumeration is 1-based for readability
        lineListFilename = listDirectory+k.TempListFileHead+'{0}'.format(index+1)+k.TempListFileExt
        firstLine = '# Line list ({0}/{1}) for spectra:{2}  Created:{3}\n'.format(index+1, len(indexList), spectraFile,  dt.datetime.now().strftime("%Y/%m/%d %H:%M:%S"))
        try:
            outfile = open(lineListFilename, 'w')
        except IOError:
            errorOutput('Unable to open:{0} for writing'.format(lineListFilename))
        outfile.write(firstLine)
        
        if index < len(sortedLines) - k.MaxLineListFileLength:
            outputList = sortedLines[lineIndex : lineIndex+k.MaxLineListFileLength]
        else:
            outputList = sortedLines[lineIndex:]
        for lineData in outputList:
            outfile.write(k.MOOGListFormat.format(lineData[0], lineData[1], lineData[2], lineData[3], lineData[4], 99.9, lineData[5])+'\n')
        outfile.close()    
        lineListFileList.append(lineListFilename)
# Return our list of list files
    return lineListFileList

def cleanFiles(fileList):
    for filename in fileList:
        try:
            os.remove(filename)
        except OSError:
            # It's aaaaaallllready gone (Eagles fans?)
            # ...or a directory. Either way, we don't care
            pass
    return
    
def parseInFilename(spectraFile):
    logFile = ''
    logDir = ''
    
    # We assume that the input spectra filename has already been changed to 
    # reflect our desired format: 
    #           NNNN_SSS_ssssss_II_OOO.fits
    # Where: 
    #           NNNN  = NGC catalog number of parent cluster 
    #                    (note: clusters not contained in the NGC catalog will be prefaced by a letter
    #                     representative of their catalog. ex: I = IC
    #           SSS    = Reference indexing system for this cluster. 
    #                    (ex: MJP = Montgomery, Janes, and Phelps)
    #           ssssss = Star index number from the aformetioned index.
    #           II     = Two letter code noting the observing telescope/instrument
    #                    reference: 
    #                               KH = Keck Hires
    #                               VU = VLT UVES
    #                               EL = Elodie
    #                               SO = Sophie
    #                               HA = HARPS
    #                               SH = Subaru HDS
    #           OOO    = Internal numbering for this spectra.
    
    try:
        pathRef, indexName, indexNo, instrument, endString = spectraFile.split('_')
        internalNum = endString.replace('.fits','')
    except:
        errorOutput("Unable to split/parse input spectra filename. Continuing.")
        return logFile, logDir
        
    # pathRef will contain the entire path of the filename, so parse out everything before the last '/'
    catRef = pathRef.split('/')[-1]
    # truncate the logfile extension for now. We will add on the .log extension and a counter, later
    logFile = spectraFile.replace('.fits','')
    
    catStr = 'NGC-'
    
    try:
        catNo = int(catRef)
    except ValueError:
        # Assume the conversion failed due to the presence of a letter, representing
        # a catalog other than NGC. If the assumption is incorrect, we're gonna fail with 
        # another ValueError...
        catStr = catRef[0]
        catNo = int(catRef[1:])

    logDir = d.EQWMLogs + '/' + catStr + '{0:04d}'.format(catNo)
    if not os.path.isdir(logDir):
        print ('Creating log directory, {0} for {1}{2}.'.format(logDir, catStr, catNo))
        os.makedirs(logDir)
        
    return logFile, logDir
    

def SpectraProcessor():
    # Returns True, as long as there might be files left to process
    
    thisFile = GetFileToProcess()
    if thisFile == '':
        return False
    # The Voigt fitting algorithms use the DB to look up lines to measure
    eqwm.eqwm(interactiveFitting = False, spectraFilename = thisFile, continuum1 = False)
    
    try:
        shutil.move(thisFile, d.ProcessedSpectra)
    except shutil.Error:
        print('Error in moving log or spectra file to destination. Process continuing.')
    
    return True

if __name__ == '__main__':

    # Parse your command line call
    temp = os.sys.argv
    argv = deque(temp)
    
    while len(argv) > 0:
        flag = argv.popleft()
        if flag == '-h':
            printHelpText()
            exit()
        elif flag == '-v':
            gVerboseMode = True
            print('Verbose Mode enabled.')
        elif str.find(flag,thisScriptName) != -1 :
            if gVerboseMode: print('Executing command.')
        elif str.find(flag,interpreterName) != -1 :
            if gVerboseMode: print('...')
        #else: ...
        # Ignore everything else, or not...

    while True:
        while SpectraProcessor():
            pass

    # Processor completed the available files, so sleep for a half hour, and check for more
        print('Sleeping for a bit...')
        time.sleep(1800)


    exit()
