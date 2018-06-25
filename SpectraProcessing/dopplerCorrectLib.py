#! /usr/bin/env python
# 
# Program: dopcorandstack.py
#
# Author: Mike Lum
#
# Usage: (Library calls) dopplerCorrectLib.py
#
# Description: Script to measure Doppler shift, and return a recommended Doppler
#               correction for the spectra
#
# Revision History:
#    Date        Vers.    Author        Description
#    04/05/2018   1.0f1    Lum        Code cleanup for github submission
#    Today        0.0a0    Lum        First checked in
#
# To Do:
#    
#

from constants import constants as k
from EqWidthMeasure import eqwm as EQW
from SpectraProcessing import SpectraData as SD

from datetime import datetime
import os
import numpy as np

from SpectraProcessing import dopcorLines as Dlines

# Functions
def makeDoplerListFile(lines):
    listFilename = k.MOOGTempListName[:-4]+datetime.now().strftime('%f')+k.MOOGTempListName[-4:]
    with open(listFilename, 'w') as listFile:
        listFile.write('# Temporary MOOG line list file: {0}'.format(datetime.now().strftime('%y-%b-%d:%H:%M:%S\n')))
        for line in lines:
            listFile.write('{0:8.3f}  {1:9.1f}  {2:6.2f}  {3:7.2f}                              99.9\n'.format(line[0], line[1], line[2], line[3]))
    return listFilename
    
def dopCor(inFile, interactive=True):

    objectName, numApertures, apertureSize, apBoundaries = SD.ReadSpectraInfo(inFile)
    hDoppler = 0.
    numHdoppler = 0
    if interactive:
        allLines = [line for line in Dlines.dopplerLines] # Generator, just in case...
        lineLogFile = makeDoplerListFile(allLines)
        EQWidths, badLines, dopplerCount, dopplerCheck = EQW.MeasureWidthsInter(inFile, lineLogFile)
        os.remove(lineLogFile)
    else:
        nonHLines = [line for line in Dlines.dopplerLines if line[1] != 1.0]
        lineLogFile = makeDoplerListFile(nonHLines)
        EQWidths, badLines, dopplerCount, dopplerCheck = EQW.MeasureWidthsAuto(inFile, lineFile=lineLogFile)
        os.remove(lineLogFile)
        
    dopplerCorrection = dopplerCheck/dopplerCount
    outputLines(dopplerCorrection, EQWidths)
    
    return dopplerCorrection, EQWidths

def outputLines(result, EQWidths):
    print('--------- Used lines: -----------')
    for line in EQWidths:
        print ('{0:7.2f}(Delta={3:+5.2f}) (ion:{2:4.1f}), weight: {1:3.1f}'.format(line[0], line[4], line[1], line[6]-line[0]))
    
    print('\n--- Overall (weighted) average: {0:+6.3f} Km/s'.format(result))
    return
    
# Once upon a time, we had this as a command line call for testing purposes.    
'''
if __name__ == '__main__':

    thisScriptName = 'N/A'      
    interpreterName = 'python'

    # Parse your command line call
    temp = os.sys.argv
    argv = deque(temp)
    inputFilename = ''
    interactive = True
    
    while len(argv) > 0:
        flag = argv.popleft()
        if flag == '-h':
            printHelpText()
            exit()
        elif flag == '-a':
            interactive = False
        elif flag == '-v':
            g.VerboseMode = True
            print('Verbose Mode enabled.')
        elif str.find(flag,thisScriptName) != -1 :
            if g.VerboseMode: print('Executing command.')
        elif str.find(flag,interpreterName) != -1 :
            if g.VerboseMode: print('...')
        else: 
            inputFilename = flag
    
    if inputFilename == '':
        print('No input spectra specified. Exiting.')
        exit()
    
    result, widths = dopCor(inputFilename, interactive)
    outputLines(result, widths)
    exit()
'''
