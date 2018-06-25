#! /usr/bin/env python
# 
# Module: MOOGInterface.py
#
# Author: Mike Lum
#
# Description: A series of functions designed to assist with calling MOOG 
#       functions. Mainly reading and writing MOOG-readable formats
#
# Revision History:
#    Date        Vers.    Author        Description
#    09/20/16    1.0f1    Lum        First checked in
#
# To Do:
#    
#

# Imports
import datetime
from collections import deque
import numpy as np
import os
import subprocess as sub
import sys
import time

from constants import constants as k
from constants import directories as dirs
from Databases import LineLookup as LL
from utilities import elements as EL
from utilities import utilities as u


# Functions

def PareMOOGLogFile(logName=k.MOOGTempLogName):
    minLines = 30 # Need 30 lines for significance in the Welch's t-test
    
    allLines = ReadMOOGLineFile(logName)
    allElements = list(set([line[1] for line in allLines]))
    
    lineData = []
    for elem in allElements:
        theseLines = [line for line in allLines if line[1] is elem]
        if len(theseLines) >= minLines:
            lineData.extend(theseLines)
    
    timestamp = datetime.datetime.now().isoformat()
    headStr = ' Line width measures pared from:{0} on: '.format(logName)+timestamp
    
    # No matter what log file name is passed, we are going to write over
    # the _temporary_ MOOG log file.
    np.savetxt(fname=k.MOOGTempLogName, X=lineData, fmt=k.MOOGLogFormat, header=headStr)
    
    return logName


def BuildMOOGParFile(parFilename=k.MOOGTempParName, modelName=k.MOOGTempModelName, logFilename=k.MOOGTempLogName, blendsDriver=False, ion=0.0):
    outFileHead = modelName[:-2]
    if blendsDriver:
        parFileStr = k.MOOGBlendParHead + k.MOOGParout1Head + outFileHead + '.out1\n' +\
                        k.MOOGParout2Head + outFileHead + '.out2\n' +\
                        k.MOOGParModelHead + modelName + '\n' +\
                        k.MOOGParLinesHead + logFilename + '\n' +\
                        k.MOOGBlendParTail.format(ion)
    else:
        parFileStr = k.MOOGParFileHead + k.MOOGParout1Head + outFileHead + '.out1\n' +\
                        k.MOOGParout2Head + outFileHead + '.out2\n' +\
                        k.MOOGParModelHead + modelName + '\n' +\
                        k.MOOGParLinesHead + logFilename + '\n'
    try:
        parFile = open(parFilename, 'w')
    except IOError:
        print('Unable to open {0} for writing. Exiting'.format(parFileName))
        sys.exit()
        
    parFile.write(parFileStr)
    parFile.close()
    
    return parFilename


def MakeMOOGScriptFile(parFilename, scriptFilename=k.MOOGScriptName):
    try:
        scriptFile = open(scriptFilename, 'w')
        scriptFile.write(parFilename+'\n')
        scriptFile.close()
    except IOError:
        print('Unable to open {0} for writing. Exiting'.format(scriptFilename))
        sys.exit()

    return scriptFilename


def MakeMOOGEQWLogFile(lines, logFilename=k.MOOGTempLogName):
    try:
        logFile = open(logFilename,'w')
        logFile.write('# Temporary MOOG log file created: {0}\n'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    except IOError:
        print('Unable to open MOOG script file: {0} for writing. Exiting.'.format(logFilename))
        sys.exit()
        
    for line in lines:
        lineInfo, lineRefs = LL.getLines(elements=[line[0]], wavelengthRange=[(line[1]-0.05, line[1]+0.05)], dataFormat=None)
        if len(lineInfo) > 0:
            try:
                logFile.write(k.MOOGListFormat.\
                format(line[1],line[0],lineInfo[0][2],lineInfo[0][3],lineInfo[0][4],line[4],'\n'))
            except ValueError:
            # Some unicode strings snuck in as floats...try typecasting
                logFile.write(k.MOOGListFormat.\
                format(float(line[1]),float(line[0]),float(lineInfo[0][2]),float(lineInfo[0][3]),float(lineInfo[0][4]),float(line[4]),'\n'))
    logFile.close()
    return logFilename
    

def ParseMOOGOut2File(out2Name):
    try:
        out2File = open(out2Name)
    except IOError:
        print('Unable to open {0} for line data read. Exiting.'.format(out2Name))
        sys.exit()
    
    fileData = out2File.read()
    out2File.close()

    lines = deque(str.split(fileData,'\n'))
    
    # MOOG can occasionally "bail" from an modeling run. Of course there's no easy
    # way to catch MOOG printing "...I QUIT!", so we manually read the second
    # line of the .out2 file: "Good" runs will have the first (comment) line
    # from the .log file. "Bad" runs will be the intermediate output product
    # starting with "element" - or nothing at all.
    waitForUser = False # For now, we're assuming automated operations:
    if len(lines)<3:
        print('MOOG unable to complete element:(unknown) for {0}.'.format(out2Name[:-4]))
        if k.isPython3:
            if waitForUser: input('Press enter to acknowledge and continue')
        else:
            if waitForUser: raw_input('Press enter to acknowledge and continue')
        return {}
        
    if lines[1].split()[0] == 'element':
        print('MOOG unable to complete element:{0} for {1}.'.format(lines[1].split()[1], out2Name[:-4]))
        if k.isPython3:
            if waitForUser: input('Press enter to acknowledge and continue')
        else:
            if waitForUser: raw_input('Press enter to acknowledge and continue')
        return {}

    # elementDict is a dictionary with the element symbol + 
    # ionization state (float) as key and an array of elements:
    # [Wavelength, Ex.Pot., logGf, eqw, logRW, abund]
    # for all the lines of that element as the value.
    # Example:   elementDict = {26.0:[[4347.23, 0.00, -5.503, 89.8, -4.68, 8.95]...]}
    elementDict = {}
    currentElement = 1.0 # HI is a flag for "No element"
    lineAbs = []
    
    while len(lines) > 1: 
        thisLine = lines.popleft().split();
        if len(thisLine) < 5:
            continue
        elif thisLine[0] == 'Abundance':
            # We have a new element. Put the last list (if any)
            # into the dictionary, and start clean
            if len(lineAbs) > 0:
                if currentElement in elementDict.keys():
                    elementDict[currentElement].extend(lineAbs)
                else:
                    elementDict[currentElement] = lineAbs

            currentElement = float(EL.symbolLU(thisLine[4])[0]) + (EL.ionLU(thisLine[5])-1.)/10.
            lineAbs = []
        elif u.is_number(thisLine[0]):
            # Note: MOOG2014 added the atom.ion ID as the second column in each line
            # So, now the input looks like:
            # [Wavelength, ID, Ex.Pot., logGf, eqw, logRW, abund, delta Avg.]
            #
            # If the abundance (or delta avg) field is 999.990, it is MOOG's
            # way of denoting one of the components of a blended line. Since 
            # the abundance calculation for this line will be listed elsewhere,
            # we can skip this line in the .out2 file
            if float(thisLine[6]) > 900.:
                continue
            
            # Abundance line, but do a touch of editing...
            lineAbs.append([float(thisLine[0]), float(thisLine[2]), float(thisLine[3]), float(thisLine[4]), float(thisLine[5]),float(thisLine[6])])

            # The following is the code for pre-2014 versions of MOOG:
            # lineAbs.append([float(x) for x in thisLine[:-1]])
            
    # Were we working on an element when the file ended?
    if len(lineAbs) > 0:
        if currentElement in elementDict.keys():
            elementDict[currentElement].extend(lineAbs)
        else:
            elementDict[currentElement] = lineAbs
   
    return elementDict


def GetMOOGAbs(logName=k.MOOGTempLogName, modelName=k.MOOGTempModelName):
# From above, the returned dictionary is keyed by element/ion, with values of 
# a list of measured lines, and corresponding abundances of the form:
# [Wavelength, Ex.Pot., logGf, eqw, logRW, abund]
# Example:   lineDict = {26.0:[[4347.23, 0.00, -5.503, 89.8, -4.68, 8.95]...]}

    parFile = BuildMOOGParFile(parFilename=modelName[:-1]+'par', logFilename=logName, modelName=modelName)
    scriptFile = MakeMOOGScriptFile(parFile, scriptFilename=modelName[:-1]+'scr')

    MOOGStr = dirs.MOOGLoc + ' < ' + scriptFile
    
    result = sub.call(MOOGStr, shell=True)
    
    lineDict = ParseMOOGOut2File(modelName[:-1] + 'out2')
    
    # Clean up the MOOG files:
    if os.path.isfile(scriptFile): os.remove(scriptFile)
    if os.path.isfile(modelName[:-1] + 'par'): os.remove(modelName[:-1] + 'par')
    if os.path.isfile(modelName[:-1] + 'out1'): os.remove(modelName[:-1] + 'out1')
    if os.path.isfile(modelName[:-1] + 'out2'): os.remove(modelName[:-1] + 'out2')
    
    return lineDict


def GetMOOGBlendedAbs(logName=k.MOOGTempLogName, modelName=k.MOOGTempModelName, ion=0.0):
# Build the .par file and run a MOOG 'blends' driver. Technically, the returned
# dictionary contains only the one entry, but we still return it as a
# dictionary, to conform to the same format as the "GetMOOGAbs" function
# above.

    parFile = BuildMOOGParFile(parFilename=modelName[:-1]+'par', logFilename=logName, modelName=modelName, blendsDriver=True, ion=ion)
    scriptFile = MakeMOOGScriptFile(parFile, scriptFilename=modelName[:-1]+'scr')

    MOOGStr = dirs.MOOGLoc + ' < ' + scriptFile
    
    result = sub.call(MOOGStr, shell=True)
    
    lineDict = ParseMOOGOut2File(modelName[:-1] + 'out2')
    
    # Clean up the .scr, .par, .out1, & .out2 files:
    if os.path.isfile(scriptFile): os.remove(scriptFile)
    if os.path.isfile(modelName[:-1] + 'par'): os.remove(modelName[:-1] + 'par')
    if os.path.isfile(modelName[:-1] + 'out1'): os.remove(modelName[:-1] + 'out1')
    if os.path.isfile(modelName[:-1] + 'out2'): os.remove(modelName[:-1] + 'out2')
    
    return lineDict

def WriteMOOGAtmModel(theModel, outfile, teff, logg, met, vel=2.0):
# Writes the model as a MOOG-readable model file
    outHeader = '''KURUCZ
#Kss72: T={0:5.0f},[g]={1:4.2f},[Fe/H]={2:4.2f},vt={3:3.1f}E+05
NTAU            {4:2d}\n'''.format(teff,logg,met,vel,len(theModel))

    outTail = '''     {0:3.1f}E+05
NATOMS           1  {1:4.2f}
      3.00     3.30
NMOL             4
     607.0     108.0     106.0     107.0'''.format(vel, met)
     
    with open(outfile, 'w') as fp:
        # Had some trouble with trying to write too fast...
        time.sleep(1.)
        fp.write(outHeader)
        time.sleep(1.)       
             
        for dataline in theModel:
            fp.write(' {0:10.8E}  {1:7.1f} {2:5.3E} {3:5.3E} {4:5.3E} {5:5.3E} {6:5.3E}\n'.format(dataline[0], dataline[1], dataline[2], dataline[3], dataline[4], dataline[5], dataline[6]))
        time.sleep(1.)            
        fp.write(outTail)


def ReadMOOGLineFile(lineFilename):
# Read through a MOOG .log or .lis file and extract the line data

    allLines = []
    with open(lineFilename, 'r') as lineFile:
        data = lineFile.readlines()
            
    if len(data) == 0:
    # There's at least one line file, with no entries.
        return allLines
        
    for dataLine in data:
        vals = dataLine.split()

        if len(vals) == 0 or vals[0][0] == '#':
        # Blank or comment line:
            continue

        if len(vals) == 8:
        # Data line contains vdw/c6, molecular dissociation energy, eqw, and reference
            allLines.append([float(vals[0]), float(vals[1]), float(vals[2]), float(vals[3]), float(vals[4]), float(vals[5]), vals[7]])
        elif len(vals) > 4:
        # Were missing 1-4 of the optional parameters. We have to assume a
        # "formatted read" type of file, and go by number position:
            try:
                vdwVal = float(dataLine[43:55])
            except ValueError:
                vdwVal = 0.
                
            try:
                deVal = float(dataLine[56:63])
            except ValueError:
                deVal = 0.
            # Note: the eqw will show up inbetween columns 63 and 70, but we don't care...
            if len(dataLine) > 70:
                comment = dataLine[71:]
                
            allLines.append([float(vals[0]), float(vals[1]), float(vals[2]), float(vals[3]), vdwVal, deVal, comment])
        else:
        # Missing all the optional fields
            allLines.append([float(vals[0]), float(vals[1]), float(vals[2]), float(vals[3]), 0., 0., ''])
    
    return allLines
