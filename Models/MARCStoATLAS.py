#! /usr/bin/env python
# 
# Program: MARCStoATLAS.py
#
# Author: Mike Lum
#
# Usage: ./MARCStoATLAS.py [-vh] [directory]
#
# Description: Converts all the MARCS model atmosphere files (*.mod) in the
#   passed directory (or current directory, if none passed), into the ATLAS
#   format.
#   The conversion functions are as those in the marcstoatlas.f code, in the 
#   GALA project. (Mucciarelli et al. 2013, ApJ, 766, 78)
#
# Revision History:
#    Date              Vers.    Author        Description
#    2015-04-02        1.0f1    Lum         First release
#    2015-03-18        0.0a0    Lum         First checked in
#
# To Do:
#    
#

from collections import deque
import os
import glob

import numpy as np

thisScriptName = 'MARCStoATLAS.py'
interpreterName = 'python'

# Constants

# Column indices (0-based) of variables in a MARCS model file
k_t_idx = 4
k_Pe_idx = 5    # Note: We currently assume the column order is t, Pe, Pg, and PRad,
k_Pg_idx = 6    #  with no intervening columns
k_PRad_idx = 7
k_xKapR_idx = 2
k_RhoX_idx = 7

# Globals
gVerboseMode = False

# Functions
def printHelpText():
    print('Program: MARCStoATLAS.py\n')
    print('Converts all the MARCS model atmosphere files (*.mod) in the')
    print('passed directory (or current directory, if none passed), into the ATLAS format.')
    print(' The conversion functions are as those in the marcstoatlas.f code, in the')
    print('GALA project. (Mucciarelli et al. 2013, ApJ, 766, 78)\n')
    # Flags and options here
    print('Usage: ./MARCStoATLAS.py [-hv] [directory]\n')
    print('Options:')
    print('-h: Print this help text.')
    print('-v: Use verbose progress and error messages')
    print('(directory): Optionally direct conversion to all files in this directory.')
    print('             Default directory is the current one.')

def readMARCS(thisFile):
    marcsData = {}
    try:
        infile = open(thisFile, 'r')
    except IOError:
        print('Unable to open {0}. Exiting'.format(thisFile))
        exit()
    
    dataLines = infile.readlines()
    infile.close()
    
    teff = round(float(dataLines[1].split()[0]))
    logg = round(np.log10(float(dataLines[3].split()[0])),2)
    vturb = round(float(dataLines[4].split()[0]),2)
    metal = round(float(dataLines[6].split()[0]),2)
    mass = round(float(dataLines[5].split()[0]),1)
    alpha = round(float(dataLines[9].split()[0]), 2)
    # burning (ha ha?) other convection params (nu, y, beta)
    xx, yy, zz = [float(x) for x in (dataLines[10].split()[0:3])]
    yfrac = yy/(4.-3.*yy)
    xfrac = 1. - yfrac      # assume zfrac is negligible
    headerInfo = [teff, logg, vturb, metal, mass, alpha, yfrac, xfrac]
    
    abundArray = np.array([])
    for line in dataLines[12:22]:
        abundArray = np.append(abundArray, [float(x) for x in line.split()])
    # Adjust here for metallicity. Note: H and He are not adjusted.
    abundArray = np.append(abundArray[:2], abundArray[2:]-metal-12.0)
    numLayers = int(dataLines[22].split()[0])
  
    structArray = []
    for line in dataLines[25:25+numLayers]:
        structArray.append([float(x) for x in line.split()[k_t_idx:k_PRad_idx+1]])
    
    for lineNum in range(26+numLayers, 26+2*numLayers):
        thisLine = dataLines[lineNum].split()
        structIdx = lineNum-26-numLayers
        structArray[structIdx].extend([float(thisLine[k_xKapR_idx]), float(thisLine[k_RhoX_idx])])
 
    # Because the acceleration due to radiation (pressure) uses the trapezoidal rule,
    # we have to special case the first and last elements of the data structure, outside
    # of the loop.
    structArray[0].append((structArray[1][3]-structArray[0][3])/(structArray[1][5]-structArray[0][5]))
    structArray[-1].append((structArray[-1][3]-structArray[-2][3])/(structArray[-1][5]-structArray[-2][5]))
    
    for structIdx in range(1, len(structArray)-1):
        structArray[structIdx].append(0.5*((structArray[structIdx][3]-structArray[structIdx-1][3])/(structArray[structIdx][5]-structArray[structIdx-1][5]) + (structArray[structIdx+1][3]-structArray[structIdx][3])/(structArray[structIdx+1][5]-structArray[structIdx][5])))

    return headerInfo, abundArray, structArray
    
def makeOutfileName(hI):
    if hI[3] < 0.:
        metalChar = 'm'
    else:
        metalChar = 'p'
        
    nameStr = 't{0:04d}g{1:02d}{2}{3:03d}k{4:02d}.dat'.format(int(hI[0]), int(hI[1]*10), metalChar, int(abs(hI[3]*100)), int(hI[2]))
    return nameStr

def writeOutfileHeader(outfileName, header):
    try:
        outfile = open(outfileName, 'w')
    except IOError:
        print('Unable to open {0} for output. Quitting'.format(outfileName))
        exit()
    
    outfile.write('TEFF {0:5.0f}.   GRAVITY {1:7.5f} LTE\n'.format(header[0], header[1]))
    outfile.write('TITLE {0} GRID  [{1:+04.2f}]   VTURB {2:03.1f} KM/S    L/H 1.25\n'.format(outfileName, header[3], header[2]))
    outfile.write(' OPACITY IFOP 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 0 0 0 0 0\n')
    outfile.write('CONVECTION ON   1.50 TURBULENCE OFF  0.00  0.00  0.00  0.00\n')
    outfile.write('ABUNDANCE SCALE   {0:7.5f}  ABUNDANCE CHANGE 1 {1:07.5f} 2 {2:07.5f}\n'.format(10**header[3], header[6], header[5]))
    outfile.close()
    return

def writeOutfileAbBlock(outfileName, abundArray):
    try:
        outfile = open(outfileName, 'a')
    except IOError:
        print('Unable to open {0} for output. Quitting'.format(outfileName))
        exit()
    atlasArray = [-20.0 if x <=-20.0 else x for x in abundArray]
    for index in range(2,len(atlasArray),6):
        outfile.write(' ABUNDANCE CHANGE {0:2d} {1:>6.2f} {2:2d} {3:>6.2f} {4:2d} {5:>6.2f} {6:2d} {7:>6.2f} {8:2d} {9:>6.2f} {10:2d} {11:>6.2f}\n'.format(index+1, atlasArray[index], index+2, atlasArray[index+1], index+3, atlasArray[index+2], index+4, atlasArray[index+3], index+5, atlasArray[index+4], index+6, atlasArray[index+5]))
    outfile.write(' ABUNDANCE CHANGE 93 -20.00 94 -20.00 95 -20.00 96 -20.00 97 -20.00 98 -20.00\n')
    outfile.write(' ABUNDANCE CHANGE 99 -20.00\n')
    outfile.close()
    return

def writeOutfileStructBlock(outfileName, struct, vel):
    try:
        outfile = open(outfileName, 'a')
    except IOError:
        print('Unable to open {0} for output. Quitting'.format(outfileName))
        exit()
    outfile.write('READ DECK6 {0:2d} RHOX,T,P,XNE,ABROSS,ACCRAD,VTURB\n'.format(len(struct)))
    for stats in struct:
        outfile.write('{0:>10.8e} {1:>8.1f} {2:>5.3E} {3:>5.3E} {4:>5.3E} {5:>5.3E} {6:>5.3E} 0.000E+00 0.000E+00\n'.format(stats[5], stats[0], stats[2], stats[1]/(1.38054e-16*stats[0]), stats[4], stats[6], vel*1E5))
    outfile.write('PRADK {0:>6.4e}\n'.format(struct[0][3]))
    outfile.write('BEGIN                    ITERATION  15 COMPLETED\n')
    outfile.close()
    return
    
def MARCStoATLAS(directory):
    fileList = glob.glob(directory+'*.mod')
    for thisFile in fileList[:5]:
        header, abund, struct = readMARCS(thisFile)
        outfileName = makeOutfileName(header)
        
        writeOutfileHeader(outfileName, header)
        writeOutfileAbBlock(outfileName, abund)
        writeOutfileStructBlock(outfileName, struct, header[2])
        
        print(outfileName, len(struct))
        
        
if __name__ == '__main__':

    # Parse the command line call
    temp = os.sys.argv
    argv = deque(temp)
    directory = ''
    
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
        elif os.path.isdir(flag):
            directory = flag
        else:
            print('Unrecognized command line string \'{0}\' - Exiting now.'.format(flag))
            exit()
        
    MARCStoATLAS(directory)
    exit()
