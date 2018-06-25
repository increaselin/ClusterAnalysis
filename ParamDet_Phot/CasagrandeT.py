#! /usr/bin/env python
# 
# Program: CasagrandeT.py
#
# Author: Mike Lum
#
# Usage: ./CasagrandeT.py [-vh]
#
# Description: Using the relations from Casagrande et al. 2010, this function 
#    returns the Teff values for all of the passed photometry.
#
# Revision History:
#    Date        Vers.    Author        Description
#    10/22/2015  1.0f1     Lum        Converted to module format
#    Today        0.0a0    Lum        First checked in
#
# To Do:
#    
#

from constants import constants as k
from constants import Photometry as kPhot
from ParamDet_Phot import RamMelT as RMT

# Functions
def CasagrandeT(magDict, metal, extinct):
    passedMags = {}
    deltas = RMT.getMagDeltas(magDict, extinct)
    Teffs = {}
    sigmas = {}
    
    goodKeys = [x for x in deltas.keys() if deltas[x] is not k.NoCoefficientData and RMT.isInRange(x, metal, deltas[x], kPhot.CasagrandeRanges)]
    goodPhot = {k:v for k,v in deltas.iteritems() if k in goodKeys and k in kPhot.CasagrandeTable4.keys()}
    for thisKey in goodPhot.keys() :
        params = kPhot.CasagrandeTable4[thisKey]
        sigmas[thisKey] = params['sigma']
        Teffs[thisKey] = 5040/(params['a0']+params['a1']*goodPhot[thisKey]+params['a2']*(goodPhot[thisKey]**2)+params['a3']*goodPhot[thisKey]*metal+params['a4']*metal+params['a5']*(metal**2))
    return Teffs, sigmas

def HuangFGKT(magDict, metal, extinct):
    passedMags = {}
    deltas = RMT.getDeltasForTable(magDict, kPhot.HuangFGKDwarfTable3, extinct=extinct)
    Teffs = {}
    
    for thisKey in deltas.keys() :
        params = kPhot.HuangFGKDwarfTable3[thisKey]
        Teffs[thisKey] = 5040/(params['a0']+params['a1']*deltas[thisKey]+params['a2']*(deltas[thisKey]**2)+params['a3']*deltas[thisKey]*metal+params['a4']*metal+params['a5']*(metal**2))
    return Teffs

def HuangGiantT(magDict, metal, extinct):
    passedMags = {}
    deltas = RMT.getDeltasForTable(magDict, kPhot.HuangGiantTable4, extinct=extinct)
    Teffs = {}
    
    for thisKey in deltas.keys() :
        params = kPhot.HuangGiantTable4[thisKey]
        Teffs[thisKey] = 5040/(params['a0']+params['a1']*deltas[thisKey]+params['a2']*(deltas[thisKey]**2)+params['a3']*deltas[thisKey]*metal+params['a4']*metal+params['a5']*(metal**2))
    return Teffs


# The following lines were from a previous incarnation as a command-line call
'''
from collections import deque
import os

thisScriptName = 'CasagrandeT.py'
interpreterName = 'python'

# Globals
verboseMode = False
 
def printHelpText():
    print('Program: CasagrandeT.py\n')
    print('This function calculates the effective surface temperature of a')
    print('star, using the relationship from Cassagrande, et al. 2010\n')
    print('Usage: ./CasagrandeT.py <photometry> -M x.xx [-E x.xx] [-hv]\n')
    print('Options:')
    print('General Options:')
    print('\t-h: Print this help text.')
    print('\t-v: Use verbose progress and error messages')
    print('Photometry (at least two required):')
    print('\t-B xx.xx: Johnson-Cousins B band magnitude')
    print('\t-V xx.xx: Johnson-Cousins V band magnitude')
    print('\t-R xx.xx: Johnson-Cousins R band magnitude')
    print('\t-I xx.xx: Johnson-Cousins I band magnitude')
    print('\t-J xx.xx: 2MASS J band magnitude')
    print('\t-H xx.xx: 2MASS H band magnitude')
    print('\t-K xx.xx: 2MASS K band magnitude')
    print('\t-BT xx.xx: Tycho B band magnitude')
    print('\t-VT xx.xx: Tycho V band magnitude')
    print('\t-b xx.xx: Stromgren b band magnitude')
    print('\t-y xx.xx: Stromgren y band magnitude')
    print('Stellar parameters:')
    print('\t-M x.xx: (required) Metallicity of the star')
    print('\t-E x.xx: Extinction for a given band')
    
if __name__ == '__main__':

    # Parse your command line call
    temp = os.sys.argv
    argv = deque(temp)
    
    Metallicity = kNoData
    Extinction = 0.0
    
    while len(argv) > 0:
        flag = argv.popleft()
        if flag == '-h':
            printHelpText()
            exit()
        elif flag == '-v':
            verboseMode = True
            showIRAFMessages = 1
            print('Verbose Mode enabled.')
        elif flag == '-M':
            try:
                Metallicity = float(argv.popleft())
            except ValueError:
                Metallicity = kNoData
        elif flag == '-E':
            try:
                Extinction = float(argv.popleft())
            except ValueError:
                Extinction = 0.0
        elif flag in ('-B','-V','-R','-I','-J','-H','-K','-BT','-VT','-b','-y'):
            temp = argv.popleft()
            try: 
                tempMag = float(temp)
            except ValueError:
                continue
            passedMags[flag[1:]] = tempMag
            if verboseMode: print('{0} Magnitude: {1:3.2f}'.format(kFilterDict[flag[1:]], tempMag))
        elif str.find(flag,thisScriptName) != -1 :
            if verboseMode: print('Executing command.')
        elif str.find(flag,interpreterName) != -1 :
            if verboseMode: print('...')
        #else: ...
        # Ignore everything else, or not...
    
    if (Metallicity == kNoData):
        print("No metallicity value passed. Use the -M flag.")
        exit()
        
    temps = CasagrandeT(passedMags, Metallicity, Extinction)
    
    for t in temps.keys():
        print("Teff = {0:5.1f} ({1})".format(temps[t], t))
    
    exit()
'''
