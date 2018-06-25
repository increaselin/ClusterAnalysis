#! /usr/bin/env python
# 
# Program: MakeAtmModel.py
#
# Author: Mike Lum
#
# API Entry: MakeModelFile(Teff, LogG, Metalicity, V_turb, outfileName='', modelPath='')
#
# Description: Library functions to interpolate between the Teff, LogG, and
# metallicity points of either the Kurucz "ATLAS9" models, as produced in 2011,
# or the 2013 MARCS models (either spherical or plane-parallel)
#
# Usage: Pretty self-explanitory, I hope. Pass your desired parameters, the output
# file name (optional*), and the directory where the ATLAS or MARCS grids are 
# located.
# *: If no output file name is passed, the output will be placed in the local directory
#    with the name: 
#                   XXXXX_TTTT_GGG_SMM_VV.m
#                   XXXXX = "MARCS" or "ATLAS" depending on the passed models
#                   TTTT = Teff
#                   GGG = LogG * 100 (integer)
#                   SMM = Metalicity * 10 (integer) with the sign as 'm' or 'p'
#                   VV = V turb * 10 (integer)
#
# Revision History:
#    Date       Vers.    Author        Description
#   2015-04-02  1.0f1    Lum        First release
#   2015-02-20  0.0a0    Lum        First checked in
#
# To Do:
#    
#

import glob
import numpy as np
import time

from utilities import utilities as u
from MOOGInterface import MOOGInterface as MI
from constants import constants as k
from Models import MARCStoATLAS as m2a


def parms_in_range(teff=-1., logg=-1., vturb=-1.):
# Returns True for teff, logg, and vturb in the acceptable
# range for our models.
    return ((teff == -1.) or u.in_range(teff, k.TEffRange[0], k.TEffRange[1])) and\
           ((logg == -1.) or u.in_range(logg, k.LogGRange[0], k.LogGRange[1])) and\
           ((vturb == -1.) or u.in_range(vturb, k.VturRange[0], k.VturRange[1]))
           
###################
# Loading functions
###################

def findDataFiles(searchPath):
    # Returns a list of all potential model files from a passed
    # directory. We assume that the directory is correctly populated
    # with the desired model files.
    # Assumes: 
    #           *.mod = MARCS models            *.dat = ATLAS models
    # if we think both are present in a directory, we opt for the larger
    # number of _files_. Note that this favors MARCS models, with a single
    # model per file.
    # We should do some more rigorous checking (ex: validating the
    # first line, using the expected format for the expected model). However,
    # we'll assume that the user is smart enough to verify that the model files 
    # are in the passed directory, before calling us.
    MARCSFileList = glob.glob(searchPath+'/*.mod')
    ATLASFileList = glob.glob(searchPath+'/*.dat')
    
    if len(MARCSFileList) >= len(ATLASFileList):
        return MARCSFileList
    else: return ATLASFileList


def loadATLASModels(filename):
# The current format (as of: Jan. 2015) for the ATLAS models subdivides
# models into files, based on metallicity and velocity. There is only
# one metallicity/velocity combination per file. We assume this, but leave the
# possibility of multiple met./vel. per file, open. 
    try:
        model = open(filename,'r')
        lines = model.readlines()
    except IOError:
        print('Error in reading: {0}'.format(filename))
        exit()
    model.close()

    modelStartLines = []
    count = 0

    for thisLine in lines:
        theseWords = thisLine.split()
        if len(theseWords) is not 0:
            if theseWords[0][-3:] == 'EFF': 
             # Note: several of the ATLAS model files are missing the first 
             # character of each line, so some files have a first model line
             # starting with "TEFF", while others have "EFF"
                modelStartLines.append(count)
        count+=1
    
    tableDict = {}
    pradkValue = {}
    
    for startLineNo in modelStartLines:
        Tstr = lines[startLineNo].split()[1]
        Gstr = lines[startLineNo].split()[3]
        Mstr = lines[startLineNo+1].split()[3].strip('[]')
        # Special case:
        if Mstr == 'SOLAR': Mstr = '0.0'
        Vstr = lines[startLineNo+1].split()[5]
        # Special case:
        if Vstr == 'WITH': Vstr = '1.5'
        TabLenstr = lines[startLineNo+22].split()[2]
        # Our table lengths need to be the same... otherwise, matrix 
        # interpolation isn't gonna work! For now, we just toss short/long
        # files. In the future, it might be better to try to convert them
        # to the same length.
        if not (u.is_number(Tstr) and u.is_number(Gstr) and u.is_number(Mstr) and u.is_number(Vstr) and u.is_number(TabLenstr)): continue
        if int(TabLenstr) != 72: continue
        if len(lines[startLineNo+24].split()) == 7:
        # Some ATLAS models omit the last two columns when they are all zeroes
            tableDict[(float(Tstr),float(Gstr))] = [[float(t) for t in s.split()]+[0.,0.] for s in lines[startLineNo+23:startLineNo+23+int(TabLenstr)]]
        else:
            tableDict[(float(Tstr),float(Gstr))] = [[float(t) for t in s.split()] for s in lines[startLineNo+23:startLineNo+23+int(TabLenstr)]]
        pradkValue[(float(Tstr),float(Gstr))] = float(lines[startLineNo+23+int(TabLenstr)].split()[1])

    theKeys = tableDict.keys()
    # Just re-use the met. & vel. values read from the final table
    met = float(Mstr)
    vel = float(Vstr)
    
    return met, vel, tableDict, pradkValue
    
def LoadModels(modelFileList):
# Assume that MARCS models end with '.mod' (and ATLAS with '.dat')
    MARCSModels = (modelFileList[0][-4:] == '.mod')

    kuruczAtms = {}
    pradks = {}

    if MARCSModels:
        # if False: print('Loading MARCS models.') 
        for model in modelFileList:
            header, abund, struct = m2a.readMARCS(model)
            # header = [teff, logg, vturb, metal, mass, alpha, yfrac, xfrac]
            # pradk = struct[0][3]
            if header[4] not in kuruczAtms.keys():
                kuruczAtms[header[4]] = {}
                pradks[header[4]] = {}
            if header[2] not in kuruczAtms[header[4]].keys():
                kuruczAtms[header[4]][header[2]] = {}
                pradks[header[4]][header[2]] = {}
            if header[3] not in kuruczAtms[header[4]][header[2]].keys():
                kuruczAtms[header[4]][header[2]][header[3]] = {}
                pradks[header[4]][header[2]][header[3]] = {}
            kuruczAtms[header[4]][header[2]][header[3]][(header[0], header[1])] = [[x[5], x[0], x[2], x[1]/(1.38054e-16*x[0]), x[4], x[6], header[2]*1E5, 0.0, 0.0] for x in struct]
            pradks[header[4]][header[2]][header[3]][(header[0], header[1])] = struct[0][3]
       
    else:
        solarAtm = {}
        solarPradk = {}
        # if False: print('Loading ATLAS models.') 
        for model in modelFileList:
            met, vel, modelDict, pradk = loadATLASModels(model)
            if vel not in solarAtm.keys():
                solarAtm[vel] = {}
                solarPradk[vel] = {}
            if met not in solarAtm[vel].keys():
                solarAtm[vel][met] = {}
                solarPradk[vel][met] = {}
            solarAtm[vel][met].update(modelDict)
            solarPradk[vel][met].update(pradk)
        # ATLAS models only have one mass: Solar 
        # However, we return the models as having both 1.0 and 0.0, due to 
        #   'weirdness' in the MARCS models (or programmer laziness...)
        kuruczAtms[0.0] = solarAtm
        kuruczAtms[1.0] = solarAtm
        pradks[0.0]= solarPradk
        pradks[1.0] = solarPradk
    return kuruczAtms, pradks



#########################
# Interpolation functions
#########################


def interpForTEff(kuruczAtms, pradks, TEff, LogG):
# Perform a linear interpolation between the closest two TEff models
    TeffLogGChoices = kuruczAtms.keys()
    TeffChoices = sorted(set([tLogGTuple[0] for tLogGTuple in TeffLogGChoices if tLogGTuple[1] == LogG]))
    loTeffIdx, hiTeffIdx = u.bracket(TEff,TeffChoices)
    loTeff = TeffChoices[loTeffIdx]
    hiTeff = TeffChoices[hiTeffIdx]

    if TEff == loTeff or hiTeff == loTeff:
        # # if False: print('\t\t\tMatched low TEff. {0:5.0f} :'.format(loTeff))
        theModel =  np.array(kuruczAtms[(loTeff, LogG)]), pradks[(loTeff, LogG)]
    elif TEff == hiTeff:
        # # if False: print('\t\t\tMatched high TEff. {0:5.0f} :'.format(hiTeff))
        theModel = np.array(kuruczAtms[(hiTeff, LogG)]), pradks[(hiTeff, LogG)]
    else:
        # # if False: print('\t\t\tLow TEff {0:5.0f} :'.format(loTeff))
        loTEffModel = np.array(kuruczAtms[(loTeff, LogG)])
        lopradk = pradks[(loTeff, LogG)]
        # # if False: print('\t\t\tHigh TEff {0:5.0f} :'.format(hiTeff))
        hiTEffModel = np.array(kuruczAtms[(hiTeff, LogG)])
        hipradk = pradks[(hiTeff, LogG)]
        linInterpFactor = (TEff-loTeff)/(hiTeff - loTeff)
        theModel = ((1-linInterpFactor)*loTEffModel + hiTEffModel*linInterpFactor), (1-linInterpFactor)*lopradk + hipradk*linInterpFactor
    
    return theModel
    

def interpForLogG(kuruczAtms, pradks, TEff, LogG):
# Perform a log interpolation between the closest two LogG models
    TeffLogGChoices = kuruczAtms.keys()
    logGChoices = sorted(set([tLogGTuple[1] for tLogGTuple in TeffLogGChoices]))
    loLogGIdx, hiLogGIdx = u.bracket(LogG,logGChoices)
    loLogG = logGChoices[loLogGIdx]
    hiLogG = logGChoices[hiLogGIdx]
    
    if len(logGChoices) == 1 or LogG == loLogG or hiLogGIdx == 0:
        # This will also handle calls where LogG is below the range
        # # if False: print('\t\tMatched low logG. {0:+4.2f} :'.format(loLogG))
        return interpForTEff(kuruczAtms, pradks, TEff, loLogG)
    elif LogG == hiLogG:
        # # if False: print('\t\tMatched high logG. {0:+4.2f} :'.format(hiLogG))
        return interpForTEff(kuruczAtms, pradks, TEff, hiLogG)
    elif LogG > hiLogG:
        # In the case of LogGs above the available range, we will extrapolate
        # using the highest two available tables. This should be the only case 
        # where the linear interpolation factor is greater than 1.
        # Note that we assume that the limit of how high we will interpolate
        # is set and enforced elsewhere. It's probably not good to
        # interpolate more than 0.5-1.0 above the actual grids' range.
        
        # Included, but not required
        # hiLogGIdx = -1; loLogGIdx = -2
        loLogG = logGChoices[-2]
        hiLogG = logGChoices[-1]
        loLogGModel, lopradk = interpForTEff(kuruczAtms, pradks, TEff, loLogG)
        hiLogGModel, hipradk = interpForTEff(kuruczAtms, pradks, TEff, hiLogG)
        # Note the interpolation is slightly different, as we are above the high value,
        # and we want to use the high value, as opposed to the "low", as our "base".
        linInterpFactor = ((10**LogG-10**hiLogG)/(10**hiLogG-10**loLogG))
        return ((1+linInterpFactor)*hiLogGModel - loLogGModel*linInterpFactor), ((1+linInterpFactor)*hipradk - lopradk*linInterpFactor)
    else:
        # # if False: print('\t\tLow logG {0:+4.2f} :'.format(loLogG))
        loLogGModel, lopradk = interpForTEff(kuruczAtms, pradks, TEff, loLogG)
        # # if False: print('\t\tHigh logG {0:+4.2f} :'.format(hiLogG))
        hiLogGModel, hipradk = interpForTEff(kuruczAtms, pradks, TEff, hiLogG)
        linInterpFactor = u.linlogfraction(loLogG, hiLogG, LogG)
        return ((1-linInterpFactor)*loLogGModel + hiLogGModel*linInterpFactor), ((1-linInterpFactor)*lopradk + hipradk*linInterpFactor)

    
    
def interpForMetal(kuruczAtms, pradks, Metal, TEff, LogG):
# Perform a log interpolation between the closest two metallicity models
    loMetIdx,hiMetIdx = u.bracket(Metal, sorted(kuruczAtms.keys()))
    loMetal = sorted(kuruczAtms.keys())[loMetIdx]
    hiMetal = sorted(kuruczAtms.keys())[hiMetIdx]
    
    if Metal == loMetal or loMetal == hiMetal:
        # if False: print('\tMatched low Met. {0:+4.2f} :'.format(loMetal))
        return interpForLogG(kuruczAtms[loMetal], pradks[loMetal], TEff, LogG)
    elif Metal == hiMetal:
        # if False: print('\tMatched high Met. {0:+4.2f} :'.format(hiMetal))
        return interpForLogG(kuruczAtms[hiMetal], pradks[hiMetal], TEff, LogG)
    else:
        # if False: print('\tLow Met. {0:+4.2f} :'.format(loMetal))
        loMetModel, lopradk = interpForLogG(kuruczAtms[loMetal], pradks[loMetal], TEff, LogG)
        # if False: print('\tHigh Met. {0:+4.2f} :'.format(hiMetal))
        hiMetModel, hipradk = interpForLogG(kuruczAtms[hiMetal], pradks[hiMetal], TEff, LogG)
        linInterpFactor = u.linlogfraction(loMetal, hiMetal, Metal)
        return ((1-linInterpFactor)*loMetModel + linInterpFactor*hiMetModel), ((1-linInterpFactor)*lopradk + linInterpFactor*hipradk)



def interpForVelocity(kuruczAtms, pradks, TEff, LogG, Metal, VTurb):
# This functions finds the set of 'bracketing models for mass, and interpolates
# between them.
    # Note that we only have ATLAS models for VTurb = 2.0...
    loVIdx,hiVIdx = u.bracket(VTurb, sorted(kuruczAtms.keys()))
    loV = sorted(kuruczAtms.keys())[loVIdx]
    hiV = sorted(kuruczAtms.keys())[hiVIdx]

    if VTurb == loV or loV == hiV:
        # if False: print('\tMatched low VTurb {0:+4.2f} :'.format(loV))
        return interpForMetal(kuruczAtms[loV], pradks[loV], Metal, TEff, LogG)
    elif VTurb == hiV:
        # if False: print('\tMatched high VTurb {0:+4.2f} :'.format(hiV))
        return interpForMetal(kuruczAtms[hiV], pradks[hiV], Metal, TEff, LogG)
    else:
        # if False: print('\tLow VTurb {0:+4.2f} :'.format(loV))
        loVModel, lopradk = interpForMetal(kuruczAtms[loV], pradks[loV], Metal, TEff, LogG)
        # if False: print('\tHigh VTurb {0:+4.2f} :'.format(hiV))
        hiVModel, hipradk = interpForMetal(kuruczAtms[hiV], pradks[hiV], Metal, TEff, LogG)
        linInterpFactor = (VTurb-loV)/(hiV - loV)
        return ((1-linInterpFactor)*loVModel + linInterpFactor*hiVModel), ((1-linInterpFactor)*lopradk + linInterpFactor*hipradk)
 
def interpForMass(kuruczAtms, pradks, TEff, LogG, Metal, VTurb, mass):
# The atmospheric models we are passed are 5-D (MARCS Spherical), 4-D
# (MARCS Plane Parallel), or 3-D (ATLAS) grids with Mass, (micro-) turbulent
# velocity, metallicity, Teff, and LogG as axes. 4-D models omit mass, and 3-D
# models also omit turbulent velocity.
# This functions finds the set of 'bracketing models for mass, and interpolates
# between them.
    # Note: 4-D MARCS, and 3-D ATLAS models assume mass=1.0, but index it as 0.0
    loMassIdx,hiMassIdx = u.bracket(mass, sorted(kuruczAtms.keys()))
    loMass = sorted(kuruczAtms.keys())[loMassIdx]
    hiMass = sorted(kuruczAtms.keys())[hiMassIdx]

    if mass == loMass or loMass == hiMass:
        # if False: print('\tMatched low mass {0:+4.2f} :'.format(loMass))
        return interpForVelocity(kuruczAtms[loMass], pradks[loMass],TEff, LogG, Metal, VTurb)
    elif mass == hiMass:
        # if False: print('\tMatched high mass {0:+4.2f} :'.format(hiMass))
        return interpForVelocity(kuruczAtms[hiMass], pradks[hiMass], TEff, LogG, Metal, VTurb)
    else:
        # if False: print('\tLow mass {0:+4.2f} :'.format(loMass))
        loMassModel, lopradk = interpForVelocity(kuruczAtms[loMass], pradks[loMass], TEff, LogG, Metal, VTurb)
        # if False: print('\tHigh mass {0:+4.2f} :'.format(hiMass))
        hiMassModel, hipradk = interpForVelocity(kuruczAtms[hiMass], pradks[hiMass], TEff, LogG, Metal, VTurb)
        linInterpFactor = (mass-loMass)/(hiMass - loMass)
        return ((1-linInterpFactor)*loMassModel + linInterpFactor*hiMassModel), ((1-linInterpFactor)*lopradk + linInterpFactor*hipradk)

    return interpForVelocity(kuruczAtms[mass], pradks[mass], TEff, LogG, Metal, VTurb)

def MakeAtmModel(kuruczAtms, pradks, TEff, LogG, Metal, VTurb, mass=0.0):
    return interpForMass(kuruczAtms, pradks, TEff, LogG, Metal, VTurb, mass)
 
###################
# Writing functions
###################

def WriteModel(theModel, outfile, pradk, teff, logg, met, vel=2.0):
# Outputs the model in ATLAS format
    outHeader = 'TEFF  {0:5.0f}.  GRAVITY {1:7.5f} LTE\nTITLE SDSC GRID  [{2:3.1f}]   VTURB {3:3.1f} KM/S    L/H 1.25\n'.format(teff, logg, met, vel)+''' OPACITY IFOP 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 0 0 0 0 0
 CONVECTION ON   1.25 TURBULENCE OFF  0.00  0.00  0.00  0.00
ABUNDANCE SCALE   1.00000 ABUNDANCE CHANGE 1 0.91100 2 0.08900
 ABUNDANCE CHANGE  3 -10.88  4 -10.89  5  -9.44  6  -3.48  7  -3.99  8  -3.11
 ABUNDANCE CHANGE  9  -7.48 10  -3.95 11  -5.71 12  -4.46 13  -5.57 14  -4.49
 ABUNDANCE CHANGE 15  -6.59 16  -4.83 17  -6.54 18  -5.48 19  -6.82 20  -5.68
 ABUNDANCE CHANGE 21  -8.94 22  -7.05 23  -8.04 24  -6.37 25  -6.65 26  -4.37
 ABUNDANCE CHANGE 27  -7.12 28  -5.79 29  -7.83 30  -7.44 31  -9.16 32  -8.63
 ABUNDANCE CHANGE 33  -9.67 34  -8.69 35  -9.41 36  -8.81 37  -9.44 38  -9.14
 ABUNDANCE CHANGE 39  -9.80 40  -9.54 41 -10.62 42 -10.12 43 -20.00 44 -10.20
 ABUNDANCE CHANGE 45 -10.92 46 -10.35 47 -11.10 48 -10.18 49 -10.58 50 -10.04
 ABUNDANCE CHANGE 51 -11.04 52  -9.80 53 -10.53 54  -9.81 55 -10.92 56  -9.91
 ABUNDANCE CHANGE 57 -10.82 58 -10.49 59 -11.33 60 -10.54 61 -20.00 62 -11.04
 ABUNDANCE CHANGE 63 -11.53 64 -10.92 65 -11.94 66 -10.94 67 -11.78 68 -11.11
 ABUNDANCE CHANGE 69 -12.04 70 -10.96 71 -11.28 72 -11.16 73 -11.91 74 -10.93
 ABUNDANCE CHANGE 75 -11.77 76 -10.59 77 -10.69 78 -10.24 79 -11.03 80 -10.95
 ABUNDANCE CHANGE 81 -11.14 82 -10.19 83 -11.33 84 -20.00 85 -20.00 86 -20.00
 ABUNDANCE CHANGE 87 -20.00 88 -20.00 89 -20.00 90 -11.92 91 -20.00 92 -12.51
 ABUNDANCE CHANGE 93 -20.00 94 -20.00 95 -20.00 96 -20.00 97 -20.00 98 -20.00
 ABUNDANCE CHANGE 99 -20.00\n'''+'READ DECK6 {0:d} RHOX,T,P,XNE,ABROSS,ACCRAD,VTURB\n'.format(len(theModel), pradk)
    
    with open(outfile, 'w') as fp:
        fp.write(outHeader)
        for dataline in theModel:
            fp.write(' {0:10.8E}  {1:7.1f} {2:5.3E} {3:5.3E} {4:5.3E} {5:5.3E} {6:5.3E} {7:5.3E} {8:5.3E}\n'.format(dataline[0], dataline[1], dataline[2], dataline[3], dataline[4], dataline[5], dataline[6], dataline[7], dataline[8]))
        fp.write('PRADK {0:6.4E}\nBEGIN                    ITERATION  15 COMPLETED\n'.format(pradk))
        
        
def WriteMOOGModel(theModel, outfile, teff, logg, met, vel=2.0):
    print('Models/MakeAtmModel.py WriteMOOGModel depricated.')
    print('Use: MOOGInterface/MOOGInterface WriteMOOGAtmModel instead')
    return MI.WriteMOOGAtmModel(theModel, outfile, teff, logg, met, vel)

def MakeModelFile(teff, logg, met, vel, outfile='', searchPath=''):
    modelFiles = findDataFiles(searchPath)
    kuruczAtms, pradks = LoadModels(modelFiles)
    thisModel, pradk = MakeAtmModel(kuruczAtms, pradks, teff, logg, met, vel)

    if outfile == '':
        modelType = 'ATLAS'
        if modelFiles[0][-4:] == '.mod':
            modelType = 'MARCS'
            
        metalChar = 'm'
        if met >= 0.0:
            metalChar = 'p'
            
        outfile = '{5}_{0:4.0f}_{1:03d}_{2}{3:02d}_{4:02d}.m'.format(int(teff), int(logg*100.), metalChar, int(abs(met*10.)), int(vel*10.), modelType)
        
    MI.WriteMOOGAtmModel(thisModel, outfile, teff, logg, met, vel)



# Note: this was originally tested using a command line call, 
# now that we're live, the following code is superfluous, but left
# in for reference.
# Also note that the syntax of most of the calls have changed, so
# this code actually won't work as-is.

# Usage: ./MakeAtmModel.py -t <Temperature> -g <LogG> -m <metallicity> [-w <output File Name>][-v <microt. vel.>][-Vh]
#

#from collections import deque
#import os

#thisScriptName = 'MakeAtmModel.py'
#interpreterName = 'python'

#def printHelpText():
#    print('Program: MakeAtmModel.py\n')
#    print('Insert your \"help text\" here.\n')
#    # Be sure to list your full list of flags and options here
#    print('./MakeAtmModel.py -t <Temperature> -g <LogG> -m <metallicity> [-w <output File Name>][-v <microt. vel.>][-Vh]\n')
#    print('Options:')
#    print('-d: Use atmospheric models located in the passed directory.')
#    print('-g: Use the passed LogG value.')
#    print('-h: Print this help text.')
#    print('-t: Use the passed Teff.')
#    print('-v: Use the passed microturbulent velocity. Note, the model grid only has values of 2.0')
#    print('-V: Use verbose progress and error messages')
#    print('-w: Write the output to the specified file. Default: aMMTTTTGG.m')


#if __name__ == '__main__':
#
#    # Parse your command line call
#    temp = os.sys.argv
#    argv = deque(temp)
#    
#    searchPath = 'ATLASModels'
#    outfile = ''
#    teff = 5000.0
#    logg = 4.2
#    met = 0.0
#    vel = 2.0
#    
#    while len(argv) > 0:
#        flag = argv.popleft()
#        if flag == '-h':
#            printHelpText()
#            exit()
#        elif flag == '-V':
#            g.VerboseMode = True
#            print('Verbose Mode enabled.')
#        elif flag == '-d':
#            searchPath = argv.popleft()
#        elif flag == '-w':
#            outfile = argv.popleft()
#        elif flag == '-t':
#            tempStr = argv.popleft()
#            if u.is_number(tempStr):
#                teff = float(tempStr)
#            else:
#                # if False: print('Unable to convert {0} to a valid temperature. Exiting.'.format3(tempStr))
#                exit()
#        elif flag == '-g':
#            loggStr = argv.popleft()
#            if u.is_number(loggStr):
#                logg = float(loggStr)
#            else:
#                # if False: print('Unable to convert {0} to a valid gravity. Exiting.'.format#(loggStr))
#                exit()
#        elif flag == '-m':
#            metStr = argv.popleft()
#            if u.is_number(metStr):
#                met = float(metStr)
#            else:
#                # if False: print('Unable to convert {0} to a valid metallicity. Exiting.'.format(metStr))
#                exit()
#        elif flag == '-v':
#            vturbStr = argv.popleft()
#            if u.is_number(vturbStr):
#                vel = float(vturbStr)
#            else:
#                # if False: print('Unable to convert {0} to a valid velocity. Exiting.'.format(vturbStr))
#                exit()
#        elif str.find(flag,thisScriptName) != -1 :
#            # if False: print('Executing command.')
#        elif str.find(flag,interpreterName) != -1 :
#            # if False: print('...')
#        else:
#            # if False: print('Unrecognized flag: {0}. Exiting.'.format(flag))
#            exit()
#    if outfile == '':
#        outfile = 'a{0:02d}{1:04d}{2:02d}.m'.format(int(abs(met)*10), int(teff), int(logg*10))
#        
#    modelFiles = findDataFiles(searchPath)
#    
#    kuruczAtms, pradks = LoadModels(modelFiles)
#    
#    thisModel, pradk = MakeAtmModel(kuruczAtms, pradks, teff, logg, met, vel)
#    
#    WriteModel(thisModel, outfile, pradk, teff, logg, met, vel)
#    exit()
