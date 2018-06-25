#! /usr/bin/env python
# 
# Module: Abundance.py
#
# Author: Mike Lum
#
# Description: Functions to calculate elemental abundances from 
#   spectroscopic data.
#
# Contents:
# CorrectNLTEAbs(ion, lineList, parmTuple) :
#       Perform NLTE corrections for the passed lines.
#   Parameters: 
#       ion - ion to correct:
#           format: (float) atno.ionState
#       lineList - List of lines to correct
#           format:  List:
#                    [[Wavelength, Ex.Pot., logGf, eqw, logRW, abund],...]
#       parmTuple - tuple with the star's atmospheric parameters:
#           format: (Teff, LogG, VTurb, [Fe/H], Mass)
#   Returns: List of the lines (same format as the passed lineList), with
#           the corrections applied.
#
#
# CalcAbs(starName, parmTuple, ionList = [], modelAtms = None, 
#         pradks = None, useDAOlines=False):
#       Calculate the abundances for the passed ions in the passed star.
#   Parameters:
#       starName: Text string containing the star's ID and cluster name.
#           format: "CCCCCC SSSSSS", where CCCC is the cluster name 
#               (any length), and SSSS is the star designation within the
#                cluster. The two strings are separated by a space.
#       parmTuple - tuple with the star's atmospheric parameters:
#           format: (Teff, LogG, VTurb, [Fe/H], Mass)
#       ionList - list with the ions to check.
#           format: [(float) atno.ionState, ...]
#       modelAtms: All of the atmospheric models for this type of star. Note:
#           This is done as an optimization courtesy - Loading the models
#           takes time, so if the caller has a lot of stars, they can opt to
#           pre-load the models.
#       pradks: All of the radiative pressure constants for the passed
#           models. Note: these can be loaded in Models/MakeAtmModel
#       useDAOlines: Boolean telling whether to use DAOSpec to measure
#           the line EQWs, or to use the pre-measured values in our 
#           database.
#   Returns: Dictionary of the form: {ion X:[[X/H], stdev, #lines]...}
#
#
# CalcAbsAndLines(starName, starParms, ionList = [], modelAtms = None, 
#       pradks = None, calcAbRange=False, useDAOlines=False):
#       As above, but also returns the full list of all measured lines,
#       and a dictionary of the minimum detectable abundance for a given line
#       based on varying the measured equivalent widths by the instrumental
#       "sigma" (1/(S/N ratio))*(spectrograph dispersion). The maximum
#       abundance dictionary contains the maximum abundance that will result
#       from the widest EQW that is still within the linear portion of the 
#       curve of growth.
#   Parameters:
#       starName: Text string containing the star's ID and cluster name.
#           format: "CCCCCC SSSSSS", where CCCC is the cluster name 
#               (any length), and SSSS is the star designation within the
#                cluster. The two strings are separated by a space.
#       starParms - tuple with the star's atmospheric parameters:
#           format: (Teff, LogG, VTurb, [Fe/H], Mass)
#       ionList - list with the ions to check.
#           format: [(float) atno.ionState, ...]
#       modelAtms: All of the atmospheric models for this type of star. Note:
#           This is done as an optimization courtesy - Loading the models
#           takes time, so if the caller has a lot of stars, they can opt to
#           pre-load the models.
#       pradks: All of the radiative pressure constants for the passed
#           models. Note: these can be loaded in Models/MakeAtmModel
#       calcAbRange: Boolean to return minimum and maximum ab limits.
#       useDAOlines: Boolean telling whether to use DAOSpec to measure
#           the line EQWs, or to use the pre-measured values in our 
#           database.
#   Returns: Abundance dictionary, line dictionary, 
#               minimum ab dict, max. ab dict:
#       format(s):
#           Abundance dictionary: 
#              {elem.ion:[abundance mean, ab. std.dev., # lines]}
#           Line dictionary:
#              {elem.ion:[[Wavelength, Ex.Pot., logGf, eqw, logRW, abund],...]}
#           Minimum abundance dictionary (only if calcAbRange is True):
#               Format: Same as Line dictionary, above
#           Maximum abundance dictionary (only if calcAbRange is True):
#               Format: Same as Line dictionary, above
#
#
# CalcAbsFromLines(allLines,  parmTuple, star='', ionList = [], 
#       modelAtms = None, pradks = None, calcAbRange=False):
#       Given the passed list of lines/eqw measures and stellar parameters,
#       return up to four abundance dictionaries, as above.
#   Parameters:
#       allLines - List of lines to use to calculate this star's abundances
#           format: List:
#               [[ion, wl, (unused), (unused), EQW]...]
#       parmTuple - tuple with the star's atmospheric parameters:
#           format: (Teff, LogG, VTurb, [Fe/H], Mass)
#       star - Text string representing the star's name. Only used to create
#           the temporary files. Suggest using a unique identifier (like a 
#           timestamp) in the name, so that there are no conflicts in 
#           multi-threaded calls.
#       ionList - List of ions to measure abundances. Empty for all. 
#       modelAtms: All of the atmospheric models for this type of star. Note:
#           This is done as an optimization courtesy - Loading the models
#           takes time, so if the caller has a lot of stars, they can opt to
#           pre-load the models.
#       pradks: All of the radiative pressure constants for the passed
#           models. Note: these can be loaded in Models/MakeAtmModel
#       calcAbRange: Boolean to return minimum and maximum ab limits.
#   Returns: Abundance dictionary, line dictionary, 
#               minimum ab dict, max. ab dict:
#       format(s):
#           Abundance dictionary: 
#              {elem.ion:[abundance mean, ab. std.dev., # lines]}
#           Line dictionary:
#              {elem.ion:[[Wavelength, Ex.Pot., logGf, eqw, logRW, abund],...]}
#           Minimum abundance dictionary (only if calcAbRange is True):
#               Format: Same as Line dictionary, above
#           Maximum abundance dictionary (only if calcAbRange is True):
#               Format: Same as Line dictionary, above
#
#
# Revision History:
#    Date        Vers.    Author        Description
#    Today        0.0a0    Lum        First checked in
#    8-1-2017     1.0f1    Lum        Added NLTE correction code
#    4-16-2018    1.1f1    Lum        New documentation and code cleanup
#
# To Do:
#    
#

# Imports
import datetime
import os

import numpy as np

from constants import constants as k
from constants import DBFormats as DB
from Databases import LineLookup as LL
from MOOGInterface import MOOGInterface as MI
from Models import MakeAtmModel as mk
from NLTE_O_Correction import getNLTE_O_EQW as NLTEO
from utilities import utilities as u

# Functions
def CorrectNLTEAbs(ion, lineList, parmTuple):
# Perform NLTE corrections for the passed lines.
# Passed lines are of the format:
# [[Wavelength, Ex.Pot., logGf, eqw, logRW, abund],...]
# The passed list is returned with the corrections (if any) applied to the
# abundance.
    (Teff, LogG, VTurb, Met, Mass) = parmTuple
    correctedList = []
    lines = lineList
    # A list with no items...really!?! OK, fine. BE that way!
    if len(lineList) == 0:
        return correctedList
        
    # A line list with one element can be problematic:
    if not u.is_list(lines[0]):
        lines = []
        lines.append(lineList)
    
    # OI triplet (7771,7774,7775) - Takeda 2003
    if ion == 8.0:
        for line in lines:
            if int(line[0]) in [7771, 7774, 7775]:
                correction = NLTEO.getNLTE_O_EQW(int(line[0]), Teff, LogG, VTurb, line[3])
                correctedList.append([line[0], line[1], line[2], line[3], line[4], line[5]+correction])
            else:
                correctedList.append(line)
    return correctedList
        

def CalcAbs(starName, parmTuple, ionList = [], modelAtms = None, pradks = None, useDAOlines=False):
# Assumes starName is of the form: "CCCCCC SSSSSS", where CCCC is the cluster
# name (any length), and SSSS is the star designation within the cluster. The
# two parameters are separated by a space.
# Returns a dictionary of the form: {ion X:[[X/H], stdev, #lines]...}
    (Teff, LogG, VTurb, Met, Mass) = parmTuple  # Why did I do this?
    abDict, unused1, unused2, unused3 = CalcAbsAndLines(starName, \
                        (Teff, LogG, VTurb, Met, Mass), ionList,\
                        modelAtms=modelAtms, pradks=pradks,\
                        useDAOlines=useDAOlines)
    return abDict


def CalcAbsAndLines(starName, starParms, ionList = [], modelAtms = None, pradks = None, calcAbRange=False, useDAOlines=False):
# As above, but also returns the full list of all measured lines
# Loading the atmospheric models takes time. So, the caller may
# optionally load (and pass) the models, in order to save time
# when calling this function for multiple stars.
#
# Data is returned as two dictionaries, detailed in CalcAbsFromLines, below.

    cluster, star = starName.split()
    # Lines from DB
    allLines = LL.getLinesForStar(cluster, star, elements=ionList, wavelengthRange=[4500, 8700], daoLinesOnly=useDAOlines)
    # getLinesForStar returns the mean of all measurements when multiple measurements
    # of a single line exist from multiple spectra.
    # "allLines" is sorted by ion and wavelength. Each line returned is of the
    # form:
    # [ion, DB wl, filename/"composite", meas WL, eqw, fwhm]
    filteredLines = [line for line in allLines if np.log10(line[4]/line[1])<k.LinearCOGLimit+3. and abs(line[3]-line[1])<k.LambdaVarianceLimit]

    return CalcAbsFromLines(filteredLines, starParms, star=star, ionList=ionList, modelAtms=modelAtms, pradks=pradks, calcAbRange=calcAbRange)


def CalcAbsFromLines(allLines,  parmTuple, star='', ionList = [], modelAtms = None, pradks = None, calcAbRange=False):
# Given the passed list of lines/eqw measures and stellar parameters, return
# two abundance dictionaries. The keys for each are element.ion
# combos (ie: 26.1 for FeII). The "abDict" abundance dictionary has entries in
# the form: {elem.ion:[abundance mean, ab. std.dev., # lines]}
# The line list has the form:
# {elem.ion:[[Wavelength, Ex.Pot., logGf, eqw, logRW, abund],...]}
#
# If the passed parameter "calcAbRange" is True, we'll also return
# minimum and maximum abundance arrays, based on varying the measured equivalent
# widths by the instrumental "sigma" (1/(S/N ratio))*(spectrograph dispersion).
    (Teff, LogG, VTurb, Met, Mass) = parmTuple
    
    if modelAtms is None or pradks is None:        
        if u.isGiant(Teff,LogG):
            modelPath = k.GiantModelPath
        else:
            modelPath = k.DwarfModelPath
        modelFiles = mk.findDataFiles(modelPath)
        modelAtms, pradks = mk.LoadModels(modelFiles)

    # If the passed ion list is empty, it's actually a notation to "measure all"
    if len(ionList) == 0:
        ionList = set([x[0] for x in allLines])
    # If a particular ion contains with HFS, we will need to run the 'blends'
    # driver for that entire element. We're building a list of elements with 
    # blending, and elements without.
    HFSIons = []
    unBlendedIons = []
    for ion in ionList:
        linesToCheck = [l for l in allLines if l[0] == ion]
        if LL.containsHFS(linesToCheck):
            HFSIons.append(ion)
        else:
            unBlendedIons.append(ion)
    
    # Create the MOOG-readable atmosphere model - good for both blended and unblended
    # ions.
    atmModel, pradk = mk.MakeAtmModel(modelAtms, pradks, Teff, LogG, Met, VTurb, Mass)
    modelFile = '{0}_{1:06d}.m'.format(star, datetime.datetime.now().microsecond)
    MI.WriteMOOGAtmModel(atmModel, modelFile, Teff, LogG, Met, VTurb)
    
    # Do our unblended lines, first, and all at once.
    unBlendedLines = \
            sorted([line for line in allLines if line[0] in unBlendedIons],\
                   key=lambda l: (l[0],l[1]))
    
    eqwLogFile = '{0}_{1:06d}.log'.format(star, datetime.datetime.now().microsecond)
    
    # Put the lines into MOOG-readable format. Note, the returned logFile name
    # is the same as the one passed.
    lineLogFile = MI.MakeMOOGEQWLogFile(unBlendedLines, logFilename=eqwLogFile)
    
    # Calc. abs for the measured lines
    MOOGLines = MI.GetMOOGAbs(eqwLogFile, modelFile)
    if len(MOOGLines.keys()) == 0 and len(HFSIons)==0:
    # MOOG failed, and we have nothing to work with, so break back out, and
    # hopefully whoever called us is able to deal with blank dictionaries.
        # Clean up the log and model files
        u.clearFiles([eqwLogFile, modelFile])
        return {}, MOOGLines, {}, {}
    
    if calcAbRange:
        # Make two copies of the line list for the min/max arrays
        minLines = [[l[0], l[1], l[2], l[3], \
            l[4]-k.EQWMeasSigma if l[4]>k.EQWMeasSigma else 0., \
            l[5]] for l in unBlendedLines]
        maxLines = [[l[0], l[1], l[2], l[3], \
            l[4]+k.EQWMeasSigma, l[5]] for l in unBlendedLines]
        minLogFile = 'Min_'+eqwLogFile
        MI.MakeMOOGEQWLogFile(minLines, logFilename=minLogFile)
        minAbs = MI.GetMOOGAbs(minLogFile, modelFile)
        u.clearFile(minLogFile)
        
        maxLogFile = minLogFile = 'Max_'+eqwLogFile
        MI.MakeMOOGEQWLogFile(maxLines, logFilename=maxLogFile)
        maxAbs = MI.GetMOOGAbs(maxLogFile, modelFile)
        u.clearFile(maxLogFile)
    else:
        minAbs = {}
        maxAbs = {}
    # Do our blended ions individually (since we need to call the blends
    # driver for each set).
    for ion in HFSIons:
        # Lose the old log file, used for unblended/non-HFS lines:
        u.clearFile(eqwLogFile)
        
        eqwLogFile = '{0}_{1:06d}.log'.format(star, datetime.datetime.now().microsecond)
        if calcAbRange:
            minLogFile = 'Min_'+eqwLogFile
            maxLogFile = 'Max_'+eqwLogFile
            
        try:
            logFile = open(eqwLogFile,'w')
            logFile.write('# Temporary MOOG log file created: {0}\n'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        except IOError:
            u.clearFile(eqwLogFile)
            break # Next ion/return
        try:
            if calcAbRange:
                minLog = open(minLogFile, 'w')
                minLog.write('# Temporary MOOG log file created: {0}\n'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
                maxLog = open(maxLogFile, 'w')
                maxLog.write('# Temporary MOOG log file created: {0}\n'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        except IOError:
            # There's a question here. Should we continue with calculating the
            # abundances, even if we can't return a range? For now, let's just fail
            # and address more robust error handling at a later date.
            u.clearFiles([minLogFile, maxLogFile])
#            calcAbRange = False
            break # Next ion/return
                        
        for line in [l for l in allLines if l[0] == ion]:
        # For loop to write out the lines for this element to the eqwLogFile
            # Need the ExPot and the vdw damping factor
            lineInfo, trash = LL.getLines(elements=[ion], wavelengthRange=[(line[1]-0.05, line[1]+0.05)], dataFormat=None)
            if len(lineInfo) == 0:
                # Couldn't get line info from the DB -
                # Skip this line
                continue
                
            blends = LL.GetHFSBlends(line[1], ion)
            if len(blends) > 0:
            # Got a series of blends
                # First one should be positive (WL), and contain the eqw for the entire blend
                logFile.write(k.MOOGListFormat.\
                    format(blends[0][0], ion, lineInfo[0][2], blends[0][1],lineInfo[0][4],line[4],'\n'))
                if calcAbRange:
                    minLog.write(k.MOOGListFormat.\
                    format(blends[0][0], ion, lineInfo[0][2], blends[0][1],lineInfo[0][4],\
                    line[4]-k.EQWMeasSigma if line[4]-k.EQWMeasSigma>0. else 0.,'\n'))
                    maxLog.write(k.MOOGListFormat.\
                    format(blends[0][0], ion, lineInfo[0][2], blends[0][1],lineInfo[0][4],line[4]+k.EQWMeasSigma,'\n'))
                # Subsequent ones should be negative, and the eqw is ignored (0.0)
                for blend in blends[1:]:
                    logFile.write(k.MOOGListFormat.format(-blend[0], ion, \
                        lineInfo[0][2], blend[1],lineInfo[0][4], 0.0,'\n'))
                    if calcAbRange:
                        minLog.write(k.MOOGListFormat.format(-blend[0], ion, \
                            lineInfo[0][2], blend[1],lineInfo[0][4], 0.0,'\n'))
                        maxLog.write(k.MOOGListFormat.format(-blend[0], ion, \
                            lineInfo[0][2], blend[1],lineInfo[0][4], 0.0,'\n'))
                    
            else:
            # Just the one line
                logFile.write(k.MOOGListFormat.\
                    format(line[1],line[0],lineInfo[0][2],lineInfo[0][3],\
                    lineInfo[0][4],line[4],'\n'))
                if calcAbRange:
                    minLog.write(k.MOOGListFormat.\
                        format(line[1],line[0],lineInfo[0][2],lineInfo[0][3],\
                        lineInfo[0][4],\
                        line[4]-k.EQWMeasSigma if line[4]-k.EQWMeasSigma>0. else 0.,'\n'))
                    maxLog.write(k.MOOGListFormat.\
                        format(line[1],line[0],lineInfo[0][2],lineInfo[0][3],\
                        lineInfo[0][4],line[4]+k.EQWMeasSigma,'\n'))
                    
        logFile.close()
        if calcAbRange:
            minLog.close()
            maxLog.close()
            
        # Calc. abs for the measured lines
        tempAbDict = MI.GetMOOGBlendedAbs(eqwLogFile, modelFile, ion)
        # MOOG can fail, and return an empty dictionary - This isn't a big deal
        # here, but we should handle that case, just in case.
        if len(tempAbDict.keys()) == 0:
            continue
        MOOGLines.update(tempAbDict)
        
        u.clearFile(eqwLogFile)

        if calcAbRange:
            minLog.close()
            minAbs.update(MI.GetMOOGBlendedAbs(minLogFile, modelFile, ion))
            maxLog.close()
            maxAbs.update(MI.GetMOOGBlendedAbs(maxLogFile, modelFile, ion))
            u.clearFiles([minLogFile, maxLogFile])

    u.clearFiles([eqwLogFile, modelFile])
    
    abDict = {}
    # Calc the abs:
    for ion in MOOGLines.keys():
        lineAbs = [x[5] for x in MOOGLines[ion]]
        abDict[ion] = (np.mean(lineAbs), np.std(lineAbs), len(lineAbs))
    
    # Clean up
    
    return abDict, MOOGLines, minAbs, maxAbs
