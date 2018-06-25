#! /usr/bin/env python
# 
# Module: ReferenceCorrections.py
#
# Author: Mike Lum
#
# Description: functions to "correct" measured lines/abundances for 
# "astrophysically determined" loggf values. Done by comparing the abundance 
# determined by the same line in a reference spectrum (Solar for dwarfs, 
# Arcturus for giants).
#
# Contents: 
#
# Revision History:
#    Date        Vers.    Author        Description
#     1/16/18    1.1f0    Lum        Rolled changes from Output.STRPlots in
#    10/6/16     1.0f0    Lum        First checked in
#
# To Do:
#    
#

# Imports
import numpy as np
from scipy import stats

from Abundances import Abundance as AB
from constants import constants as k
from constants import ptoe_abs as PA
from Databases import LineLookup as LL
from Models import MakeAtmModel as mk
from Output import psPlots as ps
from utilities import elements as el
from utilities import utilities as u

## Hack-y hard-coded limits. We probably don't need these anymore, since we have
## a better "Solar LogGf" algorithm.
#DwarfLimits = {
#            22.0:[-3.0, 99.0, -99, 5.5],
#            22.1:[-5.0, 99.0, -99, 5.5],
#            24.0:[-5.0, 99.0, -99, 6.3],
#            24.1:[-6.0, 99.0, -99, 6.3],
#            25.0:[-5.0, 99.0, -99, 6.3],
#            26.0:[-6.5, 99.0, -99, 8.0],
#            26.1:[-7.0, 99.0, -99, 8.0],
#            27.0:[-4.6, 99.0, -99, 6.0],
#            28.0:[-5.5, 99.0, -99, 7.0],
#            68.1:[-3.0, 2.0, -0.5, 2.0]}
#
#
#GiantLimits = {
#            6.0:[-12., 99.0, -99, 9.5],
#            8.0:[-11.3, 99.0, -99, 9.5],
#            20.0:[-8.0, 99.0, -99, 8.0],
#            22.0:[-4.0, -2, 4.0, 6.0],
#            22.1:[-5.0, 99.0, 4.0, 6.0],
#            24.0:[-5.5, -3.0, 4.0,  6.0],
#            24.1:[-7.0, 99.0, 4.0, 6.2],
#            25.0:[-5.5, 99.0, -99, 6.3],
#            26.0:[-8.0, -5.0, 6.7, 8.3],
#            26.1:[-9.0, 99.0, 6.7, 8.1],
#            27.0:[-6.0, 99.0, -99, 6.3],
#            28.0:[-6.5, -4.0, 5.4, 7.2],
#            68.1:[-3.0, -1.0, -0.5, 2.0]}
#

# Utility Functions
def isGiantStar(parms):
# Star is a giant if the Teff<5200 and LogG<3.5
    return u.isGiant(parms[0], parms[1])


# "SortAndFilterLines is the bread-and-butter function, here. 
def SortAndFilterLines(lines, ion, starParms, filterBlends=False, 
                       solarCorrect=False, solarLines=None, solarCorrs=None, 
                       lineWeights=None, modelAtms=None , pradks=None):
# Filter the passed line list for blended lines. Correct each line measurement
# with an "astrophysically determined" loggf (ie: figure out by how much a given
# line had to be adjusted in a Solar spectrum to yield the known Solar abundance
# for that element (Asplund, et al. 2009). Given a Solar adjustment, split the
# passed line set by how much each line had to be adjusted.

# Formats:
# lines (input):
#       [[Wavelength, Ex.Pot., logGf, eqw, logRW, abund],...]
# adjustedSet (output):
#       [[ab, line STR score, wl, "quality"], ...]
# allSet (output):
#       [[ab, line STR score, wl],...]

    sortedLines = sorted(lines, key = lambda l:l[0])

    isGiant = isGiantStar(starParms)
    
    if isGiant:
        modelPath = k.GiantModelPath
#        hardLimits = GiantLimits
    else:
        modelPath = k.DwarfModelPath
#        hardLimits = DwarfLimits

    if solarCorrect and (solarLines == None or solarCorrs==None or lineWeights==None):
    # We need all of them, otherwise, look them up again.
        if modelAtms == None or pradks == None:
        # No models passed, so look 'em up on our own
            modelFiles = mk.findDataFiles(modelPath)
            modelAtms, pradks = mk.LoadModels(modelFiles)

        if isGiant:
            modelAtms, pradks = mk.LoadModels(mk.findDataFiles(k.DwarfModelPath))
            solarCorrsD, solarLinesD, lineWeightsD = \
                    GetDwarfCorrections(ionList=[ion], 
                                        modelAtms=modelAtms, 
                                        pradks=pradks)
#                    GetGiantCorrections(ionList=[ion], 
#                                        modelAtms=modelAtms, 
#                                        pradks=pradks)
        else:
            solarCorrsD, solarLinesD, lineWeightsD = \
                    GetDwarfCorrections(ionList=[ion], 
                                        modelAtms=modelAtms, 
                                        pradks=pradks)
        if ion in solarLinesD.keys() and len(solarLinesD[ion]) >0:
            solarCorrs = solarCorrsD[ion]
            solarLines = solarLinesD[ion]
            lineWeights = lineWeightsD[ion]

    if len(solarLines) == 0:
    # No corrections found for this ion, turn 'em off
        solarCorrect = False
        solarCorrs = None
        solarLines = None
        lineWeights = None
    # We needed the solar values for filtering blends. Since we don't have
    # it, fake it with the variance of our star's equivalent lines
        solarStdev = np.std([l[5] for l in sortedLines])
    else:
        solarStdev = np.std([l[5] for l in solarLines])
    # Set a minimum for line variance
    if solarStdev < 0.03:
        solarStdev = 0.03
    
    # Throw out lines which are below our detection limit, or above the 
    # linear CoG limit (+ a small margin)
    loSlope, loInt, rv, pv, err = \
        ps.GetDetectionLimit(tuple(starParms),
                         ion, 
                         modelAtms=modelAtms, 
                         pradks=pradks)
    hiSlope, hiInt, rv, pv, err = \
        ps.GetCoGLimit(tuple(starParms), 
                       ion, 
                       modelAtms=modelAtms, 
                       pradks=pradks)
    if not(loSlope==0 and loInt==0) and not(hiSlope==0 and hiInt==0):
    # Our method for calculating limits doesn't work for HFS lines. In which
    # case, we will have zeroes returned from the detection and cog limit
    # calculation functions.
        if ion>56.1:
        # We're going to relax our minimum detection for high At. Nos:
            tempLines = [line for line in sortedLines if \
                hiSlope*line[5]+hiInt > el.STR(line[2], line[1], starParms[0])]
        else:
            tempLines = [line for line in sortedLines if \
                loSlope*line[5]+loInt < el.STR(line[2], line[1], starParms[0]) and \
                hiSlope*line[5]+hiInt+0.5 > el.STR(line[2], line[1], starParms[0])]
        sortedLines = tempLines
    
#    # Exercise our "Hard" limits:
#    if ion in hardLimits.keys():
#        STRLimitLo = hardLimits[ion][0]
#        STRLimitHi = hardLimits[ion][1]
#        AbLimitLo = hardLimits[ion][2]
#        AbLimitHi = hardLimits[ion][3]
#        tempLines = [line for line in sortedLines if line[5] < AbLimitHi and\
#                     line[5] > AbLimitLo and\
#                     el.STR(line[2], line[1], starParms[0]) > STRLimitLo and \
#                     el.STR(line[2], line[1], starParms[0]) < STRLimitHi]
#        sortedLines = tempLines

    if solarCorrect and solarLines is not None:
    # and correctDict is not None and ion in correctDict.keys():
    # - this is assumed from above.
        # Note: We will ONLY return lines which have Solar corrections!
        correctedLines = []
        # Don't look up the keys every loop step
        correctedWLs = solarCorrs.keys()
        
        # See Lum & Boesgaard (2018) for discussion as to why these quality
        # limits are selected.
        nLines = len(sortedLines)
        if nLines <= k.abQRange[0]:
            qList = [w for w in k.AbWeights] # Use a generator, since we're popping
        elif nLines >= k.abQRange[-1]:
            qList = [w for w in k.AbWeights[2:]]
        else:
            qList = [w for w in k.AbWeights[1:]]
        qList.reverse()
        while len(correctedLines) == 0 and len(qList)>0:
            correctedLines = []
            currentQ = qList.pop()
            for line in sortedLines:
                if line[0] in lineWeights.keys() and \
                   lineWeights[line[0]] >= currentQ:
                   correctedLines.append([line[0], line[1], line[2], line[3],\
                                       line[4], line[5]+solarCorrs[line[0]]])
#        # We want to end up with a minimum number of lines. If we don't hit
#        # that limit, then use both corrected and non-corrected lines.
#        if len(correctedLines)<min(len(solarLines)/2,5):
        # If we didn't get any lines, then use the uncorrected ones
        if len(correctedLines) == 0:
            for line in sortedLines:
                correctedWLs = [l[2] for l in correctedLines]
                if line[2] not in correctedWLs:
                    correctedLines.append(line)
            correctedLines = sorted(correctedLines, key=lambda l:l[2])
    else:
        correctedLines = sortedLines

    filteredLines=[]
    if filterBlends:
    # Filtered blends are defined as any (corrected) line which differs from
    # the "standard" by 0.75sigma (wide) or 0.5sigma (narrow)
        try:
            abSlope, abIntercept, rv, pv, err = \
                    stats.linregress(np.array(correctedLines)[:,1], 
                                     np.array(correctedLines)[:,5])
        except IndexError:
            pass

        for line in correctedLines:
            expectedAb = abIntercept + abSlope*line[1]
            if (line[5]-expectedAb) > k.MaxBlendSTDev*solarStdev or \
               (line[5]-solarStdev) < -k.MinBlendSTDev*solarStdev:
                continue
            else:
                filteredLines.append(line)
    else:
        filteredLines = correctedLines

    if len(lineWeights)>0:
        weightedLines = lineWeights.keys()
        
        # Lines with weights have had corrections made. The amount of correction
        # is represented by the weight, and can act as a measure of the 
        # "quality" of the line - Lines with lower adjustments are more
        # reliable indicators of abundance.
        adjustedSet = np.array([[line[5], 
                              el.STR(line[2], line[1], starParms[0]), line[0],
                              lineWeights[line[0]]] 
                                for line in \
                                    sorted(filteredLines, key=lambda l:l[0]) \
                                        if line[0] in weightedLines])
    else:
        adjustedSet = []

    # Raw points
    allSet = np.array([[line[5], 
                        el.STR(line[2], line[1], starParms[0]), line[0]] 
                            for line in 
                                sorted(sortedLines, key=lambda l:l[0])])
    return adjustedSet, allSet


# Functions to provide "Astrophysical" log(gf) corrections, based on
# measurement of the same line in a star with "known" abundances - The Sun
# for dwarf/main-sequence stars, and Arcturus for giant stars.
def GetDwarfCorrections(ionList=[], modelAtms=None,
                        pradks=None, useDAOSpec=False):
    return GetSolarCorrections(ionList=ionList, modelAtms=modelAtms,
                               pradks=pradks, useDAOSpec=useDAOSpec)


def GetSolarCorrections(ionList=[], modelAtms=None,
                        pradks=None, useDAOSpec=False):
# Returns a dictionary of lines (for the passed ion list) with their
# abundance corrections, based on our line measurement of Solar spectra,
# and the accepted abundance (from Asplund 2009).
# If the passed ion list is empty, all corrected lines are returned.
#
# The three returned dictionaries:
# solarCorrs:{ion:{wl:correction,...},...}
# solarLines:{ion:[[wl, Ex.Pot., logGf, eqw, logRW, abund],...]}
# weights:   {ion:{wl:weight,...},...}
    solarCorrs = {}
    weights = {}

    # First, check to see if we already have calculated Solar abs
    allLines = LL.LookupSolarAbundances(ionList=ionList)
    if len(allLines.keys()) > 0:
        solarLines = allLines
    else:
        # Don't have them. So, calculate, and enter them. We have the line list
        # now.
        if modelAtms is None or pradks is None:
            modelPath = k.DwarfModelPath
            modelFiles = mk.findDataFiles(modelPath)
            modelAtms, pradks = mk.LoadModels(modelFiles) 

        allLines = LL.getSolarLines(ionList, [(4500.,8750.)],\
                                    useDAOSpec=useDAOSpec)
        # We need a slightly different form if we have to calculate the abs.
        compositeLines = sorted([[l[0], l[1], 'filler', 'filler', l[2]] \
                            for l in allLines], key = lambda l:(l[1], l[0]))

        unused1, solarLines, unused2, unused3 = AB.CalcAbsFromLines(\
                        compositeLines, k.SolarRefParms, star='SolarRef',\
                        ionList=ionList, modelAtms=modelAtms, pradks=pradks)
        # MOOG can fail, in which case, our dictionary will be empty:
        if len(solarLines.keys()) == 0:
            return {}, {}, {}
        # Now, put them into our DB, so we don't have to do this again.
        newLines = []
        for ion in sorted(solarLines.keys()):
            newLines.extend([[k.SolarRefStarName, ion, l[0], l[1], l[2], l[3], l[4], l[5]] for l in solarLines[ion]])
        LL.EnterCalibrationAbundances(newLines)
        
    if len(ionList) == 0:
        ionList = solarLines.keys()

    for ion in ionList:
        solarAb = PA.ptoe[int(ion)-1][PA.abIdx]
        solarCorrs[ion] = {}
        weights[ion] = {}
        if ion not in solarLines.keys():
            continue
        for line in solarLines[ion]:
            correction = solarAb-line[5]
            solarCorrs[ion][line[0]] = correction
            weights[ion][line[0]] = GetRefLineWeight(correction)

    return solarCorrs, solarLines, weights


def GetRefLineWeight(correction):
# Returns a weight, based on the ranges and values in constants.py, for
# the passed correction (in dex), to be used for overall abundance calcs.
# The lesser the correction value, the greater the recommended weight of
# the line.
    checkValue = abs(correction)
    # Take advantage of numpy array's boolean mask as a lookup table
    try:
        returnValue = np.array(k.AbWeights)[np.array([checkValue >=i[0] and \
                               checkValue<i[1] for i in k.AbWeightRanges])][0]
    except IndexError:
        returnValue = k.NoAbWeight

    return returnValue


def GetGiantCorrections(ionList=[], modelAtms=None, 
                        pradks=None, useDAOSpec=False):
# Returns a dictionary of lines (for the passed ion list) with their
# abundance corrections, based on our line measurement of Arcturus spectra,
# and the accepted abundance (from Ramirez and Allende Prieto, 2011).
# If the passed ion list is empty, all corrected lines are returned.

    # For now, the Arcturus corrections....suck, so just do Solar corrections.
    # Not the best solution, but hopefully better than nothing.
    return GetSolarCorrections(ionList=ionList)

#    if modelAtms is None or pradks is None:
#        modelPath = k.GiantModelPath
#        modelFiles = mk.findDataFiles(modelPath)
#        modelAtms, pradks = mk.LoadModels(modelFiles) 
#    allLines = LL.getCalibrationLines(starName=k.ArcturusRefStarName, 
#                            elements=ionList, wavelengthRange=[(4500.,8750.)])
#
#    # Each entry looks like: [Wl, ion, median EQW, delta EQW]
#    # Unfortunately, line details are not stored in the calibration table,
#    # so we need to do a second DB lookup for those, s we'll need to "fake"
#    # that these lines were read in the same manner as "normal" spectral lines.
#    compositeLines = sorted([[l[0], l[1], 'filler', 'filler', l[2]] for l in allLines], key = lambda l:(l[1], l[0]))
#    
#    unused1, giantLines, unused2, unused3 = \
#                         AB.CalcAbsFromLines(compositeLines,\
#                         k.ArcturusRefParams, star='GiantRef',\
#                         ionList=ionList, modelAtms=modelAtms, pradks=pradks)
#
#    giantCorrs = {}
#    weights = {}
#
#    if len(ionList) == 0:
#        ionList = giantLines.keys()
#
#    for ion in ionList:
#        arctAb = PA.ptoe[int(ion)-1][PA.giantAbIdx]
#        giantCorrs[ion] = {}
#        weights[ion] = {}
#        if ion not in giantLines.keys():
#            continue
#        for line in giantLines[ion]:
#            correction = arctAb-line[5]
#            giantCorrs[ion][line[0]] = correction
#            weights[ion][line[0]] = GetRefLineWeight(correction)
#
#    return giantCorrs, giantLines, weights
    
