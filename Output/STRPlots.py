#! /usr/bin/env python
#
# Module: STRPlots.py
#
# Author: Mike Lum
#
# Description: Functions to create a series of abundance-line strength
#              plots to help with line selection and abundance
#              calculations.
#
# Contents:
#   Function xxx: (Description)
#
# Revision History:
#    Date       Vers.  Author   Description
#    2017-10-26 1.1f1  Lum      Functions to plot Ab to atmospheric
#                               params. 
#                               Indentations now adhere to Python
#                               Programming Style Guide:
#                             https://www.python.org/dev/peps/pep-0008/
#    2016-06-09 1.0f1  Lum      First checked in
#
# To Do:
#       Variable and function names do not conform to style guidelines
#       Functions need to be generalized for other clusters
#

# Standard Library Imports
import os

# Astro- and 3rd-party Library Imports
import numpy as np
from matplotlib import pyplot
from scipy import stats

# Local Imports
from Abundances import Abundance as AB
from Abundances import ReferenceCorrections as RC
from constants import constants as k
from constants import directories as dirs
from constants import ptoe_abs as PA
from Databases import StarDB as SD
from Models import MakeAtmModel as mk
from Output import psPlots as ps
from utilities import elements as el
from utilities import utilities as u

# Constants
kLaTexHeader1 = \
'''\\documentclass{article}

\\usepackage{pdflscape}
\\usepackage{geometry}
\\usepackage{fancyhdr}
\\setlength{\\headheight}{12.0pt}
\\pagestyle{fancyplain}
\\fancyhf{}
'''

# Physical line and measurement info, which should move to a DB access
NLTEIons = [8.0]
# Format for Synth directories:
# ion:{star:[(region, ab, sigma, # lines (weight)),...]...}
SynthAbsDict = {
6.0 :{
            'PLA-0350':[(7115.0, 8.28, 0.20, 5), (7442.0, 8.28, 0.20, 3)],
            'PLA-0356':[(7115.0, 8.23, 0.20, 5), (7442.0, 8.33, 0.20, 3)],
#            'PLA-0506':[(7115.0, 8.13, 0.20, 5), (7442.0, 8.18, 0.20, 3)],
            'PLA-0687':[(7115.0, 8.13, 0.20, 5), (7442.0, 8.18, 0.20, 3)],
            'PLA-1089':[(7115.0, 8.23, 0.20, 5), (7442.0, 8.28, 0.20, 3)],
            'PLA-1172':[(7115.0, 8.33, 0.20, 5), (7442.0, 8.23, 0.20, 3)],
            'PLA-0391':[(7115.0, 8.23, 0.20, 5)],
#            'PLA-0413':[(7115.0, 8.53, 0.20, 5)],
            'PLA-0429':[(7115.0, 8.38, 0.20, 5)],
            'PLA-0520':[(7115.0, 8.03, 0.20, 5)],
#            'PLA-0552':[(7115.0, 8.43, 0.20, 5)],
            'PLA-0575':[(7115.0, 8.43, 0.20, 5)],
            'PLA-0699':[(7115.0, 8.38, 0.20, 5)],
            'PLA-0701':[(7115.0, 8.23, 0.20, 5)],
            'PLA-0786':[(7115.0, 8.43, 0.20, 5)],
            'PLA-0790':[(7115.0, 8.53, 0.20, 5)],
            'PLA-0791':[(7115.0, 8.33, 0.20, 5)],
            'PLA-0828':[(7115.0, 8.43, 0.20, 5)],
#            'PLA-0859':[(7115.0, 8.23, 0.20, 5)],
            'PLA-0864':[(7115.0, 8.28, 0.20, 5)],
            'PLA-0889':[(7115.0, 8.13, 0.20, 5)],
            'PLA-0921':[(7115.0, 8.42, 0.20, 5)],
            'PLA-0964':[(7115.0, 8.23, 0.20, 5)],
            'PLA-0993':[(7115.0, 8.28, 0.20, 5)],
            'PLA-0999':[(7115.0, 8.23, 0.20, 5)],
            'PLA-1012':[(7115.0, 8.23, 0.20, 5)],
            'PLA-1017':[(7115.0, 8.23, 0.20, 5)],
            'PLA-1107':[(7115.0, 8.28, 0.20, 5)],
            'PLA-1270':[(7115.0, 8.23, 0.20, 5)]
            },
7.0 :{
            'PLA-0350':[(7115.0, 7.98, 0.20, 3), (7442.0, 8.12, 0.20, 5)],
            'PLA-0356':[(7115.0, 8.03, 0.20, 3), (7442.0, 8.12, 0.20, 5)],
#            'PLA-0506':[(7115.0, 8.23, 0.20, 3), (7442.0, 8.17, 0.20, 5)],
            'PLA-0687':[(7115.0, 8.23, 0.20, 3), (7442.0, 8.23, 0.20, 5)],
            'PLA-1089':[(7115.0, 7.73, 0.20, 3), (7442.0, 8.12, 0.20, 5)],
            'PLA-1172':[(7115.0, 7.93, 0.20, 3), (7442.0, 8.12, 0.20, 5)]
            }
}

#kAllIonList = [
#               6.0, 7.0, 8.0, 11.0, 12.0, 13.0, 14.0, 14.1, 16.0, 20.0,
#               21.0, 21.1, 22.0, 22.1, 23.0, 23.1, 24.0, 24.1, 25.0, 26.0,
#               26.1, 27.0, 27.1, 28.0, 29.0, 30.0, 38.1, 39.0, 39.1, 40.0,
#               40.1, 56.1, 58.1, 59.1, 60.1, 62.1, 63.1, 64.1, 66.1, 68.1]
kAllIonList = [\
            6.0, 7.0, 8.0, 11.0, 12.0, 13.0, 14.0, 14.1, 20.0, 21.0, 21.1,\
            22.0, 22.1, 23.0, 23.1, 24.0, 24.1, 25.0, 26.0, 26.1, 27.0, 28.0,\
            29.0, 30.0, 39.1, 40.0, 40.1, 56.1, 58.1, 60.1, 62.1]
            
kGoodIonList = [\
            6.0, 7.0, 8.0, 11.0, 12.0, 13.0, 14.0, 20.0, 21.1,\
            22.0, 22.1, 23.0, 23.1, 24.0, 24.1, 25.0, 26.0, 26.1, 27.0, 28.0,\
            29.0, 30.0, 39.1, 40.0, 40.1, 56.1, 58.1, 60.1, 62.1]

# Obtain physical data for stars
def GetAllStarParms(clusterName='NGC-0752', dMass=0.0, gMass=1.6, 
                    metallicity=0.0):
# Uses the database lookup function to get a star's data, and create the 
# format we want:
    allParms = GetDwarfStarParms\
                (clusterName=clusterName, mass=dMass, metallicity=metallicity)
    allParms.extend(GetGiantStarParms\
                (clusterName=clusterName, mass=gMass, metallicity=metallicity))
    return allParms
    
    
def GetDwarfStarParms(clusterName='NGC-0752', mass=0.0, metallicity=-0.03):
    dwarfInfo = SD.GetDwarfInfoForCluster(clusterName=clusterName)
    # Returned Data format:
    # [RA, DEC, Epoch, clusterName, starID, Teff, LogG, vTurb, Vmag, Bmag, s2n]
    # We want:
    # [(starID, Teff, LogG, Vturb, metallicity, mass, s2n)...]
    retData = [(s[4], s[5], s[6], s[7], metallicity, mass, s[10])\
        for s in dwarfInfo]
    return retData
    
    
def GetGiantStarParms(clusterName='NGC-0752', mass=1.7, metallicity=-0.03):
    giantInfo = SD.GetGiantInfoForCluster(clusterName=clusterName)
    # Returned Data format:
    # [RA, DEC, Epoch, clusterName, starID, Teff, LogG, vTurb, Vmag, Bmag, s2n]
    # We want:
    # [(starID, Teff, LogG, Vturb, metallicity, mass, s2n)...]
    retData = [(s[4], s[5], s[6], s[7], metallicity, mass, s[10])\
        for s in giantInfo]
    return retData
    
def GetStarFeH(starAbs):
# Get the Fe/H for this star.
# If both ions are measured, weight 75/25 for FeI/FeII
    # Calculate Fe vs. solar:
    # Note: Fe is atno 26, lut is 0-based, so index = 26-1
    solarFeH = PA.ptoe[25][PA.abIdx]
    starIons = sorted(starAbs.keys())
    if 26.0 in starIons and 26.1 in starIons:
        return 0.75*starAbs[26.0][0]\
                 + 0.25*starAbs[26.1][0]\
                 - solarFeH
    elif 26.0 in starIons:
        return starAbs[26.0][0]-solarFeH
    else:
        return k.BadAbundanceValue


def GetAbsForStar(starData, clusterName='NGC-0752', ions=kAllIonList, \
                  filterBlends=False, refCorrect=False,\
                  refLines=None, refDict=None, lineWeights=None, 
                  useDAOSpec=False, modelAtms=None, pradks=None):
# Get all the elemental abundances for the single passed star (params):
# starData: (starName, Teff, LogG, VTurb, Met, Mass)
#
# Returns a dictionary of {ion:(Ab, std, #lines, avg Quality)...}
    starName = starData[0]
    starParms = tuple(starData[1:6])

    abdict, uncorrLines, unusedMins, unusedMaxes = \
            AB.CalcAbsAndLines(clusterName+' '+starName, 
                               starParms, 
                               ionList=ions, 
                               modelAtms=modelAtms, 
                               pradks=pradks, 
                               useDAOlines=useDAOSpec)
    # uncorrLines:
    # # {elem.ion:[[Wavelength, Ex.Pot., logGf, eqw, logRW, abund],...]}
    
    starAbDict = {}
    ionLUT = sorted(list(abdict.keys()))
    ionLUT = sorted(set(ionLUT + [ion for ion in SynthAbsDict.keys()\
                        if starName in SynthAbsDict[ion].keys()]))
    if refCorrect and (refLines is None or refDict is None):
        refCorrect=False
         
    for ion in ionLUT:
        # We'll make one line per ion. Note, synth elements won't 
        # necessarily have lines in the table/dict.
        if ion not in list(uncorrLines.keys()) + list(SynthAbsDict.keys()):
            continue
        if refCorrect and ion not in refLines.keys():
            correctionsAvail = False
        else:
            correctionsAvail = True

        # Does this element have NLTE corrections available? Note: we do
        # this BEFORE Solar corrections, which assumes that the same
        # NLTE corrections are already applied to any Solar corrections
        # we use.
        if ion in NLTEIons:
            LTELines = AB.CorrectNLTEAbs(ion, uncorrLines[ion], starParms)
        else:
            if ion in uncorrLines.keys():
                LTELines = uncorrLines[ion]
            else:
                # Either synthesized lines, or none available
                LTELines = np.array([])
                
        # Now apply "Solar" or "Astrophysical Log gf" corrections, if requested
        if correctionsAvail and (refCorrect or filterBlends) and ion in uncorrLines.keys():
            tempAdj,tempAll = \
                    RC.SortAndFilterLines(LTELines, 
                                       ion, 
                                       starParms, 
                                       filterBlends=filterBlends, 
                                       solarCorrect=refCorrect, 
                                       solarLines=refLines[ion], 
                                       solarCorrs=refDict[ion], 
                                       lineWeights=lineWeights[ion])

            # correctedLines:
            #           [[ab, line STR score, wl, "quality"], ...]
            # allLines:
            #           [[ab, line STR score, wl],...]
            # We want to use np.arrays, so...
            allLines = np.array(tempAll)
            correctedLines = np.array(tempAdj)
                
        else:
            if len(LTELines)>0:
                allLines = np.array(\
                    [[l[5], el.STR(l[2], l[1], starParms[0]), l[0]] \
                                                for l in LTELines])
                correctedLines = np.array(\
                    [[l[0], l[1], l[2], k.AbWeights[-1]] \
                                                for l in allLines])
            else:
                allLines = correctedLines = np.array([])
                
        if ion in SynthAbsDict.keys() and \
             starName in SynthAbsDict[ion].keys():
            # We have a C I or N I synthesis:
            if ion == 6.0:
                nonSynthLines = np.array([l for l in correctedLines\
                               if l[2]<7110 or l[2]>7120])
            elif ion==7.0:
                # N I synthesis uses both the 7115 region (giants)
                # and the 7442 region (all). Only the 7442 region 
                # has a potentially measured line, tho.
                nonSynthLines = np.array([l for l in correctedLines \
                               if np.floor(l[2]) != 7442])
            else:
                nonSynthLines=[]
            
            synthArray = GetSynthLines(starName, ion,\
                         clusterName=clusterName)
            # Returned tuple: (synthMean, synthStd, synthCount)
            
            # Include our synthesis as the number of lines measured
            if len(nonSynthLines) == 0:
                starAbDict[ion] = [synthArray[0], synthArray[1],\
                                   synthArray[2],k.NoAbWeight]
            else:
                mean = np.mean(nonSynthLines[:,0])
                std = np.std(nonSynthLines[:,0])
                ct = len(nonSynthLines)
                
                totalCt = ct+synthArray[2]
                totalMean = ((mean*ct)+(synthArray[0]*synthArray[2]))/\
                                (totalCt)
                totalStd = np.sqrt(\
                        (ct*(totalMean-mean)**2+
                         synthArray[2]*(synthArray[0]-mean)**2)/\
                        (totalCt-1))
                starAbDict[ion] = [totalMean,totalStd,totalCt,k.NoAbWeight]
            # We have already created the abundance dictionary for this ion,
            # so skip out on the next part.
            continue
        else:
            # No synthesized lines need to be added.
            pass
        
        if refCorrect and len(correctedLines)>0:
            linesToCount = correctedLines
            avgQ = np.mean(correctedLines[:,3])
        else:
            linesToCount = allLines
            avgQ = k.NoAbWeight
        
        if len(linesToCount)>0:
            starAbDict[ion] = [np.mean(linesToCount[:,0]), 
                               np.std(linesToCount[:,0]),
                               len(linesToCount), avgQ]
        else:
            # No lines to count. Make sure we don't have an entry:
            if ion in starAbDict.keys():
                del starAbDict[ion]
    return starAbDict


def GetAbTable(clusterName='NGC-0752', starDataList=None, 
               ions=kAllIonList, filterBlends=False, referenceCorrect=False, 
               useDAOSpec=False, 
               gModelAtms=None, gPradks=None, 
               dModelAtms=None, dPradks=None):
# Returns a dictionary with each star name (key) has a list of elements
# (secondary key), with a list containing [Abundance in [X/H], variance in the
# line measurements, count of lines measured, average quality score]:
# returnDict={starName:{ion:[[X/H],stdev, count, avg.quality],...},...}
    if starDataList==None:
        starDataList = GetAllStarParms(clusterName=clusterName)
    
    # Since we're doing a bunch of lines (uh...), we'll need both giant and
    # dwarf references. Giant references are prefixed by a 'g', dwarfs by 'd'.
    if gModelAtms==None or gPradks==None:
        gModelFiles = mk.findDataFiles(k.GiantModelPath)
        gModelAtms, gPradks = mk.LoadModels(gModelFiles)

    if dModelAtms==None or dPradks==None:
        dModelFiles = mk.findDataFiles(k.DwarfModelPath)
        dModelAtms, dPradks = mk.LoadModels(dModelFiles)
        
        
    # We need iron to do relative abundances:
    if 26.0 not in ions:
        ions.append(26.0)
        
# Turns out our giant corrections are not good
#    gCorrectDict, gReferenceLines, gLineWeights = \
#           GetGiantCorrections(ionList=ions, 
#                               modelAtms=gModelAtms, 
#                               pradks=gPradks)
# So, use the Solar corrections for now
    gCorrectDict, gReferenceLines, gLineWeights = \
            RC.GetSolarCorrections(ionList=ions, 
                                modelAtms=dModelAtms, 
                                pradks=dPradks)


    dCorrectDict, dReferenceLines, dLineWeights = \
            RC.GetDwarfCorrections(ionList=ions, 
                                modelAtms=dModelAtms, 
                                pradks=dPradks)

    tableDict = {}
    
    for starData in starDataList:
        starName = starData[0]
        starParms = tuple(starData[1:6])
        if RC.isGiantStar(starParms):
            modelAtms = gModelAtms
            pradks = gPradks
            referenceLines = gReferenceLines
            referenceDict = gCorrectDict
            lineWeights = gLineWeights
        else:
            modelAtms = dModelAtms
            pradks = dPradks
            referenceLines = dReferenceLines
            referenceDict = dCorrectDict
            lineWeights = dLineWeights
            
        tableDict[starName] = GetAbsForStar(starData, ions=ions,\
                filterBlends=filterBlends, refCorrect=referenceCorrect,\
                refDict=referenceDict, refLines=referenceLines, \
                lineWeights=lineWeights, useDAOSpec=useDAOSpec, \
                modelAtms=modelAtms, pradks=pradks)

    return tableDict


def GetSynthLines(starID, ion, clusterName='NGC-0752'):
# Looks up a synthesis result for the passed star and element. If one exists,
# return a tuple: (synth. Ab result, synth. Error, num lines (weight))
    if ion in SynthAbsDict.keys() and starID in SynthAbsDict[ion].keys():
        starSynths = SynthAbsDict[ion][starID]
    else:
        return ()
    
    accumAbs = []
    synthCount = 0
    accumVar = 0.0
    for synth in starSynths:
        accumAbs.extend([synth[1]]*synth[3])
        synthCount += synth[3]
        accumVar += synth[2]**2*(synth[3]-1)
    synthMean = np.mean(accumAbs)
    synthStd = np.sqrt(accumVar/(synthCount-1))
    
    return (synthMean, synthStd, synthCount)


def GetAbsForStars(starList=None, clusterName='NGC-0752', ionList=kAllIonList,\
                   gModelAtms=None, gPradks=None,\
                   dModelAtms=None, dPradks=None,\
                   referenceCorrect=False):
# Measures the abundances for the passed stars, and returns a dictionary
# with entries like:
# {ion:[ion/H mean, std.dev., # stars measured], ...}
    if starList==None:
        starList = GetAllStarParms(clusterName=clusterName)

    # Want to use GetAbTable to take advantage of synthesis!
    starAbDict = GetAbTable(clusterName=clusterName, starDataList=starList, 
               ions=ionList, gModelAtms=gModelAtms, gPradks=gPradks, 
               dModelAtms=dModelAtms, dPradks=dPradks,
               referenceCorrect=referenceCorrect)
    # returns:
    # returnDict={starName:{ion:[[X/H],stdev, count],...},...}

    accumulator = {}
    for star in starList:
        starName = star[0]
        if starName not in starAbDict.keys(): continue
        starAbs = starAbDict[starName]
        for ion in starAbs.keys():
            if ion not in accumulator.keys():
                accumulator[ion] = [starAbs[ion][0]]
            else:
                accumulator[ion].append(starAbs[ion][0])
    
    abDict = {}
    for ion in accumulator.keys():
        abDict[ion] = [np.mean(accumulator[ion]), np.std(accumulator[ion]),
                       len(accumulator[ion])]
    return abDict
    
    
#------------------------------------------------------------------------------
# LaTex table code:
#------------------------------------------------------------------------------

# MakeAbTable creates a LaTeX table containing rows with each ion abundance
# [X/Fe] or [Fe/H] for the cluster average, dwarf and giant star averages, 
# then each individual star in the cluster.
def MakeAbTable(clusterName='NGC-0752', starDataList=None, 
                ions=kAllIonList, outFilename=k.TempAbOutFilename, 
                filterBlends=False, referenceCorrect=False, useDAOSpec=False, 
                headerLeft='', headerRight=''):

    if starDataList==None:
        starDataList = GetAllStarParms(clusterName=clusterName)

    # Lookup the abundance(s) for the passed star(s), and create a LaTeX chart

    tableDict = GetAbTable(clusterName, starDataList, ions, 
                           filterBlends, referenceCorrect, useDAOSpec)

    giantNames = np.array([s[0] for s in starDataList \
        if RC.isGiantStar(s[1:5])])
#    dwarfNames = np.array([s[0] for s in starDataList 
#        if not RC.isGiantStar(s[1],s[2])])
    
    # Iron abundance is listed as [Fe/H], all others will be [X/Fe].

    # Compile a cluster/dwarfs/giants average table entry.
    # Entries are of the form:
    # {Ion:[[star ab, star line quality], ...], ...}
    clusterDict = {}
    dwarfDict = {}
    giantDict = {}
    for star in tableDict.keys():
        # Place this stars measures into the proper category dictionary
        if star in giantNames:
            catDict = giantDict
        else:
            catDict = dwarfDict
            
        starFe = GetStarFeH(tableDict[star])

        for ion in tableDict[star].keys():
            solarIonH = PA.ptoe[int(ion)-1][PA.abIdx]
            if tableDict[star][ion][0] == 0.:
                continue

            if ion not in [26.0, 26.1]:
            # Adjust values for [X/Fe]
                tableDict[star][ion][0] = tableDict[star][ion][0]\
                                          - solarIonH\
                                          - starFe

            if ion in clusterDict.keys():
                clusterDict[ion].append([tableDict[star][ion][0],\
                                         tableDict[star][ion][3]])
            else:
                clusterDict[ion] = [[tableDict[star][ion][0],\
                                     tableDict[star][ion][3]]]
                
            if ion in catDict.keys():
                catDict[ion].append([tableDict[star][ion][0],\
                                     tableDict[star][ion][3]])
            else:
                catDict[ion] = [[tableDict[star][ion][0],\
                                 tableDict[star][ion][3]]]

    idxName = clusterName+' (all)'
    tableDict[idxName] = {}
    for ion in clusterDict.keys():
    # If each star's ab measurement for are "independent measures" of the
    # cluster abundance, we should divide the stdev by sqrt(num_measures)...
    # Calcs broken into separate lines for readability.
        cMean = np.mean([l[0] for l in clusterDict[ion]])
        cStd = np.std([l[0] for l in clusterDict[ion]])
        cScoreMean = np.mean([l[1] for l in clusterDict[ion]])
        tableDict[idxName][ion] = \
                        [cMean, cStd, len(clusterDict[ion]), cScoreMean]
        
    idxName = clusterName+' (dwarfs)'
    tableDict[idxName] = {}
    for ion in dwarfDict.keys():
        dMean = np.mean([l[0] for l in dwarfDict[ion]])
        dStd = np.std([l[0] for l in dwarfDict[ion]])
        dScoreMean = np.mean([l[1] for l in dwarfDict[ion]])
        tableDict[idxName][ion] = \
                    [dMean, dStd, len(dwarfDict[ion]), dScoreMean]
                    
    idxName = clusterName+' (giants)'
    tableDict[idxName] = {}
    for ion in giantDict.keys():
        gMean = np.mean([l[0] for l in giantDict[ion]])
        gStd = np.std([l[0] for l in giantDict[ion]])
        gScoreMean = np.mean([l[1] for l in giantDict[ion]])
        tableDict[idxName][ion] = \
                    [gMean, gStd, len(giantDict[ion]), gScoreMean]

    ions = sorted(list(set(ions)))
    # The cluster name starts with "N", so it will be first in the list
    nameList = sorted(tableDict.keys())

    outfile = open(outFilename, 'w')

    latexHeader = kLaTexHeader1 + '\\rhead{\\textbf{' + headerRight\
                  + '}}\n\\lhead{\\textbf{' + headerLeft\
                  + '}}\n\\begin{document}\n'
    outfile.write(latexHeader)
    # New page/table:
    for tableCount in range(int(len(nameList)/6) +1):
        outfile.write('\\begin{landscape}\n'\
                      + '\\hspace*{-5cm}\n'\
                      + '\\begin{tabular}{|l|l|l|l|l|l|l|')
        outfile.write('}\n\\multicolumn{1}{l}{Ion}')
        for name in nameList[tableCount*6:(tableCount+1)*6]:
            outfile.write(' & \\multicolumn{{1}}{{l}}{{{0}}}'.format(name))
        outfile.write('\\\\\n\\hline\n')
        
        # A little bit of stupidity to put FeI and FeII at the top of the table:
        printIons = sorted(ions)
        printIons.remove(26.)
        printIons.remove(26.1)
        printIons = [26.0, 26.1] + printIons
        
        for ion in printIons:
            if int(ion) == 26:
                outfile.write('[{0}/H]'.format(el.getIonName(ion)))
            else:
                outfile.write('[{0}/Fe]'.format(el.getIonName(ion)))

            for star in nameList[tableCount*6:(tableCount+1)*6]:
                if ion in tableDict[star].keys() and tableDict[star][ion][2]>0:
                    outfile.write(' & {0:3.2f}$\pm${1:3.2f} ({2:d}, {3:1.1f})'.\
                                  format(tableDict[star][ion][0], 
                                         tableDict[star][ion][1], 
                                         tableDict[star][ion][2], 
                                         tableDict[star][ion][3]))
                else:
                    outfile.write(' & \multicolumn{1}{c|}{---}')

            outfile.write('\\\\\n\\hline\n')
            # Place a double line after the Fe/H measures
            if ion==26.1:
                outfile.write('\\hline\n')

        outfile.write('\\label{{tab:752Abs-{0}}}\n'.format(tableCount+1)
                      + '\\end{tabular}\n\\end{landscape}\n'
                      + '\\clearpage\n')

    outfile.write('\\end{document}\n')
    outfile.close()


#------------------------------------------------------------------------------
# Plots and plotting code:
#------------------------------------------------------------------------------

# PlotClusterAbs creates a png plot of the mean giant and dwarf star abundance
# for all ions requested for the passed cluster.
def PlotClusterAbs(clusterName='NGC-0752', ionList=kAllIonList, 
                   fileTag='', mass=1.2, metallicity=0.0,
                   solarRelative = False, referenceCorrect=False):
    giantParms = GetGiantStarParms(clusterName=clusterName,\
                                   mass=mass, metallicity=metallicity)
    modelFiles = mk.findDataFiles(k.GiantModelPath)
    gAtms, gPradks = mk.LoadModels(modelFiles)
    modelFiles = mk.findDataFiles(k.DwarfModelPath)
    dAtms, dPradks = mk.LoadModels(modelFiles)
    
    gAbs = GetAbsForStars(giantParms,\
                          gModelAtms=gAtms, gPradks=gPradks,\
                          dModelAtms=dAtms, dPradks=dPradks,\
                          referenceCorrect=referenceCorrect)
    
    gIons = sorted([ion for ion in gAbs.keys() if ion in ionList])
    if solarRelative:
        gElems = [gAbs[i][0] for i in gIons]
    else:
        solarFeH = PA.ptoe[25][PA.abIdx]
        starIons = gAbs.keys()
        starFeH = GetStarFeH(gAbs)
        
        gElems = [gAbs[i][0]-PA.ptoe[int(i)-1][PA.abIdx]-starFeH for i in gIons]

    gVar = [gAbs[i][1] for i in gIons]
    
    dwarfParms = GetDwarfStarParms(clusterName=clusterName,\
                                   mass=mass, metallicity=metallicity)
    
    dAbs = GetAbsForStars(dwarfParms,\
                          gModelAtms=gAtms, gPradks=gPradks,\
                          dModelAtms=dAtms, dPradks=dPradks,\
                          referenceCorrect=referenceCorrect)
    
    dIons = sorted([ion for ion in dAbs.keys() if ion in ionList])
    if solarRelative:
        dElems = [dAbs[i][0] for i in dIons]
    else:
#        solarFeH = PA.ptoe[25][PA.abIdx] - calculated in giants, above
        starIons = dAbs.keys()
        starFeH = GetStarFeH(dAbs)
        
        dElems = [dAbs[i][0]-PA.ptoe[int(i)-1][PA.abIdx]-starFeH for i in dIons]
        
    dVar = [dAbs[i][1] for i in dIons]
    
    ionLUT = {}
    
    for num, ion in enumerate(sorted(set(dIons+gIons))):
        ionLUT[ion] = num
    
    ionLabels = [' -',' ']
    for ion in sorted(ionLUT.keys()): ionLabels.append(el.getIonName(ion))
    ionLabels.append(' ')
    
    fig = pyplot.figure()
    ax = fig.gca()
    
    ax.errorbar([ionLUT[ion] for ion in gIons], gElems, yerr=gVar, fmt='rs',\
                label=clusterName+' giants')
    ax.errorbar([ionLUT[ion] for ion in dIons], dElems, yerr=dVar, fmt='bo',\
                label=clusterName+' dwarfs')
                
    ax.xaxis.set_major_locator(pyplot.MultipleLocator(1))
    ax.set_xticklabels(ionLabels, rotation='vertical')
    ax.set_xlabel('Ion')
    ax.set_xlim(left=-1, right=len(ionLUT.keys()))
    if solarRelative:
        ax.set_ylabel('[X/H]')
    else:
        ax.set_ylabel('[X/Fe]')
#    ax.set_ylim(bottom=-0.4, top=0.5)
    if not solarRelative:
        ax.set_ylim(bottom=np.average(gElems+dElems)-0.5, top=np.average(gElems+dElems)+0.5)
    else:
        ax.set_ylim(bottom=0, top=9)
    ax.legend(numpoints=1)
    ax.axhline(y=0.0, xmin=-1, xmax=len(ionLUT.keys()), color='g', linestyle='--')
#    pyplot.show()
    pyplot.savefig(dirs.PlotDir+clusterName+'Abs'+fileTag+'.png')
    pyplot.close()
    return


# Excitation potential plots allow us to see if there are trends between the
# abundance calculated from a given line to the excitation potential for that
# line. If the Teff value for the atmospheric model of a given star is
# "correct", there should be no trend.
def PlotXPAbs(starData, clusterName='NGC-0752', 
            ionList=kAllIonList, fileTag='', labelPlot=True, 
            labelPoints=False, showTrendLine=False,
            modelAtms=None, pradks=None, referenceCorrect=False):
# Make XP vs. Ab for the passed star
# One element per plot.
    starName = starData[0]
    starParmTuple = tuple(starData[1:])

    isGiant = RC.isGiantStar(starParmTuple)
    if isGiant:
        modelPath = k.GiantModelPath
    else:
        modelPath = k.DwarfModelPath

    if modelAtms == None or pradks == None:
        modelFiles = mk.findDataFiles(modelPath)
        modelAtms, pradks = mk.LoadModels(modelFiles)

    abdict, uncorrLines, unusedMin, unusedMax = \
                AB.CalcAbsAndLines(clusterName+' '+starName, 
                                   tuple(starData[1:6]), ionList=ionList, 
                                   modelAtms=modelAtms, pradks=pradks)
    # uncorrLines:
    # # {elem.ion:[[Wavelength, Ex.Pot., logGf, eqw, logRW, abund],...]}

    if referenceCorrect:
        if isGiant: # Obligatory comment on bad Giant corrections, and using
                    # Solar instead.
            correctDict, referenceLines, lineWeights = \
            RC.GetSolarCorrections(ionList=ionList, 
                                modelAtms=modelAtms, 
                                pradks=pradks)
        else:
            correctDict, referenceLines, lineWeights = \
            RC.GetDwarfCorrections(ionList=ionList, 
                                modelAtms=modelAtms, 
                                pradks=pradks)

    correctionsAvailable = False
    if len(correctDict) > 0 and len(referenceLines)>0:
        correctionsAvailable = True
        
    for ion in ionList:
        if ion not in uncorrLines.keys():
            continue
        # Does this element have NLTE corrections available? Note: we do
        # this BEFORE Solar corrections, which assumes that the same
        # NLTE corrections are already applied to any Solar corrections
        # we use.
        if ion in NLTEIons:
            LTELines = AB.CorrectNLTEAbs(ion, uncorrLines[ion], 
                                         tuple(starData[1:6]))
        else:
            if ion in uncorrLines.keys():
                LTELines = uncorrLines[ion]
            else:
                # Either synthesized lines, or none available
                LTELines = np.array([])
        
        # Do we want the "reference corrected" abundances?
        if referenceCorrect and correctionsAvailable:
            tempAdj,tempAll = \
                    RC.SortAndFilterLines(LTELines, 
                                       ion, 
                                       tuple(starData[1:6]), 
                                       solarCorrect=referenceCorrect, 
                                       solarLines=referenceLines[ion], 
                                       solarCorrs=correctDict[ion], 
                                       lineWeights=lineWeights[ion])

            # correctedLines:
            #           [[ab, line STR score, wl, "quality"], ...]
            # allLines:
            #           [[ab, line STR score, wl],...]
            # We want to use np.arrays, so...
            allLines = np.array(tempAll)
            correctedLines = np.array(tempAdj)
            if len(allLines)==0 or len(correctedLines)==0:
                correctionsAvailable = False
            elif len(allLines)==1 or len(correctedLines)==1:
                print('Single Line determination:{0}'.format(starData[0]))
                print(allLines, correctedLines)
        # One plot per ion.
        if labelPlot:
            plotLabel = 'XP vs Ab for [{2}/H] in {0} {1}.'.\
                        format(clusterName, starName, el.getIonName(ion))
        else:
            plotLabel = ''
        
        if referenceCorrect and correctionsAvailable:
            tempPoints = []
            for line in uncorrLines[ion]:
                correctedAbs = [l[0] for l in correctedLines if \
                                u.in_range(l[2],line[0]-0.05, line[0]+0.05)]
                if len(correctedAbs)>0:
                    tempPoints.append([line[1],np.mean(correctedAbs),line[3],line[0]])
                else:
                    tempPoints.append([line[1],line[5],line[3],line[0]])
            XPAbPoints = np.array(tempPoints)
        else:
            XPAbPoints = np.array([[line[1],line[5],line[3],line[0]]\
                                for line in uncorrLines[ion]])
        if labelPoints:
        # Label the points with the wavelength
            pointLabels = ['{0:2.3f}'.format(point[3]) \
                           for point in XPAbPoints]
        else:
            pointLabels = None

        ps.XPAbPlot(XPAbPoints, starName, ion, fileTag=fileTag+'XPAb', 
                    plotTitle=plotLabel, pointLabels=pointLabels,
                    showTrendLine=showTrendLine)

# Helper functions for NGC 752:
def PlotGiantXPAbs(clusterName='NGC-0752', mass=1.6, metallicity=0.00,
                   fileTag='', labelPlot=True ,labelPoints=False, 
                   ionList=kAllIonList, showTrendLine=False, 
                   referenceCorrect=False):
    modelFiles = mk.findDataFiles(k.GiantModelPath)
    modelAtms, pradks = mk.LoadModels(modelFiles)
    giantData =  GetGiantStarParms(clusterName=clusterName,\
                                   mass=mass, metallicity=metallicity)
    for star in giantData:
        PlotXPAbs(star, ionList=ionList, fileTag=fileTag, 
                labelPlot=labelPlot, labelPoints=labelPoints,
                showTrendLine=showTrendLine,
                modelAtms=modelAtms, pradks=pradks,
                referenceCorrect=referenceCorrect)


def PlotDwarfXPAbs(clusterName='NGC-0752', mass=0.0, metallicity=0.00,
                   fileTag='', labelPlot=True ,labelPoints=False, 
                   ionList=kAllIonList, showTrendLine=False, 
                   referenceCorrect=False):
    modelFiles = mk.findDataFiles(k.DwarfModelPath)
    modelAtms, pradks = mk.LoadModels(modelFiles)
    dwarfData = GetDwarfStarParms(clusterName=clusterName,\
                                  mass=mass, metallicity=metallicity)
    for star in dwarfData:
        PlotXPAbs(star, ionList=ionList, fileTag=fileTag, 
                labelPlot=labelPlot, labelPoints=labelPoints, 
                showTrendLine=showTrendLine,
                modelAtms=modelAtms, pradks=pradks,
                referenceCorrect=referenceCorrect)


def PlotAllXPAbs(clusterName='NGC-0752', mass=0.0, metallicity=0.00,
                 fileTag='', labelPlot=True ,labelPoints=False, 
                 ionList=kAllIonList, showTrendLine=False,
                 referenceCorrect=False):
    PlotGiantXPAbs(fileTag=fileTag, labelPlot=labelPlot,
                   labelPoints=labelPoints, ionList=ionList,
                   showTrendLine=showTrendLine,
                   referenceCorrect=referenceCorrect)
    PlotDwarfXPAbs(fileTag=fileTag, labelPlot=labelPlot,
                   labelPoints=labelPoints, ionList=ionList,
                   showTrendLine=showTrendLine,
                    referenceCorrect=referenceCorrect)

                        
# Curve-of-Growth plots let us look at individual elements for a given star. 
# Comparing the calculated abundance vs. excitation potential or, more 
# accurately, a "STR" value calculated from Ex.Pot., log(gf), and Teff, as per
# Sneden et al. (can't recall the date), lets us determine whether a given
# line measurement is too weak to be measured for a given star, too strong
# and therefore in the non-linear region of the curve-of-growth, or too broad
# and probably contaminated by another blended line, or (least likely) too
# narrow, and probably a mis-measurement.
def PlotGiantCOGs(clusterName='NGC-0752', mass=1.7, metallicity=-0.10,
                  fileTag='', labelPlot=True ,labelPoints=False, 
                  ionList=kAllIonList, makeQualityPlot=False):
    modelFiles = mk.findDataFiles(k.GiantModelPath)
    modelAtms, pradks = mk.LoadModels(modelFiles)
    giantData =  GetGiantStarParms(clusterName=clusterName,\
                                   mass=mass, metallicity=metallicity)
    for star in giantData:
        COGPlot(clusterName, star, ionList, fileTag=fileTag, 
                labelPlot=labelPlot, labelPoints=labelPoints, 
                modelAtms=modelAtms, pradks=pradks,
                makeQualityPlot=makeQualityPlot)


def PlotDwarfCOGs(clusterName='NGC-0752', mass=0.0, metallicity=-0.10,
                  fileTag='', labelPlot=True ,labelPoints=False, 
                  ionList=kAllIonList, makeQualityPlot=False):
    modelFiles = mk.findDataFiles(k.DwarfModelPath)
    modelAtms, pradks = mk.LoadModels(modelFiles)
    dwarfData = GetDwarfStarParms(clusterName=clusterName,\
                                  mass=mass, metallicity=metallicity)
    for star in dwarfData:
        COGPlot(clusterName, star, ionList, fileTag=fileTag, 
                labelPlot=labelPlot, labelPoints=labelPoints, 
                modelAtms=modelAtms, pradks=pradks,
                makeQualityPlot=makeQualityPlot)


def PlotAllCOGs(clusterName='NGC-0752', dMass=0.0, gMass=1.7, 
                metallicity=-0.10, fileTag='', labelPlot=True, 
                labelPoints=False, ionList=kAllIonList, makeQualityPlot=False):
    PlotGiantCOGs(clusterName=clusterName, mass=gMass,
                  metallicity=metallicity, fileTag=fileTag, 
                  labelPlot=labelPlot, labelPoints=labelPoints, 
                  ionList=ionList, makeQualityPlot=makeQualityPlot)
    PlotDwarfCOGs(clusterName=clusterName, mass=dMass,
                  metallicity=metallicity, fileTag=fileTag, 
                  labelPlot=labelPlot, labelPoints=labelPoints, 
                  ionList=ionList, makeQualityPlot=makeQualityPlot)

def COGPlot(clusterName, starData, ions, fileTag='', labelPlot=True, 
            labelPoints=False, modelAtms=None, pradks=None, 
            makeQualityPlot=False):
# Make Curve of Growth (COG) plots for the passed star
# One element per plot.
# If the makeQualityPlot flag is set, a second plot will be made, with the
# points shaded relative to how well they correspond to the linear fit of
# measurements with similar excitation potentials.
    starName = starData[0]
    starParmTuple = tuple(starData[1:])

    isGiant = RC.isGiantStar(starParmTuple)
    if isGiant:
        modelPath = k.GiantModelPath
    else:
        modelPath = k.DwarfModelPath

    if modelAtms == None or pradks == None:
        modelFiles = mk.findDataFiles(modelPath)
        modelAtms, pradks = mk.LoadModels(modelFiles)

    abdict, uncorrLines, unusedMin, unusedMax = \
                AB.CalcAbsAndLines(clusterName+' '+starName, 
                                   tuple(starData[1:6]), ionList=ions, 
                                   modelAtms=modelAtms, pradks=pradks)

    for ion in ions:
    # One plot per ion.

        if ion not in uncorrLines.keys():
            print('No {0:2.1f} lines for {1}.'.format(ion, starName))
            continue

        if labelPlot:
            plotLabel = 'Curve of Growth for [{2}/H] in {0} {1}.'.\
                        format(clusterName, starName, el.getIonName(ion))
        else:
            plotLabel = ''

        if labelPoints:
        # Label the points with the excitation potential
            pointLabels = ['{0:2.3f}'.format(point[1]) \
                           for point in np.array(uncorrLines[ion])]
        else:
            pointLabels = None

        ps.COGPlot(np.array(uncorrLines[ion]), starName, ion, fileTag=fileTag,
                            plotTitle=plotLabel, pointLabels=pointLabels)

        if makeQualityPlot:
            if labelPlot:
                plotLabel = 'Measurement quality for [{2}/H] in {0} {1}.'.\
                            format(clusterName, starName, el.getIonName(ion))
            else:
                plotLabel = ''
            ps.QualityPlot(np.array(uncorrLines[ion]), starName, ion, 
                           fileTag=fileTag+'QP', plotTitle=plotLabel, 
                           pointLabels=pointLabels)


def PlotAbs(clusterName, starData, ions, fileTag='', plotTitle='', 
            filterBlends=False, referenceCorrect=False, labelPoints=False, 
            useDAOSpec=False, modelAtms=None, pradks=None):
    pointLabels = None
# Lookup the abundance(s) for the passed star, and plot them 
# (into a .png file). One element per plot.
    starName = starData[0]
    starParmTuple = tuple(starData[1:6])

    isGiant = RC.isGiantStar(starParmTuple)

    if isGiant:
        modelPath = k.GiantModelPath
    else:
        modelPath = k.DwarfModelPath

    if modelAtms == None or pradks == None:
        modelFiles = mk.findDataFiles(modelPath)
        modelAtms, pradks = mk.LoadModels(modelFiles)
    
    # uncorrLines format: 
    #       {elem.ion:[[Wavelength, Ex.Pot., logGf, eqw, logRW, abund],...]}
    # abDict format:
    #       {elem.ion:[abundance mean, ab. std.dev., # lines]}
    abdict, uncorrLines, unusedMin, unusedMax = \
                AB.CalcAbsAndLines(clusterName+' '+starName, 
                                   starParmTuple, 
                                   ionList=ions, 
                                   modelAtms=modelAtms, 
                                   pradks=pradks,
                                   useDAOlines=useDAOSpec)

    if isGiant:
        # Obtain the reference corrections for a giant star - Note: this is
        # badly broken! So, we're just going to use the Solar corrections for
        # now:
        correctDict, referenceLines, lineWeights = \
                RC.GetDwarfCorrections(ionList=ions, 
                                    modelAtms=modelAtms, 
                                    pradks=pradks, 
                                    useDAOSpec=useDAOSpec)
#        correctDict, referenceLines, lineWeights = \
#                GetGiantCorrections(ionList=ions, 
#                                    modelAtms=modelAtms, 
#                                    pradks=pradks)
    else:
        # ...or for a dwarf star.
        correctDict, referenceLines, lineWeights = \
                RC.GetDwarfCorrections(ionList=ions, 
                                    modelAtms=modelAtms, 
                                    pradks=pradks, 
                                    useDAOSpec=useDAOSpec)

    for ion in ions:

    # We'll make one plot per ion.
        pointLabels = []
        redData = greenData = blueData = []
        if ion not in uncorrLines.keys():
            print('No {0:2.1f} lines for {1}.'.format(ion, starName))
            continue
        if (referenceCorrect or filterBlends) and \
        ion in referenceLines.keys() and ion in correctDict.keys() and \
        ion in lineWeights.keys():
            adjustedLines, dataPoints = \
                    RC.SortAndFilterLines(uncorrLines[ion], 
                                       ion, 
                                       starParmTuple, 
                                       filterBlends=filterBlends, 
                                       solarCorrect=referenceCorrect, 
                                       solarLines=referenceLines[ion], 
                                       solarCorrs=correctDict[ion], 
                                       lineWeights=lineWeights[ion])

            tempData = []
            for line in uncorrLines[ion]:
                corrLine = u.find(lambda l: l[0]==line[0], dataPoints)
                if corrLine is not None:
                    tempData.append(corrLine)
                else:
                    tempData.append([line[5], el.STR(line[2], line[1], \
                                     starParmTuple[0]),line[0]])
            redData = np.array(tempData)
        else:

            dataPoints = np.array([[line[5], 
                                    el.STR(line[2], 
                                    line[1], 
                                    starParmTuple[0]),
                                    line[0]] \
                                    for line in uncorrLines[ion]])

        if labelPoints:
            pointLabels = ['{0:4.3f}'.format(line[2]) for line in dataPoints]
            pointLabels.extend(['{0:4.3f}'.\
                                format(line[2]) for line in redData])
            pointLabels.extend(['{0:4.3f}'.\
                                format(line[2]) for line in greenData])
            pointLabels.extend(['{0:4.3f}'.\
                                format(line[2]) for line in blueData])

        loSlope, loIntercept, rv, pv, err = \
                ps.GetDetectionLimit(starParmTuple,
                                     ion, 
                                     modelAtms=modelAtms, 
                                     pradks=pradks)
        hiSlope, hiIntercept, rv, pv, err = \
                ps.GetCoGLimit(starParmTuple, 
                               ion, 
                               modelAtms=modelAtms, 
                               pradks=pradks)

        ps.AbSTRPlot(starName, ion, dataPoints, redSet=redData, 
                     greenSet=greenData, blueSet=blueData, 
                     lowLimit=(loSlope, loIntercept), 
                     hiLimit=(hiSlope, hiIntercept), 
                     fileTag=fileTag, plotTitle=plotTitle, 
                     labelPoints=pointLabels)


# "Quick" Plot Functions - Not so sure these are necessary - For a
# generalized output package, we should have access to these as an
# external call, with a giant/main-sequence flag, if that is the
# desired functionality.
def qPlotGiants(clusterName='NGC-0752', mass=1.6, metallicity=0.0,
                fileTag='', plotTitle='', filterBlends=False, 
                referenceCorrect=False, labelPoints=False, useDAOSpec=False,
                ionList=kAllIonList):
    # Function produces a T vs. [X/H] plot for all giant stars in NGC 752
    modelFiles = mk.findDataFiles(k.GiantModelPath)
    modelAtms, pradks = mk.LoadModels(modelFiles)
    giantData =  GetGiantStarParms(clusterName=clusterName,\
                                   mass=mass, metallicity=metallicity)
    for star in giantData:
        PlotAbs(clusterName, star, ionList, fileTag=fileTag,
                plotTitle=plotTitle, filterBlends=filterBlends,
                referenceCorrect=referenceCorrect, labelPoints=labelPoints,
                useDAOSpec=useDAOSpec, modelAtms=modelAtms, pradks=pradks)


def qPlotDwarfs(clusterName='NGC-0752', mass=0.0, metallicity=-0.10,
                fileTag='', plotTitle='', filterBlends=False, 
                referenceCorrect=False, labelPoints=False, useDAOSpec=False,
                ionList=kAllIonList):
    # Function produces a T vs. [X/H] plot for all dwarf stars in NGC 752
    modelFiles = mk.findDataFiles(k.DwarfModelPath)
    modelAtms, pradks = mk.LoadModels(modelFiles)
    dwarfData =  GetDwarfStarParms(clusterName=clusterName,\
                                   mass=mass, metallicity=metallicity)
    for star in dwarfData:
        PlotAbs(clusterName, star, ionList, fileTag=fileTag,
                plotTitle=plotTitle, filterBlends=filterBlends,
                referenceCorrect=referenceCorrect, labelPoints=labelPoints,
                useDAOSpec=useDAOSpec, modelAtms=modelAtms, pradks=pradks)


def qPlotAll(clusterName='NGC-0752', dMass=0.0, gMass=1.7, metallicity=-0.10,
        fileTag='', plotTitle='', filterBlends=False, referenceCorrect=False,
        labelPoints=False, useDAOSpec=False, ionList=kAllIonList):
    # Function produces a T vs. [X/H] plot for all stars in NGC 752
    qPlotGiants(clusterName=clusterName, mass=gMass, metallicity=metallicity,
                fileTag=fileTag, plotTitle=plotTitle,
                filterBlends=filterBlends, referenceCorrect=referenceCorrect,
                labelPoints=labelPoints, useDAOSpec=useDAOSpec,
                ionList=ionList)
    qPlotDwarfs(clusterName=clusterName, mass=dMass, metallicity=metallicity,
                fileTag=fileTag, plotTitle=plotTitle,
                filterBlends=filterBlends, referenceCorrect=referenceCorrect,
                labelPoints=labelPoints, useDAOSpec=useDAOSpec,
                ionList=ionList)

