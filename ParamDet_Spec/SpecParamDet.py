#! /usr/bin/env python
# 
# Module: SpecParamDet.py
#
# Author: Mike Lum
#
# Description: Functions used to determine the atmospheric parameters
#   of a star with the passed spectra.
#
# Contents: 
#   Function xxx: (Description)
#
# Revision History:
#    Date        Vers.    Author        Description
#    2015-06-24    0.0a0    Lum        First checked in
#
# To Do:
#    
#
import os
import glob
import datetime
from collections import deque
import shutil
import sqlite3
import sys
import time

from scipy import stats
import numpy as np
from matplotlib import pyplot as PyP

from Abundances import ReferenceCorrections as RC
from constants import constants as k
from constants import DBFormats as DB
from Databases import LineLookup as LL
from Isochrones import IsoProb as IP
from MOOGInterface import MOOGInterface as MI
from Models import MakeAtmModel as mk
from utilities import elements as el
from utilities import utilities as u

# Functions
def isGiant(Teff, Logg):
    print('SpecParamDet.py: isGiant depricated. Use utilities.py:isGiant, instead.')
    return u.isGiant(Teff, Logg)
        
        
def getGuesses(logFile, guessFilename):
    tGuess = k.SolarTeff
    gGuess = k.SolarLogG
    vGuess = k.SolarVturb

    if guessFilename is not None:
    # Assume the log file has the correct format of: NNNN_SSS_ssssss_II_OOO.fits
    # Where: 
    #           NNNN  = NGC catalog number of parent cluster 
    #                    (note: clusters not contained in the NGC catalog will be prefaced by a letter
    #                     representative of their catalog. ex: I = IC
    #           SSS    = Reference indexing system for this cluster. 
    #                    (ex: MJP = Montgomery, Janes, and Phelps)
    #           ssssss = Star index number from the aformetioned index.
    # ...which means that we can look up the guess by the star name.
        logFileWords = os.path.basename(logFile).split('_')
#        starName = '_'.join(logFileWords[0:3])
        starName = logFileWords[0]
        print('Looking for star: %s' % starName)
        try:
            guessFile = open(guessFilename, 'r')
            guessLines = guessFile.readlines()
        except:
            u.errorOut('Error in parameter determination:\nUnable to open file:%s' % guessFilename, True)

        try:
            guessStr = [guessLine for guessLine in guessLines if guessLine.find(starName)==0][0]
        except IndexError:
            return tGuess, gGuess, vGuess

        guesses = guessStr.split()
        if len(guesses) > 0 and u.is_number(guesses[1]):
            tGuess = float(guesses[1])
        if len(guesses) > 1 and  u.is_number(guesses[2]):
            gGuess = float(guesses[2])
        if len(guesses) > 2 and  u.is_number(guesses[3]):
            vGuess = float(guesses[3])
                    
    return tGuess, gGuess, vGuess
    
def getMOOGLisLines(FeILines, FeIILines, TiILines, TiIILines):
    # Put all the passed lines into a MOON-readable form.
    lines = []
    
    if u.is_list(FeILines):
    # Normal list case
        lines = [[x[0], 26.0, x[1], x[2], 2.0, x[3]] for x in FeILines]
    elif len(FeIILines) > 0:
        # One line - We're actually pretty hosed, statistically, if we only
        # have one line, but that's for someone else to decide.
        lines = [[FeILines[0], 26.0, FeILines[1], FeILines[2], 2.0, FeILines[3]]]
    else:
        # No FeI Lines - I don't care who you are. We're giving up now.
        return lines
    
    if u.is_list(FeIILines):
        lines = np.concatenate((lines, [[x[0], 26.1, x[1], x[2], 2.0, x[3]] for x in FeIILines]))
    elif len(FeIILines) > 0:
        # One line
        lines = np.concatenate((lines, [[FeIILines[0], 26.1, FeIILines[1], FeIILines[2], 2.0, FeIILines[3]]]))
    #else:
        # No FeII Lines
        #pass
    if u.is_list(TiILines):
        lines = np.concatenate((lines, [[x[0], 22.0, x[1], x[2], 2.5, x[3]] for x in TiILines]))
    elif len(TiILines) > 0:
        lines = np.concatenate((lines, [[TiILines[0], 22.0, TiILines[1], TiILines[2], 2.5, TiILines[3]]]))
        
    if u.is_list(TiIILines):
        lines = np.concatenate((lines, [[x[0], 22.1, x[1], x[2], 2.5, x[3]] for x in TiIILines]))
    elif len(TiIILines) > 0:
        lines = np.concatenate((lines, [[TiIILines[0], 22.1, TiIILines[1], TiIILines[2], 2.5, TiIILines[3]]]))
        
    return lines


def selectLinesForElement(lines, ion, modelAtms, pradks, starParms, varLimit = np.finfo(np.float64).max, enforceMin=False, refLines=None, refCorrs=None, lineWeights=None):
    # From the passed list, select all the lines of the passed ion
    # (in XX.x form), with abundances within +/-varLimit of the mean.
    # An additional selection criteria is that the line eqw falls on
    # the linear portion of the curve of growth.
    # If enforceMin = True, the selection will be limited to the lesser
    # of the total number of lines for this element, or the constant
    # MinNumStatLines.
    # To save time in looped, or multiple calls, solar (or giant) corrections
    # can be passed in refLines and refCorrs. This is basically the output
    # of RC.getGiantCorrections or RC.getSolarCorrections for the passed ion.

    corrLines, allLines = RC.SortAndFilterLines(lines[ion], ion, starParms, filterBlends=True, refCorrect=True, modelAtms=modelAtms , pradks=pradks, refLines=refLines, refCorrs=refCorrs, lineWeights=lineWeights)
    
    # If we did not receive enough (any) corrected lines, so return a subset of
    # all the lines.
    if (corrLines==None or len(corrLines) < k.MinNumStatLines):
        if len(allLines) <= k.MinNumStatLines:
            # Well, we can't make more lines, so just return what we have
            return sorted(allLines, key=lambda line: line[0])
        else:
            meanAb = np.mean(np.array(allLines)[:,5])
            allSorted = sorted(allLines, key=lambda line: abs(line[5]-meanAb))[:k.MinNumStatLines]
            return sorted(allSorted, key=lambda line: line[0])
    elif enforceMin:
        # We have enough lines to consider filtering for the variance limit.
        meanAb = np.mean(np.array(corrLines)[:,5])
        limitedLines = [line for line in corrLines if abs(line[5]-meanAb)<=varLimit]
        if len(limitedLines) < k.MinNumStatLines:
        # Variance filtering was a bit much, return the best we can
            return sorted(sorted(corrLines, key=lambda line: abs(line[5]-meanAb))[:k.MinNumStatLines], key=lambda l: l[0])
        else:
            return sorted(limitedLines, key=lambda line: line[0])
    else:
        # enforceMin not set, so filter away!
        meanAb = np.mean(np.array(corrLines)[:,5])
        rtnLines = [x for x in corrLines if abs(x[5]-meanAb) < abs(varLimit)]
        return sorted(rtnLines, key=lambda line: line[0])



def selectParamLines(lineList, ion):
    # Returns a sub-list of the passed list, consisting of only those lines
    # which have the "ForParams" flag set in the line database
    newList = []
#    print(lineList)
    with sqlite3.connect(DB.ClusterDBName) as DBConnection:
        DBCursor = DBConnection.cursor()

        for line in lineList:
            DBCursor.execute('''SELECT {0}, {1} FROM {2} WHERE ion={3:2.1f} AND Wavelength BETWEEN {4:8.3f} AND {5:8.3f} ORDER BY Wavelength;'''.format(DB.LineTableFieldNames[0], DB.LineTableFieldNames[8], DB.LineTableName, ion, line[0]-0.05,  line[0]+0.05))
            DBLines = DBCursor.fetchall()
            # Only check the first result
            for retdLine in DBLines:
                if int(retdLine[1]) == 1:
                    newList.append(line)
#    print(newList)
#    raw_input("Pausing...")
    
    return newList
    
def selectLines(lines, modelAtms, pradks, starParms, varLimit = np.finfo(np.float64).max, enforceMin = False, refDict=None):
    # Select lines where the equivalent width measures are still within the 
    # linear region of the curve of growth, and where the measure abundance 
    # does not vary by more than a (passed) limit.
    # To streamline multiple calls to MOOG for abundance calculations, a
    # reference line dictionary can be passed. The dictionary would be of the 
    # form:
    # {ion:(refLines, refCorrs)}
    FeILines, FeIILines, TiILines, TiIILines = [],[],[],[]

    if refDict is not None:
        try:
            # FeI line selection
            ion = 26.0
            FeILines = selectLinesForElement(lines, ion, modelAtms, pradks,\
                    starParms, varLimit = varLimit, enforceMin = enforceMin,\
                    refLines=refDict[ion][0], refCorrs=refDict[ion][1],\
                    lineWeights=refDict[ion][2])
                
            # FeII abundances for comparison
            ion=26.1
            FeIILines = selectLinesForElement(lines, ion, modelAtms, pradks,\
                    starParms, varLimit = varLimit, enforceMin = enforceMin,\
                    refLines=refDict[ion][0], refCorrs=refDict[ion][1],\
                    lineWeights=refDict[ion][2])
                    
            # TiI lines - We have to accept a wider range of variances with Ti lines
            ion = 22.0
            TiILines =selectLinesForElement(lines, ion, modelAtms, pradks,\
                    starParms, varLimit = varLimit*2, enforceMin = enforceMin,\
                    refLines=refDict[ion][0], refCorrs=refDict[ion][1],\
                    lineWeights=refDict[ion][2])
            
            # TiII
            ion = 22.1
            TiIILines = selectLinesForElement(lines, ion, modelAtms, pradks,\
                    starParms, varLimit = varLimit*2, enforceMin = enforceMin,\
                    refLines=refDict[ion][0], refCorrs=refDict[ion][1],\
                    lineWeights=refDict[ion][2])
                    
        except KeyError:
        # Occurs when MOOG bombs out, or we could not find lines for one of the element/ionization pairs.
        # Just return what we have so far
            pass
    else:
        try:
            # FeI line selection
            FeILines = selectLinesForElement(lines, 26.0, modelAtms, pradks,\
                    starParms, varLimit = varLimit, enforceMin = enforceMin)
                
            # FeII abundances for comparison
            FeIILines = selectLinesForElement(lines, 26.1, modelAtms, pradks,\
                    starParms, varLimit = varLimit, enforceMin = enforceMin)
                    
            # TiI lines - We have to accept a wider range of variances with Ti lines
            TiILines =selectLinesForElement(lines, 22.0, modelAtms, pradks,\
                    starParms, varLimit = varLimit*2, enforceMin = enforceMin)
            
            # TiII
            TiIILines = selectLinesForElement(lines, 22.1, modelAtms, pradks,\
                    starParms, varLimit = varLimit*2, enforceMin = enforceMin)

        except KeyError:
        # Occurs when MOOG bombs out, or we could not find lines for one of the element/ionization pairs.
        # Just return what we have so far
            pass

    return FeILines, FeIILines, TiILines, TiIILines
        
def evaluateParams(FeILines, FeIILines, TiILines, TiIILines):
    # The evaluation score is the product of the p-values for matching
    # I and II ionization state abundances, and having a slope of 0.0
    # (ie: no correlation between) excitation potential and measured 
    # abundance for a given line.
    #
    # Note that we are being VERY loose with our definition of "probabilities"
    # Specifically, we are using the P-score of the linear regression and of
    # the T-Test as a probability. This is technically incorrect, as the 
    # P-score should really only be used for rejection of a null hypothesis 
    # that the abundance vs. excitation potential slope is zero, or that the
    # two populations have the same mean.
    # Since we are really only interested in getting a relative scoring for the
    # quality of the linear fit to a zero slope, and the quality of the match
    # of the abundances of the two population, the use of a P-score as a 
    # scoring mechanism is valid.
    if len(FeILines) < k.MinNumStatLines:
        print('Warning: Only {0:2d} FeI Lines!'.format(len(FeILines)))
        # If the number is _really_ small, we realistically can't do anything...
        if len(FeILines) <  k.MinNumStatLines/5 or not u.is_list(FeILines): return 0., 0.

    FeIExPotList = np.array(FeILines)[:,1].tolist()
    FeIAbList = np.array(FeILines)[:,5].tolist()
    FeSlope, intercept, r_value, FeP_value, stderr = stats.linregress(FeIExPotList,FeIAbList)

    if len(FeIILines) < k.MinNumStatLines/2:
    # Not enough FeII lines, so we have to assume that the Fe I abundance is definitive
        FeAbTProb = k.BadProbability
    else:
        FeIIAbList = np.array(FeIILines)[:,5].tolist()
        FeAbTScore, FeAbTProb = stats.ttest_ind(FeIAbList, FeIIAbList)

    if len(TiILines) < k.MinNumStatLines/2:
    # Not enough TiI lines, so we will just use the Fe I line slope
        TiP_value = k.BadProbability
        TiAbTProb = k.BadProbability
    else:
        TiIExPotList = np.array(TiILines)[:,1].tolist()
        TiIAbList = np.array(TiILines)[:,5].tolist()
        TiSlope, intercept, r_value, TiP_value, stderr = stats.linregress(TiIExPotList,TiIAbList)

        if len(TiIILines) < k.MinNumStatLines/4:
        # Not enough TiII lines, so we will just use the Fe I/II comparison
            TiAbTProb = k.BadProbability
        else:
            TiIIAbList = np.array(TiIILines)[:,5].tolist()
            TiAbTScore, TiAbTProb = stats.ttest_ind(TiIAbList, TiIIAbList)
        
        if not u.is_number(TiAbTProb):
            TiAbTProb = k.BadProbability

    # Returns slope probability, ionization balance probability
    # Note that python floating point calculations should treat
    # any value of around 10**-16 as zero, so we force calculations
    # with values of < 10**-14 to 0.0
    # Testing with Fe Only:
#    slopeProb = FeP_value*TiP_value
    slopeProb = FeP_value
    if slopeProb < k.BadProbability/100.:
        slopeProb = 0.
    # I'm not sure I like using the Ti I/I abundance balance (NLTE effects?)
    # Weight Fe by 2x over Ti:
    balanceProb = (FeAbTProb**2*TiAbTProb)**(1./3.)
#    balanceProb = np.sqrt(FeAbTProb*TiAbTProb)
#    balanceProb = FeAbTProb
    if balanceProb < k.BadProbability/100.:
        balanceProb = 0.
    return slopeProb, balanceProb


def getProbKDE(clustMetal, clustAge, starTeff):
# Get a kernel density estimator for the isochrone matching the
# passed age and metallicity
    clustFullIso = IP.GetPARSECIsochrone(clustAge, clustMetal)
    if starTeff > 5100:   # MS/turn-off?
        # Convert LogTeff to Teff, and only take data Sub-Giant stage and earlier
        TLogGPoints = np.array([[10**x[5], x[6]] for x in clustFullIso if x[17] < 3])
    else: # (Sub-) giant
         # Convert LogTeff to Teff, and only take data Sub-Giant stage and later
        TLogGPoints = np.array([[10**x[5], x[6]] for x in clustFullIso if x[17] < 5 and x[17] > 2 ])
    
    TLogGIso = IP.EvenSpacedIso(TLogGPoints, 500)
    probKDE = IP.BuildProbArray(TLogGIso)
    
#    # Display our KDE
#    T_flat = np.r_[TLogGIso[:,0].min():TLogGIso[:,0].max():200j]
#    LogG_flat = np.r_[TLogGIso[:,1].min():TLogGIso[:,1].max():200j]
#    
#    t,g = np.meshgrid(T_flat, LogG_flat)
#    gridCoords = np.append(t.reshape(-1,1), g.reshape(-1,1), axis=1)
#    Z = probKDE(gridCoords.T)
#    Z = Z.reshape(200,200)
#    
#    PyP.scatter(TLogGIso[:,0], TLogGIso[:,1])
#    PyP.show()
#    PyP.imshow(Z)
#    PyP.show()
        
    return probKDE
    
def evaluatePoints(paramScores, metallicities, guesses, currentM):
# Compare the "probabilities"* of the passed list of guesses, and return the
# best option. If there is a score for all three probabilities, we simply
# select the option with the highest score (product).
# In the case where one or more are zero, we examine the sub-scores, and 
# select Teff, LogG, based on best abundance vs. excitation potential, and
# x.0 and x.1 ionization state abundance comparison. If _both_ of those
# abundances are zero, we select the result with the closest metallicity.
#
# The paramScores triplet tuple is of the form: 
#                   (isoProb, slopeProb, balanceProb, photometricProb)
#

# *: See the discussion above in re: P-scores vs. probabilities.
    nonZeroScores = [[x if x>10**-10 else 10**-10 for x in i] for i in paramScores]

#    scores = [np.product(i)*len(i) if len(i)>0 else 0. for i in nonZeroScores]
    scores = [np.product(i) for i in nonZeroScores]
    
    bestStep = np.argmax(scores)
    if scores[bestStep] == 0.:
    # All of the probabilities on every step are zero, so choose based on
    # P-scores, and isochrone/photometric probabilities.
        bestStep = -1
        
        # If the ionization balance (or ab/ex. pot. slope) score is  
        # zero, we still want to have a score, but multiply in a 
        # "minimization" factor to give priority to scores with a non-zero
        # (but still low) value.
        minimization = 10e-5
        
        # The ab/ex.pot. slope is most affected by the Teff
        bestTStep = np.argmax([i[1] for i in paramScores])
        
        # Must be better by at least 1%
        if paramScores[bestTStep][1] == 0. or abs(paramScores[bestTStep][1]-paramScores[0][1])/paramScores[bestTStep][1] < 0.01:
            # Our best step wasn't great (or was still 0), so choose the best of the rest.
            # We'll want to compare the best of scores with, and without ionization
            # balance P-scores, but give priority to ones with.
            bestGStep = np.argmax([i[2] for i in paramScores])
            # Using ionization balance for logG
            if paramScores[bestGStep][2] == 0. or abs(paramScores[bestGStep][2]-paramScores[0][2])/paramScores[bestGStep][2] < 0.01:
                # Our best G step wasn't great (or was still 0), so choose the best of the rest
                TScoreArray = [i[0]*i[2]*i[3] if i[2]>0. else i[0]*i[3]*minimization for i in paramScores]
                GScoreArray = [i[0]*i[1]*i[3] if i[1]>0. else i[0]*i[3]*minimization  for i in paramScores]
                if max(TScoreArray) >= max(GScoreArray):
                    if max(TScoreArray) == 0.:
                        # We're seeing a bunch of "0.00" scores, so pick the
                        # one with the metallicity closest to the current guess.
                        bestStep = np.argmin([abs(x-currentM) for x in metallicities])
                        # If our metallicity really isn't going to change much, just
                        # return our current guess, and hope that we have more line
                        # selection or other criteria to make the decision at a 
                        # later time.
                        if abs(metallicities[bestStep] - currentM) < k.MetStep*2:
                            bestStep = 0
                    else:
                        bestStep = np.argmax(TScoreArray)
                else:
                    bestStep = np.argmax(GScoreArray)
            else:
                bestStep = bestGStep
        else:
            # We have a good "slope" P-score - Go with it
            bestStep = bestTStep

    (bestT, bestG, bestV) = guesses[bestStep]
    bestM = metallicities[bestStep]

    return bestT, bestG, bestV, bestM, bestStep


def determineParms(guessTuple, kuruczAtms, pradks, mass, clustMetal = 0.0, clustAge = 4500, logFile=None, starname=None, photoGuess = True, refDict=None):
# External entry point to determine the atmospheric parameters (Teff,
# LogG, VTurb, [Fe/H]) for a star, given a set of absorption line
# measures, and a starting "guess" - usually made
# from photometry. If the photoGuess flag is set, we operate on
# the assumption that the passed guess is based on tangible
# data, and will use a Gaussian distribution around the Teff
# value (sigma=200K) as a (quasi-) Baysian prior. Otherwise
# we just use the guess as a starting point and have no
# prior for the Teff value.
#
# Params are determined for the star passed, either by name 
# (In the form "NGC-XXXX IIII-YYYY"), or by a "log" file, containing
# the lines for the star.
    # Needless to say, we need at least one of them to make a determination
    assert logFile != None or starname != None
    
    # Testing parameter visibility error
    paramScoreList = []
    metallicityList = np.array([])
    
    (tGuess, logGGuess, vTurbGuess) = guessTuple
    
    bestT = tGuess
    bestG = logGGuess
    bestV = vTurbGuess
    bestM = clustMetal
    
    # We're going to give some credit to the original guess. Let's try
    # a normal distribution around the guess with a sigma of 200K
    if photoGuess:
        photoDist = stats.norm(tGuess,k.AtmTGuessSigma)
        photoGuessPDF = photoDist.pdf(tGuess) # For normalization
    else:
        vTurbDist = stats.gamma(2.5)
        vTurbDistMax = max([vTurbDist.pdf(i) for i in np.arange(0., 6., 0.04)])
#        photoProb = 1.0
        
    probKDEs = {}
    probKDE = getProbKDE(clustMetal, clustAge, tGuess)
    probKDEs[clustMetal] = probKDE
    
    # Something (probably MOOG - FORTRAN is dumb!) is choking on long file path names
    if logFile!=None and len(logFile) > 32:
        shortLogFilename = k.MOOGTempLogName
        if os.path.isfile(shortLogFilename):
            if not os.path.samefile(shortLogFilename, logFile):
                os.remove(shortLogFilename)
        else:
            shutil.copy2(logFile, shortLogFilename)
    elif logFile == None:
    # No log file, so look up the lines in the DB, and make one.
        clusterID = starname.split()[0]
        starID = starname.split()[1]
        lines = LL.getLinesForStar(clusterID, starID, elements=[22.0, 22.1, 26.0, 26.1], dataFormat='MOOG')
        shortLogFilename = k.MOOGTempLogName
        with open(shortLogFilename, 'w') as logFile:
            fileHeader = '# Temporary line log file for {0} {1}\n'.format(clusterID, starID)
            logFile.write(fileHeader)
            for line in lines:
                logFile.write(line+'\n')
    else:
        shortLogFilename = logFile
                
    # Copy the list by slicing, since we will be using "pop"
    varLimitList = k.EQWAbVariances[:]
    
    # Any line with a starting abundance two orders off from the rest
    # of the lines needs to be tossed immediately.
    #currentVarLimit = 2.0
    # Initial read and prune of the line log
    tempFileHead = 'temp_{0:05d}_{1:06d}'.format(os.getpid(),datetime.datetime.now().microsecond)
    modelFilename=tempFileHead+'.m'
    
    with open("mylogFile.txt", 'a') as myLogFile:
        myLogFile.write('New star: {0} {1}\n'.format(clusterID, starID))
        for currentVarLimit in varLimitList:
            # Calculate starting abundances, based on the current model
            myLogFile.write("    Starting new abundance limit: % 1.3f\n" % currentVarLimit)

            thisModel, temppradk = mk.MakeAtmModel(kuruczAtms, pradks, bestT, bestG, bestM, bestV, mass)
            MI.WriteMOOGAtmModel(thisModel, modelFilename, bestT, bestG, bestM, bestV)
#            mk.WriteMOOGModel(thisModel, modelFilename, bestT, bestG, bestM, bestV)
            lines = MI.GetMOOGAbs(shortLogFilename, modelFilename)
            FeILines, FeIILines, TiILines, TiIILines = selectLines(lines, kuruczAtms, pradks, (bestT, bestG), varLimit=currentVarLimit, enforceMin=True, refDict=refDict)
            
            if len(FeILines)+len(FeIILines)+len(TiILines)+len(TiIILines) <1:
                myLogFile.write("    Excessive abundance {1:1.3f} filtering: {0:3d} lines.\n".format(len(lines), currentVarLimit))
                break
                #currentVarLimit = 0.
                #continue

            # Clean up old temp files
            u.clearFiles(glob.glob(tempFileHead[:-7]+'*'))

            # Once we pare the line list, don't deal with the other lines, anymore
            linesToLog = getMOOGLisLines(FeILines, FeIILines, TiILines, TiIILines)
            
            headStr = '# Line width measures pared from:{0} on: '.format(shortLogFilename)+datetime.datetime.now().isoformat()
            # Had an error in the output at one point. This code just left in
            # "just in case".
            try:
                np.savetxt(fname=shortLogFilename, X=linesToLog, fmt=k.MOOGLogFormat, header=headStr)
            except ValueError:
                print ("--------- {0:1.3f} --------------".format(currentVarLimit))
                print (linesToLog)
                print ("Error in lines parameter")
                exit(0)
                
            # Take big steps first, then try smaller ones
            for stepSize in [10., 5., 2., 1.]:
                myLogFile.write('        New Step: T:+/-{0:3.0f};  LogG: +/-{1:1.3f}; VTurb:+/-{2:1.3f}\n'.format(stepSize*k.TeffStep, stepSize*k.LogGStep, stepSize*k.VTurbStep))
                slopeProb, balanceProb = evaluateParams(FeILines, FeIILines, TiILines, TiIILines)
                
                # Need to verify that our isochrone is at least close to our 
                # current guess' metallicity.
                isoMets = probKDEs.keys()
                closestIsoMet = min(range(len(isoMets)), key=lambda i: abs(isoMets[i]-bestM))
                probKDE = probKDEs[isoMets[closestIsoMet]]
                
                if abs(isoMets[closestIsoMet]-bestM) > 0.50:
                # Need a new isochrone:
                    if bestM > k.MaxParsecMetal:
                        if k.MaxParsecMetal in isoMets:
                            probKDE = probKDEs[k.MaxParsecMetal]
                        else:
                            myLogFile.write('***** New Isochrone: T:{0:3.0f}; M:{1:+1.3f}({2:+1.3f})\n'.format(bestT, bestM, k.MaxParsecMetal))
                            probKDE = getProbKDE(k.MaxParsecMetal, clustAge, bestT)
                            probKDEs[k.MaxParsecMetal] = probKDE
                    else:
                        myLogFile.write('***** New Isochrone: T:{0:3.0f}; M:{1:+1.3f}\n'.format(bestT, bestM))
                        probKDE = getProbKDE(bestM, clustAge, bestT)
                        probKDEs[bestM] = probKDE
                        
                isoProb = IP.GetProb((bestT, bestG), probKDE)
                if photoGuess:
                    photoProb = photoDist.pdf(bestT)/photoGuessPDF
                else: 
                    photoProb = (vTurbDist.pdf(bestV)/vTurbDistMax)**3
                
                scores = [slopeProb, balanceProb, isoProb, photoProb]
                myLogFile.write('        Current guess: T:{0:4.0f} - G:{1:1.3f} - V:{2:1.3f} - M:{3:+1.3f}: {4:1.4e} * {5:1.4e} * {6:1.3e} * {7:1.3e} = {8:1.3e}\n'.format(bestT, bestG, bestV, bestM, slopeProb, balanceProb, isoProb, photoProb, np.product([x for x in scores if not x==0.])))
                
                bestStep = -1
                lastScore = 0.
                lastParms = (bestT, bestG, bestV, bestM)
                while bestStep != 0:
                   
                    # When we first enter, the best step is always the first, signified by 
                    # a "bestStep" value of "-1", which we need to change to "0". We can't
                    # start the loop with the actual index, for obvious reasons.
                    if bestStep < 0: bestStep = 0
                    paramScoreList = []
                    metallicityList = np.array([])
                    
                    # Get the scores for our current guess and 
                    # the seven potential steps around it.
                    guessArray = np.array([(bestT+stepSize*item[0], bestG+stepSize*item[1], bestV+stepSize*item[2]) for item in k.ParmStepArray])
                    
                    for (thisT, thisG, thisV) in guessArray:
                        # New step in the range of our models?
                        if not mk.parms_in_range(teff=thisT, logg=thisG, vturb=thisV):
                            myLogFile.write('    *** Parmeter out of range: Continuing.\n')
                            paramScoreList.append((-1., -1., -1.))
                            metallicityList = np.append(metallicityList, thisM)
                            continue

                        tempFileHead = 'temp_{0:05d}_{1:06d}'.format(os.getpid(),datetime.datetime.now().microsecond)
                        modelFilename=tempFileHead+'.m'
                        thisModel, temppradk = mk.MakeAtmModel(kuruczAtms, pradks, thisT, thisG, bestM, thisV)    
                        MI.WriteMOOGAtmModel(thisModel, modelFilename, thisT, thisG, bestM, thisV)
                        # We re-read (but do not re-prune!) the lines for every guess, since abundances change
                        # Note: the data format is different from above.
                        lines = MI.GetMOOGAbs(shortLogFilename, modelFilename)
                        FeILines, FeIILines, TiILines, TiIILines = selectLines(lines, kuruczAtms, pradks, (thisT, thisG), varLimit=currentVarLimit, enforceMin=True, refDict=refDict)
                        
                        # Clean up our temp files
                        u.clearFiles(glob.glob(tempFileHead[:-7]+'*'))
                        if len(FeILines) < 7:
                            if len(FeILines) == 0 or not u.is_list(FeILines):
                            # Most frequently, we will have zero FeI lines as a
                            # result of a failed MOOG step. We'll just record a
                            # "bad" probability and continue with our guesses.
                                myLogFile.write('        *** 0/1 FeI Lines: Continuing.\n')
                                paramScoreList.append((-1., -1, -1.))
                                metallicityList = np.append(metallicityList, thisM)
                                continue
                            elif u.is_list(FeILines):
                                myLogFile.write('        *** {0:1d} FeI Lines: FeI[0]:{1}\n'.format(len(FeILines), repr(FeILines[0])))
                            else:
                                myLogFile.write('        *** {0:1d} FeI Lines: Continuing.\n'.format(len(FeILines)))
                                
                        # We were getting a bug where the FeILines list 
                        # (of lists) only had one member, and wasn't showing as
                        # a list of a single list of lines. Rather it was just
                        # showing up as a single list of line parameters. This
                        # ugliness was desinged to figure out what was going on
                        # and can probably be removed...eventually, once I am
                        # convinced that the original issue is fixed.
                        try:
                            thisM = np.mean(np.array(FeILines)[:,5])-k.SolarFeH
                        except (IndexError, TypeError):
                            myLogFile.write('    Error, in FeILines indexing\n')
                            myLogFile.write('    --------------------\n')
                            for line in FeILines:
                                print (line)
                            print('--------------------')
                            print(FeIILines)
                            print('--------------------')
                            print(TiILines)
                            print('--------------------')
                            print(TiIILines)
                            #raw_input('Hit return to continue')
                            continue
                        # If the new params create a different metallicity value, recalculate
                        # with a "compromise" value
                        if abs(thisM-bestM) > 0.03:
                            newM = bestM + (thisM-bestM)/2
                            tempFileHead = 'temp_{0:05d}_{1:06d}'.format(os.getpid(),datetime.datetime.now().microsecond)
                            modelFilename=tempFileHead+'.m'
                            thisModel, temppradk = mk.MakeAtmModel(kuruczAtms, pradks, thisT, thisG, newM, thisV)    
                            MI.WriteMOOGAtmModel(thisModel, modelFilename, thisT, thisG, newM, thisV)
#                            mk.WriteMOOGModel(thisModel, modelFilename, thisT, thisG, newM, thisV)            
                            lines = MI.GetMOOGAbs(shortLogFilename, modelFilename)
                     
                            u.clearFiles(glob.glob(tempFileHead[:-7]+'*'))

                            FeILines, FeIILines, TiILines, TiIILines = selectLines(lines, kuruczAtms, pradks, (thisT, thisG), varLimit=currentVarLimit, enforceMin=True, refDict=refDict)
                                                    
                        slopeProb, balanceProb = evaluateParams(FeILines, FeIILines, TiILines, TiIILines)
                        isoProb = IP.GetProb((thisT, thisG), probKDE)
                        
                        if photoGuess:
                            photoProb = photoDist.pdf(thisT)/photoGuessPDF
                        else: 
                            photoProb = (vTurbDist.pdf(thisV)/vTurbDistMax)**3
                        
                        scores = [slopeProb, balanceProb, isoProb, photoProb]
                        
                        myLogFile.write('             Evaluating: T:{0:4.0f} - G:{1:1.3f} - V:{2:1.3f} - M:{3:+1.3f}: {4:1.4e} * {5:1.4e} * {6:1.3e} * {7:1.3e} = {8:1.3e}\n'.format(thisT, thisG, thisV, thisM, slopeProb, balanceProb, isoProb, photoProb, np.product([x for x in scores if not x==0.])))
                        # Place a tuple of "probabilities" (insert P-score discussion, here)
                        # for later evaluation
                        paramScoreList.append((isoProb, slopeProb, balanceProb, photoProb))
                        
                        metallicityList = np.append(metallicityList, thisM)
                    
                    # Select the guess with the highest score. 
                    # If it's our current point, then the loop will continue, 
                    # tighten the abundance variance restriction, and re-measure.
                    newT, newG, newV, newM, newStep = evaluatePoints(paramScoreList, metallicityList, guessArray, bestM)
                    if (newT, newG, newV, lastParms[3]) == lastParms:
                        # We're oscillating between two steps. Don't
                        bestStep = 0
                    else:
                        lastParms = (bestT, bestG, bestV, bestM)
                        bestT, bestG, bestV, bestM, bestStep = newT, newG, newV, newM, newStep
                        
                    try:
                        time.sleep(0.1)
                        myLogFile.write('        Best guess: T:{0:4.0f} - G:{1:1.3f} - V:{2:1.3f} M:{3:+1.3f}\n'.format(bestT, bestG, bestV, bestM))
                        time.sleep(0.1)
                        myLogFile.flush()
                    except IOError:
                        # We can actually try flushing/writing before the last write has
                        # completed. If this is the case, don't sweat it. We'll just get
                        # it next time.
                        # just "Pass-ing" doesn't seem to be ignoring the exception, so...
                        print("Error - Handling")
#                        pass

                    # If our current point has re-evaluated to be a lower score
                    # than our previous, go back!
                    newScore = paramScoreList[newStep][0]*paramScoreList[newStep][1]*paramScoreList[newStep][2]*paramScoreList[newStep][3]
                    if newScore <= lastScore and newScore > 0.:
                    # Restore the last, best parms, exit this iteration and
                    # either force the next paring, or end the run
                        (bestT, bestG, bestV, bestM) = lastParms
                        try:
                            myLogFile.write('        Going back!: T:{0:4.0f} - G:{1:1.3f} - V:{2:1.3f} - M:{3:+1.3f} - ({4:1.3e} - {5:1.3e})\n'.format(bestT, bestG, bestV, bestM, newScore, lastScore))
                        except IOError:
                            time.sleep(0.5)
                            myLogFile.flush()
                            myLogFile.write('        Going back!: T:{0:4.0f} - G:{1:1.3f} - V:{2:1.3f} - M:{3:+1.3f} - ({4:1.3e} - {5:1.3e})\n'.format(bestT, bestG, bestV, bestM, newScore, lastScore))
                            pass
                        bestStep = 0
                        lastScore = 0.
                        continue
                    else:
                        time.sleep(0.5)
                        lastScore = newScore
                        
             
             # This worked most of the time, but occasionally would spawn an IOError, so changed it to a for loop
 #           # Get the next Ab. limit, and carry on.
 #           if currentVarLimit > 0.0: 
 #               currentVarLimit = varLimitList.pop(0)
       
    return bestT, bestG, bestV, bestM
        
# Command line calls for testing only - Broken now, anyway.
'''
if __name__ == '__main__':
    moduleName = 'SpecParamDet.py'
    resultsFilename = 'StellarParams_'+datetime.datetime.now().strftime('%Y%m%d_%H%M%S')+'.txt'
    # Parse your command line call
    temp = os.sys.argv
    argv = deque(temp)
    
    inputDirectory = ''
    guessFile = None
    
    while len(argv) > 0:
        flag = argv.popleft()
        if flag == '-v':
            g.VerboseMode = True
            print('Verbose Mode enabled.')
        elif flag.split('/')[-1] == moduleName:
            pass
        else:
            if inputDirectory == '':
                inputDirectory = flag
            else:
                guessFile = flag
                
    # We'll be running multiple passes, but only want to load models once. Otherwise,
    # we could just call mk.MakeModelFile.
    # To implement giant atmosphere calculation, we would need to use the 
    # "spherical" MARCS models. In case it's not obvious: d... = dwarf model file, 
    # g... = giant model file.
    u.errorOut('Reading plane-parallel model files:\n')
#    dModelPath = '/netdisks/galileo18/mikelum/Dropbox/CodeCloud/MyTools/ClusterAnalysis/Models/MARCS_PlaneParallel'
    dModelPath = '/home/mikelum/Dropbox/CodeCloud/MyTools/ClusterAnalysis/Models/MARCS_PlaneParallel'
    dModelFiles = mk.findDataFiles(dModelPath)
    dKuruczAtms, dPradks = mk.LoadModels(dModelFiles)
    dRefCorrD, dRefLineD, dWeightD = RC.getSolarCorrections(ionList=[26.0, 26.1, 22.0, 22.1], modelAtms=dKuruczAtms, pradks=dPradks)
    dRefDict = {}
    for ion in dRefCorrD.keys():
        dRefDict[ion] = (dRefLineD[ion], dRefCorrD[ion], dWeightD[ion])
    
    u.errorOut('Reading spherical model files:\n')
#    gModelPath = '/netdisks/galileo18/mikelum/Dropbox/CodeCloud/MyTools/ClusterAnalysis/Models/MARCS_Spherical'
    gModelPath = '/home/mikelum/Dropbox/CodeCloud/MyTools/ClusterAnalysis/Models/MARCS_Spherical'
    gModelFiles = mk.findDataFiles(gModelPath)
    gKuruczAtms, gPradks = mk.LoadModels(gModelFiles)
    gRefCorrD, gRefLineD, gWeightD = RC.getGiantCorrections(ionList=[26.0, 26.1, 22.0, 22.1], modelAtms=gKuruczAtms, pradks=gPradks)
    gRefDict = {}
    for ion in gRefCorrD.keys():
        gRefDict[ion] = (gRefLineD[ion], gRefCorrD[ion], gWeightD[ion])
        
    starList = [\
        ['NGC-0752 PLA-0300', 5578., 4.25, 0.76],\
        ['NGC-0752 PLA-0350', 4848., 2.58, 2.83],\
        ['NGC-0752 PLA-0356', 4880., 2.65, 2.76],\
        ['NGC-0752 PLA-0361', 5464., 4.32, 0.82],\
        ['NGC-0752 PLA-0391', 5437., 4.55, 0.80],\
        ['NGC-0752 PLA-0413', 6220., 4.37, 1.67],\
        ['NGC-0752 PLA-0429', 5124., 4.62, 0.80],\
        ['NGC-0752 PLA-0475', 5918., 4.50, 1.18],\
        ['NGC-0752 PLA-0506', 4898., 2.65, 2.77],\
        ['NGC-0752 PLA-0520', 5944., 4.28, 1.42],\
        ['NGC-0752 PLA-0552', 5977., 4.50, 1.23],\
        ['NGC-0752 PLA-0575', 5373., 4.57, 0.80],\
        ['NGC-0752 PLA-0687', 4833., 2.62, 1.36],\
        ['NGC-0752 PLA-0699', 5745., 4.40, 1.12],\
        ['NGC-0752 PLA-0701', 5655., 4.39, 0.94],\
        ['NGC-0752 PLA-0786', 5600., 4.40, 0.80],\
        ['NGC-0752 PLA-0790', 6199., 4.08, 1.61],\
        ['NGC-0752 PLA-0791', 6114., 4.15, 1.55],\
        ['NGC-0752 PLA-0828', 5243., 4.54, 0.59],\
        ['NGC-0752 PLA-0859', 5429., 4.37, 1.01],\
        ['NGC-0752 PLA-0864', 6068., 4.29, 1.39],\
        ['NGC-0752 PLA-0889', 5815., 4.23, 1.51],\
        ['NGC-0752 PLA-0921', 6092., 4.36, 1.50],\
        ['NGC-0752 PLA-0964', 6042., 4.30, 1.24],\
        ['NGC-0752 PLA-0993', 5546., 4.54, 0.87],\
        ['NGC-0752 PLA-0999', 5644., 4.52, 1.01],\
        ['NGC-0752 PLA-1000', 6448., 4.10, 1.50],\
        ['NGC-0752 PLA-1012', 6154., 4.19, 1.55],\
        ['NGC-0752 PLA-1017', 5753., 4.51, 1.05],\
        ['NGC-0752 PLA-1089', 4991., 2.85, 1.47],\
        ['NGC-0752 PLA-1107', 5619., 4.51, 1.05],\
        ['NGC-0752 PLA-1172', 4820., 2.70, 1.45],\
        ['NGC-0752 PLA-1270', 5437., 4.54, 0.90],\
        ['NGC-0752 PLA-1365', 5591., 4.54, 0.86],\
        ]
    with open(resultsFilename,'w') as resultFile:
#        for logFile in glob.glob(inputDirectory+'/*.log'):
#            # Check for dwarf vs. giant:
#            tGuess, logGGuess, vTurbGuess = getGuesses(logFile, guessFile)
        for starData in starList:
            name = starData[0]
            tGuess = starData[1]
            logGGuess = starData[2]
            vTurbGuess = starData[3]
            if logGGuess < k.dwarfGiantLimit:
                kuruczAtms = gKuruczAtms
                pradks = gPradks
                mass = 1.15
                refDict = gRefDict
            else:
                kuruczAtms = dKuruczAtms
                pradks = dPradks
                mass = 0.0
                refDict = dRefDict
            tResult, gResult, vResult, mResult = determineParms((tGuess, logGGuess, vTurbGuess), kuruczAtms, pradks, mass, clustMetal = -0.05, clustAge = 1350,  logFile=None, starname=name, refDict=refDict)
#            tResult, gResult, vResult, mResult = determineParms((tGuess, logGGuess, vTurbGuess), kuruczAtms, pradks, mass, clustMetal = -0.05, clustAge = 1350,  logFile=logFile, starname=None)
#            starName = (logFile.split('/')[-1])[:-4]
            resultFile.write('{0:15s} {1:4.0f} {2:1.3f} {3:1.3f} {4:1.3f}\n'.format(name,tResult, gResult, vResult, mResult))
            resultFile.flush()
    exit()
'''      
