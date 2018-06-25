#! /usr/bin/env python
# 
# Module: SpecParmPlots.py
#
# Author: Mike Lum
#
# Description: Three routines to create <parameter> vs. <spectroscopic parameter>
#              as an aid in determining stellar Teff, LogG, and Vturb.
#              The assumption is that certain characteristics will be minimized
#              when a parameter is "correct" for a star.
#
# Contents:
#   MakeVPlots: Plots the microturbulent velocity parameter against the variance
#               in the Ca I abundance measurements from different lines. 
#               The "correct" Vturb parameter is where the variance is
#               minimized.
#
#   MakeTeffPlots: Plots Teff vs. the difference between the abundance measured
#               by two ionization states of Ti and Fe. It's assumed that the
#               delta will be minimized (zero) at the correct temperature.
#               Note that the spectroscopic temperature determination is 
#               conventionally determined by minimizing the excitation potential
#               vs. measured [Fe/H] trend from Fe I lines.
#
#   MakeLogGPlots: Plots LogG vs. the difference between the abundance measured
#               by two ionization states of Ti and Fe. It's assumed that the
#               delta will be minimized (zero) at the correct LogG value.
#
# Revision History:
#    Date        Vers.    Author        Description
#    12-14-2017  1.0f1    Lum        Code separated from test/sample code directory
#
# To Do:
#    
#

# Imports
# Astro- and 3rd-party Library Imports
import numpy as np
from matplotlib import pyplot
import scipy as sp

# Local Imports
from Abundances import Abundance as AB
from Abundances import ReferenceCorrections as RC
from constants import constants as k
from Models import MakeAtmModel as mk
from Output import STRPlots as STP
from utilities import elements as el
from utilities import utilities as u

# Functions

# We can use abundance vs. parameter plots to refine atmospheric parameters
# (Teff, logG, V_turb), which are originally determined from photometric
# values.
def MakeVPlots(clusterName='NGC-0752', starDataList=None, fileTag='', 
               filterBlends=False, referenceCorrect=False, useDAOSpec=False, ions=[20.0]):
# Function plots V_turb vs. sigma [Ca/H] (sigma is the line-to-line
# standard deviation of the measured Ca lines), for a range of V_turb values. 
# Intended to provide a spectroscopic confirmation/adjustment
# for photometrically-determined parameters.
    if starDataList is None:
        starDataList = STP.GetAllStarParms(clusterName=clusterName)
        
    dModelFiles = mk.findDataFiles(k.DwarfModelPath)
    dModelAtms, dPradks = mk.LoadModels(dModelFiles)
    dCorrectDict, dReferenceLines, dLineWeights = \
            RC.GetDwarfCorrections(ionList=ions, 
                                modelAtms=dModelAtms, 
                                pradks=dPradks)
    
    gModelFiles = mk.findDataFiles(k.GiantModelPath)
    gModelAtms, gPradks = mk.LoadModels(gModelFiles)
    # Note: We're using dwarf corrections for giant stars...
    # because they're better than giant corrections for the 
    # giant stars. This needs to be fixed. :(
    gCorrectDict, gReferenceLines, gLineWeights = \
            RC.GetSolarCorrections(ionList=ions, 
                                modelAtms=dModelAtms, 
                                pradks=dPradks)    
# For now, CaI appears to work best for all stars. 
# VI, Ti I/II, Fe I/II, Cr I/II also work, but not as well - (See Reddy, et al. 2012).
#    ions = [20.0]
    abDict = {}
    
    colorArray=['r','orange','gold','darkgreen','navy',\
                'm','saddlebrown','skyblue','hotpink']    
    
    for starParms in starDataList:
        if RC.isGiantStar(starParms[1:]):
            modelAtms = gModelAtms
            pradks = gPradks
            refLines = gReferenceLines
            refDict = gCorrectDict
            lineWeights = gLineWeights
        else:
            modelAtms = dModelAtms
            pradks = dPradks
            refLines = dReferenceLines
            refDict = dCorrectDict
            lineWeights = dLineWeights
            
        starPoints = {}
        for v in np.linspace(0.5, 3.0, 26):
            parmTuple = (starParms[0], starParms[1], starParms[2], v, \
                         starParms[4], starParms[5], starParms[6])
                         
            abDict[v] = STP.GetAbsForStar(parmTuple, ions=ions,\
                filterBlends=filterBlends, refCorrect=referenceCorrect,\
                refDict=refDict, refLines=refLines, lineWeights=lineWeights,\
                useDAOSpec=useDAOSpec, modelAtms=modelAtms, pradks=pradks)
                
            for ion in abDict[v].keys():
                if ion not in starPoints.keys():
                    starPoints[ion]={}
                if abDict[v][ion][0] == 0.:
                    continue
                starPoints[ion][v] = abDict[v][ion][1]

        colorIdx=0
        
        ymaxes = []
        for ion in starPoints.keys():
            Xs = sorted(np.array(list(starPoints[ion].keys())))
            if len(Xs)==0 or ion in STP.SynthAbsDict.keys():
                continue
            Ys = np.array([starPoints[ion][x] for x in Xs])
            if len(Xs) < 2 or len(Ys) < 2:
                continue
            Ys = Ys/min(Ys)
            ymaxes.append(max(Ys))
            try:
                minV = Xs[np.where(Ys==min(Ys))[0][0]]
            except IndexError:
                continue
            pyplot.plot(Xs, Ys, label=r'$({1}) V_{{turb}}={0:3.1f}km/sec$'.format(minV,el.getIonName(ion)), color=colorArray[colorIdx])
            pyplot.scatter(Xs, Ys, color=colorArray[colorIdx])
            pyplot.axvline(minV, linestyle=':', color=colorArray[colorIdx])
            colorIdx += 1
            if colorIdx == len(colorArray): colorIdx=0
        ax = pyplot.gca()
        ax.set_xlabel(r'$V_{{turb}} km/s$')
        ax.set_ylabel('Normalized [X/H] variance')
        ax.set_ylim((0., min(5.0, max(ymaxes))))
        pyplot.legend(fontsize=8)
        pyplot.savefig(k.ParmPlotDir+'Vturb/'+starParms[0]+fileTag+'.png')
        pyplot.close()


def MakeLogGPlots(clusterName='NGC-0752', starDataList=None, 
                  fileTag='', filterBlends=False, referenceCorrect=False, 
                  useDAOSpec=False):
# Function plots logG vs. [TiI/TiII] for a range of LogG values
# Intended to provide a spectroscopic confirmation/adjustment
# for photometrically-determined parameters.

    if starDataList is None:
        starDataList = STP.GetAllStarParms(clusterName=clusterName)
        
    dModelFiles = mk.findDataFiles(k.DwarfModelPath)
    dModelAtms, dPradks = mk.LoadModels(dModelFiles)
    gModelFiles = mk.findDataFiles(k.GiantModelPath)
    gModelAtms, gPradks = mk.LoadModels(gModelFiles)
 # Want to watch TiI/TiII
    ions = [22.0, 22.1]
    abDict = {}
    dLogGRange = np.linspace(3.75, 5.00, 26)
    gLogGRange = np.linspace(1.75, 3.00, 26)

    for s in starDataList:
        if u.isGiant(s[1], s[2]):
            logGRange = gLogGRange
        else:
            logGRange = dLogGRange
        
        abDeltaPts = []
        logGPts = []
        for logG in logGRange:
            tempParms = [(s[0], s[1], logG, s[3], s[4], s[5], s[6])]
            abDict = STP.GetAbTable(clusterName=clusterName,
                                      starDataList=tempParms, ions=ions, 
                                      filterBlends=filterBlends, 
                                      referenceCorrect=referenceCorrect, 
                                      useDAOSpec=useDAOSpec, 
                                      gModelAtms=gModelAtms, gPradks=gPradks, 
                                      dPradks=dPradks, dModelAtms=dModelAtms)
            star = list(abDict.keys())[0]
            
            # TiI/TiII points:
            if 22.0 in abDict[star].keys() and \
               22.1 in abDict[star].keys():
                logGPts.append(logG)
                abDeltaPts.append(abDict[star][22.0][0] - abDict[star][22.1][0])

        MakeOneLogGPlot(logGPts, abDeltaPts, starName=s[0], fileTag=fileTag)


def MakeOneLogGPlot(logGPts, abPts, starName='', fileTag=''):
# Draw a logG plot for a single star.
    # We want to use a couple of np.array tricks, so...
    if len(logGPts)<2 or len(abPts)<2: return
    Xs = np.array(logGPts)
    Ys = np.array(abPts)
    
    ax = pyplot.gca()
    pyplot.scatter(Xs, Ys)
    pyplot.axhline(0.0,linestyle=':')
    interpGs = np.interp([0.1,0.,-0.1],-Ys,Xs)
    pyplot.plot(Xs, Ys, 
                label=r'${0}({1:4.2f}^{{+{2:4.2f}}}_{{-{3:4.2f}}}$)'.\
                format(el.getIonName(22.0),interpGs[1],\
                       interpGs[1]-interpGs[2],\
                       interpGs[0]-interpGs[1]))
    ax.set_xlabel(r'$Log(G)$')
    ax.set_ylabel('[TiI/TiII]')
    pyplot.legend()
    pyplot.savefig(k.ParmPlotDir+'Logg/'+starName+fileTag+'_logG.png')
    pyplot.close()


def MakeTeffPlots(clusterName='NGC-0752', starDataList=None, 
                 fileTag='', referenceCorrect=False):
# Function plots Teff vs. slope of (XP vs. [Fe/H]) for a range of temperatures.
# Intended to provide a spectroscopic confirmation/adjustment
# for photometrically-determined parameters. We assume that at the "correct"
# temperature, the slope of XP vs. [Fe/H] would be zero.

    if starDataList is None:
        starDataList = STP.GetAllStarParms(clusterName=clusterName)
        
    dModelFiles = mk.findDataFiles(k.DwarfModelPath)
    dModelAtms, dPradks = mk.LoadModels(dModelFiles)
    gModelFiles = mk.findDataFiles(k.GiantModelPath)
    gModelAtms, gPradks = mk.LoadModels(gModelFiles)
 # Map Fe I Ab vs. Ex pot.
    ions = [26.0]
    abDict = {}
    dTeffRange = np.linspace(5000, 7000, 81)
    gTeffRange = np.linspace(4000, 6000, 81)
    
    for star in starDataList:
        starName = star[0]

        isGiant = RC.isGiantStar(star[1:])
        if isGiant:
            modelAtms = gModelAtms
            pradks = gPradks
            tRange = gTeffRange
            if referenceCorrect:
                correctDict, referenceLines, lineWeights = \
                        RC.GetGiantCorrections(ionList=ions, 
                                               modelAtms=modelAtms, 
                                               pradks=pradks)
        else:
            modelAtms = dModelAtms
            pradks = dPradks
            tRange = dTeffRange
            if referenceCorrect:
                correctDict, referenceLines, lineWeights = \
                        RC.GetDwarfCorrections(ionList=ions, 
                                               modelAtms=modelAtms, 
                                               pradks=pradks)
        logG = star[2]
        vTurb = star[3]
        met = star[4]
        
        slopes = []
        for Teff in tRange:
            unusedDict, uncorrLines, unusedMin, unusedMax = \
                AB.CalcAbsAndLines(clusterName+' '+starName, 
                           tuple([Teff, star[2], star[3], star[4], star[5]]),\
                           ionList=ions, modelAtms=modelAtms, pradks=pradks)
            XPs = []
            Abs = []
            if referenceCorrect:
                adjLines, allLines = RC.SortAndFilterLines(uncorrLines[26.0],
                           26.0,
                           tuple([Teff, star[2], star[3], star[4], star[5]]),
                           solarCorrect=True, 
                           solarLines=referenceLines[26.0], 
                           solarCorrs=correctDict[26.0], 
                           lineWeights=lineWeights[26.0])
                if len(adjLines) > 10:
                # Assume 10 Fe I lines needed for nice plots
                    for line in adjLines:
                        # XP isn't returned, so have to do a lookup into the
                        # original list.
                        XPs.append(u.find(lambda l: l[0]==line[2],
                                                    uncorrLines[26.0])[1])
                        Abs.append(line[0])

            if len(XPs) == 0 or len(Abs)==0:
            # Either no corrections, or corrections failed
                XPs = [line[1] for line in uncorrLines[26.0]]
                Abs = [line[5] for line in uncorrLines[26.0]]

            slope, intercept, rVal, pVal, err = sp.stats.linregress(XPs, Abs)
            
            slopes.append(slope)
        
        # We have to flip because interp expects the xvalues to be 
        # monotomically increasing
        interpTs = np.interp([-0.02,0.,0.02],slopes[::-1],tRange[::-1])
        if interpTs[0]>interpTs[1]:
            minusR = interpTs[1]-interpTs[2]
            plusR = interpTs[0]-interpTs[1]
        else:
            minusR = interpTs[1]-interpTs[0]
            plusR = interpTs[2]-interpTs[1]
            
        fig = pyplot.figure()
        TeffLabel = r'$T_{{eff}} = {0:4.0f}^{{+{1:4.0f}}}_{{-{2:4.0f}}}$)'.\
                    format(interpTs[1],minusR,plusR)
        pyplot.scatter(tRange, slopes, label=TeffLabel)
        pyplot.axhline(0.0,linestyle=':')
        
        pyplot.legend()
        pyplot.savefig(k.ParmPlotDir+'Teff/'+star[0]+fileTag+'_FeSlope.png')
        pyplot.close()

def MakeTeffIonPlots(clusterName='NGC-0752', starDataList=None, 
                  fileTag='', filterBlends=False, referenceCorrect=False, 
                  useDAOSpec=False):
# Function plots Teff vs. [TiI/TiII] for a range of temperatures.
# Intended to provide a spectroscopic confirmation/adjustment
# for photometrically-determined parameters.

    if starDataList is None:
        starDataList = STP.GetAllStarParms(clusterName=clusterName)
        
    dModelFiles = mk.findDataFiles(k.DwarfModelPath)
    dModelAtms, dPradks = mk.LoadModels(dModelFiles)
    gModelFiles = mk.findDataFiles(k.GiantModelPath)
    gModelAtms, gPradks = mk.LoadModels(gModelFiles)
 # Want to watch TiI/TiII
    ions = [22.0, 22.1]
    abDict = {}
    dTeffRange = np.linspace(4500, 7000, 101)
    gTeffRange = np.linspace(3500, 6000, 101)

    for TeffIdx in range(len(dTeffRange)):
        tempParmList = []
        for s in starDataList:
            if s[1]<k.giantTeffLimit and s[2]<k.giantLogGLimit:
                Teff = gTeffRange[TeffIdx]
            else:
                Teff = dTeffRange[TeffIdx]
            tempParmList.append((s[0], Teff, s[2], s[3], s[4], s[5], s[6]))
        abDict[Teff] = STP.GetAbTable(clusterName=clusterName,
                                  starDataList=tempParmList, ions=ions, 
                                  filterBlends=filterBlends, 
                                  referenceCorrect=referenceCorrect, 
                                  useDAOSpec=useDAOSpec, 
                                  gModelAtms=gModelAtms, gPradks=gPradks, 
                                  dPradks=dPradks, dModelAtms=dModelAtms)
    starPoints = {}
    for Teff in abDict.keys():
        for star in abDict[Teff].keys():
            if star not in starPoints.keys():
                starPoints[star]={}
            # TiI/TiII points:
            if 22.0 in abDict[Teff][star].keys() and \
               22.1 in abDict[Teff][star].keys():
                if 22.0 not in starPoints[star].keys():
                    starPoints[star][22.0]={}
                starPoints[star][22.0][Teff] =\
                    abDict[Teff][star][22.0][0] - abDict[Teff][star][22.1][0]
    
    for star in starPoints.keys():
        fig = pyplot.figure()
        
        starData = [s for s in starDataList if s[0] == star]
        isGiant = False
        if len(starData)>0: isGiant = STP.isGiantStar(starData[0][1:])
        for ion in starPoints[star].keys():

            Xs = np.array(sorted(starPoints[star][ion].keys()))
            Ys = np.array([starPoints[star][ion][x] for x in Xs])
            if isGiant:
                Xs = Xs-1000.

            interpTs = np.interp([-0.1,0.,0.1],Ys,Xs)
            pyplot.scatter(Xs, Ys, 
                        label=r'${0}({1:4.0f}^{{+{2:3.0f}}}_{{-{3:3.0f}}}$)'.\
                        format(el.getIonName(ion),interpTs[1],\
                               interpTs[1]-interpTs[0],\
                               interpTs[2]-interpTs[1]), marker='.')
            pyplot.axhline(0.0,linestyle=':')
        ax = pyplot.gca()
        ax.set_xlabel(r'$T_{{eff}}$')
        ax.set_ylabel('[TiI/TiII]')
        pyplot.legend()
        pyplot.savefig(k.ParmPlotDir+'Teff/'+star+fileTag+'_TiRel.png')
        pyplot.close()
