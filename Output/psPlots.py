#! /usr/bin/env python
# 
# Module: psPlots.py
#
# Author: Mike Lum
#
# Description: Functions to create a variety of plots from the cluster data.
#               plots can be displayed on screen, saved as .ps files, or both.
#
#              These plots are mostly used to observe trends with line
#               characteristics and abundance measurements, in an effort to
#               analyze trends to improve our error correction. As our analysis
#               has improved, we don't use these plots as much. Therefore, 
#               the functions herein are to be considered:
#               "use at your own risk"
#
# Contents: 
#   Function xxx: (Description)
#
# Revision History:
#    Date        Vers.    Author        Description
#    2018-03-31  1.1f1    Lum           "Use at your own risk"
#    2016-07-12  1.0f1    Lum           First checked in
#
# To Do:
#    
#

# Imports
from matplotlib import pyplot
import numpy as np
import os
import scipy as sp
import sys

from Abundances import Abundance as AB
from Abundances import ReferenceCorrections as RC
from constants import constants as k
from constants import directories as dirs
from constants import ptoe_abs as PA
from Databases import StarDB as sdb
from Models import MakeAtmModel as mk
from utilities import elements as el
from utilities import utilities as u

# Functions
def GetCoGLimit(starData, ion=26.0, modelAtms = None, pradks = None):
# Returns the parameters for a line (using the same return parameters as 
# scipy.stats.linregress), representing the maximum limit for an EQW measurement
# assuming that the non-linear portion of the Curve of Growth (CoG) begins
# around a value of -4.60, as measured by Log(EQW/Wavelength).
# The return parameters are (in order): 
#                               slope, intercept, r-value, p-value, stderr
    starName = k.MaxDetFakeStarName
    clustName = k.CalibClusterName
    abdict, lines, unusedMins, unusedMaxes = AB.CalcAbsAndLines(clustName+' '+starName, starData, ionList=[ion], modelAtms = modelAtms, pradks = pradks)
    if ion not in lines.keys():
        # Hack for N
        if ion==7.0:
            return 0,0,0,0,0
        else:
            print('Ion:{0:2.1f} will not use CoG limits (possibly HFS).'.format(ion))
            return 0,0,0,0,0
    Xs = [line[5] for line in sorted(lines[ion], key=lambda l:l[0])]
    Ys = [el.STR(line[2], line[1], starData[0]) for line in sorted(lines[ion], key=lambda l:l[0])]
    
    return sp.stats.linregress(Xs, Ys)


def GetDetectionLimit(starData, ion=26.0, modelAtms = None, pradks = None):
# Returns the parameters for a line (using the same return parameters as 
# scipy.stats.linregress), representing the minimum detection limit for the
# passed ion (default: FeI). Measurements below, and to the left of the line 
# are below the calculated detection limit of the spectra, and should be 
# ignored.
# Since detection limit is affected by parameters, passing a stellar parameter
# tuple of the form: (Teff, LogG, Vturb, [Fe/H], Giant mass (=0.0 for dwarfs))
# is needed.
# The return parameters are (in order): 
#                               slope, intercept, r-value, p-value, stderr

    starName = k.MinDetFakeStarName
    clustName = k.CalibClusterName
    abdict, lines, unusedMin, unusedMax = AB.CalcAbsAndLines(clustName+' '+starName, starData, ionList=[ion], modelAtms = modelAtms, pradks = pradks)

    if ion not in lines.keys():
        # Hack for N
        if ion==7.0:
            return 0,0,0,0,0
        else:
            print('Ion:{0:2.1f} will not use detection limits (possibly HFS).'.format(ion))
            return 0,0,0,0,0
    Xs = [line[5] for line in sorted(lines[ion], key=lambda l:l[0])]
    Ys = [el.STR(line[2], line[1], starData[0]) for line in sorted(lines[ion], key=lambda l:l[0])]
    
    return sp.stats.linregress(Xs, Ys)        

def AbSTRPlot(starName, ion, points, redSet = np.array([]), greenSet =np.array([]), blueSet =np.array([]), lowLimit=None, hiLimit=None, fileTag='', plotTitle='', labelPoints=None):
    # For now, just bailing if no points were passed
    if len(points) == 0:
        return
        
    assert len(points)>0

    # Labels?
    blackLabels = []
    redLabels = []
    greenLabels = []
    blueLabels = []
# Plot the passed (np.array) point array in an external file
    abundances = points[:,0]
    STRs = points[:,1]
#    Xmax = max(abundances)
#    Xmin = min(abundances)
    Ymin = min(STRs)
    Ymax = max(STRs)
    
#    # Use an arbitrary Solar+1.0/-1.5 in [X/H] space for the plot limits
#    solarAB = el.atnoLU(int(ion))[3]
#    pyplot.axis([np.round(solarAB-0.5,1), np.round(solarAB+1.0,1), np.round(Ymin,1)-0.5, np.round(Ymax,1)+0.5])
    
    axes = pyplot.gca()
    
    if len(points) < 500:
        pointSize = 5
        if labelPoints is not None:
            # Must have one label for each point (even if it is '')
            if len(labelPoints) < len(abundances):
                # Must have labels for all the "black" (data) points
                labelPoints = None
            else:
            # Note that array slicing will automatically return empty members
            # for indices beyond the sliced array length
                blackLabels = labelPoints[:len(abundances)]
                redLabels = labelPoints[len(abundances):len(abundances)+len(redSet)]
                greenLabels = labelPoints[len(abundances)+len(redSet):len(abundances)+len(redSet)+len(greenSet)]
                blueLabels = labelPoints[len(abundances)+len(redSet)+len(greenSet):len(abundances)+len(redSet)+len(greenSet)+len(blueSet)]
    else:
        labelPoints = None # Seriously, don't label more than 30 points
        pointSize = 3
    pyplot.scatter(abundances, STRs, marker='.', s=pointSize+2)
    
    # Write in a (tiny) wavelength reference next to each point
    if len(blackLabels)>0:
        for coord in zip(abundances.tolist(), STRs.tolist(), blackLabels):
            # Place the label below and right of the point
            pyplot.text(coord[0]+0.02, coord[1]-0.03, coord[2], color='k', fontsize='xx-small')
    
    if len(redSet) > 0:
        if len(redSet) < 30:
            pointSize = 5
        else:
            pointSize = 3
        redAbs = redSet[:,0]
        redSTRs = redSet[:,1]
        Ymin = min(Ymin, min(redAbs))
        Ymax = max(Ymax, max(redSTRs))
        pyplot.scatter(redAbs, redSTRs, marker='x', color='red', s=pointSize)
        for i in range(len(redLabels)):
            pyplot.text(redSet[i,0]+0.02, redSet[i,1]-0.03, redLabels[i], color='r', fontsize='xx-small')
            
    if len(greenSet) > 0:
        if len(greenSet) < 30:
            pointSize = 5
        else:
            pointSize = 3
        greenAbs = greenSet[:,0]
        greenSTRs = greenSet[:,1]
        Ymin = min(Ymin, min(greenAbs))
        Ymax = max(Ymax, max(greenSTRs))
        pyplot.scatter(greenAbs, greenSTRs, marker='x', color='green', s=pointSize)
        for i in range(len(greenLabels)):
            pyplot.text(greenSet[i,0]+0.02, greenSet[i,1]-0.03, greenLabels[i], color='g', fontsize='xx-small')
    if len(blueSet) > 0:
        if len(blueSet) < 30:
            pointSize = 5
        else:
            pointSize = 3
        blueAbs = blueSet[:,0]
        blueSTRs = blueSet[:,1]
        Ymin = min(Ymin, min(blueAbs))
        Ymax = max(Ymax, max(blueSTRs))
        pyplot.scatter(blueAbs, blueSTRs, marker='x', color='blue', s=pointSize)
        for i in range(len(blueLabels)):
            pyplot.text(blueSet[i,0]+0.02, blueSet[i,1]-0.03, blueLabels[i], color='b', fontsize='xx-small')

# DEBUGGING
#    axes.set_ylim([min(STRs)-0.1, max(STRs)+0.1])
    axes.set_xlim([min(abundances)-0.1, max(abundances)+0.1])
    axes.set_ylim([np.round(Ymin,1)-0.2, np.round(Ymax,1)+0.2])
    axes.set_xlabel('[{0}/H]'.format(el.getIonName(ion)))
    axes.set_ylabel('STR. value')
    if plotTitle is not '':
        axes.set_title(plotTitle)

    # Detection Limit?
    if lowLimit is not None:
        detXs = np.array([min(abundances)-1.0, max(abundances)+1.0])
        detYs = lowLimit[0]*detXs + lowLimit[1]
        pyplot.plot(detXs.tolist(), detYs.tolist(), 'k--')

    # CoG Limit?
    if hiLimit is not None:
        detYs = hiLimit[0]*detXs + hiLimit[1]
        pyplot.plot(detXs.tolist(), detYs.tolist(), 'k-.')

    # Plot the best-fit abundance line
    # Note: We're trying to fit a (hopefully) vertical line, which is ummmm...
    # ...bad for typical linear regression methods. So, we swap X and Y, and
    # pretend we're fitting a horizontal line (which also give us the ability to
    # use p-values to claim "vertical-ness"
    allPoints = zip(STRs.tolist(), abundances.tolist())
#    abFitPoints = np.array([point for point in allPoints if point[0]>(lowLimit[0]*point[1] + lowLimit[1])])
    abFitPoints = np.array(allPoints)
#    avgExPot = np.mean(abFitPoints[:,0])

    # Red = Low corrections
#    weightArray = [10.0*(4.0-abs(avgExPot-x)) for x in abFitPoints[:,0]]
#    abSlope, abIntercept = sp.polyfit(abFitPoints[:,0], abFitPoints[:,1], 1, w=weightArray)
#    pyplot.text(0.80, 0.9, 'slope={0:5.1f}\n[{1}/H]={2:1.2f}+/-{3:1.2f}'.format(1./abSlope, el.getIonName(ion), np.mean(abundances), np.std(abundances)), transform=axes.transAxes, color='r', fontsize='smaller')
#    detYs = np.array([np.round(Ymin,1)-0.5, np.round(Ymax,1)+0.5])
#    detXs = abSlope*detYs + abIntercept
#    pyplot.plot(detXs.tolist(), detYs.tolist(), 'r:')

    # Green = unweighted fit
#    abSlope, abIntercept, rv, pv, err = sp.stats.linregress(abFitPoints[:,0], abFitPoints[:,1])
#    pyplot.text(0.80, 0.9, 'slope={0:5.1f}\nR^2-value={1:1.4f}\n[{2}/H]={3:1.2f}+/-{4:1.2f}({5:3d})'.format(1./abSlope, rv**2, el.getIonName(ion), np.mean(abundances), np.std(abundances),len(abundances)), transform=axes.transAxes, color='g', fontsize='smaller')
    if len(redSet) > 0:
        pyplot.text(0.75, 0.95, '[{0}/H]={1:1.2f}+/-{2:1.2f}({3:3d})'.format(el.getIonName(ion), np.mean(redAbs), np.std(redAbs),len(redSet)), transform=axes.transAxes, color='r', fontsize='smaller')
#    else:
#        pyplot.text(0.75, 0.95, '(No data)', transform=axes.transAxes, color='r', fontsize='smaller')
        
    if len(greenSet) > 0:
        pyplot.text(0.75, 0.90, '[{0}/H]={1:1.2f}+/-{2:1.2f}({3:3d})'.format(el.getIonName(ion), np.mean(greenAbs), np.std(greenAbs),len(greenSet)), transform=axes.transAxes, color='g', fontsize='smaller')
#    else:
#        pyplot.text(0.75, 0.90, '(No data)', transform=axes.transAxes, color='g', fontsize='smaller')
    
    if len(blueSet)>0:
        pyplot.text(0.75, 0.85, '[{0}/H]={1:1.2f}+/-{2:1.2f}({3:3d})'.format(el.getIonName(ion), np.mean(blueAbs), np.std(blueAbs),len(blueSet)), transform=axes.transAxes, color='b', fontsize='smaller')
#    else:
#        pyplot.text(0.75, 0.85, '(No Data)', transform=axes.transAxes, color='b', fontsize='smaller')
        
    pyplot.text(0.75, 0.80, '[{0}/H]={1:1.2f}+/-{2:1.2f}({3:3d})'.format(el.getIonName(ion), np.mean(abundances), np.std(abundances),len(abundances)), transform=axes.transAxes, color='k', fontsize='smaller')
    
#    detYs = np.array([np.round(Ymin,1)-0.5, np.round(Ymax,1)+0.5])
#    detXs = abSlope*detYs + abIntercept
#    pyplot.plot(detXs.tolist(), detYs.tolist(), 'g:')
   
    # Do the same for just the points in the linear region of the CoG
#    linearPoints = np.array([point for point in allPoints if point[0]<(hiLimit[0]*point[1] + hiLimit[1])])
#    abSlope, abIntercept, rv, pv, err = sp.stats.linregress(linearPoints[:,0], linearPoints[:,1])

#    pyplot.text(0.50, 0.9, 'slope={0:5.1f}\nR^2-value={1:1.4f}\nAb.@0.0={2:1.2f}'.format(1./abSlope, rv**2, abIntercept), transform=axes.transAxes, color='g', fontsize='smaller')
#    detXs = abSlope*detYs + abIntercept
#    pyplot.plot(detXs.tolist(), detYs.tolist(), 'g:')
    
#    pyplot.xlabel('[{0}/H]'.format(el.getIonName(ion)))
#    pyplot.ylabel('Ex.Pot.')
    saveDirName = dirs.XPAbPlotsDir+el.getIonName(ion)+'/'
    try:
        os.stat(saveDirName)
    except:
        os.mkdir(saveDirName[:-1])
    pyplot.savefig('{0}{1}_{2}{3}.png'.format(saveDirName, starName, el.getIonName(ion), fileTag))
    pyplot.close()


def COGPlot(points, starName, ion, fileTag='', plotTitle='', pointLabels=None, fitXPLines=False):
# Plot a curve of growth (COG)
# The linear portion of the CoG is modeled by least squares (numpy polyfid for
# a 1-d polynomial.
# The passed data, in the array 'points' is of the form: 
# [wl, Ex.Pot., Log(gf), eqw, Log(eqw/wl), [X/H]]

    # For now, just bailing if no points were passed
    if len(points) == 0:
        return
        
    assert len(points)>0

    # We'll weight our linear COG fits by the individual points' [X/H] variance
    abMed = np.median(points[:,5])
#    abVar = np.std(points[:,5])

## Quick experiment - what does the plot look like if we use the median abundance
## instead of the measured?
## We'll color the points with ex.pot. for this one.
#    abStrs = np.array([abMed+line[2]+np.log10(line[0]) for line in points])
    abStrs = np.array([line[5]+line[2]+np.log10(line[0]) for line in points])
    
    # We're rounding 'cause...Python.
    XPs = np.around(points[:,1],4)
    minXP = min(XPs)
    maxXP = max(XPs)
    uniqueXPs = sorted(set(XPs))
    deltaXP = maxXP - minXP
    
    # Arrange points to analyze: 
    # I just want you to be confused by the different indices
    pointData = np.array(zip(abStrs, points[:,4], XPs, points[:,5], points[:,3], points[:,0]))
    
    linFits={}
    plotPoints={}
    
    XPBins = [uniqueXPs[0]]
    for xp in uniqueXPs:
        if xp>XPBins[-1]+k.evBinSize:
            XPBins.append(xp)

    minminErr = 10.
    maxmaxErr = 0.
    for xp in XPBins:
        fitPoints = np.array([p for p in pointData if p[2]>=xp and p[2]<=xp+k.evBinSize])
        if fitPoints.shape[0]<3:
            linFits[xp] = None
            plotPoints[xp]=None
            continue
        # Weight = distance to group median - max = 100 (ie: 1/sigma <=100)
        deltas = [(abs(ab-abMed))**2 for ab in fitPoints[:,3]]
        pointWeights = [1./(500.*d) if d > 0.002 else 1. for d in deltas]
        
        # The linfits dictionary entry will now contain:
        # [[Linear fit coefficients], [parameter estimated errors], [covariance matrix], 
        #  [estimated error in x], [estimated error in y]]
        linFits[xp], coeffErrs, covariance, xErrs, yErrs = u.odRegress(fitPoints[:,0], fitPoints[:,1], weights=pointWeights)
        pointErrs = [np.sqrt(p[0]**2+p[1]**2) for p in zip(xErrs, yErrs)]
        plotPoints[xp] = np.array(zip(fitPoints[:,0],fitPoints[:,1], pointErrs, fitPoints[:,5]))

        minErr = min(pointErrs)
        maxErr = max(pointErrs)
        if maxErr>maxmaxErr: maxmaxErr = maxErr
        if minErr<minminErr: minminErr = minErr

    # Normalize the errors for plotting
    deltaErr = maxmaxErr-minminErr
    for xp in XPBins:
        if plotPoints[xp] == None:
            continue
        pointErrs = plotPoints[xp][:,2]
        normedErrs = np.array([(pe-minminErr)/deltaErr for pe in pointErrs])
        plotPoints[xp] = np.array(zip(plotPoints[xp][:,0], plotPoints[xp][:,1], normedErrs, plotPoints[xp][:,3]))

    fig = pyplot.figure()
    axes = pyplot.gca()
    # Y axis is negative
    Xmin = min(pointData[:,0])*1.1
    Xmax = max(pointData[:,0])*0.9
    Xmed = np.median(pointData[:,0])
    Ymin = min(pointData[:,1])*0.9
    Ymax = max(pointData[:,1])*1.1
    
    legendCount=0
#    fp =open('/home/mikelum/Dropbox/CodeCloud/MyTools/ClusterAnalysis/Plots/slopes.txt','a')
    for xp in XPBins:
        if plotPoints[xp] is None:
            continue
#        fp.write('---------------\n{0} {1:1.2f}eV ({2:2d}): {3:2.3f}\n    '.format(el.getIonName(ion), xp, len(plotPoints[xp]), linFits[xp][0]))
        # Plot the fit
        minX = min(plotPoints[xp][:,0])
        minY = linFits[xp][0]*minX + linFits[xp][1]
        maxX = max(plotPoints[xp][:,0])
        maxY = linFits[xp][0]*maxX + linFits[xp][1]
        normedXP = ((xp-minXP)/deltaXP)
        plt = axes.plot([minX, maxX],[minY,maxY], color=pyplot.get_cmap('CMRmap')(normedXP))
        
        # Plot the points:
        # Color by distance to fit:
#        axes.scatter(plotPoints[xp][:,0], plotPoints[xp][:,1], figure=fig, marker='o', s=8, c=plotPoints[xp][:,2], cmap=pyplot.get_cmap('CMRmap'), edgecolors='None')
        # Color by Ex. Pot.
        axes.scatter(plotPoints[xp][:,0], plotPoints[xp][:,1], figure=fig, marker='o', s=8, color=pyplot.get_cmap('CMRmap')(normedXP), edgecolors='None')
        # Legend Key
        axes.text(0.75, 0.35-legendCount*0.05, r'{0:1.2f}eV: {1:2.3f}'.format(xp, linFits[xp][0]), transform=axes.transAxes, color=pyplot.get_cmap('CMRmap')(normedXP), fontsize='smaller')
        legendCount += 1
        
#        slopes = np.array([[(p[3],(p[1]-a[1])/(p[0]-a[0])) if (p[0]-a[0]) != 0. else (p[3],float('nan')) for p in plotPoints[xp]] for a in plotPoints[xp]])
        
#        diagIndex = 0
#        for p in slopes:
#            var = abs(linFits[xp][0]/np.nanmean(p[:,1]))
#            if var > 2. or var< 0.5: 
#                fp.write(' -**')
#            fp.write('{0:4.3f}:{1:2.3f}'.format(p[diagIndex,0], np.nanmean(p[:,1])))
#            if var > 2. or var< 0.5: 
#                fp.write('**- ')
#            else:
#                fp.write('   ')
#            diagIndex += 1
#        fp.write('\n')
#    fp.close()
    # Plot the linear CoG limit
    axes.plot([4,14],[k.cogPlotLimit,k.cogPlotLimit], linestyle=':')
    
#    cb = pyplot.colorbar()
#    cb.set_label('Excitation Potential (eV)')
#    cb.set_ticks([0., 0.5, 1.0])
#    cb.set_ticklabels(['{0:1.2f}'.format(minXP), '{0:1.2f}'.format(minXP+(maxXP-minXP)/2.), '{0:1.2f}'.format(maxXP)])

    
    # Write in a (tiny) label next to each point
    if pointLabels is not None:
        for coord in zip(abStrs.tolist(), XPs.tolist(), pointLabels):
            # Place the label below and right of the point
            axes.text(coord[0]+0.02*(Xmax-Xmin), coord[1]-0.03*(Ymax-Ymin), coord[2], color='k', fontsize='xx-small')

    # Axes labels and ranges
    axes.set_ylim([Ymin, Ymax])
    axes.set_xlim([Xmin, Xmax])
    axes.set_xlabel(r'$Log(A*gf*\lambda)$')
    axes.set_ylabel(r'$Log(EQW/\lambda)$')
    
    if plotTitle is not '':
        axes.set_title(plotTitle)

    saveDirName = dirs.CoGPlotsDir+el.getIonName(ion)+'/'
    try:
        os.stat(saveDirName)
    except:
        os.mkdir(saveDirName[:-1])
    pyplot.savefig('{0}{1}_{2}{3}.png'.format(saveDirName, starName, el.getIonName(ion), fileTag), figure=fig)
    pyplot.close()


def QualityPlot(points, starName, ion, fileTag='', plotTitle='', pointLabels=None):
# Plot an [X/H]vs.line strength plot, with the points colored by "quality" - 
# or how well they correspond to the linear fit (along the linear section of
# the COG) of other measurements with similar excitation potential values.
# The linear portion of the CoG is modeled by an Orthogonal Data Reduction (ODR)
# using a linear function (from sp.odr). 
# The quality of a point is measured as an inverse of its distance to the linear
# fit as measured by the ODR.
# The passed data, in the array 'points' is of the form: 
# [wl, Ex.Pot., Log(gf), eqw, Log(eqw/wl), [X/H]]
    # For now, just bailing if no points were passed
    nPoints = len(points)
    if nPoints == 0:
        return
#    assert nPoints>0
    # This abSTR measure is basically the log10 of the gravity*abundance*temp
    abStrs = np.array([line[5]+line[2]+np.log10(line[0]) for line in points])
    # We're rounding 'cause...Python.
    XPs = np.around(points[:,1],4)
    minXP = min(XPs)
    maxXP = max(XPs)
    uniqueXPs = sorted(set(XPs))
    XPDelta = maxXP - minXP
    normedXPs = [(xp-minXP)/XPDelta for xp in XPs]
    
    # Arrange points to analyze: 
    # I just want you to be confused by the different indices
    pointData = np.array(zip(abStrs, points[:,4], XPs, points[:,5], points[:,3]))
    
    # Bin our points by nearby XPs
    XPBins = [uniqueXPs[0]]
    for xp in uniqueXPs:
        if xp>XPBins[-1]+k.evBinSize:
            XPBins.append(xp)

    # We'll weight our linear COG fits by the individual points' [X/H] variance
    # Note: we'll change the median once we 'weed out' measurements outside of
    # the linear portion of the CoG.
    abMed = np.median(points[:,5])
#    abVar = np.std(points[:,5])

    linFits = {}
    plotPoints = {}
    minminErr = 10.
    maxmaxErr = 0.
    for xp in XPBins:
        fitPoints = np.array([p for p in pointData if p[2]>=xp and p[2]<=xp+k.evBinSize and p[1]<k.cogPlotLimit-xp*k.cogCorrFactor])
#        fitPoints = np.array([p for p in pointData if p[2]>=xp and p[2]<=xp+k.evBinSize])
        if fitPoints.shape[0]<2:
            linFits[xp] = []
            plotPoints[xp]=[]
            continue
        # Weight = distance to group median - max = 100 (ie: 1/sigma <=100)
        deltas = [abs(ab-abMed) for ab in fitPoints[:,3]]
        pointWeights = u.normalize([1./(500.*d) if d > 0.002 else 1. for d in deltas])
        
        # The linfits dictionary entry will now contain:
        # [[Linear fit coefficients], [parameter estimated errors], [covariance matrix], 
        #  [estimated error in x], [estimated error in y]]
        linFits[xp], coeffErrs, covariance, xErrs, yErrs = u.odRegress(fitPoints[:,0], fitPoints[:,1], weights=pointWeights)
        pointErrs = [np.sqrt(p[0]**2+p[1]**2) for p in zip(xErrs, yErrs)]
        plotPoints[xp] = np.array(zip(fitPoints[:,3],fitPoints[:,1],pointErrs,fitPoints[:,2],pointWeights))

        minErr = min(pointErrs)
        maxErr = max(pointErrs)
        if maxErr>maxmaxErr: maxmaxErr = maxErr
        if minErr<minminErr: minminErr = minErr

    # Calculate the point 'score' based on priors
    deltaErr = maxmaxErr-minminErr
    for xp in XPBins:
        if len(plotPoints[xp]) == 0:
            continue
        # Normalize the errors for the first prior
        normedErrs = np.array([(pe-minminErr)/deltaErr for pe in plotPoints[xp][:,2]])
        # Normalized XP for second prior, however, we want 0.0eV to be highest
        # values, and the highest XP to be valued no lower than 0.5
        normedXPs = np.array([1.-((px-minXP)/XPDelta)/2. for px in plotPoints[xp][:,3]])
        
        combPriors = u.normalize(np.array([(ps[0]*(ps[1]**2)*ps[2])**0.25 for ps in zip(normedErrs, normedXPs, plotPoints[xp][:,4])]))
        plotPoints[xp] = np.array(zip(plotPoints[xp][:,0], plotPoints[xp][:,1], combPriors))

        
    fig = pyplot.figure()
    axes = pyplot.gca()
    # Y axis is negative
    Ymin = min(pointData[:,1])*1.05
    Ymax = max(pointData[:,1])*0.95
    Xmin = min(pointData[:,3])*0.95
    Xmax = max(pointData[:,3])*1.05
        
    if len(points) < 200:
        pointSize = 9
        if pointLabels is not None:
            # Must have one label for each point (even if it is '')
            if len(pointLabels) < len(abStrs):
                # Must have labels for all the points
                pointLabels = None
    else:
        pointLabels = None # Seriously, don't label more than 100 points
        pointSize = 5
    
    allPlotPts = [] # For median/variance calcs.
    for xp in XPBins:
        if len(plotPoints[xp]) != 0:
            pyplot.scatter(plotPoints[xp][:,0], plotPoints[xp][:,1], figure=fig, marker='o', s=pointSize, c=1.-plotPoints[xp][:,2], cmap=pyplot.get_cmap('CMRmap'), edgecolors='None')
            allPlotPts.extend(plotPoints[xp])
    
    # It's possible that we've weeded a bit too heavily, in which case, 
    # we're done. 
    if len(allPlotPts) < 2:
        pyplot.close()
        return
        
    abMed = np.median(np.array(allPlotPts)[:,0])
    abVar = np.std(np.array(allPlotPts)[:,0])
    # Plot the median abundance value
    axes.plot([abMed,abMed],[Ymin, Ymax], linestyle='--', color='r')
    axes.plot([abMed-abVar, abMed-abVar],[Ymin, Ymax], linestyle=':', color='c')
    axes.plot([abMed+abVar, abMed+abVar],[Ymin, Ymax], linestyle=':', color='c')

    pyplot.text(0.65, 0.25, r'{0}({1:3d}): {2:1.3f}+/-{3:1.3f}'.format(el.getIonName(ion), nPoints, abMed, abVar), transform=axes.transAxes, fontsize='smaller')

    cb = pyplot.colorbar()
    cb.set_label('Point \"Quality\" (Score)')
    cb.set_ticks([0., 0.5, 1.0])
#    cb.set_ticklabels(['{0:1.2f}'.format(minminErr), '{0:1.2f}'.format(minminErr+(maxmaxErr-minminErr)/2.), '{0:1.2f}'.format(maxmaxErr)])
    cb.set_ticklabels(['1.00', '0.50', '0.00'])
    
    # Write in a (tiny) label next to each point
    if pointLabels is not None:
        for coord in zip(abStrs.tolist(), wls.tolist(), pointLabels):
            # Place the label below and right of the point
            pyplot.text(coord[0]+0.02*(Xmax-Xmin), coord[1]-0.03*(Ymax-Ymin), coord[2], color='k', fontsize='xx-small')

    # Axes labels and ranges
    axes.set_ylim([Ymin, Ymax])
    axes.set_xlim([Xmin, Xmax])
    axes.set_xlabel('[{0}/H]'.format(el.getIonName(ion)))
    axes.set_ylabel(r'$Log(EQW/\lambda)$')
    
    if plotTitle is not '':
        axes.set_title(plotTitle)

    saveDirName = dirs.CoGPlotsDir+el.getIonName(ion)+'/'
    try:
        os.stat(saveDirName)
    except:
        os.mkdir(saveDirName[:-1])
    pyplot.savefig('{0}{1}_{2}{3}.png'.format(saveDirName, starName, el.getIonName(ion), fileTag), figure=fig)
    pyplot.close()


def XPAbPlot(points, starName, ion, fileTag='', plotTitle='', pointLabels=None, showTrendLine=False):
# Plot the Excitation potential vs abundance for the passed points.
# The passed data, in the array 'points' is of the form: 
# [[Ex. Pot, [X/H], EQW (mA), wl](A), ...]
    # For now, just bailing if no points were passed
    if len(points) == 0:
        return
    assert len(points)>0

    # We're rounding 'cause...Python.
    XPs = np.around(points[:,0],4)
    Abs = points[:,1]

    fig = pyplot.figure()
    axes = pyplot.gca()

    Xmin = min(XPs)*0.9
    if Xmin==0.:
        Xmin = -0.25
    Xmax = max(XPs)+0.25

    Ymin = min(Abs)-0.2
    Ymax = max(Abs)+0.2
    
    pyplot.scatter(XPs, Abs, figure=fig, marker='o', s=9)

    # Write in a (tiny) label next to each point
    if pointLabels is not None:
        for coord in zip(XPs.tolist(), Abs.tolist(), pointLabels):
            # Place the label below and right of the point
            pyplot.text(coord[0]+np.random.random()*0.04*(Xmax-Xmin), coord[1]-np.random.random()*0.06*(Ymax-Ymin), coord[2], color='k', fontsize='xx-small')
    
    # Fit a line, if desired
    if showTrendLine:
        slope, intercept, rVal, pVal, err = sp.stats.linregress(XPs, Abs)
        pyplot.plot([Xmin,Xmax],[intercept+slope*Xmin, intercept+slope*Xmax], 'g:')
        axes.text(0.95, 0.01, 'Slope:{0:4.2f}, R-squared:{1:4.2f}'.format\
                (slope,rVal**2), fontsize=9, verticalalignment='bottom',\
                horizontalalignment='right', transform=axes.transAxes,)
    # Show mean and standard deviation range:
    mean = np.mean(Abs)
    std = np.std(Abs)
    axes.axhline(y=mean, ls='dashed', color='r', label=r'$[{0}/H]:{1:4.2f}\pm{2:4.2f}$'.format(el.getIonName(ion), mean, std))
    axes.axhline(y=mean-std, ls='dotted', color='r')
    axes.axhline(y=mean+std, ls='dotted', color='r')
    pyplot.plot([Xmin,Xmax],[intercept+slope*Xmin, intercept+slope*Xmax], 'g:')
    pyplot.plot([Xmin,Xmax],[intercept+slope*Xmin, intercept+slope*Xmax], 'g:')

    # Axes labels and ranges
    axes.set_ylim([Ymin, Ymax])
    axes.set_xlim([Xmin, Xmax])
    axes.set_xlabel('Excitation Potential (eV)')
    axes.set_ylabel('Abundance [X/H]')
    axes.legend()
    
    if plotTitle is not '':
        axes.set_title(plotTitle)

    saveDirName = dirs.XPAbPlotsDir+el.getIonName(ion)+'/'
    try:
        os.stat(saveDirName)
    except:
        os.mkdir(saveDirName[:-1])
    pyplot.savefig('{0}{1}_{2}{3}.png'.format(saveDirName, starName, el.getIonName(ion), fileTag), figure=fig)
    pyplot.close()


# Obsolete function to plot Teff vs. Ab for all stars in the cluster. Now 
# superceded by: TAbPlots.MakeTAbPlots
'''
def ClusterAbPlots(clusterName, ions, outFilename='temp_clustPlot.ps'):
    saveDirName = dirs.PlotDir

    clustMet = -0.10

    # Load giant models
    gModelFiles = mk.findDataFiles(k.GiantModelPath)
    gModelAtms, gPradks = mk.LoadModels(gModelFiles)
    giantMass = 1.15
    giantMet = 0.00

    # Load Dwarf models
    dModelFiles = mk.findDataFiles(k.DwarfModelPath)
    dModelAtms, dPradks = mk.LoadModels(dModelFiles)
    dwarfMass = 0.00
    
    starNames = sdb.GetStarsForCluster(clusterName)
    parmDict = {}
    abDict = {}
    lineDict = {}
    
    for starName in starNames:
        parmDict[starName] = sdb.GetStarParms(starName)
        starParms = (parmDict[starName][0], parmDict[starName][1], parmDict[starName][2], clustMet, dwarfMass)
        starAbs, starLines, unusedMin, unusedMax = AB.CalcAbsAndLines(clusterName+' '+starName, starParms, modelAtms=dModelAtms, pradks=dPradks)
        # For now, we're just going to use a straight mean. However, we know
        # that individual lines can (and will) be removed due to blends and 
        # (weak) strength
        abDict[starName] = starAbs
        lineDict[starName] = starLines
        
    for ion in sorted(ions):
        plotPoints = []

        for star in starNames:
            starIndex = starNames.index(star)
            if ion in abDict[star].keys():
                # We can start our line filtering here:
                # shortDict[ion] = [line for line in abDict[star][ion] if line[4]<-4.70]
#                print lineDict[star][ion]
#                elemAvg = np.mean([line[5] for line in lineDict[star][ion]])
#                elemStd = np.std([line[5] for line in lineDict[star][ion]])

                plotPoints.append([starIndex, abDict[star][ion][0], abDict[star][ion][1]])
        
        Xs = [x[0] for x in plotPoints]
        Ys = [x[1] for x in plotPoints]
        Yerrors = [x[2] for x in plotPoints]
        
        fig = pyplot.figure()
        axes = pyplot.gca()
        
        pyplot.xlabel('Stars from {0}'.format(clusterName))
        pyplot.xticks(np.arange(-1,len(Xs)+1,1))
        pyplot.ylabel('[{0}/H]'.format(el.getIonName(ion)))

        axes.set_xlim([0, len(Xs)])
#        Ymin = min(Ys)
#        Ymax = max(Ys)
        Ymean = np.mean(Ys)
        axes.set_ylim([round(Ymean-1.), round(Ymean+1.)])
        axes.set_xticklabels(['']+starNames+[''], rotation=90, fontsize=8)
        
        axes.errorbar(Xs, Ys, yerr=Yerrors, fmt='o')
                
        pyplot.savefig('{0}{1}_{2}{3}.png'.format(saveDirName, clusterName, el.getIonName(ion), plotTag), figure=fig)
        pyplot.close()
'''
