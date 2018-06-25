#! /usr/bin/env python
# 
# Module: Baseline.py
#
# Author: Mike Lum
#
# Description: Functions to calculate (spline-) baselines for
#               a spectral region
#
# Contents: 
#
# Revision History:
#    Date        Vers.    Author        Description
#    2018-04-03  1.1f1    Lum        No more pyspeckit
#    2016-04-14  1.0f1    Lum        Functions moved from EqWidthMeasure.eqwm
#
# To Do:
#    
#

# Imports
import numpy as np
import scipy.signal as sig
import scipy.interpolate as interp

# Functions:

def FitSplineBaseline(loPix, hiPix, smoothSpec):
# Create a polynomial (spline) baseline/continuum fit. We assume that
# the most important region for the fit will lie in the middle third 
# of the passed window.
#
# smoothSpec is assumed to be a list of two equal-length lists, wl and flux
# Note: This once assumed that the sommthed spectra was a pyspeckit object
# We found it easier to just deal with the actual wl and flux arrays.

    # Our (spline) baseline is a "connect the dots" of all the local
    # maxima, which are within 2 standard deviations of the mean of
    # the local maxes. This method works great where the spectra is
    # relatively "flat" and "quiet", but can get a bit hairy when it
    # is noisy (eg: <400nm)
    maxes = sig.argrelmax(smoothSpec[1][loPix:hiPix])[0]
    yArr = [smoothSpec[1][pixel+loPix] for pixel in maxes]
    pairs = zip(maxes, yArr)

    highMean = np.mean(yArr)
    highSTD = np.std(yArr)
    filteredMean = np.mean([x for x in yArr if abs(x-highMean)<2*highSTD])
    
    first3rdIdx = (hiPix-loPix)/3
    second3rdIdx = 2*(hiPix-loPix)/3
    # Select points in the first 1/3 of the window
    first3rdPts = [(smoothSpec[0][x[0]+loPix], x[1]) for x in pairs if x[0] <= first3rdIdx and abs(x[1]-filteredMean)<2*highSTD]
    # Select points in the 2nd 1/3 of the window
    second3rdPts = [(smoothSpec[0][x[0]+loPix], x[1]) for x in pairs if first3rdIdx <= x[0] <= second3rdIdx and abs(x[1]-filteredMean)<2*highSTD]
    # Select points in the final 1/3 of the window
    final3rdPts = [(smoothSpec[0][x[0]+loPix], x[1]) for x in pairs if x[0] >= second3rdIdx and abs(x[1]-filteredMean)<2*highSTD]

    # 
    # We need an average of the peaks in the first and last thirds of the
    # "window" so that the "ends" of our spline are well-defined. In the 
    # case where we have no local maxes in these regions, we use the 
    # overall average of the spectral "highs".
    highPoints = []
    if len(first3rdPts) > 0:
        highPoints = sorted(first3rdPts, key=lambda x:x[1])[-3:]
        firstAvg = np.mean([x[1] for x in first3rdPts])
    else:
        firstAvg = highMean
    if len(second3rdPts) > 0:
        highPoints.extend(sorted(second3rdPts, key=lambda x:x[1])[-3:])
    if len(final3rdPts) > 0:
        highPoints.extend(sorted(final3rdPts, key=lambda x:x[1])[-3:])
        lastAvg = np.mean([x[1] for x in final3rdPts])
    else:
        lastAvg = highMean
    
    highPoints.extend([(smoothSpec[0][loPix], firstAvg),(smoothSpec[0][hiPix], lastAvg)])
    highPoints.sort(key=lambda x:x[0])
    
    plotX = [x[0] for x in highPoints]
    plotY = [x[1] for x in highPoints]
    if len(plotX) < 4:
        polyOrder = len(plotX) -1
    else:
        polyOrder = 3

    pointWeights = 2.*np.ones(len(plotX))/np.std(plotY)
    smoothingFactor = len(pointWeights) + int(np.sqrt(len(pointWeights)))-2
    
    spline = interp.UnivariateSpline(plotX, plotY, k=polyOrder, w=pointWeights, s=smoothingFactor)
    return spline   



# Obsolete function(s)
'''      
def findLocalBaseline(sp, wlLower, wlUpper):

    minIdx = sp.xarr.x_to_pix(wlLower)
    maxIdx = sp.xarr.x_to_pix(wlUpper)
    localSpectra = sp.data[minIdx:maxIdx]
    pointMask = localSpectra == localSpectra
    
    # Two iterations is sort of a minimum. I believe sdas IRAF continuum uses 3 as a default
    # Each iteration keeps points between 0.5s.d below and 3.5s.d above the mean.
    mean = np.mean(localSpectra)
    std = np.std(localSpectra)
    trimmedData = [x for x in localSpectra if ((x > mean-0.5*std) and (x < mean+3.5*std))]
    mean = np.mean(trimmedData)
    std = np.std(trimmedData)
    trimmedData = [x for x in trimmedData if ((x > mean-0.5*std) and (x < mean+3.5*std))]
    # Scaling the continuum by a factor to account for line blanketing.
    # Note: The proper compensation for line blanketing should be based on much more
    # than trial-and-error. It should take into account metallicity, atmospheric conditions,
    # etc. We should also have a wavelength dependence. Blanketing is an issue at UV wls,
    # but not so much as you get closer to IR.
    mean = np.mean(trimmedData) * k.ContinuumBlanketFactor
    
    # if our data is sketchy (ie: lots of CRs, etc.), we might be way out of range
    # so, put us above a minimum for a continuum-ed spectra
    if mean < k.MinBaseline:
        mean = k.MinBaseline
    
    return mean
'''


