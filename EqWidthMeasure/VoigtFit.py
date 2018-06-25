#! /usr/bin/env python
# 
# Module: gaussianFitTest.py
#
# Author: Mike Lum
#
# Description: Creating my own line fitting functions
#
# Contents: List externally-callable functions here
#   Function xxx: (Description)
#
# Revision History:
#    Date        Vers.    Author        Description
#    Today        0.0a0    Lum        First checked in
#
# To Do:
#    
#

# Imports
from matplotlib import pyplot
import numpy as np
from scipy.optimize import curve_fit
import scipy.signal as sig
import scipy.interpolate as interp
from scipy.integrate import quad
#import pyfits as pf
import astropy.io.fits as pf
import pylab as pl


from utilities import utilities as u
from constants import constants as k
from EqWidthMeasure import Baseline as BL
from Databases import LineLookup as LL
from SpectraProcessing import SpectraData as SD

# Global (ugh! - but required for scipy's curve_fit)
global gBaselineSpline

# Functions
def trueVoigt(x, scale, offset, gamma):
    val = gaussian(x, scale, offset, sigma) * lorenzian(x, 1.0, offset, gamma)
    return val
       
def continuum(x):
    global gBaselineSpline
    try:
        return gBaselineSpline(x)
    except:
        print("No spline")
        return 1.01
    
def gaussian(x, scale, offset, sigma):
    # Gaussian distribution:
    # scale * (exp(-(x-b)**2/((2*sigma)**2))/(sigma * sqrt(2*pi))
    # 1/sqrt(2*pi) = 0.3989422
    val = 0.3989422*scale * np.exp(-(x-offset)**2/(2*sigma**2))/(sigma)
    return continuum(x) - val

def lorenzian(x, scale, offset, gamma):
    # Lorenzian distribution:
    # scale * (gamma/(pi*((x-offset)**2+gamma**2))
    val = scale/np.pi  * (gamma/((x-offset)**2+gamma**2))
    return continuum(x) - val

def lorenzWings(xArray, scale, offset, gamma, sigma, coreFraction):
    # A "made-up" line profile which uses a "pure" Gaussian line core,
    # "pure" Lorenzian wings, and a linear interpolation between the two 
    # models for intermediate values. The transition point is passed as a fit 
    # parameter.
    
    # Note: This isn't what the code below actually does...Hmm...
    
    weight = 1. if coreFraction > 1. else 0. if coreFraction < 0. else coreFraction
    wingModel = weight*np.array(gaussian(xArray, scale, offset, sigma)) + (1.-weight)*np.array(lorenzian(xArray, scale, offset, gamma))

    return wingModel

def estVoigtFWHM(gamma, sigma):
    # Estimation from Thompson et al. 1987
    # fsubg = 2.*sigma*sqrt(2*ln(2))      
    fsubg = 2.354820*sigma
    fsubl = 2.*gamma     
    try:
        fwhm = (fsubg**5 + 2.69269*fsubg**4*fsubl + 2.42843*fsubg**3*fsubl**2 + 4.47163*fsubg**2*fsubl**3 + 0.07842*fsubg*fsubl**4 + fsubl**5)**0.20
    except Warning:
        print("Warning: %.3f  %.3f  %.3f" % (gamma, sigma, fwhm))
        pass
    return fwhm
    
def fitCurve(xRange, yRange, eRange, guesses):
    # Parms for lorenzWings = scale, offset, gamma, sigma, coreFraction
    bounds = ((0., -1., 0., 0., 0.),(np.inf, 1., 1., 1., 1.))
    try:
        return curve_fit(lorenzWings, xRange, yRange, sigma = eRange, p0=guesses, bounds=bounds)
    except ValueError:
    # This occurs when the passed points cannot be fit...
    # There are various higher-level handlers which respond, based on which
    # parpmeters were trying to be fit.
        raise RuntimeError('Unable to fit passed ranges')

def calcEQW(popt, pcov, x):
    global gBaselineSpline # Referenced for the eqm calculation
    
    # Determine the x range to fit - 2.5*sigma, based on the gaussian fit
    # lineCenter = popt[1]
    # sigma = popt[3]
    minX = u.closest(popt[1]-2.5*abs(popt[3]),x)
    maxX = u.closest(popt[1]+2.5*abs(popt[3]),x)

    # scipy.integrate.quad - integrate the fn (area under the plot)
    contEQWResult, contError = u.intSimpsons(gBaselineSpline, minX, maxX)
    
    # Typical eqw fitting assumes that the continuum is a constant line,
    # normalized to 1.0. We have a spline, which may, or may not be a line,
    # and probably doesn't have a constant value of 1.0 - So, we'll need
    # to normalize our eqw result to an assumed continuum of 1.0.
    normFactor = contEQWResult/(maxX-minX)
    
    fitEQWResult, fitError = u.intSimpsons(lorenzWings, minX, maxX, popt[0], popt[1], popt[2], popt[3], popt[4])
    
    eqw = (contEQWResult - fitEQWResult)*1000./normFactor
    eqwError = (abs(contError) + abs(fitError))*1000./normFactor
    return eqw, eqwError
    
    
def showLineFit(popt, pcov, x, y, e, centerwl):
    global gBaselineSpline
            
#    pyplot.errorbar(x+centerwl, y, yerr=e, linewidth=1, color='black', fmt='none')
    pyplot.plot(x+centerwl, y, 'b+')
    
    xm = np.linspace(x[0], x[-1], 150)
    
    pyplot.plot(xm+centerwl, lorenzWings(xm, popt[0], popt[1], popt[2], popt[3], popt[4]), color='red')
    pyplot.plot(xm+centerwl, gBaselineSpline(xm), color='blue')
#    if popt[4] > 0. and popt[4] < 1.:
    # Only show the components if we have a mix.
#        pyplot.plot(xm+centerwl, lorenzian(xm, popt[0], popt[1], popt[2]), color='blue')
#        pyplot.plot(xm+centerwl, gaussian(xm, popt[0], popt[1], popt[3]), color='green')

    pyplot.show()


def measureVoigtEQWs(specFilename, lineData=None):
    global gBaselineSpline

    objectName, numApertures, apertureSize, apBoundaries = SD.ReadSpectraInfo(specFilename)
    specList = SD.LoadAndSmoothSpectra(specFilename)
    measuredLines = []
    badLines = []
    dopplerCheck = 0.
    dopplerCount = 0.

    for smoothSpec in specList:
        wls = smoothSpec[0]
        fluxes = smoothSpec[1]
        if not SD.IsGoodAperture(smoothSpec):
            continue
        if lineData==None:
            lineData, unusedRefs = LL.getLines(wavelengthRange=[(wls[0], wls[-1])], dataFormat=None)   
            
        for line in lineData:    
            # lineData entries are variable length, depending on how many 
            # optional parameters are included. So, we have to manually assign
            # the ones we need.
            wl = line[0]
            ion = line[1]
            ep = line[2]
            gf = line[3]
            if len(line) > 4:
            # This is the only optional parameter we might need.
                c6 = line[4]
            else:
                c6 = ''
                  
            wlUpper = wl + k.HalfWindowWidth
            wlLower = wl - k.HalfWindowWidth

            loApPix = u.closestIndex(wls, wlLower)
            centerPix = u.closestIndex(wls, wl)
            hiApPix = u.closestIndex(wls, wlUpper)

            # The algorithm for curve_fit needs values "close" to one. So,
            # we normalize our wavelength-calibrated spectrum window to be
            # centered at "0", and range from -2 to +2. We also assume that
            # our continuum normalized spectrum has amplitudes near 1. 
            x = wls[loApPix:hiApPix+1] - wls[centerPix]
            y = fluxes[loApPix:hiApPix+1]
            # At some point, I suppose we should use the actual S/N of the spectra
            # for the error array...
            e = 0.01*np.ones(len(x))
            
            if len(x) < 10:
                # Line is too close to the order boundary, so skip it
                badLines.append([wl, ion, k.badSpectralRegion, 'Too close to order boundary'])
                continue
                
            try:
            # Because of the normalization needed for curve fitting, we need to perform
            # some shenanigans on the baseline spline
                rawSpline = BL.FitSplineBaseline(loApPix, hiApPix, smoothSpec)
                gBaselineSpline = interp.UnivariateSpline(x, rawSpline(x+wls[centerPix]))
            except:
                # We've had errors when trying to fit bad regions - just skip to the next one
                print('Fitting error for {0:4.3f}, ({1:2.1f}).'.format(wl, ion))
                continue
            
            maxIdxs = [mx for mx in sig.argrelmax(y)[0].tolist() if y[mx]<2.5]
            maxIdxs.append(0)
            maxIdxs.sort()
            maxIdxs.append(-1)

            minIdxs = [mn for mn in sig.argrelmin(y)[0].tolist() if y[mn]> 0.05]
            minIdxs.append(0)
            minIdxs.sort()
            minIdxs.append(-1)
            mins = [x[idx] for idx in minIdxs]
            
            closestMin = mins.index(u.closest(0, mins))
            
            # Find the two local maxes which bracket the central wavelength (offset=0.)
            maxRange = u.bracket(0., [x[idx] for idx in maxIdxs])
            fitRange = [maxIdxs[maxRange[0]], maxIdxs[maxRange[1]]+1]
            
            # guess that the scale is the difference between the local max and mins,
            # although we could opt to use the difference between the continuum and the min...
            scaleGuess = (y[maxIdxs[maxRange[0]]]+y[maxIdxs[maxRange[1]]])/2. - y[minIdxs[closestMin]]
            # The line center should occur at, or really close to the offset
            offsetGuess = 0.0
            # Gamma and sigma are Lorenzian and Gaussian parms...need a better idea for guesses
            gammaGuess = 0.10
            sigmaGuess = 0.10
            # The core fraction is basically how much of the line is "pure" Gaussian
            coreFractionGuess = 0.65
            guessParms = [scaleGuess, offsetGuess, gammaGuess, sigmaGuess, coreFractionGuess]
            curveFitted = False
            
            if len(x[fitRange[0]:fitRange[1]]) < 5:
            # Way too few points to even do an interpolation. Try extending to the 
            # two "maxes" on either side of the current
                loRange = 0 if maxRange[0]<2 else maxRange[0]-1
                hiRange = -1 if maxRange[1]>len(maxIdxs)-2 or maxRange[1]==-1 else maxRange[1]+1

                fitRange = [maxIdxs[loRange], maxIdxs[hiRange]+1 if hiRange>0 else -1]
                xToFit = x[fitRange[0]:fitRange[1]]
                yToFit = y[fitRange[0]:fitRange[1]]
                eToFit = e[fitRange[0]:fitRange[1]]
                                
            if len(x[fitRange[0]:fitRange[1]]) < 5:
            # STILL too small of a range - Bail, with a bad line.
                badLines.append([wl, ion, k.badSpectralRegion, 'Line region too small'])
                continue
                
            if len(x[fitRange[0]:fitRange[1]]) < 20:
            # Too few points?
            # Try adding more (interpolated) data points, along a spline
                if fitRange[1]==-1:
                    ySpline = interp.interp1d(x[fitRange[0]:], y[fitRange[0]:], kind='cubic')
                    eSpline = interp.interp1d(x[fitRange[0]:], e[fitRange[0]:], kind='cubic')

                    xToFit = np.linspace(x[fitRange[0]], x[-1], num=100)
                    yToFit = ySpline(xToFit)
                    eToFit = eSpline(xToFit)
                else:
                    ySpline = interp.interp1d(x[fitRange[0]:fitRange[1]+1], y[fitRange[0]:fitRange[1]+1], kind='cubic')
                    eSpline = interp.interp1d(x[fitRange[0]:fitRange[1]+1], e[fitRange[0]:fitRange[1]+1], kind='cubic')
                    xToFit = np.linspace(x[fitRange[0]], x[fitRange[1]], num=100)
                    yToFit = ySpline(xToFit)
                    eToFit = eSpline(xToFit)
            else:
                xToFit = x[fitRange[0]:fitRange[1]]
                yToFit = y[fitRange[0]:fitRange[1]]
                eToFit = e[fitRange[0]:fitRange[1]]
                    
            try:
                popt, pcov = fitCurve(xToFit, yToFit, eToFit, guessParms)
                if popt[4] > 0. and popt[4] < 1.:
                # Even though we have a fit, if it's "pure" Lorenzian or Gaussian
                # pretend we can do better with some data manipulation
                    curveFitted = True
#                else:
#                    curveFitted = False                
            except RuntimeError:
            # Occurs when curve_fit cannot converge (Line not fit)
            # We need to catch it in order to try our own data manipulation
            # for a good fit.
#                print('Initial fit failure for line at: %.3f' % wl)
                pass
            
            if not curveFitted:
            # Try to refit, by addressing typical issue(s)
                # Assume not enough of the line is visible (no wings)
                # Add ~10 points along the continuum, in the "wings"
                loBound = fitRange[0] - 5
                if loBound<0: loBound=0
                
                if fitRange[1] == -1:
                    hiBound = -1
                    newY = np.concatenate((gBaselineSpline(x[loBound:fitRange[0]]), y[fitRange[0]:hiBound]))
                else:
                    hiBound = fitRange[1]+5
                    if hiBound > len(x)-1: hiBound=-1
                    newY = np.concatenate((gBaselineSpline(x[loBound:fitRange[0]]), y[fitRange[0]:fitRange[1]+1], gBaselineSpline(x[fitRange[1]+1:hiBound])))

                newX = x[loBound:hiBound]
                newE = e[loBound:hiBound]
                
                if len(newX)!=len(newY) or len(newX)!=len(newE):
                # Probably just have an indexing problem in the code immediately above, 
                # but just skip the line for now.
                    badLines.append([wl, ion, k.unknownReason, k.unknownStr])
                    continue
                    
                try:
                    popt, pcov = fitCurve(newX, newY, newE, guessParms)
                    # As long as something fits, here - we'll take it.
                    curveFitted = True
                except RuntimeError:
                # Well, and truly screwed?
    #                pyplot.errorbar(x+smoothSpec.xarr[centerPix], y, yerr=e, linewidth=1, color='black', fmt='none')
    #                pyplot.show()
                    badLines.append([wl, ion, k.unknownReason, k.unknownStr])
                    continue

            eqw, eqwError = calcEQW(popt, pcov, x)
            
#            showLineFit(popt, pcov, x, y, e, smoothSpec.xarr[centerPix].value)
            # FWHM of the Gaussian core is sqrt(8*ln(2))*sigma = 2.355*sigma
            # The factor of 1000 is to convert to miliangstroms.
            measuredLines.append([wl, ion, ep, gf, c6, eqw, wl+popt[1], 2355.*popt[3]])
            thisDoppler = ((popt[1])/wl) * k.SpeedofLight
            dopplerCheck += thisDoppler
            dopplerCount += 1
   
    #        if popt[4] <= 0.:
    #            print('Best fit is a pure Lorenzian')
    #        elif popt[4] >= 1.:
    #            print('Best fit is a pure Gaussian')
                

#    nans = len([x for x in measuredLines if np.isnan(x[5])])
    badLines.extend([[x[0], x[1], k.badEQWMValue, k.nanStr] for x in measuredLines if np.isnan(x[5])])
    badLines.extend([[x[0], x[1], k.emissionLine, "Emission strength: %.3f"%(abs(x[5]))] for x in measuredLines if x[5]<0.])
    badLines.extend([[x[0], x[1], k.centerWLOff, "Line center delta: %.3f"%(x[0]-x[6])] for x in measuredLines if abs(x[0]-x[6]) > 0.25])
#    print("nan eqw: %d  -  emission: %d  - bad offset: %d" % (nans, emiss, bads))
#    
#    for line in sorted(measuredLines, key=lambda x: x[0]):
#        if not np.isnan(line[5]) and line[5]>0. and abs(line[0]-line[6]) < 0.25:
#            print("%2.1f :%4.3f(%4.3f)  =  %.1f (%.1f)" % (line[1], line[0], line[6], line[5], line[7]))
    # Measured WL within 250mA of expected, non-emission lines, only
    filteredLines = [l for l in measuredLines if not np.isnan(l[5]) and l[5]>0. and abs(l[0]-l[6]) < 0.25]
     # Want ascending wls within ascending ions
    filteredLines.sort(key=lambda x: x[0])
    filteredLines.sort(key=lambda x: x[1])
#    print "*****************"
#    print "*******{0} bad *********".format(len(badLines))
#    print "*******{0} measured *********".format(len(measuredLines))
#    print "*******{0} filtered **********".format(len(filteredLines))
   
    return filteredLines, badLines, dopplerCount, dopplerCheck
    
