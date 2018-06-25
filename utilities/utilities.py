#! /usr/bin/env python
# 
# Module: utilities.py
#
# Author: Mike Lum
#
# Functions:
#       Boolean in_range(num, lo, hi) : Returns whether num is between lo and hi, inclusive
#       Boolean is_number(testString) : Returns whether the string can be cast as a number
#       (int, int) bracket(value, sortedRange): Returns indices of the two values in the passed list
#                                           which bracket the passed value.
#
#
# Description: Utilities for the SAAT functions
#
# Revision History:
#    Date        Vers.    Author        Description
#    2015-06-05  1.0f1     Lum       Added utility to merge improperly split strings from .csv files
#    2014-07-31  0.0a0    Lum        First checked in
#
# To Do:
#    
#
import numpy as np
from constants import constants as k
import scipy.odr as odr
import scipy.optimize as opt
import scipy.stats as stats
import matplotlib.pyplot as pyplot
import os


# Math
def Cartesian(*arrays):
# Returns the cartesian product of the passed arrays. Blatantly stolen from:
# https://stackoverflow.com/questions/1208118/using-numpy-to-build-an-array-of-all-combinations-of-two-arrays
    mesh = np.meshgrid(*arrays)  # standard numpy meshgrid
    dim = len(mesh)  # number of dimensions
    elements = mesh[0].size  # number of elements, any index will do
    flat = np.concatenate(mesh).ravel()  # flatten the whole meshgrid
    reshape = np.reshape(flat, (dim, elements)).T  # reshape and transpose
    return reshape

def Gaussian(x,peak=1.0, mu=0.0, sigma=1.0):
# Returns P(x) for a Gaussian distribution with mean=mu and variance=sigma^2
# (This is probably defined in scipy, or numpy, but I'm too lazy to go looking)
    k = 1.0/np.sqrt(2.*np.pi)
    return peak*(k/sigma)*np.e**-((x-mu)**2/(2*sigma**2))

def Linear(coefficients, x):
# Returns a linear conversion of the passed value (x), calculated from the 
# passed coefficient array. The coefficients are in descending order, ie:
#                       y = c[0]*x + c[1] 
# This function is mostly used as the model for ODR calculations.
    # Must be an array with the two coefficients
    assert len(coefficients) == 2
    return (coefficients[0]*x+coefficients[1])
    
    
def odRegress(x, y, weights=[], linreg=None):
# Produce Orthogonal Data Reduction (fit) for the passed x,y pairings.
# Returns: [fit parameters] , [parm. errors], [[covariance matrix]],
#           [estimated error in x vars], [estimated error in y vars]
# Code modified from: Robin Wilson, see the site for copyright notice:
# http://blog.rtwilson.com/orthogonal-distance-regression-in-python/
#
    # x,y both arrays of same length
    assert len(x) == len(y)
    
    if linreg is None:
        linreg=stats.linregress(x,y)
        
    # Linear function defined above
    model = odr.Model(Linear)
    if len(weights) != len(x):
        data = odr.Data(x, y=y)
    else:
        data = odr.Data(x, y=y, we=weights)
        
    # Bad variable naming, I know...so caps to avoid(?) confusion. :P
    ODR = odr.ODR(data, model, beta0=linreg[0:2])
    
    out = ODR.run()
    return out.beta, out.sd_beta, out.cov_beta, out.delta, out.eps
    
def intSimpsons(function, minX, maxX, *args):
# Function: intSimpsons
#   Description: Integrates a passed function, using Simpson's
#                Rule
#   Inputs: function - a function to integrate of the form:
#                      y = function(x, args)
#           xMin, xMax - The range to integrate (must be finite)
#           args - Additional arguments to pass to the integrand
#
#   Outputs: Total - float value of the integrated function.
#            Error - A remainder/error estimate, based on:
#                   Rn = 1/90*deltaX**5*Ybar
#                       where we are using the average function
#                       value (Ybar)at a given point as a proxy for
#                       the fourth derivative.
#
    if minX-maxX == 0.:
        return 0., 0.
    # We'll use 100 points for now. We may want to use
    # fewer for speed, later?
    numPoints = 100
    deltaX = (float(maxX)-float(minX))/float(numPoints)
    x = np.arange(float(minX), maxX+deltaX, deltaX)
    
    y = function(x, *args)
    try:
        simpsons = [deltaX*(y[index]+4*y[index+1]+y[index+2])/3. for index in range(0,numPoints,2)]
    except:
        print (minX, maxX, deltaX)
        print (x)
        print (y)
    return sum(simpsons), (deltaX**5)*np.mean(simpsons)/90.
    
def normalize(xArray):
    minX = min(xArray)
    maxX = max(xArray)
    deltaX = maxX-minX
    if deltaX == 0.:
        return np.array([1.0 for x in xArray])
    else:
        return np.array([(x-minX)/deltaX for x in xArray])
    

def boxcar(x, size=3):
    # Performs a boxcar smoothing of the passed array. Note the resulting array
    # will be shortened by 2*(size-2) members.
    cs = np.cumsum(np.insert(x, 0, 0)) 
    return (cs[size:] - cs[:-size]) / float(size)
    
# Utilities
def errorOut(errorStr, quit=False, verbose=False):
    if verbose: print (errorStr)
    if quit: os.exit()

def clearFiles(fileList):
# Deletes the passed files
    for thisFile in fileList:
        clearFile(thisFile)
    return

def clearFile(fileName):
    try:
        if os.path.isfile(fileName): os.remove(fileName)
    except OSError:
    # "No such file or directory"
        pass # That's fine - We just won't delete it, then
    

def in_range(num, lo, hi):
# Returns whether num is in the range of (lo, hi)
    return (num >= lo and num <= hi)

def is_list(testList):
# We have run into problems where a list of lists (of floats) with only one
# item turn into a single list of (floats). This tests to see if we have a
# list of floats, or something else.
    try:
        if testList is None or len(testList) == 0:
            return False
    except TypeError:
        # You can't take the len of a float (or int):
        return False
        
    try:
        temp = testList[0]
        return True
    except TypeError:
        # You can't index into a non-list
        return False

def is_number(testString):
    # A little trick for nan: (nan == nan) = False, but we're not using that here...
    # even though it's probably faster than two attempted float conversions...
    
    if testString is None:
        return False
# Note: "real" nan values are actually floats, and won't be passed as strings.
# Basically, this means that EVERY call to this function (with a string parameter)
# will hit the following exception.  
    try:
        if np.isnan(testString):
            return False
    except TypeError:
        pass
 
    try:
        float(testString)
    except ValueError:
        return False
  
    if np.isnan(float(testString)):
        return False
        
    return True
    
    
def closest(value, valueList):
    return min(valueList, key=lambda x:abs(x-value))

def count(item, itemList):
    return len([x for x in itemList if x==item])

def closestIndex(passedRange, value):
    # As index, but will return the index of the nearest value to the passed.
    # Just in case we get a numpy array:
    sortedRange = list(passedRange)
    if value >= sortedRange[-1]: return -1
    if value <= sortedRange[0]: return 0
    return sortedRange.index(closest(value, sortedRange))
    

def bracket(value, sortedRange):
    closestIdx = closestIndex(sortedRange, value)
    if closestIdx==-1: 
        return -1, -1
    elif closestIdx==0: 
        return 0, 0
    elif sortedRange[closestIdx] >= value:
        # We know that the index is not zero, from above
        return closestIdx-1, closestIdx
    else:
        # ...also that it is not the last index.
        return closestIdx, closestIdx+1


def linlogfraction(lolog, hilog, desired):
    return((10**desired-10**lolog)/(10**hilog-10**lolog))
    

def find(f, seq):
# find function from Ryan Tomayko (http://tomayko.com/cleanest-python-find-in-list-function)
# Note: Link is now broken (2-2017)
# This function is preferable to filter or reduce, as it will only traverse the list
# until it finds the item.
#
# The call requires the use of a passed function, most commonly a "lambda"
# expression. For example, if I have a list of "people" data, which looks like:
# [['Smith', '123-45-6789', '123 Main St.',...],...], I can search for the first
# "person" in the list on "Main St." by the following call:
# find(lambda p: 'Main St.' in p[2], peopleList)

    for item in seq:
        if f(item):
            return item
    # Traversed the list, no successes so return None
    return None

def findAll(f, seq):
# Same function as above, but returns a list of _all_ items in the list which 
# match the passed condition

    returnList = []
    for item in seq:
        if f(item):
            returnList.append(item)
    # Traversed the list, no successes so return None
    return returnList


def hrs2degs(hourString):
# Converts a (RA) hh:mm:ss.ss coordinate into a degree value
    try:
        hh, mm, ss = [float(x) for x in hourString.replace(':',' ').split()]
        assert in_range(hh, 0., 24.) and in_range(mm, 0., 60.) and in_range(ss, 0., 60.)
        dd = hh*15. + mm/4. + ss/240.
        return dd
    except ValueError:
        errorOut('Hour conversion failed (incorrect number of elements).')
        return -999.999
    except AssertionError:
        raise ValueError('Hour conversion failed (element out of range).')
    
def degs2digit(degString):
# Converts a dd:mm:ss.ss coordinate into a digital degree value
    try:
        dd, mm, ss = [float(x) for x in degString.replace(':',' ').split()]
        absd = abs(dd)
        sign = dd/absd
        assert in_range(absd, 0., 360.) and in_range(mm, 0., 60.) and in_range(ss, 0., 60.)
        digd = dd + (sign*mm/60.) + (sign*ss/3600.)
        return digd
    except ValueError:
        errorOut('Degree conversion failed (incorrect number of elements).')
        return -999.999
    except AssertionError:
        raise ValueError('Degree conversion failed (element out of range).')
    

def getObsInfo(obs):
    if obs == k.KeckTelescope:
        obsLat = k.Keck1Latitude
        obsLong = k.Keck1Longitude
        obsAlt = k.Keck1Altitude
        obsTZ = k.KeckTimeZone
    else:
        obsLat = k.OHPLatitude
        obsLong = k.OHPLongitude
        obsAlt = k.OHPAltitude
        obsTZ = k.OHPTimeZone
       
    return obsLat, obsLong, obsAlt, obsTZ

def dictMerge(d1, d2, merge_fn = lambda x,y: x+y):
# Function to merge two dictionaries, non-destructively
# Taken from: 
# stackoverflow.com/questions/38987/how-can-i-merge-two-dictionaries-in-a-single-expression
# Note that the default behavior is to add the two values, in the case of a collision
    result =dict(d1)
    for k,v in d2.iteritems():
        if k in result:
            result[k] = merge_fn(result[k], v)
        else:
            result[k] = v
    return result
    

def mergeSplitStrings(badSplitList):
# Function takes a passed list of split strings, from a .csv file, and re-splits 
# them under the assumption that there was a bad splitting, due to quoted commas.
# ie. Given a data line from a .csv file of: 
#   7771.947,8.0,9.140,0.3598,2.50,"2,6,11,12","Loggf: (2:0.25, 6:0.369,11:0.46,12:0.36)",,,,1408.167
# str.split(',') will yield:
#   [7771.947, 8.0, 9.140, 0.3598, 2.50, '"2', 6, 11, '12"', '"Loggf: (2:0.25', '6:0.369', '11:0.46', '12:0.36)"', '', '', '', 1408.167]
# (note how the string "2,6,11,12" was split into 4 entries)
    foundStringStart = False
    newList = []
    newElement = ''
    for elem in badSplitList:
        if foundStringStart:
            newElement += ',{0}'.format(elem)
            if elem.strip()[-1] == '"':
                foundStringStart = False
                newList.append(newElement[:-1])
                newElement = ''
        else:
            if len(elem) > 0 and elem.strip()[0] == '"':
                foundStringStart = True
                newElement += '{0}'.format(elem[1:])
            else:
                newList.append(elem)
                
#    print("Merged: {0}".format(newList))
    return newList
    
    
def tupleList2ListList(tupleList):
# The Python database design API requires that the "fetchall" function returns a list of tuples.
# This is a utility to convert a passed list of tuples to a list of lists.
# This is O(cN) complexity, where c will be the size of each tuple. This can become ugly if the 
# tuples are large!
    return [[element for element in tupleItem] for tupleItem in tupleList]
    
def EdvardssonVTurb(Teff, LogG):
# This uses the relationship from Edvardsson (A&A 275, p101-152, 1993) to 
# calculate a microturbulent velocity given Teff, and LogG. 
# Metallicity, (expressed as [Fe/H] relative to solar) is not a direct factor,
# but limits the range where this calculation is valid.
# Technically, the valid ranges for this calculation are:
# 5550 < Teff < 6800
# 3.50 < LogG < 4.5
# -1.10 < [Fe/H] < 0.3
    VTurb = 1.25+8*(10**-4)*(Teff-6000) - 1.3*(LogG-4.5)
    if VTurb < 0.:
        VTurb = 0.
    return VTurb

def isGiant(Teff, Logg):
# Eventually, this function should do a photometric determination...
    if Teff < k.giantTeffLimit and Logg < k.giantLogGLimit:
        return True
    else:
        return False
        

def angsep(ra1deg,dec1deg,ra2deg,dec2deg):
# Determine separation in degrees between two celestial objects 
# arguments are RA and Dec in decimal degrees. 
# From angsep.py Written by Enno Middelberg 2001
# Blatantly "Google-grabbed" from: http://www.stsci.edu/~ferguson/software/pygoodsdist/pygoods/angsep.py

    ra1rad=ra1deg*np.pi/180.
    dec1rad=dec1deg*np.pi/180.
    ra2rad=ra2deg*np.pi/180.
    dec2rad=dec2deg*np.pi/180.
    
    # calculate scalar product for determination
    # of angular separation
    
    x = np.cos(ra1rad) * np.cos(dec1rad) * np.cos(ra2rad) * np.cos(dec2rad)
    y = np.sin(ra1rad) * np.cos(dec1rad) * np.sin(ra2rad) * np.cos(dec2rad)
    z = np.sin(dec1rad) * np.sin(dec2rad)
    
    rad = np.arccos(x+y+z) # Sometimes gives warnings when coords match
    
    # use Pythargoras approximation if rad < 1 arcsec
    sep = np.choose( rad<0.000004848 , (
        np.sqrt((np.cos(dec1rad)*(ra1rad-ra2rad))**2+(dec1rad-dec2rad)**2),rad))
        
    # Angular separation
    sep = sep*180./np.pi

    return sep


def make2dKDE(Xs, Ys, bw=None):
# Returns a normalized KDE for the passed array
    
    # No Errs passed, so just use the 'standard' scipy KDE
    values = np.vstack([Xs, Ys])
    kde = stats.gaussian_kde(values, bw_method=bw)
    preSamp = np.reshape(kde(values).T, np.array(Xs).shape)
    samples = preSamp.flatten()
    
    # Don't want to risk re-calculating this every time, so...
    maxSamp = max(samples)
    minSamp = min(samples)
    
    normKDE = lambda x: (kde(x)-minSamp)/(maxSamp-minSamp)
#    plot2dKDE(normKDE, min(Xs)-1.0, max(Xs)+1.0, min(Ys)-1., max(Ys)+1.)
    

    return normKDE


def plot2dKDE(kde, xmin, xmax, ymin, ymax, title='',\
        xaxislabel='', yaxislabel='', outfilename=None, overplotPoints=None):
    # 100x100 should be fine (for now?)
    try:
        Xs, Ys = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    except AttributeError:
        print(xmin,xmax,ymin,ymax)
        input()
    pointPos = np.vstack([Xs.ravel(), Ys.ravel()])
    kdearray = np.array([kde(point) for point in pointPos.T])

    Zs = np.reshape(kdearray, Xs.shape)
    
    fig = pyplot.figure()
    ax = fig.gca()
    
    # Doesn't matter if we set the labels to empty strings
    ax.set_xlabel(xaxislabel, fontsize=18)
    ax.set_ylabel(yaxislabel, fontsize=18)
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.set_title(title)
    

    if xmin is not None and xmax is not None:
        ax.set_xlim(xmin, xmax)
    if ymin is not None and ymax is not None:
        ax.set_ylim(ymin, ymax)
    
    colMap = ax.contourf(Xs, Ys, Zs, cmap='Blues')
    contours = ax.contour(Xs, Ys, Zs, colors='k')
    ax.clabel(contours, inline=1, fontsize=10)
    if overplotPoints != None:
        ax.plot(overplotPoints[0], overplotPoints[1], 'r.')
    if outfilename is not None:
        pyplot.savefig(outfilename, figure=fig)
    else:
        pyplot.show()
    
    pyplot.close()
    return
    
    
def make1dKDE(Xs, errs=None):
# As above, but for an array (1-d)
    errarray = errs
    if errarray is not None and len(errarray)!=len(Xs):
        print('Error! Cannot make KDE with len(X errors)!=len(X array)')
        errarray = None

    if errarray is None:
        kde = stats.gaussian_kde(Xs)
        sig = np.std(Xs)
        testVals = np.arange(min(Xs)-4*sig,max(Xs)+4.*sig, (max(Xs)-min(Xs))/1000.)
        samples = np.array([kde(s)[0] for s in testVals])
        # Set these up as constants, so there's no risk of continual re-calculation
        maxSamp = max(samples)
        minSamp = min(samples)
        returnKDE = lambda x: (kde(x)[0]-minSamp)/(maxSamp-minSamp)
    else:
    # Warning: This is horribly inefficient, like O(n^2), because we are
    # creating a Gaussian at every point, and then summing all the results
        pointArray = zip(Xs, errarray)    # Separated for readability
        kde = lambda x: np.sum([stats.norm(p[0], p[1]).pdf(x) for p in zip(Xs, errarray)])
        testVals = np.arange(min(Xs)-4.*max(errs),max(Xs)+4.*max(errs), (max(Xs)-min(Xs))/1000.)
        samples = np.array([kde(s) for s in testVals])
        # Set these up as constants, so there's no risk of continual re-calculation
        maxSamp = max(samples)
        minSamp = min(samples)
        normSamples = np.array([(kde(x)-minSamp)/(maxSamp-minSamp) for x in testVals])

        returnKDE = lambda x: np.interp(x, testVals, normSamples, left=0., right=0.)
    
    return returnKDE

def plot1dKDE(kde, xmin, xmax, title='', xaxislabel='', outfilename=None, plotGaussian=False, parmsOnly=False):
    if xmin >= xmax:
        return
        
    unfilteredXs = np.linspace(xmin, xmax, 100)
    unfilteredYs = np.array([kde(x) for x in unfilteredXs])

    plotXs = [p0 for (p0,p1) in zip(unfilteredXs, unfilteredYs) if not np.isinf(p1) and not np.isnan(p1)]
    plotYs = [p1 for (p0,p1) in zip(unfilteredXs, unfilteredYs) if not np.isinf(p1) and not np.isnan(p1)]
    fittedParms = [0,0,0]

    if plotGaussian:
        gaussian = lambda x, peak, mu, sig: Gaussian(x, peak=peak, mu=mu, sigma=sig)
        initVals = [max(plotYs), np.mean(plotXs), np.std(plotXs)]
        try:
            fittedParms, covar = opt.curve_fit(gaussian, plotXs, plotYs, p0=initVals)
        except (RuntimeError, ValueError):
        # RuntimeError: non-Gaussian-able plots, so just give up
        # ValueError: nans in the arrays
            return fittedParms
            
        gaussYs = gaussian(plotXs, fittedParms[0], fittedParms[1], fittedParms[2])
        
    fig = pyplot.figure()
    ax = fig.gca()
    
    # Doesn't matter if we set the labels to empty strings
    ax.set_xlabel(xaxislabel, fontsize=12)
    ax.set_ylabel('Score', fontsize=15)
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.set_title(title)
    
    # The default (I think), but meh...
    if xmin is not None and xmax is not None:
        ax.set_xlim(xmin, xmax)
    pyplot.plot(plotXs, plotYs, c='k')
    if plotGaussian:
        pyplot.plot(plotXs, gaussYs, 'r--')
        pyplot.axvline(x=fittedParms[1], color='r', linestyle=':', label='Median={0:3.1f}'.format(fittedParms[1]))
        pyplot.axvline(x=fittedParms[1]+k.sigmaFor90*fittedParms[2], color='b', linestyle=':', label=r'$10^{{th}}/90^{{th}}$ percentile ({0:3.1f}  {1:3.1f})'.\
                format(fittedParms[1]-k.sigmaFor90*fittedParms[2], \
                       fittedParms[1]+k.sigmaFor90*fittedParms[2]))
        pyplot.axvline(x=fittedParms[1]-k.sigmaFor90*fittedParms[2], color='b', linestyle=':')
        pyplot.legend()
        
    if outfilename is not None:
        pyplot.savefig(outfilename, figure=fig)
    elif not parmsOnly:
        pyplot.show()
    
    pyplot.close()
    
    return fittedParms
    
