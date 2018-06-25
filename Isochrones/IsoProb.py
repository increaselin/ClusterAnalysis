#! /usr/bin/env python
# 
# Program: IsoProb.py
#
# Author: Mike Lum
#
# Usage: N/A (API-based)
#
# Description: Functions to read in a set of isochrones, and calculate the 
#       probability of a Teff, LogG pair lying in a given set.
#       Currently, the only isochrone data returned is in the Teff, LogG space.
#
# Revision History:
#    Date              Vers.    Author        Description
#    2015/04/21        0.0a0    Lum           First checked in
#    2017/07/24        1.0f1    Lum           CMD changed data formats
#
# To Do:
#    
#

from bs4 import BeautifulSoup as BS
import glob

import numpy as np
from scipy import stats, signal
import time
import urllib
from constants import constants as k
from utilities import utilities as u

# Py2->3
if k.isPython3:
    import mechanicalsoup as mech
    from urllib.request import urlopen
else:
    import mechanize as mech
    import urllib2

# Functions

def ReadIsos(inputDirectory='.'):
# Read in a set of isochrones from a passed directory.
# The data is returned in a {Age:[[Teff,LogG],...]} format
    IsoDict = {}
    IsoFiles = glob.glob(inputDirectory+'/*.dat.gz')
    for thisFile in IsoFiles:
        
        try:
            age = float(thisFile.split('/')[-1][2:-9])
        except ValueError:
            # The filename is not of the expected format.
            continue
        data = np.genfromtxt(thisFile, usecols=(5,6), converters={5:lambda x: 10**(float(x or 0.0))})
        IsoDict[age] = data
    
    return IsoDict
    
def BuildProbArray(dataArray):
# Take the passed set of data points, and build a Kernel-Density Estimate
# Note: This should only be a temporary substitute for a real probability
#       space estimator. A KDE observes the frequency of points in
#       a given area, as opposed to a true(-er) estimate of probability.
#       ie: A very "curvey" area of the isochrone (ex: RGB) will have 
#       more data points than a straight region (ex: MS), even though the
#       probability of a star lying on the straight portion may be higher.
#
#       I suspect that the proper way to fix this would be to map a spline
#       to the isochrone, and then re-parameterize it. Then create a new
#       evenly-spaced point map based on arc length.
    kernel = stats.gaussian_kde(dataArray.T)
    return kernel
    
def GetProb(point, probArray):
# Returns a (non-log) probability of the passed point lying on the isochrone(s)
# represented in the passed probability array. Of course, the probability of 
# any single point lying in a probability space is zero, so we are actually 
# returning a sum of the pdf in a space bounded by the passed Teff and LogG,
# based on their units, and the number of points in the KDE space.
# Note: This is currently using a KDE - see above comment.
    teffRange = 10.
    loggRange = 0.02
    prob = np.sqrt(probArray.integrate_box([point[0]-teffRange, point[1]-loggRange],[point[0]+teffRange, point[1]+loggRange]))
#    prob = np.sqrt(probArray.evaluate((point[0], point[1]))[0])
    # Zero is zero, regardless of what Python fp math says.
    if abs(prob) < 10**-15:
        prob = 0.
    return prob
    
    
def GetPARSECZ(FeH):
# The PARSEC form needs metallicity in terms of Z, as opposed to [Fe/H].
# Basic calculations: 
#            X+Y+Z=1.0
#            Z=((1.0-Y)*10^([Fe/H]*A))/((Xsolar/Zsolar)+10^([Fe/H]*A))
#                   Xsolar = 0.73626, Ysolar = 0.2485, Zsolar=0.01524
# PARSEC uses a Y=0.2485+1.78Z relation:
    # A is the ratio: [M/H]/[Fe/H]
    if FeH > 0.3:
    # Maximum acceptable Z for the Parsec system is 0.06
        return 0.06
    else:
        # This relation has no basis in any model, it just fits the 
        # relation "A is between 0.9 and 1.0, with A lower for high [Fe/H]
        A = 1.0 - 0.08*(2.**FeH)

    expFactor = 10**(FeH*A)
    
    # A bit of algebra, and inserting the solar values:
    # Note: Our resulting Z values differ slightly from those in Table 4 of
    # Bressan et al. 2012 - but not so much that we're going to worry about it
    Z = (0.7515*expFactor)/(48.311+2.78*expFactor)
    
    return Z

    
def GetPARSECIsochrone(age, metal):
# Goes to the CMD isochrone generator at: http://stev.oapd.inaf.it/cgi-bin/cmd
# and gets the points for an isochrone with the passed age and metallicity.
# The returned data is a series of points with the format:
#     Z, log(age/yr), M_ini, M_act, logL/Lo, logTe, logG, mag_bolo, U, B, V, R, I, J, H, K, evolutionary stage
# ...which corresponds to the old data format from the CMD website. The new 
# WEBSITE data format (which is converted to the desired format) is:
# Zini Age Mini  Mass   logL    logTe  logg  label   McoreTP C_O  period0 period1 pmode  Mloss  tau1m   X   Y   Xc  Xn  Xo  Cexcess  Z 	 mbolmag  Umag    Bmag    Vmag    Rmag    Imag    Jmag    Hmag    Kmag    Stage


    PARSECSiteURL = 'http://stev.oapd.inaf.it'
    PARSECFormURL = '/cgi-bin/cmd'
    cmdURLHead = PARSECSiteURL + PARSECFormURL
    
    ZStr = '{0:1.4f}'.format(GetPARSECZ(metal))
    AgeStr = '{0:1.3f}e9'.format(age/1000.)
    
    #Py2->3 no more mechanize code:
    if k.isPython3:
        browser = mech.StatefulBrowser()
    else:
        browser = mech.Browser()
        
    browser.set_handle_robots(False)
    browser.open(cmdURLHead)
    
    browser.select_form(nr=0)
    browser['isoc_age'] = AgeStr
    browser['isoc_zeta'] = ZStr
    content = browser.submit()
    
    try:
        resp = content.read()
        # Use beautifulsoup to parse out the data file
        soup = BS(resp, 'html.parser')
        anchors = soup.find_all('a')
        isoDataFileURL = PARSECSiteURL + anchors[0].get('href').split('..')[1]
    except IndexError:
        # We didn't get the result we expected. Wait 10 sec, and try again
        time.sleep(10)
        try:
            resp = content.read()
            # Use beautifulsoup to parse out the data file
            soup = BS(resp, 'html.parser')
            anchors = soup.find_all('a')
            isoDataFileURL = PARSECSiteURL + anchors[0].get('href').split('..')[1]
        except IndexError:
            # Still failed, we're hosed
            print ("Fatal error in loadng isochrone. Exiting.")
            exit()
# Py2-3            
#    # urllib & urllib2 code:
#    urlReq = urllib2.Request(isoDataFileURL)
#    urlRes = urllib2.urlopen(urlReq)
    urlRes = urlopen(isoDataFileURL)
    
    data = urlRes.read()
    isoPoints = [[float(item) for item in line.split() if item != ''] for line in data.split('\n') if len(line) > 0 and line[0] != '#']
    
    oldFormat = [i[0:7]+i[22:32] for i in isoPoints]
    
    return oldFormat

def ExtractPARSECMags(isoPoints):
# Function takes a passed PARSEC isochrone data set, and returns an array of
# magnitude 'dictionaries', one for each data point. Each dictionary has the
# form:
# {'U':<UMag value>, 'B':<BMag value>, ...}
    return [{'UJ':i[8] , 'B':i[9] , 'V':i[10] , 'RJ':i[11] , 'IJ':i[12] , 'JJ':i[13] , 'HJ':i[14] , 'KJ':i[15] } for i in isoPoints]

def GetPARSECStages(isoPoints):
# Analyze a PARSEC isochrone, and return a dictionary with the range(s) for
# Type V (Main sequence), IV (Sub-giant), and III (giant) stars. The directory
# is in the form:
#   {'III':(first idx, last idx), 
#    'IV':(first idx, last idx), 
#    'V':(first idx, last idx)}
# with the indexes feferring to the passes isochrone data. It is assumed that
# the passed isochrone starts at the 'base' (Late M-dwarfs) of the MS, and 
# moves through the end of the appropriate giant phase for this isochrone.
# The isochrone is assumed to be in the same form as GetPARSECIsochrone:
#     Z, log(age/yr), M_ini, M_act, logL/Lo, logTe, logG, mag_bolo, U, B, V, R, I, J, H, K, evolutionary stage
# This should be easy, if the 'stage' field were accurate (as of 7/2017 
# it is not). So, we have to figure it out.
    BminusVs = [i[9]-i[10] for i in isoPoints]
    BoloMags = [i[7] for i in isoPoints]
    stageDict = {}
    # Find the ends of the main sequence. Since we assume that the isochrone
    # starts at the 'bottom', the first index is easy...
    Vfirst = 0
    # We assume that the B-V values are decreasing throughout the MS. This is
    # generally ok, although some 'irregularity' can occur near the MSTO. 
    # Therefore, we can assume that the MSTO is at the minimum B-V value.
    Vlast = np.argmin(BminusVs)
    stageDict['V'] = (Vfirst, Vlast)

    # This is also the start of the subgiant branch.
    IVfirst = Vlast
    # The He core 'flash' signifies the end of the subgiant branch, which 
    # is seen as a local max for B-V, and a local min for bolometric Mag.
    # Note: we use 'minimum' and 'maximum' for bolometric magnitudes in terms
    # of actual brightness, not the numeric value, which is reversed.
    BoloMaxes = signal.find_peaks_cwt([0.-m for m in BoloMags],np.arange(1,10))
    BmVMaxes = signal.find_peaks_cwt(BminusVs, np.arange(1,10))

    # The first max/min combo, after the MSTO is the one we want:
    candidates = [idx for idx in BoloMaxes if idx in BmVMaxes and idx>IVfirst]
    IVlast = candidates[0]
    stageDict['IV'] = (IVfirst, IVlast)

    # While technically, the giant branch starts with the He flash, the 
    # transition is so fast, we're going to start our RGB at the 'cluster'
    # point, which is the first bolometric minimum after the 'flash'
    BoloMins = signal.find_peaks_cwt(BoloMags,np.arange(1,10))
    IIIfirst = [idx for idx in BoloMins if idx>IVlast][0]
    
    # Even though the horizontal branch is part of the giant (III) phase, we're
    # going to say that our RGB ends when the tip is hit.
    IIIlast = [idx for idx in BoloMaxes if idx>IIIfirst][0]
    stageDict['III'] = (IIIfirst, IIIlast)
    
    return stageDict
    
def EvenSpacedIso(IsoPoints, numPoints=200):
# The PARSEC isochrones have points placed along the isochrone in such a way as
# to produce an exact curve as possible to the model. This means that there are
# a higher concentration of points along areas of the curve that are changing 
# rapidly. For probability calculations, we need points that are equidistantly 
# spaced along the curve (in terms of arc length). This function takes a passed
# set of isochrone data and returns data set with the passed number of points
# spaced along the isochrone. Note: for future reference, it would probably be
# more physically appropriate to concentrate more points along the curve, based
# on the (IMF) stellar population.
#
    x = IsoPoints[:,0]
    xmin = min(x)
    x = x-xmin
    xnorm = max(x)
    x = [i/xnorm for i in x]
    
    y = IsoPoints[:,1]
    ymin = min(y)
    y = y-ymin
    ynorm = max(y)
    y = [i/ynorm for i in y]
    
    Dx = np.diff(x)
    Dy = np.diff(y)
    Dl = np.sqrt(Dx**2+Dy**2)
    
    # L contains the arclength coordinate for our given X-Y system
    L = np.cumsum(Dl)
    totalArclen = L[-1] # also = np.sum(Dl)
    
    # S contains the interpolated, evenly spaced arclength coordinates
    S = np.linspace(0., totalArclen, numPoints)
    
    xS = [x[0]]
    yS = [y[0]]
    
    # Remember that x and y are not monotonic, so the usual interpolation
    # won't work...
    for sLen in S[1:-1]:
        loCoord, hiCoord = u.bracket(sLen, L.tolist())
        if (L[hiCoord] - L[loCoord]) != 0.:
            delta = (sLen - L[loCoord])/(L[hiCoord] - L[loCoord])
        elif loCoord == 0 and sLen<L[0]:
            # Want a point before we started measuring
            hiCoord = 1
            delta = (sLen - L[0])/(L[1] - L[0]) # delta will be negative
        else:
            delta = 0.
        loCoord = loCoord+1
        hiCoord = hiCoord+1
        xS.append(x[loCoord]+delta*(x[hiCoord]-x[loCoord]))
        yS.append(y[loCoord]+delta*(y[hiCoord]-y[loCoord]))         
    xS.append(x[-1])
    yS.append(y[-1])
    
    xN = [i*xnorm + xmin for i in xS]
    yN = [i*ynorm + ymin for i in yS]
    
    evenIso = np.array(zip(xN, yN))
    return evenIso
    
# I originally tried this handy function from stackoverflow.com...but it didn't do exactly what I wanted...
## This code was taken from the following stackoverflow post:
## http://stackoverflow.com/questions/19117660/how-to-generate-equispaced-interpolating-values
#    x = IsoPoints[:,0]
#    y = IsoPoints[:,1]
#    xd =np.diff(x)
#    yd = np.diff(y)
#    dist = np.sqrt(xd**2+yd**2)
#    u = np.cumsum(dist)
#    u = np.hstack([[0],u])
#
#    t = np.linspace(0,u.max(),numPoints)
## This is where the stackoverflow function breaks for non-monotomic functions
## You can't do a straight interpolation over a coordinate for the arclength.
#    xn = np.interp(t, u, x)
#    yn = np.interp(t, u, y)
#    
#    newIso = np.array(zip(xn, yn))
#    return newIso

    
# Command line calls for testing only
if __name__ == '__main__':
    # Age in My
    testAge = 2500
    testMetal = -0.10
    
    GetPARSECIsochrone(testAge, testMetal)
