#! /usr/bin/env python
# 
# Program: RamMelT.py
#
# Author: Mike Lum
#
#
# Description: Using the relations from Ramirez & Martinez 2005, 
#    function returns the Teff values for the passed photometry.
#    Note: The Teff values are for Giants. For Dwarf stars, use 
#    the function in CasagrandeT
#
# Revision History:
#    Date        Vers.    Author        Description
#    10/22/2015  1.0f1     Lum        First checked in
#
# To Do:
#    
#


from constants import constants as k
from constants import Photometry as kPhot
import numpy as np


# Functions
def GetColorsForTable(photoTable):
# Given one of our photometric tables, return the set of filters used in 
# determining the parameters using that table
    return sorted(set(np.array([i.split('m') for i in photoTable.keys()]).flatten()))

def getDeltasForTable(colorDict, tableDict, extinct = 0.):
# Given a passed table of parameters, indexed on a color difference, find all
# the matching differences, given the passed set of color measurements
    deltaDict = {}
    passedColors = colorDict.keys() # Assuming a unique set of colors
    
    for diffKey in tableDict.keys():
        first, second = diffKey.split('m')
        if first in passedColors and second in passedColors:
            deltaDict[diffKey] = colorDict[first]-colorDict[second]
            # Adjust for extinction
            deltaDict[diffKey] -= (kPhot.RieLebTable3[first]-kPhot.RieLebTable3[second])*extinct
            
#            continue
    return deltaDict
    
def getMagDeltas(magDict, extinct):
# This is kinda 'dumb' programming...need to replace with more generalized solution
    deltas = {}
    try:
        # Both Casagrande 2010 and R&M '05 use this calculation
        deltas['BmV'] = magDict['B'] - magDict['V'] - \
        (kPhot.RieLebTable3['B'] - kPhot.RieLebTable3['V'])*extinct
    except (KeyError, IndexError):
        deltas['BmV'] = k.NoCoefficientData
    try:
        # Only R&M '05, but used for dwarfs in C10
        deltas['YVmVV'] = magDict['YV'] - magDict['VV'] - \
        (kPhot.RieLebTable3['YV'] - kPhot.RieLebTable3['VV'])*extinct
    except (KeyError, IndexError):
        deltas['YVmVV'] = k.NoCoefficientData
    try:
        # Only R&M '05, but used for dwarfs in C10
        deltas['VVmSV'] = magDict['VV'] - magDict['SV'] - \
        (kPhot.RieLebTable3['VV'] - kPhot.RieLebTable3['SV'])*extinct
    except (KeyError, IndexError):
        deltas['VVmSV'] = k.NoCoefficientData
    try:
        # Only R&M '05
        deltas['BGmVG'] = magDict['BG'] - magDict['VG'] - \
        (kPhot.RieLebTable3['BG'] - kPhot.RieLebTable3['VG'])*extinct
    except (KeyError, IndexError):
        deltas['BGmVG'] = k.NoCoefficientData
    try:
        # Only R&M '05
        deltas['BGmGG'] = magDict['BG'] - magDict['GG'] - \
        (kPhot.RieLebTable3['BG'] - kPhot.RieLebTable3['GG'])*extinct
    except (KeyError, IndexError):
        deltas['BGmGG'] = k.NoCoefficientData
    try:
        # C. '10 and R&M '05
        deltas['VmR'] = magDict['V'] - magDict['R'] - \
        (kPhot.RieLebTable3['V'] - kPhot.RieLebTable3['R'])*extinct
    except (KeyError, IndexError):
        deltas['VmR'] = k.NoCoefficientData
    try:
        # C. '10 and R&M '05
        deltas['VmI'] = magDict['V'] - magDict['I'] - \
        (kPhot.RieLebTable3['V'] - kPhot.RieLebTable3['I'])*extinct
    except (KeyError, IndexError):
        deltas['VmI'] = k.NoCoefficientData
    try:
        # C. '10 and R&M '05
        deltas['RmI'] = magDict[R''] - magDict['I'] - \
        (kPhot.RieLebTable3['R'] - kPhot.RieLebTable3['I'])*extinct
    except (KeyError, IndexError):
        deltas['RmI'] = k.NoCoefficientData
    try:
        # C. '10 and R&M '05
        deltas['BTmVT'] = magDict['BT'] - magDict['VT'] - \
        (kPhot.RieLebTable3['BT'] - kPhot.RieLebTable3['VT'])*extinct
    except (KeyError, IndexError):
        deltas['BTmVT'] = k.NoCoefficientData
    try:
        # C. '10 and R&M '05
        deltas['VmJ'] = magDict['V'] - magDict['J'] - \
        (kPhot.RieLebTable3['V'] - kPhot.RieLebTable3['J'])*extinct
    except (KeyError, IndexError):
        deltas['VmJ'] = k.NoCoefficientData
    try:
        # C. '10 and R&M '05
        deltas['VmH'] = magDict['V'] - magDict['H'] - \
        (kPhot.RieLebTable3['V'] - kPhot.RieLebTable3['H'])*extinct
    except (KeyError, IndexError):
        deltas['VmH'] = k.NoCoefficientData
    try:
        # C. '10 and R&M '05
        deltas['VmK'] = magDict['V'] - magDict['K'] - \
        (kPhot.RieLebTable3['V'] - kPhot.RieLebTable3['K'])*extinct
    except (KeyError, IndexError):
        deltas['VmK'] = k.NoCoefficientData
    try:
        # C. '10 and R&M '05
        deltas['VTmK'] = magDict['VT'] - magDict['K'] - \
        (kPhot.RieLebTable3['VT'] - kPhot.RieLebTable3['K'])*extinct
    except (KeyError, IndexError):
        deltas['VTmK'] = k.NoCoefficientData
    try:
        # C. '10 and R&M '05
        deltas['bmy'] = magDict['b'] - magDict['y'] - \
        (kPhot.RieLebTable3['b'] - kPhot.RieLebTable3['y'])*extinct
    except (KeyError, IndexError):
        deltas['bmy'] = k.NoCoefficientData
    try:
        # Only Casagrande 2010
        deltas['JmK'] = magDict['J'] - magDict['K'] - \
        (kPhot.RieLebTable3['J'] - kPhot.RieLebTable3['K'])*extinct
    except (KeyError, IndexError):
        deltas['JmK'] = k.NoCoefficientData
    try:
        # Only Casagrande 2010
        deltas['VTmJ'] = magDict['VT'] - magDict['J'] - \
        (kPhot.RieLebTable3['VT'] - kPhot.RieLebTable3['J'])*extinct
    except (KeyError, IndexError):
        deltas['VTmJ'] = k.NoCoefficientData
    try:
        # Only Casagrande 2010
        deltas['VTmH'] = magDict['VT'] - magDict['H'] - \
        (kPhot.RieLebTable3['VT'] - kPhot.RieLebTable3['H'])*extinct
    except (KeyError, IndexError):
        deltas['VTmH'] = k.NoCoefficientData

    return deltas


def isInRange(photKey, metal, delta, rangeList):
    try:
        polyRange = [x for x in rangeList[photKey] if metal > x[0][0] and metal < x[0][1]]
    except (KeyError, IndexError):
        return False

# The metallicity should only fall into one range, so...?    
#    if len(polyRange) > 1: return False

    return delta > polyRange[0][1][0] and delta < polyRange[0][1][1]

def getRamMelPolyCorrect(thisKey, metal, polyTable):
# Note: We've already verified that we are in the proper metallicity
# range, so we know that we are above -4.0 for the "else" option.
    if metal>-0.5:
        metalIdx = 0
    elif metal>-1.5:
        metalIdx = 1
    elif metal>-2.5:
        metalIdx = 2
    else:
        metalIdx = 3
        
    try:
        coeffs = polyTable[thisKey][metalIdx]
    except (KeyError, IndexError):
        coeffs = kPhot.EmptyList
    
    return coeffs

def RamMelGiantT(magDict, metal, extinct):
    passedMags = {}
    Teffs = {}
    sigmas ={}
    
    deltas = getMagDeltas(magDict, extinct)
    
    goodKeys = [x for x in deltas.keys() if deltas[x] is not k.NoCoefficientData and isInRange(x, metal, deltas[x], kPhot.RamMelGiantRanges)]
    goodPhot = {k:v for k,v in deltas.iteritems() if k in goodKeys and k in kPhot.RamMelTable3.keys()}
    for thisKey in goodPhot.keys() :
        polyCoeffs = getRamMelPolyCorrect(thisKey, metal, kPhot.RamMelTable5)
        polySum = 0.
        for power, coeff in enumerate(polyCoeffs):
            polySum += coeff*(goodPhot[thisKey])**power
            
        params = kPhot.RamMelTable3[thisKey]
        sigmas[thisKey] = params['sigma']
        
        Teffs[thisKey] = 5040/(params['a0']+params['a1']*goodPhot[thisKey]+params['a2']*(goodPhot[thisKey]**2)+params['a3']*goodPhot[thisKey]*metal+params['a4']*metal+params['a5']*(metal**2)) + polySum
    return Teffs, sigmas

def RamMelDwarfT(magDict, metal, extinct):
    passedMags = {}
    Teffs = {}
    sigmas ={}
    
    deltas = getMagDeltas(magDict, extinct)
    
    goodKeys = [x for x in deltas.keys() if deltas[x] is not k.NoCoefficientData and isInRange(x, metal, deltas[x], kPhot.RamMelDwarfRanges)]
    goodPhot = {k:v for k,v in deltas.iteritems() if k in goodKeys and k in kPhot.RamMelTable2.keys()}
    for thisKey in goodPhot.keys() :
        polyCoeffs = getRamMelPolyCorrect(thisKey, metal, kPhot.RamMelTable4)
        polySum = 0.
        for power, coeff in enumerate(polyCoeffs):
            polySum += coeff*(goodPhot[thisKey])**power
            
        params = kPhot.RamMelTable2[thisKey]
        sigmas[thisKey] = params['sigma']
        
        Teffs[thisKey] = 5040/(params['a0']+params['a1']*goodPhot[thisKey]+params['a2']*(goodPhot[thisKey]**2)+params['a3']*goodPhot[thisKey]*metal+params['a4']*metal+params['a5']*(metal**2)) + polySum
    return Teffs, sigmas


