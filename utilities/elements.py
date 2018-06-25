#! /usr/bin/env python
# 
# Program: elements.py
#
# Author: Mike Lum
#
# Usage: (module calls only)
#
# Description: Basic utilities to look up various elemental data
#
# Revision History:
#    Date          Vers.    Author        Description
#    10-10-2014    0.0a0    Lum        First checked in
#
# To Do:
#    
#

from constants import ptoe_abs as p

# Functions
def STR(loggf, exPot, Teff):
# Calculate the expected "strength" of a line using the relation:
# STR = log(gf) - 5040/Teff * Ex.Pot (most recently, Lawler 2015)
    return loggf - (5040./Teff) * exPot

def symbolLU(symbolStr):
    data = [x for x in p.ptoe if x[p.symbolIdx] == symbolStr]
    if len(data) > 0:
        return data[0]
    else:
        return None

def atnoLU(atNo):
    try:
    # Note: 1-based periodic table, 0-based python array
        return p.ptoe[atNo-1]
    except KeyError:
        return None

def nameLU(nameStr):
    data = [x for x in p.ptoe if x[p.nameIdx].lower() == nameStr.lower()]
    if len(data) > 0:
        return data[0]
    else:
        return None
def ionLU(ionStr):
    try:
        ionization = p.RomanNumbers.index(ionStr)+1
    except KeyError:
        ionization = 0
    return ionization

def getIonName(ionNum):
# Return a string with the element "name" in the form: "XxYY"
# where Xx is the standard symbol for the atom, and YY is the
# (Roman numeral) ionization state
    # Note: You _MUST_ round here for python funkiness
    # ie: for SiII, without rounding, ionIdx = 0.9999999999999964
    ionIdx = round((ionNum-int(ionNum)) * 10.)
    if ionIdx-int(ionIdx*10.) > 0.:
    # Ionization state > 11. Note: this mechanism does not differentiate
    # between ionization state 2 (xx.1) and 11 (xx.10)
        ionIdx = round((ionNum-int(ionNum)) * 100.)
    return p.ptoe[int(ionNum)-1][p.symbolIdx] + p.RomanNumbers[int(ionIdx)]

def getIonState(ionStr):
# Passed a string like: FeXIII, return a tuple with ('Fe','XIII',13)
    # First, assume the first character is part of the At. symbol:
    symbolStr = ionStr[0]
    if ionStr[1].islower():
        symbolStr += ionStr[1]
        ionState = ionStr[2:].strip()
    else:
        ionState = ionStr[1:]
        
    ionization = ionLU(ionState)
    
    return (symbolStr, ionState, ionization)

def ionStrToFloat(ionStr):
# Converts a string like: FeXIII into a float like: 26.13
    symbolStr = ionStr[0]
    if ionStr[1].islower():
        symbolStr += ionStr[1]
        ionState = ionStr[2:].strip()
    else:
        ionState = ionStr[1:]
    retVal = symbolLU(symbolStr)[p.atnoIdx]
    ionVal = ionLU(ionState)-1
    if ionVal > 9:
        retVal += ionVal/100.
    else:
        retVal += ionVal/10.
    
    return retVal
