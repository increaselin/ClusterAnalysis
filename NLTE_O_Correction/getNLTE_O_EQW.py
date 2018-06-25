#! /usr/bin/env python
import sys
import os
import sqlite3
import math
from collections import deque

from constants import constants as k

# Function getNLTE_O_EQW
#
# Accepts: (float: line, float: teff, float: logg, float: micro_v, float: eqw)
#		line : The O line to correct (eg: 7771.94 Angstrom)
#		teff : Model effective temperature (eg: 6000.0 Angstrom)
#		logg : Model log(g) value (eg: 4.40)
#		micro_v : Model microturbulent velocity (eg:  1.50 km/s)
#		eqw: Measured line width to adjust (eg: 100 mAngstrom)
#
# Requires: NLTEOxyCorrectionDB - a sqlite3 database, built from the Tanaka 2003 data
#
# Returns: (real: eqw)
#		eqw : The eqw value, correcte for NLTE effects
#

def getNLTE_O_EQW(line, teff, logg, micro_v, eqw) :
	connection = sqlite3.connect(k.NLTEODB)
	cursor = connection.cursor()
	
	if math.floor(line) == 7771:
		lookupLine = 7771.94
	elif math.floor(line) == 7774:
		lookupLine = 7774.17
	elif math.floor(line) == 7775:
		lookupLine = 7775.39
	else:
		return(eqw)

	if teff < 5000:
		lowT = 4500
		hiT = 5000
	elif teff > 6000:
		lowT = 6000
		hiT = 6500
	elif teff > 5500:
		lowT = 5500
		hiT = 6000
	else:
		lowT = 5000
		hiT = 5500
	
	if logg < 2.0:
		lowLogg = 1.0
		hiLogg = 2.0
	elif logg < 3.0:
		lowLogg = 2.0
		hiLogg = 3.0
	elif logg < 4.0:
		lowLogg = 3.0
		hiLogg = 4.0
	else:
		lowLogg = 4.0
		hiLogg = 5.0
	
	
	if micro_v < 2.0:
		lowMicroV = 1
		hiMicroV = 2
	else:
		lowMicroV = 2
		hiMicroV = 3
	
	# Interpolation function. We need to interpolate between 8 values in the table.
	# The interpolation points are denoted by point_XXX, where 'X' is one of: {H L C}.
	# H, L and C represent High, Low and Correct (interpolated) values for the quantity
	# denoted by the position of the letter. The positions go: Temperature, Log(g), Vel.
	# So, point_HCL would be the high temperature, correct (interpolated for) Log(g) and Low micro. vel.
	
	# Low T, Low Log(g) Interpolate V
	cursor.execute('''select * from NLTECorrectionsTable where line=? and teff=? and logg=?''', (lookupLine,lowT,lowLogg))
	matchedRecord = cursor.fetchall()
	thisPoint = matchedRecord[0]
	point_LLL = [thisPoint[4*lowMicroV-1], thisPoint[4*lowMicroV]]
	point_LLH = [thisPoint[4*hiMicroV-1], thisPoint[4*hiMicroV]]
	deltaA = point_LLH[0] - point_LLL[0]
	deltaB = point_LLH[1] - point_LLL[1]
	point_LLC = [point_LLL[0]+deltaA*(micro_v-float(lowMicroV)), point_LLL[1]+deltaB*(micro_v-float(lowMicroV))]
		
	# Low T, High Log(g) Interpolate V
	cursor.execute('''select * from NLTECorrectionsTable where line=? and teff=? and logg=?''', (lookupLine,lowT,hiLogg))
	matchedRecord = cursor.fetchall()
	thisPoint = matchedRecord[0]
	point_LHL = [thisPoint[4*lowMicroV-1], thisPoint[4*lowMicroV]]
	point_LHH = [thisPoint[4*hiMicroV-1], thisPoint[4*hiMicroV]]
	deltaA = point_LHH[0] - point_LHL[0]
	deltaB = point_LHH[1] - point_LHL[1]
	point_LHC = [point_LHL[0]+deltaA*(micro_v-float(lowMicroV)), point_LHL[1]+deltaB*(micro_v-float(lowMicroV))]

	# Low T, Interpolate Log(g)) (interpolated V)
	deltaA = point_LHC[0] - point_LLC[0]
	deltaB = point_LHC[1] - point_LLC[1]
	point_LCC = [point_LLC[0]+deltaA*(logg-lowLogg), point_LLC[1]+deltaB*(logg-lowLogg)]

	# High T, Low Log(g), Interpolate V
	cursor.execute('''select * from NLTECorrectionsTable where line=? and teff=? and logg=?''', (lookupLine,hiT,lowLogg))
	matchedRecord = cursor.fetchall()
	thisPoint = matchedRecord[0]
	point_HLL = [thisPoint[4*lowMicroV-1], thisPoint[4*lowMicroV]]
	point_HLH = [thisPoint[4*hiMicroV-1], thisPoint[4*hiMicroV]]
	deltaA = point_HLH[0] - point_HLL[0]
	deltaB = point_HLH[1] - point_HLL[1]
	point_HLC = [point_HLL[0]+deltaA*(micro_v-float(lowMicroV)), point_HLL[1]+deltaB*(micro_v-float(lowMicroV))]
	
	# High T, High Log(g), Interpolate V
	cursor.execute('''select * from NLTECorrectionsTable where line=? and teff=? and logg=?''', (lookupLine,hiT,hiLogg))
	matchedRecord = cursor.fetchall()
	thisPoint = matchedRecord[0]
	point_HHL = [thisPoint[4*lowMicroV-1], thisPoint[4*lowMicroV]]
	point_HHH = [thisPoint[4*hiMicroV-1], thisPoint[4*hiMicroV]]
	deltaA = point_HHH[0] - point_HHL[0]
	deltaB = point_HHH[1] - point_HHL[1]
	point_HHC = [point_HHL[0]+deltaA*(micro_v-float(lowMicroV)), point_HHL[1]+deltaB*(micro_v-float(lowMicroV))]

	# High T, Interpolate Log(g) (interpolated V)
	deltaA = point_HHC[0] - point_HLC[0]
	deltaB = point_HHC[1] - point_HLC[1]
	point_HCC = [point_HLC[0]+deltaA*(logg-lowLogg), point_HLC[1]+deltaB*(logg-lowLogg)]

	# Interpolate T, (interpolated Log(g)) (interpolated V)
	deltaA = (point_HCC[0] - point_LCC[0])/(hiT-lowT)
	deltaB = (point_HCC[1] - point_LCC[1])/(hiT-lowT)
	point_CCC = [point_LCC[0]+deltaA*(teff-lowT), point_LCC[1]+deltaB*(teff-lowT)]

#	I think the following is wrong - for two reasons:
#	One: This should be an Abundance correction, not a eqw correction.
# 	Two: The correction formula is: (Delta) = A(NLTE) - A(LTE) = a*10**(b*eqw)
#	this_eqw = eqw + point_CCC[0] * (10**point_CCC[1]) * eqw
	abundCorr = point_CCC[0] * (10**(point_CCC[1] * eqw))
	
	return(abundCorr)
