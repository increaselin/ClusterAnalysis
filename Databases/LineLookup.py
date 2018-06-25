#! /usr/bin/env python
# 
# Module: LineLookup.py
#
# Creator: Mike Lum
#
# Description: Functions to search the atomic line database
#
# Contents: List externally-callable functions here
#   containsHFS: (Description)
#
#   getRefData: (Description)
#
#   renumRefs: (Description)
#
#   addLines: (Description)
#
#   getLines: (Description)
#
#   getSolarLines: (Description)
#
#   EnterCalibrationAbundances: (Description)
#
#   LookupSolarAbundances: (Description)
#
#   getCalibrationLines: (Description)
#
#   measLinesToMOOG: (Description)
#
#   getLinesForSpectrum: (Description)
#
#   getLinesForStar: (Description)
#
#   NumLinesOfIon: (Description)
#
# Revision History:
#    Date        Vers.    Author        Description
#    2018-04-05  1.0f1    Lum       github code cleanup.
#    Today       0.0a0    Lum       First checked in
#
# To Do:
#    
#

from itertools import chain
import numpy as np
import sqlite3
import re

from constants import constants as k
from constants import DBFormats as DB
from utilities import utilities as u
    
# Functions
def containsHFS(lineList):
# Checks the passed line list against our blended line (HFS) DB,
# and returns True/False, based on the result.
# The passed linelist is assumed to have the form:
# [ion, expected wavelength, ...]
    with sqlite3.connect(DB.ClusterDBName) as DBConnection:
        DBcursor = DBConnection.cursor()
        
        for line in lineList:
            DBcursor.execute(\
            '''SELECT * FROM {0} WHERE Wavelength IS {1:4.3f} AND Ion IS {2:2.1f};'''.\
                format(DB.LineTableName, line[1], line[0]))
            lineData = DBcursor.fetchall()
            # We should have receive either no or one result, since the 
            # combination of Ion and WL should be unique
            if len(lineData) > 0:
                if lineData[0][10] == 1:
                    return True
    return False

def GetHFSBlends(BlendCoG, ion, MOOGFormat=False):
# Looks for the blend components (assumed HFS, only for now) for the 
# passed center of gravity. Results are returned in the format:
# [(Blend WL, Log(gf)), ...]
# Currently, the MOOGFormat parameter is not used, but the intent
# is to allow a return string format, suitable for directly writing
# to a MOOG line list file.
    with sqlite3.connect(DB.ClusterDBName) as DBConnection:
        DBCursor = DBConnection.cursor()
        
        DBCursor.execute('''SELECT BlendWL, Loggf FROM {0} WHERE Wavelength IS {1:4.3f} AND Ion IS {2:2.1f} ORDER BY BlendWL;'''.format(DB.HFSTableName, BlendCoG, ion))
        tempData = DBCursor.fetchall()
        lineData = [[float(line[0]), float(line[1])] for line in tempData]
    return lineData


   

def getRefData(refList, DBcursor):
# Returns the full list of reference data,
# in an alphabetized (1st Author+year) list
    # Get the full list of references from the passed list.
    # Resolve secondary (in the form "XX(YY,ZZ,...)" references, as well.
    keyStrList = set(chain.from_iterable([item.replace(')','').\
                replace('(',',').split(',') for item in refList]))
    keyList = sorted([int(item) for item in keyStrList])
    
    if len(keyList) > 0:
        DBcursor.execute(\
        '''SELECT * FROM {0} r WHERE r.RefNo IN ({1}) ORDER BY r.RefNo;'''.\
            format(DB.LineRefTableName, ','.join([str(key) \
                                                for key in keyList])))
        tempList = DBcursor.fetchall()
        return sorted(tempList, key=lambda refdata: refdata[1]+refdata[2])
    else:
        return {}

def renumRefs(lineData, refLUT):
# Function takes a passed mapping (LUT) of old#->new# for references,
# and renumbers the references field in the passed line reference list.
# A copy of the list, with the renumbered references is returned.
    # COPY the old list
    newList = [line for line in lineData]

    for line in newList:
        oldRefs = re.split("([\(\),])",line[5])
        newRefs = []
        for item in oldRefs:
            try:
                newRefs.append('{0:d}'.format(refLUT[int(item)]))
            except ValueError:
                newRefs.append(item)

        line[5] = ''.join(newRefs)
    return newList
    

def addLines(newLines):
# Inserts the passed list of lines into the "Measured_Lines" table.
# newLines list format is assumed to be:
# ['Ion', 'Wavelength', 'Filename', 'MeasWL', 'EQM', 'FWHM']
# Note that it is assumed that the line WL/Ion combo is already entered in the
# Line data table, and the filename is in the "Spectra" table
    tablename = DB.MeasuredLinesTableName
    tableFormat = DB.MeasuredLinesGenFormat
    
    with sqlite3.connect(DB.ClusterDBName) as DBconnection:
        DBcursor = DBconnection.cursor()        

        DBcursor.executemany(\
            '''INSERT INTO \"{0}\" values{1}'''.\
                    format(tablename, tableFormat), newLines)
        DBconnection.commit()
        DBcursor.close()


def getLines(elements=[], wavelengthRange=[], dataFormat='MOOG', 
                                                minRefs=1, maxEP=0):
# Function returns a list of lines (and parameters) which match the passed
# selection criteria. Setting the 'dataformat' parameter to anything but 'MOOG'
# will cause the line data to be returned in the raw "fetchall" tuple format.
# The 'MOOG' data format will be text lines, suitable for use as a MOOG line
# list.
# Note: An empty element list will return lines for all elements, and an empty
# wavelength range list will return lines in the full range (2000, 10000).

# The Python database design API requires that a "fetchall" function return a
# list of tuples. We honor that requirement in our API, and leave it to the 
# consumer to manipulate the tuples in whatever way they need.

    ReturnLines = []

    with sqlite3.connect(DB.ClusterDBName) as DBConnection:
        DBcursor = DBConnection.cursor()
        lineData = []

        if len(wavelengthRange) == 0:
            wavelengthRange = [(2000.,10000.)]
        for lo, hi in wavelengthRange:
            if len(elements) == 0:
                DBcursor.execute(\
                '''SELECT * FROM {0} l WHERE l.Wavelength BETWEEN {1:8.3f} AND {2:8.3f} ORDER BY l.Ion, l.Wavelength;'''.\
                    format(DB.LineTableName, lo, hi))
                lineData.extend(DBcursor.fetchall())
            else:
                DBcursor.execute(\
                '''SELECT * FROM {0} l WHERE l.Ion IN ({1}) AND l.Wavelength BETWEEN {2:8.3f} AND {3:8.3f} ORDER BY l.Ion, l.Wavelength;'''.\
                    format(DB.LineTableName, ','.join([str(elem) \
                                            for elem in elements]), lo, hi))
                lineData.extend(DBcursor.fetchall())
        # Produce a dictionary matching all the references used in our line
        # list to the "LineRef" table.
        # Note, we also re-number the references from 1..N
        ReturnRefs = sorted(set(chain.from_iterable([line[5].split(',') for line in lineData])))
        refList = getRefData(ReturnRefs, DBcursor)
        
        DBcursor.close()
        if dataFormat=='MOOG':
            ReturnLines = [k.MOOGListFormat.format(line[0], line[1], line[2],\
                        line[3], line[4], 99.9, line[5]) for line in lineData]
        else:
        # Raw data
            ReturnLines = lineData

    return ReturnLines, refList

def getSolarLines(elements=[], wavelengthRange=[], useDAOSpec=False):
    if useDAOSpec:
        return getCalibrationLines(starName=k.DAOSolarRefStarName, \
                elements=elements, wavelengthRange=wavelengthRange)
    else:
        return getCalibrationLines(starName=k.SolarRefStarName, \
                elements=elements, wavelengthRange=wavelengthRange)


def EnterCalibrationAbundances(calibLines):
# Enter the passed reference spectrum lines into our DB
# calibLines is assumed to be of the format:
# [[starName, Ion, WL, XP, Loggf, EQW, LogRW, Ab],...]
    with sqlite3.connect(DB.ClusterDBName) as DBConnection:
        DBcursor = DBConnection.cursor()
    
        DBcursor.executemany('''INSERT INTO \"{0}\" values{1}'''.format\
                    (DB.RefCorrTableName, DB.RefCorrGenFormat), calibLines)
            
        DBConnection.commit()
        DBcursor.close()

def LookupSolarAbundances(ionList=[]):
# We have potentially pre-calculated abundances for a Solar reference spectrum
# If they exist, then return them.
    calibLines = {}

    with sqlite3.connect(DB.ClusterDBName) as DBConnection:
        DBcursor = DBConnection.cursor()
    
        if len(ionList) > 0:
        # We want lines from particular ions
            for ion in ionList:
                DBcursor.execute(\
                '''SELECT Wavelength, ExPot, Loggf, EQW, LogRW, Abundance FROM {0} WHERE ReferenceStar is \"{1}\" AND Ion is {2:2.1f} ORDER BY Wavelength;'''.\
                    format(DB.RefCorrTableName, k.SolarRefStarName, ion))
                lineData = DBcursor.fetchall()
                if len(lineData) == 0:
                    continue
                # If we receive only a single line back, then we need to
                # make this a list of the one item:
                if u.is_list(lineData[0]):
                    calibLines[ion] = lineData
                else:
                    calibLines[ion] = [lineData]
        else:
        # We want all lines
            DBcursor.execute(\
            '''SELECT Ion, Wavelength, ExPot, Loggf, EQW, LogRW, Abundance FROM {0} WHERE ReferenceStar is \"{1}\" ORDER BY Ion, Wavelength;'''.\
                format(DB.RefCorrTableName, k.SolarRefStarName))
            dbData = DBcursor.fetchall()
            if len(dbData) == 0:
                return calibLines
            # Conceivably, we could receive only one line back, which would
            # not be the list of lists we are expecting...It's probably better
            # to NOT handle that case, and let the crash occur, since there's
            # probably something wrong already.
            lineData = np.array(dbData)
            ions = sorted(set(lineData[:,0]))
            for ion in ions:
                calibLines[ion] = [l[1:] for l in lineData if l[0]==ion]

        DBcursor.close()
        
    return calibLines


def getCalibrationLines(starName=k.SolarRefStarName, elements=[],\
                        wavelengthRange=[]):
# As "GetLines", but for the requested reference star. Since multiple spectra
# go into the EQWs measured for a single line, we return the median EQW
# measurement as the line width.
    lineData = []

    with sqlite3.connect(DB.ClusterDBName) as DBConnection:
        DBcursor = DBConnection.cursor()

        if len(wavelengthRange) == 0:
            wavelengthRange = [(2000.,10000.)]
        for lo, hi in wavelengthRange:
            if len(elements) == 0:
                DBcursor.execute(\
                '''SELECT * FROM {0} l WHERE l.Star is \"{3}\"  AND l.Wavelength BETWEEN {1:8.3f} AND {2:8.3f} ORDER BY l.Ion, l.Wavelength;'''.\
                    format(DB.CalibrationLinesTableName, lo, hi, starName))
                lineData.extend(DBcursor.fetchall())
            else:
                DBcursor.execute(\
                '''SELECT * FROM {0} l WHERE l.Ion IN ({1}) AND l.Star is \"{4}\" AND l.Wavelength BETWEEN {2:8.3f} AND {3:8.3f} ORDER BY l.Ion, l.Wavelength;'''.\
                    format(DB.CalibrationLinesTableName, ','.\
                        join([str(elem) for elem in elements]),\
                                lo, hi, starName))
                lineData.extend(DBcursor.fetchall())
                        
        DBcursor.close()

    return [[l[1], l[0], np.median(eval(l[2])), l[3]] for l in lineData]


def measLinesToMOOG(measLines):
# Convert raw measured line data from the sqlite DB into a 
# list of MOOG-formatted lines.
    lineData = []
    with sqlite3.connect(DB.ClusterDBName) as DBConnection:
        DBcursor = DBConnection.cursor()
        for line in measLines:
            DBcursor.execute(\
            '''SELECT ExPot, Loggf, vdwDamping  FROM {0} WHERE Ion is \"{1}\" and Wavelength is \"{2}\";'''.\
                format(DB.LineTableName, line[0], line[1]))
            DBData = DBcursor.fetchall()
            if len(DBData) == 0:
            # Couldn't find the line
                continue
            # Note: If we have multiple matching lines (we shouldn't), 
            # then just pick the first.
            # We also need to switch the order of the ion, wavelength fields, here
            lineData.append([line[1], line[0], DBData[0][0], DBData[0][1],\
                            DBData[0][2], line[4]])
    # Format is a string with: 
    #               Wavelength; Ion; Ex. Pot.; LogGF; vdw factor; EQW; notes
    # and appropriate spacing for a MOOG line list file
    return [k.MOOGListFormat.format(line[0], line[1], line[2], line[3],\
            line[4], line[5], '') for line in lineData]



def getLinesForSpectrum(spectraFile, elements=[], \
                        wavelengthRange=[], dataFormat=None):
# Function to fetch all the measured lines for a given spectrum file which 
# match the passed criteria (see getLines(....) for more detail).

    # Build the query string
    lineData = []
    shortName = spectraFile.split('/')[-1]
    sqlString = '''SELECT * FROM {0} WHERE Filename IS \"{1}\"'''.\
                    format(DB.MeasuredLinesTableName, shortName)
    # We'll pass the requested ions as arguments to a (?,...)-form query.    
    argList = []
    if len(elements) > 0:
        sqlString = '''{0} AND Ion IN ({1})'''.\
                    format(sqlString, ','.join(['?']*len(elements)))
        argList.extend(elements)

    if len(wavelengthRange) > 0:
        sqlString = '''{0} AND Wavelength BETWEEN {1} AND {2}'''.\
                    format(sqlString, wavelengthRange[0], wavelengthRange[1])
    
    # Query-terminating semicolon is here.
    sqlString = sqlString + ' ORDER BY Ion, Wavelength;'
    
    # Minimize DB usage time. Not that race conditions are an issue, here.
    with sqlite3.connect(DB.ClusterDBName) as DBConnection:
        DBcursor = DBConnection.cursor()
        DBcursor.execute(sqlString, argList)
        tempData = DBcursor.fetchall()
        
    # There may be issues if only one line is found. Not sure if 'fetchall'
    # will return a list of one (line) item, or just a list of the single
    # line's parameters. Well...if we break here with a IndexError, we 
    # know why...
    lineData.extend([[float(line[0]), float(line[1]), line[2], float(line[3]),\
                        float(line[4]), float(line[5])] for line in tempData])
    # lineData format is: [ion, wavelength, filename, meas. wl, eqm, fwhm]

    if dataFormat == 'MOOG':
        # Convert line data to MOOG-readable string:
        
        return measLinesToMOOG(lineData)
    else:
        return lineData
    

def getLinesForStar(cluster, starID, elements=[], wavelengthRange=[],\
                    dataFormat=None, daoLinesOnly=False):
# Function to fetch all the measured lines for a given star which match the
# passed criteria (see getLines(....) for more detail).
# Measurements of the same line from different spectra are averaged
#
# Returns: list with Strings: "WL ion exPot loggf vdw eqw" for "MOOG' format
#          list with lists: [ion, WL, spectrum filename, meas WL, eqw, 
    spectraFiles = []
    with sqlite3.connect(DB.ClusterDBName) as DBConnection:
        tempFiles = []
        DBcursor = DBConnection.cursor()
        DBcursor.execute(\
        '''SELECT Filename FROM {0} WHERE Cluster is \"{1}\" and StarID is \"{2}\";'''.\
            format(DB.SpectraTableName, cluster, starID))
        if daoLinesOnly:
            tempFiles.extend(DBcursor.fetchall())
            spectraFiles = [t for t in tempFiles \
                    if '0752_PLAD_' in t[0] or '0001_Sund' in t[0]]
        else:
            tempFiles.extend(DBcursor.fetchall())
            spectraFiles = [t for t in tempFiles \
                    if '0752_PLAD_' not in t[0] and '0001_Sund' not in t[0]]
#        spectraFiles.extend(DBcursor.fetchall())

    allLines = []    
    for spectrumFile in spectraFiles:
        allLines.extend(getLinesForSpectrum(spectrumFile[0], elements,\
                        wavelengthRange))

    # We now need to merge measurements of the same line from multiple spectra
    sortedLines = []
    measuredIons = set([x[0] for x in allLines])
    for thisIon in measuredIons:
        theWLs = set([x[1] for x in allLines if x[0] == thisIon])
        for thisWL in theWLs:
            # Some quick pre-processing here:
            # Lines must be absorption lines with EQW > 2.9mAngstroms and
            # must not be farther than 0.2Angstroms from the expected position.
            try:
                lines = [x for x in allLines 
                        if x[0] == thisIon and x[1] == thisWL \
                        and x[4] > k.MinEQWMeasurement \
                        and abs(float(x[1])-float(x[3])) < 0.2]
            except TypeError:
                # Dunno why we had non numeric values in the DB, but we did...
                # (emphasis on "did")
                print("Type Error: String in numeric field. WL:{0:4.3f} ".\
                      format(thisWL))
                continue
            numLines = len(lines)
            if numLines == 0: 
                continue
            elif numLines == 1:
                sortedLines.append(lines[0])
            else:
                # Some floats made it into the DB as unicode strings...
                # this the forced conversion
                eqws = [float(x[4]) for x in lines if float(x[4])>0.0]
                eqwMean = np.mean(eqws)
                eqwStd = np.std(eqws)
                goodLines = [x for x in lines 
                            if float(x[4])>0.0 \
                            and abs(float(x[4])-eqwMean) < 1.5*eqwStd]
                if len(goodLines) == 0:
                    continue
                sortedLines.append([\
                goodLines[0][0], \
                goodLines[0][1], \
                k.MutipleFiles, \
                np.mean([float(x[3]) for x in goodLines]), \
                np.mean([x[4] for x in goodLines]), \
                np.mean([x[5] for x in goodLines])])
                # List element format is: 
                # [ion, DB wl, filename/"composite", avg. meas WL, 
                #                              avg. eqw, avg fwhm]
    sortedLines.sort(key=lambda x: x[1])
    sortedLines.sort(key=lambda x: x[0])
    
    if dataFormat == 'MOOG':
        return measLinesToMOOG(sortedLines)
    else:
    # Raw data - don't need the ExPot, etc. values
        return sortedLines    
    

def NumLinesOfIon(lineList, ion):
# Simple utility to count the number of lines for a particular ion
# in the passed list. An ion number of 0 returns a total count.
    if ion < 1.:
        return len(lineList)
        
    return len([line for line in lineList if line[0] == ion])

