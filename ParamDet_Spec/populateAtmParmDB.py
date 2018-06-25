#! /usr/bin/env python
# 
# Module: populateAtmParmDB.py
#
# Author: Mike Lum
#
# Description: Functions to search the cluster and spectra DB, and calculate
#           spectroscopically-determined atmospheric parameters for those spectra which
#           do not have them.
#
# Contents:
#   Function xxx: (Description)
#
# Revision History:
#    Date        Vers.    Author        Description
#    09-10-2015  1.0a0    Lum         First checked in
#
# To Do:
#    
#

# Imports
import datetime
import math
import sqlite3
import numpy as np

from constants import constants as k
from constants import DBFormats as DB
from Models import MakeAtmModel as mk
from ParamDet_Spec import SpecParamDet as spd
from utilities import utilities as u

# Utility functions (move to some sort of DB utility file?)
def hasAtmParms(Filename, parmType=''):
    DBconnection = sqlite3.connect(DB.ClusterDBName)
    files = []
    with DBconnection:
        DBcursor = DBconnection.cursor()
        if parmType != '':
            DBcursor.execute('''SELECT * FROM {0} l WHERE l.Filename IS \"{1}\" AND l.ParmSource IS \"{2}\";'''.format(DB.AtmParTableName, Filename, parmType))
        else:
            DBcursor.execute('''SELECT * FROM {0} l WHERE l.Filename IS \"{1}\";'''.format(DB.AtmParTableName, Filename))
            
        files = DBcursor.fetchall()
    
    if len(files) > 0:
        return True
    else:
        return False
        
        
def filterLines(lines, DBcursor):
# Function to read through the spectral lines as measured, and filter
# them, based on delta(lambda)/lambda and difference between measured
# and expected wavelength.
# Filtered lines are returned with Loggf and excitation potential, designed
# for MOOG output
    dbLines = []
    for line in lines:
        # Omit lines with measured centroids too far from the expected WL.
        if abs(line[1] - line[3]) > 0.500:
            continue
        # Omit emission lines, and lines too small to measure:
        if line[4] < 1.0:
            continue
        # Omit lines with log(EQW/WL) > Linear COG limit
        # Note: WL in Angstroms, Delta in mAngstroms
        if math.log(line[4]/line[1]) > k.LinearCOGLimit:
            continue
        
        # Get the line data
        DBcursor.execute('''SELECT ExPot, Loggf, vdwDamping FROM {0} l WHERE l.Wavelength IS {1} AND l.Ion IS {2};'''.format(DB.LineTableName, line[1], line[0]))
        lineData = DBcursor.fetchall()[0]
        
        try:
            dbLines.append([line[1], line[0], lineData[0], lineData[1], lineData[2], line[4]])
        except IndexError:
            continue
            
    lines2return = []
    # For multiple measurements of the same line, use the median value
    lineWLs = set([x[0] for x in dbLines])
    for lineWL in lineWLs:
        measLines = [x for x in dbLines if x[0]==lineWL]
        medianEQW = np.median([x[5] for x in measLines])
        lines2return.append([measLines[0][0], measLines[0][1], measLines[0][2], measLines[0][3], measLines[0][4], medianEQW])
        
    return lines2return

def makeMOOGLogFile(lines, filename = k.MOOGTempLogName):
    with open(filename, 'w') as logFile:
        firstLine = '# temporary MOOG log file built: {0}\n'.format(datetime.datetime.now().isoformat())
        logFile.write(firstLine)
        for line in lines:
            tempStr = k.MOOGLogFormat % tuple(line)
            logFile.write(tempStr+'\n')
    
    return filename
    
    
# Functions
giantStars = ['PLA-350','PLA-356', 'PLA-506','PLA-687','PLA-858','PLA-1089','PLA-1172']
def populateAtmTable():
    # Loading the atmospheric models (takes a while, but we're gonna need 'em)
    u.errorOut('Reading plane-parallel model files:\n')
    dModelPath = '/netdisks/galileo18/mikelum/Dropbox/CodeCloud/MyTools/ClusterAnalysis/Models/MARCS_PlaneParallel'
    dModelFiles = mk.findDataFiles(dModelPath)
    dKuruczAtms, dPradks = mk.LoadModels(dModelFiles)
    
    u.errorOut('Reading spherical model files:\n')
    gModelPath = '/netdisks/galileo18/mikelum/Dropbox/CodeCloud/MyTools/ClusterAnalysis/Models/MARCS_Spherical'
    gModelFiles = mk.findDataFiles(gModelPath)
    gKuruczAtms, gPradks = mk.LoadModels(gModelFiles)
    
    DBconnection = sqlite3.connect(DB.ClusterDBName)
    files = []
    with DBconnection:
        DBcursor = DBconnection.cursor()
        DBcursor.execute('''SELECT DISTINCT Filename FROM {0};'''.format(DB.SpectraTableName))
        filenames = [x[0] for x in DBcursor.fetchall()]
        spectraNames = []
        spectraLines = {}
        
        for filename in filenames:
            DBcursor.execute('''SELECT Cluster, StarID FROM {0} s WHERE s.Filename IS \"{1}\";'''.format(DB.SpectraTableName, filename))
            clusterID = DBcursor.fetchall()[0]
            
            if not hasAtmParms(filename, parmType=k.SpectroscopicParm):
                DBcursor.execute('''SELECT Teff, Loggf, Vturb FROM {0} s WHERE s.Cluster IS \"{1}\" AND s.StarID IS \"{2}\";'''.format(DB.StarTableName, clusterID[0], clusterID[1]))
                photParms = DBcursor.fetchall()[0]
                if len(photParms) > 0 and photParms[0] != 0.:
                    newRecord = [filename, photParms[0], photParms[2], photParms[2], k.PhotometricParm]
                    DBcursor.execute('''INSERT INTO \"{0}\" {1} VALUES{2}'''.format(DB.atmospheric_parms, DB.AtmParFieldNames, DB.AtmParGenFormat), newRecord)
                else:
                    # We'll have to insert our own calculation from B-V values
                    pass

            if not hasAtmParms(filename, parmType=k.SpectroscopicParm):
                # Get the Fe and Ti lines for MOOG
                DBcursor.execute('''SELECT * FROM {0} s WHERE s.ion IN ({1}) AND Filename IS\"{2}\";'''.format(DB.MeasuredLinesTableName, '22.0, 22.1, 26.0, 26.1', filename))
                lines = DBcursor.fetchall()
            
            if (clusterID[0],clusterID[1]) not in spectraNames:
                spectraNames.append((clusterID[0],clusterID[1]))
                
            if clusterID[1] in spectraLines.keys():
                tempList = spectraLines[clusterID[1]]
                tempList.extend(lines)
                spectraLines[clusterID[1]] = tempList
            else:
                spectraLines[clusterID[1]] = lines
            

        for (cluster, star) in spectraNames:
                starLines = sorted(filterLines(spectraLines[star], DBcursor), key=lambda x: x[0])
                starLines = sorted(starLines, key=lambda x: x[1])
                
                if len(starLines) < 20:
                    print('Insufficient lines from: {0} ({1}).'.format(filename, len(starLines)))
                    continue
                
                logFile = makeMOOGLogFile(starLines, k.MOOGTempLogName)
                
                if star in giantStars:
                    (Teff, Logg, VmTurb, Metal) = spd.determineParms(logFile, (k.AtmGiantTGuess, k.AtmGiantGGuess, k.AtmGiantvGuess), gKuruczAtms, gPradks, k.AtmGiantMGuess)
                else:
                    (Teff, Logg, VmTurb, Metal) = spd.determineParms(logFile, (k.AtmDwarfTGuess, k.AtmDwarfGGuess, k.AtmDwarfvGuess), dKuruczAtms, dPradks, k.AtmDwarfMGuess)
                with open("Trial_Run_09_17_stepsize.log", "a+") as thisLog:
                    thisLog.write('{0} - T={1:4.1f} : Logg={2:1.2f} : v={3:1.2f} : m={4:+1.3f}\n'.format(star, Teff, Logg, VmTurb, Metal))
                
                
# Command line calls for testing only
if __name__ == '__main__':
    populateAtmTable()
    exit()
