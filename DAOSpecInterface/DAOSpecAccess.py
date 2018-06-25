#! /usr/bin/env python
# 
# Library: DAOSpecAccess.py
#
# Creator: Mike Lum
#
# Functions:
#       makeDAOSpecScript:
#
#       callDAOSpec:
#
#       readDAOAbs
#
#       readDAORV
#       
#       getDAOSpecEQWs
#
# Description: Functions to access the DAOSpec software
#
# Revision History:
#    Date          Vers.    Author        Description
#    2018-04-05    1.0f2    Lum         Format and comment cleanup for github
#    2017-02-05    1.0f1    Lum         First checked in
#
#
#

import datetime
import os
import shutil
import subprocess as sub

from DAOSpecInterface import DAOSpecConsts as kDAO
from SpectraProcessing import SpectraData as SD
from utilities import utilities as u

def makeDAOSpecScript(spectraFile, logfile = 'logfile', addlOpts = '', WLRange=None):
    if WLRange is not None:
        wlOpts = 'SH={0:4.2f}\nLO={1:4.2f}\n'.format(WLRange[0], WLRange[1])
    else:
        wlOpts = ''
    allOpts = kDAO.DAOSpecOpts+addlOpts+wlOpts
    return kDAO.DAOSpecScriptHead + kDAO.DAOSpecCommand + ' << DONE >{0}{1}\n{2}\n\nDONE\n'.format(logfile, allOpts, spectraFile)
    

def callDAOSpec(spectraFile, logfile = 'logfile', WLRange=None):
# Note: relies on laboratory.dat line list in the local directory
    DAOScript = makeDAOSpecScript(spectraFile, logfile=logfile, WLRange=WLRange)
    sub.call(DAOScript, shell=True)
    
def readDAOAbs(daoLogfile):

    infile = open(daoLogfile, 'r')
    data = infile.readlines()
    infile.close()
    lines = [line.split() for line in data]
    meas = [line for line in lines[2:] if len(line)>6]
    lineData = sorted([[float(line[6]), float(line[5]), daoLogfile[:-8]+'.fits', float(line[1]), float(line[2])] for line in meas], key=lambda l: l[1])
    return lineData


def readDAORV(daoLogfile):
# Reads the calculated RV from the first line of a log file.
# We're expecting the first line to read:
# "Radial velocity =    <RV value>  dispersion =..."
    rvStr='NaN'
    infile = open(daoLogfile, 'r')
    data = infile.readlines()
    infile.close()
    for line in data:
        words = line.split()
        if len(words)>0:
            if words[0] == "Final":
                rvStr = words[3]
                break
    if u.is_number(rvStr):
        rv = float(rvStr)
    else:
        rv = -999.99
    return rv
    
    
def getDAOSpecEQWs(spectraFile, outlog=None, frame=0.0):
    if outlog is not None:
        outfile = open(outlog,'a')
    objectName, numAps, apSize, apBounds = SD.ReadSpectraInfo(spectraFile)
    # DAOSpec can only handle 1d spectra, so if there are more than 1 apertures, bail
    if numAps > 1:
        print('2d spectra detected, DAOSpec failure.')
        return []
    usedTempfile = False
    # DAOSpec is Fortran, which means that it blows up if passed a filename,
    # including path, of longer than ~16 characters. If we have a long filename
    # we copy the spectra to the local directory, and nuke it when we're done.
    if len(spectraFile) > 16:
        usedTempfile = True
        shortFilename = 'tempSpec'+datetime.datetime.now().strftime('%f')+'.fits'
        shutil.copy(spectraFile, shortFilename)
    else:
        shortFilename = spectraFile
    if outlog is not None:
        logfile = outlog
    else:
        logfile = 'daolog'+datetime.datetime.now().strftime('%f')
        
    callDAOSpec(shortFilename, logfile = logfile, WLRange=(apBounds[0][0]+frame, apBounds[0][1]-frame))
    if usedTempfile:
        os.remove(shortFilename)
    if outlog is None:
        os.remove(logfile)

    lineData = readDAOAbs(shortFilename[:-4]+'daospec')
    os.remove(shortFilename[:-4]+'daospec')
    
    if len(lineData) == 0:
        if outlog == None:
            print('Abundance calculation failed for:{0}'.format(spectraFile))
        else:
            outfile.write('{0}\n'.format(spectraFile))
            outfile.close()
    return lineData
