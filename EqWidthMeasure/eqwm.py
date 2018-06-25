#! /usr/bin/env python
# 
# Program: eqwm.py
#
# Author: Mike Lum
#
# Usage: ./eqwm.py <fits spectra file> <line list file> [-ahov]
#
# Description: Program measures the equivalent widths of the lines from the passed line list
#                in the passed .fits spectra. 
#                Interactive line fitting is done using the IRAF function: 'splot'.
#                Initial .fits spectra processing requires IRAF functions:
#                 'imgets' and 'listpixels'
#
# Revision History:
#    Date        Vers.    Author        Description
#    2015-07-14    1.2b1    Lum     Added return codes for bad line measurements
#    2014-09-08    1.1f1    Lum     Re-organization of the file into a project format.
#                                   Common functions moved to other sub-modules. Globals
#                                   removed, and constants moved to a constant module.
#    2014-05-14    1.0f1    Lum        Basic functionality of Hai Fu's "ewmeasure" function
#                                implemented in a stand-alone Python program.
#    2014-05-12    0.0a0    Lum        First checked in
#
# To Do:
#    *
#
import os
import sys

from constants import constants as k
from constants import DBFormats as DB
from Databases import LineLookup as LL
#from EqWidthMeasure import Baseline as BL
from EqWidthMeasure import VoigtFit as VF
from MOOGInterface import MOOGInterface as MI
from utilities import utilities as u
from SpectraProcessing import SpectraData as SD

#from EqWidthMeasure import eqwmHelp as help

#from collections import deque
from pyraf import iraf
import numpy as np
#import pyfits as pf
#import astropy.io.fits as pf
#import pylab as pl
#from matplotlib import pyplot
#from matplotlib import figure
#import pyspeckit
#import scipy.signal as sig
#import scipy.interpolate as interp
import sqlite3

# Functions
def eqwm(interactiveFitting = True,
        includeOutputHeader = False, 
        pauseBetweenLines = False, 
        recordSkippedLines = True, 
        outputFilename = '', 
        spectraFilename = '', 
        lineFilename = '',
        lineEvalFilename = '',
        continuum1 = False):

    spectraFilename, outputFilename, lineEvalFilename, lineFilename = \
    DoFileSetup(spectraFilename, outputFilename, recordSkippedLines, lineEvalFilename, lineFilename)

    if interactiveFitting:
        EQWidths, badLines, dopplerCount, dopplerCheck = MeasureWidthsInter(spectraFilename, lineFilename)
    else:
        EQWidths, badLines, dopplerCount, dopplerCheck = MeasureWidthsAuto(spectraFilename, lineFilename, continuum1, pauseBetweenLines)
    
    EnterNewData(EQWidths, badLines, spectraFilename)
    
    # The following are legacy from when we used text files for data storage:
#    OutputLineMeasures(EQWidths, outputFilename, objectName, includeOutputHeader)
#    OutputLinesEvaluation(EQWidths, badLines, lineEvalFilename)
#    OutputDopplerVerification(dopplerCount, dopplerCheck)
    
    DoCleanup()

    return


def DoFileSetup(spectraFilename, outputFilename, recordSkippedLines, lineEvalFilename, lineFilename):
#   Function: DoFileSetup(spectraFilename, outputFilename, recordSkippedLines, lineEvalFilename, lineFilename)
#
#   Verifies the existence of the input spectra and line list files. Verifies that
#   the measured line output, skipped line log, and splot log files can be written.
    
    # Verify the spectra file - We kinda have to have this one. ;)
    try:
        # Need to have the full path:
        spectraFilename = os.path.abspath(spectraFilename)
        tempFile = open(spectraFilename,'r')
        tempFile.close()
    except IOError:
        print('Cannot access spectra file: {0}. Exiting.'.format(spectraFilename))
        exit()
    
    # Verify the output file (obsoleted by DB, but capability remains)
    if outputFilename != '':
#        outputFilename = spectraFilename+'.log' # Formerly used as a default outfile name
        if os.path.isfile(outputFilename):
            response = raw_input('File: {0} already exists. Overwrite (Yn)?'.format(outputFilename))
            if response.upper() == 'N':
                exit()
        try:
            tempFile = open(outputFilename,'w')
            tempFile.close()
        except IOError:
            print('Cannot write to output log file: {0}. Exiting.'.format(outputFilename))
            exit()
        
#    Check the line evaluation file (obsoleted by DB, but capability remains)
    if lineEvalFilename != '':
#        lineEvalFilename = spectraFilename + k.LineEvalExtension   # Old default name
        try:
            tempFile = open(lineEvalFilename,'w')
            tempFile.close()
        except IOError:
            print('Cannot write to line evaluation log: {0}. Exiting.'.format(lineEvalFilename))
            exit()

    # Clean up prior runs and logs:
    if os.path.isfile(k.SplotLogFilename):
        try:
            os.remove(k.SplotLogFilename)
        except IOError:
            print('Unable to delete old log file: {0}. Exiting.'.format(k.SplotLogFilename))
            exit()
    if lineFilename != '':
    # Not passing a line file is o.k. - we'll just look up lines in the DB
        # Verify the line list file
        try:
            # Need to have the full path:
            lineFilename = os.path.abspath(lineFilename)
            lineFile = open(lineFilename,'r')
            lineFile.close()
        except IOError:
            print('Cannot access line list: {0}. Exiting.'.format(lineFilename))
            exit()

    return spectraFilename, outputFilename, lineEvalFilename, lineFilename
    
    
def MeasureWidthsInter(spectraFilename, lineFilename):
#   Function:     MeasureWidthsInter(spectraFilename, lineData, apBoundaries)
#
#   Inputs: spectraFilename : The absolute path of the .fits file containing the spectrum to measure
#           lineFilename: A file containing information on lines to measure

#   Outputs: EQWidths - A list of the measured widths, each tuple of the form:
#               (wl, ion, ep, gf, c6, eqw, meas. wl, fwhm (gaussian)). 

#
#   Desription: This function calls the iraf 'splot' task, and asks the user to manually
#           measure each line in the passed list. The splot log file is read and the actual
#           line equivalent width is assumed to be the last measurement in the log file.
    objectName, numApertures, apertureSize, apBoundaries = SD.ReadSpectraInfo(spectraFilename)
    lineData, unused = getLinesForApertures(lineFilename, apBoundaries)
    
    EQWidths =[]
    badLines = []
    dopplerCount = 0
    dopplerCheck = 0.0
    for [wl, ion, ep, gf, c6, solarEqw, comment] in lineData:
        contentAps = [apBoundaries.index(x) for x in apBoundaries if x[0]<=wl<=x[1]]
        if len(contentAps) == 0:
            print('Wavelength ({0:8.3f}) for {1:4.1f} line not covered by spectra. Continuing.'.format(wl, ion))
            continue
        wlUpper = wl + k.HalfWindowWidth
        wlLower = wl - k.HalfWindowWidth
        print('Measuring {0:4.1f} line: {1:8.3f} (displaying first of {2:1d} orders)'.format(ion, wl,len(contentAps)))
        iraf.onedspec()
        iraf.splot(spectraFilename, line=contentAps[0]+1, band=1, xmin=wlLower+1., xmax=wlUpper+1., save_file=k.SplotLogFilename, options='nosysid, auto')
#        iraf.splot(spectraFilename, line=contentAps[0]+1, band=1, xmin=wlLower+1., xmax=wlUpper+1., ymin=k.SpecFluxMin, ymax=k.SpecFluxMax, save_file=k.SplotLogFilename, options='nosysid, auto')
        try:
            with open(k.SplotLogFilename) as filePointer:
                # Assume that the desired fit is the _last_ fit made.
                    # It is important for the blended fit ('d' option) that the
                    # actual line is the last one selected in the blend
                dataLine = filePointer.readlines()[-1].split()
                # splot data line is: line center, cursor loc. (cont. level), flux(?), 
                #                       eqw, core depth from cursor, gfwhm, lfwhm
                lineCenter = float(dataLine[0])
                print('Measured line: {0:8.3f} at: {1:5.1f} mAngstroms (FWHM = {2:5.4f})'.format(wl, float(dataLine[3])*1000.0, float(dataLine[5])))
                EQWidths.append([wl, ion, ep, gf, c6, float(dataLine[3])*1000.0, lineCenter, float(dataLine[5])])
                thisDoppler = ((lineCenter - wl)/wl) * k.SpeedofLight
                dopplerCheck += thisDoppler
                dopplerCount += 1
            os.remove(k.SplotLogFilename)
        except (IOError, ValueError):
            # IOError - No line is fit (no log was created)
            # ValueError - Log output from splot is in unrecognizeable format
            badLines.append([wl,ion, k.unknownReason, k.unknownStr])
            print('No line measured for wl: {0:8.3f}.'.format(wl))
         
    return EQWidths, badLines, dopplerCount, dopplerCheck
    
    
def isGoodFit(guesses, results):
#   Function:     isGoodFit(guesses, results)
#
#   Inputs: guesses, results: arrays of n line guesses and fits to be evaluated for similarity.
#               Each list is n*3 items long, with each "guess" consisting of:
#               line depth (above continuum is pos.), line center, line FWHM
#
#   Outputs: Boolean whether the guesses match the results.
#
#   Desription: For now, we just evaluate the results based on the expected centers and
#               FWHM.

    if len(results) != len(guesses):
        return False
        
    measCenters = sorted([results[i] for i in range(1,len(results),3)])
    guessCenters = sorted([guesses[i] for i in range(1,len(guesses),3)])
    deltaCenters = [abs(i-j) for i,j in zip(measCenters, guessCenters)]
    
    if max(deltaCenters) > k.WLRegionTolerance:
#        for i,j in zip(measCenters, guessCenters):
#            print i,j
        return False
    
    measFWHM = [results[i] for i in range(2,len(results),3)]
    guessFWHM = [guesses[i] for i in range(2,len(guesses),3)]
    
    if max(measFWHM) > k.BroadLineFWHMLimit:
        return False
    
    return True


def MeasureWidthsAuto(spectraFilename, lineFile=None, continuum1=False, pauseBetweenLines=False):
#   Function:     MeasureWidthsAuto(spectraFilename, lineData, apBoundaries)
#
#   Inputs: spectraFilename : The absolute path of the .fits file containing the spectrum to measure. We will measure all lines for this spectrum, and leave it up to the caller to parse it.
#
#           continuum1   : Should we assume a continuum of 1.00?
#           pauseBetweenLines : Mainly for debugging purposes.
#
#   Outputs: EQWidths - A list of the measured widths, each tuple of the form:
#               (wl, ion, ep, gf, c6, eqw, meas. wl., fwhm (gaussian)).
#            badLines : A list of lines we were unable to measure. Each element is a
#                       tuple of the form: (wl, ion)
#            dopplerCount : A count of the 'good' lines measured with a Doppler shift
#            dopplerCheck : A running total of the Doppler shifts of the good lines
#
#   Desription: This function automatically measures the equivalent width of each line
#               in the passed list
    
    # Un-comment below to use the home-grown Voigt line fitting function
    if lineFile is not None:
        unused1, unused2, unused3, apBoundaries = SD.ReadSpectraInfo(spectraFilename)
        lineData, unused4 = getLinesForApertures(lineFile, apBoundaries)
        
    EQWidths, badLines, dopplerCount, dopplerCheck = VF.measureVoigtEQWs(spectraFilename, lineData=lineData)
    
    # Un-comment below to use the PySpecKit fitting functionality 
    # - Currently (6-2017) unsupported!
#    EQWidths, badLines, dopplerCount, dopplerCheck = doPyspecFit(spectraFilename, lineFilename, continuum1, pauseBetweenLines)
    return EQWidths, badLines, dopplerCount, dopplerCheck

def EnterNewData(EQWidths, badLines, sourceFilename):
    DBconnection = sqlite3.connect(DB.ClusterDBName)
    DBcursor = DBconnection.cursor()
    
    shortFilename = sourceFilename.split('/')[-1]
    # Note: We assume that this is the newest measurement of the line, and should
    # supercede any previous ones for this spectrum.
    
    # For 1-d spectra, where we have split out individual orders, we need to 
    # change the filename from CCCC_XXX_NNNN_II_SS_OO.fits to omit the "Order"
    # ("OO") numbers:
    if len(shortFilename.split('_'))>5:
        # We need the first four "sections", and to add a zero to the fifth:
        newName =  '_'.join(shortFilename.split('_')[:4]) + '_0'+shortFilename.split('_')[-2]+'.fits'
        shortFilename = newName
#        sourceFilename = '/'.join(sourceFilename.split('/')[:-1])+'/'+newName
    
    SD.EnterStarName(sourceFilename)
    # Note: the star will only be entered if it's not already in our DB
    
    for line in EQWidths:
        DBcursor.execute('''SELECT * FROM {0} WHERE Filename=\"{1}\" AND Ion={2:2.1f} AND Wavelength={3:4.3f}'''.format(DB.MeasuredLinesTableName, str(shortFilename), line[1], line[0]))
        linesMatched = DBcursor.fetchall()
        if len(linesMatched) > 0:
        # Overwrite the old measured value
            DBcursor.execute('''UPDATE {0} SET EQM={1:3.1f}, MeasWL={2:4.3f}, FWHM={3:3.3f} WHERE Filename=\"{4}\" AND Ion={5:2.1f} AND wavelength={6:4.3f}'''.format(DB.MeasuredLinesTableName, line[5], line[6], line[7], str(shortFilename), line[1], line[0]))
        else:
            newRecord = [line[1], line[0], shortFilename, line[6], line[5], line[7]]
            DBcursor.execute('''INSERT INTO \"{0}\" VALUES {2}'''.format(DB.MeasuredLinesTableName, DB.MeasuredLinesFieldNames, DB.MeasuredLinesGenFormat), newRecord)

    for line in badLines:
        DBcursor.execute('''SELECT * FROM {0} WHERE Filename=\"{1}\" AND Ion={2:2.1f} AND Wavelength={3:4.3f}'''.format(DB.MeasuredLinesTableName, str(shortFilename), line[1], line[0]))
        linesMatched = DBcursor.fetchall()
        if len(linesMatched) > 0:
            DBcursor.execute('''UPDATE {0} SET EQM={1:3.1f}, MeasWL={2:4.3f}, FWHM={3:3.3f} WHERE Filename=\"{4}\" AND Ion={5:2.1f} AND wavelength={6:4.3f}'''.format(DB.MeasuredLinesTableName, k.BadWLValue, k.BadWLValue, k.BadWLValue, str(shortFilename), line[1], line[0]))
        else:
            newRecord = [line[1], line[0], shortFilename, k.BadWLValue, k.BadWLValue, k.BadWLValue]
            DBcursor.execute('''INSERT INTO \"{0}\" VALUES{2}'''.format(DB.MeasuredLinesTableName, DB.MeasuredLinesFieldNames, DB.MeasuredLinesGenFormat), newRecord)
    
    DBconnection.commit()
    DBcursor.close()
    return


def getLinesForApertures(lineFilename, apBoundaries=None):
#   Function:     getLinesForApertures(lineFilename, apBoundaries):
#
#   Inputs: lineFilename: A file containing information on lines to measure - can be "None"
#           apBoundaries: A list of a spectra's aperture boundaries - used to determine
#                         the spectral range for line lookup into the DB
#
#   Outputs: EQWidths - A list of the measured widths, each tuple of the form:
#               (wl, ion, ep, gf, c6, eqw, meas. wl., fwhm (gaussian)).
#            badLines : A list of lines we were unable to measure. Each element is a
#                       tuple of the form: (wl, ion)
#            dopplerCount : A count of the 'good' lines measured with a Doppler shift
#            dopplerCheck : A running total of the Doppler shifts of the good lines
#
#   Desription: This function automatically measures the equivalent width of each line
#               in the passed list.
    if lineFilename != None:
        lineData = MI.ReadMOOGLineFile(lineFilename)
        refs = []  # Note: It is POSSIBLE to do line lookups into our DB to find
                   # references for the lines in the passed file, but we have no
                   # way of insuring that the data in the file is actually from our
                   # DB's references....Thus the empty reference list when passed 
                   # a line file.
    elif apBoundaries != None:
        minWL = min([x[0] for x in apBoundaries])
        maxWL = max([x[1] for x in apBoundaries])
        lineData, refs = LL.getLines(elements=[], wavelengthRange=[(minWL, maxWL)], dataFormat=None)
    else:
        lineData = []
        refs = []
    
    return lineData, refs
    
    
def OutputLineMeasures(EQWidths, outputFilename, objectName, includeOutputHeader):
    
    if len(EQWidths) != 0:
        try:
            outfile = open(outputFilename,'w')
        except IOError:
            # Note: We need an error outputting function
            print('Unable to open output file: {0}. Exiting.'.format(outputFilename))
            exit()
        outfile.write('# '+ objectName+'\n')
        if includeOutputHeader: outfile.write(k.OutputFileHeader)
        for line in EQWidths:
            outfile.write('{0:8.3f}{1:11.1f}{2:11.5f}{3:11.6f}{4:11.1f}{5:17.1f}\n'.
                format(line[0],line[1],line[2],line[3],line[4],line[5]))
        outfile.close()
    return
    
def OutputLinesEvaluation(goodLines, badLines, lineFilename):

    try:
        outfile = open(lineFilename,'w')
    except IOError:
        # Note: We need an error outputting function
        print('Unable to open output file: {0}. Exiting.'.format(lineFilename))
        exit()
    outfile.write('# '+ lineFilename[:-16] + ' - Line evaluation\nIon | Expected WL | Measured WL | FWHM | Reject Reason | Notes:\n--------|------\n')
    
    outTextLines = ['{0:2.1f} {1:4.3f}    nan      nan    {2:3d}   {3}\n'.format(line[1],line[0], line[2], line[3]) for line in badLines]

    outTextLines.extend(['{0:2.1f} {1:4.3f} {2:4.3f} {3:4.3f}\n'.format(line[1],line[0], line[6], line[7]) for line in goodLines])
    
    # An alphabetical sort actually works as a numerical sort, since we have
    # our lines in ion, expected wl...etc. form.
    outTextLines.sort()
    
    for line in outTextLines:
        outfile.write(line)
    outfile.close()
    return

def OutputDopplerVerification(dopplerCount, dopplerCheck):
    # Verify our Doppler checking:
    # (**** Should this be implemented as an option? Command line activated?***)
    if dopplerCount > 0:
        avgDoppler = dopplerCheck/dopplerCount
        # Real error out
        print(' Measured average line Doppler of: {0:6.2f} km/s'.format(avgDoppler))
        if abs(avgDoppler) > k.MinDopplerCorrection:
            print('*** Warning: Average measured line shift (={0:6.2f} km/s) may require a Doppler adjustment to the input spectrum.'.format(avgDoppler))
    return

def DoCleanup():
    return
            
# This should be called as a library function, not command line:
# Besides, the command line stuff is probably broken by now.

'''

if __name__ == '__main__':

# Parse the command line arguments
    temp = os.sys.argv
    argv = deque(temp)
    
    outputFilename = ''
    spectraFilename = ''
    lineFilename = ''
    lineEvalFilename = ''
    
    interactiveFitting = True
    includeOutputHeader = False
    pauseBetweenLines = False
    recordSkippedLines = False
    continuum1 = False

    while len(argv) > 0:
        flag = argv.popleft()
        if flag == '-h':
            help.printHelpText()
            exit()
        elif flag == '-1':
            continuum1 = True
            if g.VerboseMode: print('Continuum will be fixed to 1.0.')
        elif flag == '-a':
            interactiveFitting = False
            if g.VerboseMode: print('Automated fitting enabled.')
        elif flag == '-c':
            includeOutputHeader = True
            if g.VerboseMode: print('Output file will include column headers.')
        elif flag == '-d':
            pauseBetweenLines = True
            if g.VerboseMode: print('Automated measurement will pause between lines.')
        elif flag == '-s':
            recordSkippedLines = True
            if g.VerboseMode: print('Recording skipped line measurements.')
        elif flag == '-v':
            g.VerboseMode = True
            showIRAFMessages = 1
            print('Verbose Mode enabled.')
        elif flag == '-o':
            outputFilename = argv.popleft()
        elif str.find(flag,k.EQWMModuleName) != -1 :
            if g.VerboseMode: print('Executing command.')
        elif str.find(flag,k.InterpreterName) != -1 :
            if g.VerboseMode: print('...')
        else:
            if spectraFilename == '':
            # The first non-flag is the spectra file name
                spectraFilename = flag
                if g.VerboseMode: print('Input spectra file: {0}.'.format(spectraFilename))
            elif lineFilename == '':
            # The second is the line list name
                lineFilename = flag
                if g.VerboseMode: print('Input line list: {0}.'.format(lineFilename))
            else:
                print('Unrecognized option: {0}. Exiting.'.format(flag))
                exit()

    eqwm(interactiveFitting, includeOutputHeader, pauseBetweenLines, 
        recordSkippedLines, outputFilename, spectraFilename, lineFilename, skippedLineFilename, continuum1)
    exit()
'''    
    
### Been having problems with pyspeckit...so stopped using it:
'''
def doPyspecFit(spectraFilename, lineData, continuum1=False, pauseBetweenLines=False):
#   Function:     doPyspecFit(spectraFilename, lineFilename, continuum1, pauseBetweenLines)
#
#   Inputs: spectraFilename : The absolute path of the .fits file containing the spectrum to measure
#           lineFilename: A file containing information on lines to measure
#           continuum1   : Should we assume a continuum of 1.00?
#           pauseBetweenLines : Mainly for debugging purposes.
#
#   Outputs: EQWidths - A list of the measured widths, each tuple of the form:
#               (wl, ion, ep, gf, c6, eqw, meas. wl., fwhm (gaussian)).
#            badLines : A list of lines we were unable to measure. Each element is a
#                       tuple of the form: (wl, ion)
#            dopplerCount : A count of the 'good' lines measured with a Doppler shift
#            dopplerCheck : A running total of the Doppler shifts of the good lines
#
#   Desription: This function automatically measures the equivalent width of each line
#               in the passed list, using the PySpecKit fitting routines
    objectName, numApertures, apertureSize, apBoundaries = SD.ReadSpectraInfo(spectraFilename)

    specList = SD.LoadAndSmoothSpectra(specFilename)
    
    lineData, refs = getLinesForApertures(lineFilename, apBoundaries)
    
    # Must have matplotlib in interactive mode
    EQWidths =[]
    badLines = []
    dopplerCount = 0
    dopplerCheck = 0.0
    
    pl.ion()
    plotFigure = pyplot.figure(num='Line Analysis')
    twopixels = 2*abs(apBoundaries[0][0]-apBoundaries[0][-1])/apertureSize
    
    for [wl, ion, ep, gf, c6, solarEqw, comment] in lineData:    
        contentAps = [apBoundaries.index(x) for x in apBoundaries if x[0]<=wl<=x[1]]

        if len(contentAps) == 0:
            if g.VerboseMode: print('Wavelength ({0:8.3f}) for {1:4.1f} line not covered by spectra. Continuing.'.format(wl,ion))
            continue
        wlUpper = wl + k.HalfWindowWidth
        wlLower = wl - k.HalfWindowWidth
        # Only looking at the first aperture for now. There may be a time where we
        # will want to consider that a line appears in two different apertures, though.
        newApIdx = contentAps[0]
        smoothSpec = specList[newApIdx]
        if not SD.IsGoodAperture(smoothSpec):
            badLines.append([wl, ion, k.badSpectralRegion, 'Bad aperture # %d'%(newApIdx)])
            continue
        
        # Plot for debugging.
        smoothSpec.plotter(xmin=wlLower, xmax=wlUpper, ymin=k.SpecFluxMin, ymax=k.SpecFluxMax, figure=plotFigure)
        
        loApPix = smoothSpec.xarr.x_to_pix(wlLower)
        hiApPix = smoothSpec.xarr.x_to_pix(wlUpper)

        maxes = [smoothSpec.xarr[x+loApPix].value for x in sig.argrelmax(smoothSpec.data[loApPix:hiApPix])[0]]
        maxes.extend([smoothSpec.xarr[loApPix].value, smoothSpec.xarr[hiApPix].value])
        maxes.sort()

        mins = [smoothSpec.xarr[x+loApPix].value for x in sig.argrelmin(smoothSpec.data[loApPix:hiApPix])[0]]
        mins.extend([smoothSpec.xarr[loApPix].value, smoothSpec.xarr[hiApPix].value])
        mins.sort()
        
        if continuum1: 
        # Forced baseline of 1. We 'cheat' by creating a spline of order 1, and coefficient of 1
            xArr = [smoothSpec.xarr[x].value for x in range(loApPix, hiApPix)]
            yArr = np.ones(len(xArr))
            spline = interp.UnivariateSpline(xArr, yArr, k=1)
        else:
        # or... a 'smarter' approach at a better local fit.
        # pyspeckit's baseline function only excludes lines which have been fit already.
        # Unfortunately, we don't know the line fit yet (it's our goal). 
        # Create a spline fit of the local maximums as a baseline.
            #
            # This is not unprecidented, Ramirez & Allende Prieto (arXiv:1109.4425v1, 2011)
            # used it when measuring abundances in Arcturus, for example.
            spline = BL.FitSplineBaseline(loApPix, hiApPix, smoothSpec)

        # "highs" are points we think that pyspec kit will recognize as the maximum
        # extent of line wings. For this to happen, the data needs to be (nearly) at,
        # or above the baseline (spine) fit for that point.
        # Note that there should always be at least one point which is above the fit
        # spline...
        highs = [x for x in maxes if smoothSpec.data[smoothSpec.xarr.x_to_pix(x)] > 0.99*spline(x)]

        # We can receive (0,0) or (-1,-1) results from u.bracket
        if len(highs) < 1:
            # It is possible that we were unable to fit a spline, so skip the line.
            # This will occur more frequently with a constant baseline (=1.0), in 
            # areas where the spectra lies entirely below 1.0.
            badLines.append([wl, ion, k.unknownReason, k.unknownStr])
            continue
        else:
            loIdx, hiIdx = u.bracket(wl, highs)
        
        # loIdx and hiIdx are the indices of the highs which can be easily fit
        # by pyspeckit. They will bracket the desired line, plus zero or more 
        # additional lines. In order to get a good fit, at least one additional
        # line must be fit on either side of the desired line.
        
        minIdx, maxIdx = u.bracket(wl, maxes) # The two local maxima which bracket our line
        if maxIdx-minIdx > 8: 
        # PySpeckit can't reasonably fit more than 8 lines in a multi-line fit
            # So, we take the "middle 8" lines from the window, 
            # and hope we get the right ones
            delta = (maxIdx-minIdx-7)/2 
            # This will actually round up, which is o.k., since we go one wider
            #  on each end, below.
            maxIdx -= delta
            minIdx += delta
            
        if minIdx > 0: minIdx -= 1
        if maxIdx < len(maxes)-1 and maxIdx >-1: maxIdx += 1
        # minIdx and maxIdx now include one line on either side of the line, to the extent
        # of our window
        
        # This will write over the lo and hi indices with the WL range to fit
        if highs[loIdx] > maxes[minIdx]:
            loPix = smoothSpec.xarr.x_to_pix(maxes[minIdx])
        else:
            loPix = smoothSpec.xarr.x_to_pix(highs[loIdx])

        if highs[hiIdx] < maxes[maxIdx]:
            hiPix = smoothSpec.xarr.x_to_pix(maxes[maxIdx])
        else:
            hiPix = smoothSpec.xarr.x_to_pix(highs[hiIdx])
                      
        loWL = smoothSpec.xarr[loPix].value
        hiWL = smoothSpec.xarr[hiPix].value
        
        # Plot the locations of the mins and maxes on the pyspeckit plot
        for pixel in [x for x in mins if loWL<= x <= hiWL] :
            pl.axvline(x=pixel, ymin=0., ymax=1., linewidth=2, color='g', linestyle=':')
        for pixel in [x for x in maxes if loWL<= x <= hiWL]:
            pl.axvline(x=pixel, ymin=0., ymax=1., linewidth=2, color='r', linestyle=':')
            
        lineGuesses = []
        strongGuesses = 0
        for thisCenter in [x for x in mins if loWL<x<hiWL]:
        # Our line 'guesses' are centered where the spectra hits a local minimum, with the 
        # depth of the minimum being the amplitude at that point and the FWHM being half
        # the distance between the two local maximums which bracket the local minimum.
            guessLeftIdx, guessRightIdx = u.bracket(thisCenter, maxes)
            guessLeft = smoothSpec.xarr[smoothSpec.xarr.x_to_pix(maxes[guessLeftIdx])].value
            guessRight = smoothSpec.xarr[smoothSpec.xarr.x_to_pix(maxes[guessRightIdx])].value

            guessDepth = smoothSpec.data[smoothSpec.xarr.x_to_pix(thisCenter)] - spline([thisCenter])[0]
            # If we have a large emission line as a guess, skip it
            if guessDepth > 0.2: continue
            # A small emission or absorption line isn't so bad, but let's really be sure that it is
            # an emission line, by forcing the fitter to 'guess' at it's size.
            if guessDepth > -0.03: guessDepth = -0.03
            else: strongGuesses += 1
            lineGuesses.extend([guessDepth, thisCenter, (guessRight-guessLeft)/2.])
        
        # We need at least one prominent guess...
        if strongGuesses < 1:
            badLines.append([wl, ion, k.badSpectralRegion, k.noLines])
            continue

        smoothSpec.baseline(subtract=False, order=0, prefit_continuum = spline, reset_selection=False, exclude=[smoothSpec.xarr[loApPix].value, smoothSpec.xarr[loPix].value, smoothSpec.xarr[hiPix].value, smoothSpec.xarr[hiApPix].value], save=False, highlight_fitregion=False)

        smoothSpec.specfit(fittype='gaussian',guesses=lineGuesses, quiet=True, save=False, show_components=False, use_window_limits=True, vheight=True, clear=True, reset_selection = True, exclude=[smoothSpec.xarr[loApPix].value, smoothSpec.xarr[loPix].value, smoothSpec.xarr[hiPix].value, smoothSpec.xarr[hiApPix].value])


        lineResults = smoothSpec.specfit.modelpars[:]
        
        if not isGoodFit(lineGuesses, lineResults):
            # Try widening our guess window:
            # This will write over the lo and hi indices with the WL range to fit
            if highs[loIdx] > maxes[minIdx]:
                if minIdx > 0: minIdx -=1
                loPix = smoothSpec.xarr.x_to_pix(maxes[minIdx])
            else:
                if loIdx > 0: loIdx -=1
                loPix = smoothSpec.xarr.x_to_pix(highs[loIdx])

            if highs[hiIdx] < maxes[maxIdx]:
                if maxIdx < len(maxes)-1 and maxIdx >-1: maxIdx += 1
                hiPix = smoothSpec.xarr.x_to_pix(maxes[maxIdx])
            else:
                if hiIdx < len(highs)-1 and hiIdx >-1: hiIdx += 1
                hiPix = smoothSpec.xarr.x_to_pix(highs[hiIdx])
                          
            loWL = smoothSpec.xarr[loPix].value
            hiWL = smoothSpec.xarr[hiPix].value
            
            # Plot the locations of the mins and maxes on the pyspeckit plot
            for pixel in [x for x in mins if loWL<= x <= hiWL] :
                pl.axvline(x=pixel, ymin=0., ymax=1., linewidth=2, color='g', linestyle=':')
            for pixel in [x for x in maxes if loWL<= x <= hiWL]:
                pl.axvline(x=pixel, ymin=0., ymax=1., linewidth=2, color='r', linestyle=':')
                
            lineGuesses = []
            strongGuesses = 0
            for thisCenter in [x for x in mins if loWL<x<hiWL]:
            # Our line 'guesses' are centered where the spectra hits a local minimum, with the 
            # depth of the minimum being the amplitude at that point and the FWHM being half
            # the distance between the two local maximums which bracket the local minimum.
                guessLeftIdx, guessRightIdx = u.bracket(thisCenter, maxes)
                guessLeft = smoothSpec.xarr[smoothSpec.xarr.x_to_pix(maxes[guessLeftIdx])].value
                guessRight = smoothSpec.xarr[smoothSpec.xarr.x_to_pix(maxes[guessRightIdx])].value

                guessDepth = smoothSpec.data[smoothSpec.xarr.x_to_pix(thisCenter)] - spline([thisCenter])[0]
                # If we have a large emission line as a guess, skip it
                if guessDepth > 0.2: continue
                # A small emission or absorption line isn't so bad, but let's really be sure that it is
                # an emission line, by forcing the fitter to 'guess' at it's size.
                if guessDepth > -0.03: guessDepth = -0.03
                else: strongGuesses += 1
                lineGuesses.extend([guessDepth, thisCenter, (guessRight-guessLeft)/2.])
            
#            temp=raw_input("Bad guesses, trying again.")
            smoothSpec.specfit(fittype='gaussian',guesses=lineGuesses, quiet=True, save=False, show_components=False, use_window_limits=True, vheight=True, exclude=[smoothSpec.xarr[loApPix].value, smoothSpec.xarr[loPix].value, smoothSpec.xarr[hiPix].value, smoothSpec.xarr[hiApPix].value])

            lineResults = smoothSpec.specfit.modelpars[:]
        
        
        # The chi-squared value gives us a way of measuring the quality of the fit. 
        # We should probably use it, but I haven't decided how, yet.
#        print("Guess array:")
#        for idx in range(0, len(lineGuesses), 3):
#            print("            {0:4.3f} - {1:2.2f} - {2:1.3f}".format(lineGuesses[idx+1], lineGuesses[idx], lineGuesses[idx+2]))
#        temp=raw_input('{0:2.1f} line ({1:4.3f}) fit with chi^2 of:{2:4.3f}'.format(ion, wl, smoothSpec.specfit.chi2))
        
        if len(lineResults) <= 3:
            lineAmplitude, lineCenter, lineFWHM = lineResults
        else:
            lineCenter = u.closest(wl, [lineResults[x] for x in range(1, len(lineResults),3)])
            bestIndex = lineResults.index(lineCenter)
            lineAmplitude  = lineResults[bestIndex-1]
            lineFWHM = lineResults[bestIndex+1]
            
        if lineFWHM < twopixels: 
            lineFWHM = twopixels
        elif lineFWHM > k.WideFWHMLimit:
            EQWResult = 0.0
            badLines.append([wl, ion, k.lineTooBroad, '{0:2.1f}:{1:4.3f} FWHM:{2:3.3f}'.format(ion, wl, lineFWHM)])
            continue
        xmin=lineCenter-2.*lineFWHM
        if xmin<loWL: xmin=loWL+0.25
        xmax=lineCenter+2.*lineFWHM
        if xmax<hiWL: xmax=hiWL-0.25
        try:
            # Set plot=True for visual representation of the EQW fit as a green box. 
            # Note: the plotting window will still be drawn if plot=False
            EQWResult = smoothSpec.specfit.EQW(plot=True, xmin=xmin, xmax=xmax, xunits='Angstrom', midpt_location='fitted')*1000.     
        except ValueError:
            EQWResult = 0.0
            badLines.append([wl, ion, k.lineTooBroad, k.unknownStr])
            continue
        except OverflowError:
            reply = raw_input('OverflowError. Spectra: {0}, Line:{1:2.1f}:{2:4.3f}.\n      Spec. mean:{3:2.3f}, std:{4:3.3f}.\n Continue (Y/n)'.format(spectraFilename,ion,wl, np.mean(smoothSpec.data), np.std(smoothSpec.data)))
            if reply == 'n':
                exit()
            else:
                badLines.append([wl, ion, k.unknownReason, k.unknownStr])
                continue

        if pauseBetweenLines: 
            reply = raw_input('Hit <return> to accept, or any \'n\'+<return> to skip: ')
            if reply == 'n':
                if g.VerboseMode: 
                    print('Ignoring {0:3.1f} line at: {1:8.3f}.'.format(ion, wl))
                    print('specfit parameters:')
                    print('specfit.modelpars: ',smoothSpec.specfit.modelpars)
                    print('specfit.modelerrs: ',smoothSpec.specfit.modelerrs)
                badLines.append([wl, ion, k.unknownReason, k.unknownStr])
                continue

        if g.VerboseMode: print('{3:2.1f} line at: {0:4.3f} : {1:3.1f}mA (FWHM = {2:2.4f}, chi^2={4:2.4f})'.format(wl, EQWResult, lineFWHM, ion, smoothSpec.specfit.chi2))

        # We're going to keep on eye on all the differences between measured
        # and expected line centers.
        # If we end up with a large systematic delta, we will recommend 
        # adjusting the spectra for an 
        # additional Doppler adjustment (and recommend re-running the valuation)
        thisDoppler = ((lineCenter - wl)/wl) * k.SpeedofLight
        dopplerCheck += thisDoppler
        dopplerCount += 1

        EQWidths.append([wl, ion, ep, gf, c6, EQWResult, lineCenter, lineFWHM*1000])
   
    pyplot.close(plotFigure)
    return EQWidths, badLines, dopplerCount, dopplerCheck
'''       
        

### Once upon a time, we had a far more complex fitting algorithm. 
### It was ugly, but included here for posterity.

'''            localBL = findLocalBaseline(sp, wlLower, wlUpper)
            if (not utilities.is_number(localBL)) or (localBL < 0.10):
            # The local spectrum is so bad, we couldn't even fit a baseline...
                badLines.append([wl, ion, k.poorBaselineFit, str(localBL)])
                return EQWidths, badLines, dopplerCount, dopplerCheck
                
            if g.VerboseMode: print('Modified baseline = {0:5.3f}'.format(localBL))

        sp.baseline(subtract=False, order=0, prefit_continuum = localBL)
        
        # Comment in fitters.py (fitters.py line 177) says:
        # guesses must have the format: 
        #    [height, amplitude, center, width]
        # 'height' refers to the background level (not particularly applicable for
        # optical spectra) and 'amplitude' is the height/depth of the line
        # However, guesses in that form do not work.
        # [amplitude, center, width] does seem to work.

        LineGuess=[-0.5, wl, twopixels]
        sp.specfit(fittype='gaussian',guesses=LineGuess, quiet=True)
        lineAmplitude, lineCenter, lineFWHM = sp.specfit.modelpars
        
        # The following are indicators of bad fits:
        #   EQW is small (FWHM < k.NarrowLineFWHMLimit)
        #   Emission line (amplitude > 0
        #   Too broad (usually a result of a heavily populated neighborhood,
        #      a weak line, bad continuum fit, or a combination)
        #   Line measured is far from the expected center.
        badFitResolutionCount = 0
        if not IsAGoodFit(wl, lineCenter, lineFWHM, lineAmplitude):
            guessAmplitude = -0.5
            guessCenter = wl
            guessWidth = lineFWHM
            while badFitResolutionCount <= 3:
                badFitResolutionCount += 1
                if lineFWHM < k.NarrowLineFWHMLimit or lineAmplitude>0.:
                    # Try a smaller sample (works for close or small lines)
                    # Note: this might benefit from a re-assessment of the local continuum
                    if guessWidth == 0:
                        guessWidth = k.HalfWindowWidth/20.
                    else:
                        guessWidth = guessWidth/2
                elif lineFWHM > k.BroadLineFWHMLimit:
                    # Broad lines can come from a variety of sources, so go use some "smarts"
                    badFitResolutionCount, EQWResult, lineCenter, lineFWHM = FitBroadLine(sp, wl, twopixels, pauseBetweenLines)
                    if badFitResolutionCount == 0:
                    # Good re-fit, ship it!
                        if g.VerboseMode:
                            print('Re-fit successful, wl:{0:8.3f}  eqw:{1:5.1f} (FWHM:{2:6.4f}).'.format(wl, EQWResult, lineFWHM))
                        EQWidths.append([wl, ion, ep, gf, c6, EQWResult, lineCenter, lineFWHM*1000])
                        # Since we've already fit, force a bail
                        badFitResolutionCount = -1
                        break
                    else:
                    # Aww. We need more smarts!
                        if g.VerboseMode:
                            print('Re-fit unsuccessful, wl:{0:8.3f}'.format(wl))
                        badLines.append([wl, ion, k.lineTooBroad, '{0:4.3f}'.format(-lineFWHM)])
                        break
                elif (abs(lineCenter-wl)>k.WLRegionTolerance/2.):
                    # Look for a really small line, closer
                    guessAmplitude = lineAmplitude/2.
                    guessWidth = guessWidth/4.
                else:
                    # Made it through, so everything's hunky dorey!
                    badFitResolutionCount = 0
                    lineAmplitude, lineCenter, lineFWHM = sp.specfit.modelpars
                    break
                if guessWidth > k.HalfWindowWidth/20:
                    guessWidth = k.HalfWindowWidth/20
                elif guessWidth < twopixels:
                    guessWidth = twopixels
                LineGuess=[guessAmplitude, guessCenter, guessWidth]
                sp.specfit(fittype='gaussian',guesses=LineGuess, quiet=True)
                lineAmplitude, lineCenter, lineFWHM = sp.specfit.modelpars
                continue
                
            if badFitResolutionCount > 0:
                # We still have a problem, so bail
                reasonCode = k.unknownReason
                reasonStr = 'Unknown'
                if lineFWHM < k.NarrowLineFWHMLimit:
                    # Still nothing...
                    if g.VerboseMode: print('Estimated line width is too small. Unable to measure this line.')
                    reasonCode = k.lineTooNarrow
                    reasonStr = '{0:4.3f}'.format(-lineFWHM)
                elif lineAmplitude >0.:
                    # Still nothing...
                    if g.VerboseMode: print('Emission line measured. Unable to measure this line.')
                    reasonCode = k.emissionLine
                    reasonStr = '{0:4.3f}'.format(lineFWHM)
                elif lineFWHM > k.BroadLineFWHMLimit:
                    # Still suspiciously broad
                    if g.VerboseMode: print('Estimated line width is unusually broad. Unable to measure this line.')
                    reasonCode = k.lineTooBroad
                    reasonStr = '{0:4.3f}'.format(-lineFWHM)
                elif (abs(lineCenter-wl)> k.WLRegionTolerance):
                    # Still off by...
                    if g.VerboseMode: print('Estimated line center far from expected center (Delta={0:6.1f} mA). Unable to measure this line.'.format((lineCenter-wl)*1000.))
                    reasonCode = k.centerWLOff
                    reasonStr = '{0:4.3f}'.format(lineCenter)
                    
                badLines.append([wl, ion, reasonCode, reasonStr])
                if pauseBetweenLines: raw_input('Hit <return> to continue.')
                continue
            elif badFitResolutionCount < 0:
                # Fit was resolved and recorded elsewhere
                if pauseBetweenLines: raw_input('Hit <return> to continue.')
                continue

        if lineFWHM < twopixels:
            lineFWHM = twopixels
        
        try:
#            EQWResult = sp.specfit.EQW(plot=g.VerboseMode, xmin=lineCenter-2.*lineFWHM, xmax=lineCenter+2.*lineFWHM, xunits='Angstrom')*1000.     
            EQWResult = sp.specfit.EQW(plot=False, xmin=lineCenter-2.*lineFWHM, xmax=lineCenter+2.*lineFWHM, xunits='Angstrom')*1000.     
'''
### This function is also no longer needed:
'''
def FitBroadLine(sp, wl, minWidth, pauseBetweenLines=False):
# Broad line measurements are, um, bad, mmmmmk?
# To make them better, we need to try several options
    if minWidth < k.NarrowLineFWHMLimit:
        narrowWidth = k.NarrowLineFWHMLimit
    else:
        narrowWidth = minWidth
        
    lineAmplitude, lineCenter, lineFWHM = sp.specfit.modelpars
    if lineFWHM > k.BadBaselineLimit:
    # REALLY broad - it's probably a baseline issue.
    # Let's re-fit, and see if the line width changes.
        sp.baseline(subtract=False, order=0)
        LineGuess=[-0.5, wl, narrowWidth*2.]
        sp.specfit(fittype='gaussian',guesses=LineGuess, quiet=True)
        newAmp, newCenter, newFWHM = sp.specfit.modelpars
        # Sometimes, these fits are so bad that the center is outside of the
        # current spectral order. Don't do that...
        # Usually, lines this close to the edge of an order are 'sketchy' -
        # perhaps we should just bail, and fail on this line...
        if newCenter <  sp.xarr.min().value:
            newCenter = sp.xarr.min().value + 2.*newFWHM
        elif newCenter > sp.xarr.max().value:
             newCenter = sp.xarr.max().value - 2.*newFWHM
           
        if newFWHM < k.BroadLineFWHMLimit:
        # O.K. looks good, let's estimate this line as the measured line,
        # plus a rectangular area, between the new baseline and continuum,
        # with a width of twice the FWHM. An alternative with a better 
        # difference estimation could involve fitting points in the measured
        # Gaussian into a new Gaussian, with an amplitude based on the
        # continuum.
            newBaseline = sp.baseline.baselinepars[0]
                # Assuming an order=0 fit, the constant baseline value is the
                # first, last, and only value in the fit polynomial.
            if newFWHM < narrowWidth:
                newFWHM = narrowWidth
            rectangularArea = 2*(1.0 - newBaseline)*newFWHM
            try:
                EQWMeasure = sp.specfit.EQW(plot=False, xmin=newCenter-2.*newFWHM, xmax=newCenter+2.*newFWHM, xunits='Angstrom')
            except ValueError:
            # Has happened when the line center is outside of the spectral range
                print 'Value Error:', newCenter, newFWHM, sp.xarr.dxarr[0], sp.xarr.dxarr[-1]
                exit()
            return 0, (EQWMeasure+rectangularArea)*1000., sp.specfit.modelpars[1], newFWHM
        #else:
        # No good. Do something else! Right now, we just fall through 
        # to the return 4,0,0,0 error case
    else:
    # Just kinda broad - Let's guess that this is a blended line:
    # We are going to make several assumptions in our fitting:
    #   - The lines to fit all occur within the edges of the FWHM
    #       boundaries
    #   - There are a limited number of lines that can fit here
    #       (Currently set at 7, max.) set here:
        maxLinesToFit = 7
    #       Note: High values of lines to fit make better measurements,
    #       at a geometrically growing hit to performance.
    #   - The multiplet will be centered near the desired wl. Note:
    #       This really means that the fit will be centered there...
        centerWL = wl
        fitWidth = lineFWHM
        leftWL = centerWL - fitWidth/2.
        # Start with a multiple line fit:
        lineCounter = 2
        while lineCounter <= maxLinesToFit:
            allGuesses=[[lineAmplitude, leftWL+x, narrowWidth*2.] for x in np.arange(lineCounter)*(fitWidth/(lineCounter-1))]
            # We need to flatten the sublists to make a guess array
            LineGuesses = [x for guess in allGuesses for x in guess]
            sp.specfit(fittype='gaussian',guesses=LineGuesses, quiet=True)
            # Find the fit line, closest to the desired wl:
            closestIndex = sp.specfit.modelpars.index(min(sp.specfit.modelpars, key=lambda x:abs(x-wl)))
            newAmp = sp.specfit.modelpars[closestIndex-1]
            newCenter = sp.specfit.modelpars[closestIndex]
            newFWHM = sp.specfit.modelpars[closestIndex+1]
            if newFWHM < narrowWidth: newFWHM = narrowWidth
            if pauseBetweenLines:
                print('Guessing multi-fit Amplitude: {0:5.3f}   Center: {1:8.3f}   FWHM: {2:6.4f}'.format(newAmp, newCenter, newFWHM))
                raw_input('Hit <return> to continue.')            
            if IsAGoodFit(wl, newCenter, newFWHM, newAmp):
                # O.K. looks good.
#                EQWMeasure = sp.specfit.EQW(plot=g.VerboseMode, xmin=newCenter-2.*newFWHM, xmax=newCenter+2.*newFWHM, xunits='Angstrom')
                EQWMeasure = sp.specfit.EQW(plot=False, xmin=newCenter-2.*newFWHM, xmax=newCenter+2.*newFWHM, xunits='Angstrom')
                # Re-fitting occasionally produces an emission line. We're done at that point.
                if EQWMeasure > 0:
                    return 0, (EQWMeasure)*1000., sp.specfit.modelpars[1], newFWHM
            # else: Keep subdividing the fit
            lineCounter += 1
            
    # No luck. Return the obviously bad result
    return 4, 0., 0., 0.
'''
### ...and this one:
'''
def IsAGoodFit(wl, lineCenter, lineFWHM, lineAmplitude):
    # The following are indicators of bad fits:
    #   EQW is small (FWHM < k.NarrowLineFWHMLimit)
    #   Emission line (amplitude > 0
    #   Too broad (usually a result of a heavily populated neighborhood,
    #      a weak line, bad continuum fit, or a combination)
    #   Line measured is far from the expected center.
    return (lineFWHM > k.NarrowLineFWHMLimit and
    lineFWHM < k.BroadLineFWHMLimit and
    lineAmplitude < 0. and
    abs(lineCenter-wl) < k.WLRegionTolerance/2.)

'''
