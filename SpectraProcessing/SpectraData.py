#! /usr/bin/env python
# 
# Module: SpectraData.py
#
# Author: Mike Lum
#
# Description: Functions to get and store data on spectra/spectra files
#
# Contents:
#   Function xxx: (Description)
#
# Revision History:
#    Date        Vers.    Author        Description
#    03-01-2016  0.0a0    Lum        First checked in
#
# To Do:
#    
#

# Imports
from constants import constants as k
from constants import DBFormats as DB
from utilities import utilities as u

from astropy.io import fits as pyfits
import datetime as DT
from pyraf import iraf
import numpy as np
import sqlite3
# import pyspeckit
import random

# Functions
def IsGoodAperture(aperture):
#   Function: IsGoodAperture
#
#   Input: spectra (or spectra order)
#   Output: Boolean representing whether the mean/std of the spectra data lie within
#           acceptable tolerances.
#    specMean = np.mean(aperture.data)
#    specStd = np.std(aperture.data)
#    if specMean < 0.2 or specMean>5.0 or specStd > 10.0:
#        return False
#    else:
#        return True
# Ugh! This only works for continuumed spectra
    return True
    
def GetInstrument(fullPathName):
# Try to determine the spectrograph from fits header keywords
    try:
        instrument = pyfits.getval(fullPathName, k.FitsKwdInstrument).split()
    except KeyError:
    # So far, this has only occurred with SOPHIE spectra, so...
        instrument = [k.SOPHIEInstrument]
    except IOError:
        print('Error in reading fits header for file: {0}. Exiting.'.format(fullPathName))
        exit()
    return instrument

def GetSpectraDatetime(fullPathName):
#   Function: GetSpectraDatetime
#
#   Input: Absolute path to a filename containing a .fits spectra
#   Output: String containing the date and time of the observation.
#       Note: The format of the date/time will vary with the observatory

    dateStr = ''
    timeStr = ''
    instrument = GetInstrument(fullPathName)
    shortInst = instrument[0]
    
    try:
        thisHDU = pyfits.open(fullPathName)
    except:
        # Should use some sort of error output here.
        print ("Unable to open {0} for header information. Continuing...".format(fullPathName))
        return dateTimeStr

    if shortInst == k.HiresInstrument1 or shortInst == k.HiresInstrument2 or shortInst == k.HDSInstrument:
    # Keck is pretty consistant with the date keyword (Subaru spectra are untested):
        dateStr = thisHDU[0].header[k.FitsKwdUTDate]
        try:
            timeStr = thisHDU[0].header[k.FitsKwdUTTime]
        except KeyError:
        # Two possibilities for the time keyword...
            try:
                timeStr = thisHDU[0].header[k.FitsKwdUTTime2]
            except KeyError:
            # No go for either, assign a random one
                timeStr=(DT.datetime(1900,1,1,0,0,0)+DT.timedelta(seconds=random.randrange(86400))).strftime('%H:%M:%S')
    else:
    # Assume OHP. (and any other?)
        try:
            dateStr = thisHDU[0].header[k.FitsKwdSODate]
        except KeyError:
        # Gonna bail for now. Maybe use file info in the future?
            pass
        # Assigning a random time for now. OHP does put a timestamp in, I'm just
        # too lazy to look at the header and figure it out
        timeStr=(DT.datetime(1900,1,1,0,0,0)+DT.timedelta(seconds=random.randrange(86400))).strftime('%H:%M:%S')
    
    return dateStr + ' ' + timeStr


def GetSpectraRADec(fullPathName):
#   Function: GetSpectraRADec
#
#   Input: Absolute path to a filename containing a .fits spectra
#   Output: Strings containing the RA and DEC for the observation.
#       Note: Coordinate formats may vary with the observatory

    RAStr = ''
    DECStr = ''
    
    try:
        thisHDU = pyfits.open(fullPathName)
    except:
        # Should use some sort of error output here.
        print("Unable to open {0} for header information. Continuing...".format(fullPathName))
        return RAStr, DECStr
    try:
    # Keck uses HH:MM:SS.S +/-DD:MM:SS.S
    # HARPS uses DDD.DD +/-DD.DD
        RAStr = '{0}'.format(thisHDU[0].header[k.FitsKwdRA])
        DECStr = '{0}'.format(thisHDU[0].header[k.FitsKwdDEC])
    except KeyError:
        try:
        # OHP uses a different field name, and format, of course
            RARes = thisHDU[0].header[k.FitsKwdSORA]
            if u.is_number(RARes):
                OHPRAStr = "%.2f"%RARes
            else:
                OHPRAStr = RARes
            DECRes = thisHDU[0].header[k.FitsKwdSODEC]
            if u.is_number(DECRes):
                OHPDECStr = "%.2f"%DECRes
            else:
                OHPDECStr = DECRes
            RAStr = '{0}:{1}:{2}'.format(OHPRAStr[:2],OHPRAStr[2:4],OHPRAStr[4:])
            DECStr = '{0}:{1}:{2}'.format(OHPDECStr[:3],OHPDECStr[3:5],OHPDECStr[5:])
        except KeyError:
        # Gonna bail for now.
            pass
    
    return RAStr, DECStr


def GetSpectraSN(fullPathName):
#   Function: GetSpectraSN
#
#   Input: Absolute path to a filename containing a .fits spectra
#   Output: Value for signal-to-noise ratio as reported by the spectra file.
#       Note: Coordinate formats may vary with the observatory

    SN = 0.
    
    try:
        thisHDU = pyfits.open(fullPathName)
    except:
        # Should use some sort of error output here.
        print("Unable to open {0} for header information. Continuing...".format(fullPathName))
        return SN

    instrument = GetInstrument(fullPathName)
        
    if k.HiresInstrument1 in instrument or k.HiresInstrument2 in instrument or k.HDSInstrument in instrument:
        # Keck uses Float
        SN = float(thisHDU[0].header[k.FitsKwdHIRESSN])
    elif k.HarpsInstrument in instrument:
        # HARPS has one for each order. Average them
        SN = np.mean([float(thisHDU[0].header[i]) for i in thisHDU[0].header.keys() if k.HarpsSNPrefix in i])
    else:
        SN = 0.
    return SN


def GetSpectraDispersion(fullPathName):
#   Function: GetSpectraDispersion
#
#   Input: Absolute path to a filename containing a .fits spectra
#   Output: Float of the Angstrom/pixel dispersion of the spectra.

    dispStr = '0.0'
    try:
        thisHDU = pyfits.open(fullPathName)
    except:
        # Should use some sort of error output here.
        print("Unable to open {0} for header information. Continuing...".format(fullPathName))
        return disp
    try:
        dispStr = thisHDU[0].header[k.FitsKwdDisp]
    except KeyError:
        try:
        # Keck can use a different field 
        # (this is the dispersion for the first order, only)
            dispStr = thisHDU[0].header[k.FitsKwdKHOrderDisp]
        except KeyError:
        # Gonna bail for now. Maybe use file info in the future?
            pass
    try:
        return float(dispStr)
    except ValueError:
        return 0.0

def GetCCDNum(fullPathName):
#   Function: GetCCDNum
#
#   Input: Absolute path to a filename containing a .fits spectra
#   Output: CCD number on the HIRES (Keck) Spectroscope

    ccdNo = 1
    try:
        thisHDU = pyfits.open(fullPathName)
    except:
        # Should use some sort of error output here.
        print("Unable to open {0} for header information. Continuing...".format(fullPathName))
        return ccdNo
    try:
        ccdNoStr = thisHDU[0].header[k.FitsKwdKHCCDLoc]
        ccdNo = int(ccdNoStr)
    except (KeyError, ValueError):
    # About to return "1", anyway...
        pass
        
    return ccdNo


def GetWLsAndFlux(filename):
    # Opens a fits spectrum, and returns a list containing the wavelength -
    # flux list pairs for each order.
    # [[wl[0], flux[0]], [wl[1],flux[1]]...]
    hdu = pyfits.open(filename)
    header = hdu[0].header
    data = hdu[0].data
    
    # Is this a 1-D spectrum?
    if header[k.FitsKwdNumAxes] == 1:
        # This is easy! First: Produce a wl solution
        baseWL = header[k.FitsKwdRVal]
        delta = header[k.FitsKwdDisp]
        numPoints = int(header[k.FitsKwdLen1DAxis])
    
        endWL = baseWL + delta*numPoints
        wls = np.linspace(baseWL, endWL, num=numPoints)
        
        # The data's just...the data. :P
        fluxes = data
        return [[wls, fluxes]]
    # 2-D (multiple order) spectrum
    else:
    # Not yet supported. Need to get the wl solutions from the WV_0_N cards
        return [[None, None]]
        
    
def LoadAndSmoothSpectra(specFilename):
    # Loads a (multi-order) spectrum, and returns a list
    # of spectra, one item for each order (or a single item for a 1-D spectrum)
    # We have removed much of the functionality for multi-order spectra, as it
    # used pyspeckit. Since our, continuum/Doppler correction function produces
    # a 1-d spec, we're ok for now...
    
    instrument = GetInstrument(specFilename)
        
    if k.HiresInstrument1 in instrument or k.HiresInstrument2 in instrument or k.HDSInstrument in instrument:
        specType = pyfits.getval(specFilename, k.FitsKwdCtype).split()
        if specType[0] == k.HiresArchiveReduction1 or specType[0] == k.HiresArchiveReduction2:
            specList = GetWLsAndFlux(specFilename)
        else:
            print('Unrecognized Spectrum Type: {0}'.fotmat(specType[0]))
#            specList=pyspeckit.wrappers.load_IRAF_multispec(specFilename)
    elif k.HarpsInstrument in instrument \
        or k.ElodieInstrument in instrument \
        or k.SOPHIEInstrument in instrument \
        or k.UVESInstrument in instrument \
        or k.NSOInstrument in instrument \
        or k.NOAOInstrument in instrument:
        sspecList = GetWLsAndFlux(specFilename)
    else:
        print('Unrecognized instrument: {0}. Exiting.'.format(instrument))
        exit()
    
    # Smooth all the apertures
    smoothList = []
    for wls, fluxes in specList:
    # Removed for un-continuumed spectra:
#        # Kill CR hits
#        for idx in range(len(roughSpec.data)):
#            if roughSpec.data[idx] > 10.:
#                roughSpec.data[idx] = 1.0
#            elif roughSpec.data[idx] < -3.:
#                roughSpec.data[idx] = 0.0
        for idx in range(len(fluxes)):
            if fluxes[idx] < 0.:
                fluxes[idx] = 0.0
        smoothList.append([wls[1:-1], u.boxcar(fluxes, size=3)])
    return smoothList


def GetObjectName(spectraFilename):
#   Function: GetObjectName
#
#   Input: Absolute path to a filename containing a .fits spectra
#   Output: objectName
#               objectName: Name as reported by the spectra in its .fits header
#                           or a generic instrument object name, if we can't 
#                           find a specific object name.
    objectName = ''
    try:
        testName = pyfits.getval(spectraFilename, k.FitsKwdObjectName)
        objectName = pyfits.getval(spectraFilename, k.FitsKwdTargetName)
    except:
    # We will have to assign the object name based on the instrument, below
        try:
            instName = pyfits.getval(spectraFilename, k.FitsKwdInstrument)
            if instName == k.SOPHIEInstrument:
                objectName = k.SOPHIEObjectName
            elif instName == k.ElodieInstrument:
                objectName = k.ElodieObjectName
            elif instName == k.UVESInstrument:
                objectName = k.UVESObjectName
            elif instName == k.HarpsInstrument:
            # Note: HARPS is a ESO instrument, and uses the same header format as
            # FEROS and UVES, so there shouldn't be much difference between HARPS 
            # and UVES
                objectName = k.HarpsObjectName
        except:
        # No instrument name (or Keck HIRES), just use the object as assigned
            pass
    
    if objectName is '':
        # So far, this has only happened with SOPHIE spectra
        objectName = k.SOPHIEObjectName
    return objectName
    
    
def ReadSpectraInfo(spectraFilename):
#   Function: ReadSpectraInfo
#
#   Input: Absolute path to a filename containing a .fits spectra
#   Output: (objectName, numApertures, apertureSize, apBoundaries)
#               objectName: Name as reported by GetObjectName, above.
#               numApertures: Number of apertures in a 2-D spectrum
#               apertureSize: The size of an aperture (assumed constant for multiple-aperture spectra)
#               apBoundaries: List of tuples containing (startWL, endWL) for each aperture

    numApertures = 0
    apertureSize = 0
    apBoundaries = []
    
    workingFilename = spectraFilename
    
    # What are we looking at?
    objectName = GetObjectName(spectraFilename)

    # Get the aperture information from the fits spectra
#    numApertures = int(pyfits.getval(workingFilename, k.FitsKwdNumApertures))
#    apertureSize = int(pyfits.getval(workingFilename, k.FitsKwdApertureSize))
    numApertures = 0
    # Some fits headers will tell us that they're only 1-d
    try:
        numApertures = int(pyfits.getval(workingFilename, k.FitsKwdNumAxes))
    except KeyError:
        # Odd...since that keyword is kinda required by the fits spec...
        # but...ok.
        pass
    if numApertures == 0:
        numApertures = int(pyfits.getval(workingFilename, k.FitsKwdNum2DAps))
    apertureSize = int(pyfits.getval(workingFilename, k.FitsKwdLen1DAxis))

    if objectName == k.SOPHIEObjectName or objectName == k.UVESObjectName or objectName == k.HarpsObjectName or numApertures == 1:
        # All these sources are 1-d, including 1-d HIRES spectra
        temp = iraf.listpixels(workingFilename, wcs='world', formats='%10.3f', Stdout=1)
        wlMin = float(temp[0].split()[0])
        wlMax = float(temp[-1].split()[0])
        apBoundaries.append((wlMin,wlMax))
    else:
        for aperture in range(1, numApertures+1):
        # Note: Keck Hires images have lowest wls in higher # apertures
            wlMin = float(pyfits.getval(workingFilename, '{0}{1:02d}'.format(k.FitsKwdKHOrderBase, aperture)))
            winWidth = float(pyfits.getval(workingFilename, '{0}{1:02d}'.format(k.FitsKwdKHOrderDisp, aperture)))
            wlMax = wlMin + winWidth*apertureSize
#            temp = iraf.listpixels(workingFilename+'[1,{0:d}]'.format(aperture),wcs='world',formats='%10.3f',Stdout=1)
#            print len(temp), float(temp[0].split()[0]), float(temp[-1].split()[0])
#            wlMin = float(temp[0].split()[0])
#            temp = iraf.listpixels(workingFilename+'[{0:d},{1:d}]'.format(apertureSize-1,aperture),wcs='world',formats='%10.3f',Stdout=1)
#            print len(temp)
#            wlMax = float(temp[0].split()[0])
            apBoundaries.append((wlMin,wlMax))

    return objectName, numApertures, apertureSize, apBoundaries
    # Should probably create a data object...


def ParseSpectraName(spectraFilename):
#   Function: ParseSpectraName
#
#   Input: Absolute path to a filename containing a .fits spectra
#   Output: (Cluster num., Cat. name, Cat. Number, Inst., Spectra series #)
#               Cluster num.: Number (generally NGC-###), as reported by the 
#               Catalog name: Author, or other name (ex: HIP) for star index
#               Catalog num.: Number within the above catalog for this star
#               Instrument  : Code for the instrument used to take the spectrum
#               Spectra series # : An internal index of the count for this
#                               spectrum. Used to differentiate between
#                               multiple spectra of the same star from the 
#                               same instrument.
    # We assume that the input spectra filename is in the desired format: 
    #           NNNN_SSS_ssssss_II_OOO.fits
    # Where: 
    #           NNNN  = NGC catalog number of parent cluster 
    #                    (note: clusters not contained in the NGC catalog will
    #                    be prefaced by a letter representative of their
    #                    catalog. ex: I = IC)
    #           SSS    = Reference indexing system for this cluster. 
    #                    (ex: MJP = Montgomery, Janes, and Phelps)
    #           ssssss = Star index number from the aformetioned index.
    #           II     = Two letter code for the observing telescope/instrument
    #                    reference: 
    #                               KH = Keck Hires
    #                               VU = VLT UVES
    #                               EL = Elodie
    #                               SO = Sophie
    #                               HA = HARPS
    #                               SH = Subaru HDS
    #           OOO    = Internal numbering for this spectra.
    
    try:
        nameParts = spectraFilename.split('_')
        if len(nameParts) == 5:
        # 1-d or 2-d spectrum
            clusterNo, catName, catNum, inst, endString = nameParts
        else:
        # Single orders from a 2-d spectrum
            clusterNo, catName, catNum, inst, orderNo, endString = nameParts
         
        specCount = endString.replace('.fits','')
    except:
        print("Unable to split/parse input spectra filename:{0}.\n Continuing.".format(spectraFilename))
        return k.BadInputNumberString, k.BadInputNameString, k.BadInputNumberString, k.BadInputNameString, k.BadInputNumberString
    
    return clusterNo, catName, catNum, inst, specCount


def EnterStarName(fullPathName):
#   Function: EnterStarName
#
#   Input: Checks if the passed spectra is in our database. Enters it if not.
#   Output: <N/A>
#
    spectraFilename = fullPathName.split('/')[-1]
#    print ('Entering {0}.'.format(spectraFilename))
    
    with sqlite3.connect(DB.ClusterDBName) as DBConnection:
        DBCursor = DBConnection.cursor()
        DBCursor.execute('''SELECT * FROM {0} WHERE Filename IS \"{1}\";'''.format(DB.SpectraTableName, spectraFilename))
        if len(DBCursor.fetchall()) > 0:
        # It's already there, so we're done.
            print("Already entered")
            return
        else:
        # Time to enter it...
            clusterNo, catName, catNum, inst, specCount = ParseSpectraName(spectraFilename)
            clusterStr = 'NGC-'+clusterNo
            catStr = catName+'-'+catNum
            target, numAps, apSize, apRanges = ReadSpectraInfo(fullPathName)
            minWL = min([x[0] for x in apRanges])
            maxWL = max([x[1] for x in apRanges])
            wlRange = '{0:4.3f}-{1:4.3f}'.format(minWL, maxWL)
            dateTimeStr = GetSpectraDatetime(fullPathName)
            RAStr, DECStr = GetSpectraRADec(fullPathName)
            disp = GetSpectraDispersion(fullPathName)
            if inst == 'KH':
                ccdNo = GetCCDNum(fullPathName)
            else:
                ccdNo = 1
                
            newStar = [RAStr, DECStr, 2000, clusterStr, catStr, 0.,0.,0., 0., 0.]
            try:
                DBCursor.execute('''INSERT INTO \"{0}\" VALUES {2}'''.format(DB.StarTableName, DB.StarFieldNames, DB.StarGenFormat), newStar)
            except sqlite3.IntegrityError:
            # It's actually in there already - Integrity error occurs when
            # unique ID fields are being duplicated
                print('Exists in {0} table.'.format(DB.StarTableName))
                pass
            newSpectra = [dateTimeStr, k.InstDict[inst], ccdNo, spectraFilename, disp, wlRange, clusterStr, catStr]
            try:
                DBCursor.execute('''INSERT INTO \"{0}\" VALUES {2}'''.format(DB.SpectraTableName, DB.SpectraFieldNames, DB.SpectraGenFormat), newSpectra)
            except sqlite3.IntegrityError:
            # In the DB already
                print(newSpectra)
                print('Exists in {0} table.'.format(DB.SpectraTableName))
                pass
                
