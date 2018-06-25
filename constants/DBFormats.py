#! /usr/bin/env python
# 
# Module: DBFormats.py
#
# Author: Mike Lum
#
# Description: Constants and table entry formats for the 
#       Atomic_Line, Cluster_Stars, Spectra, and Measured_Lines
#       databases
#
# Contents: List externally-callable functions here
#   Function xxx: (Description)
#
# Revision History:
#    Date        Vers.    Author        Description
#    2015-07-20  1.0a0    Lum           First checked in
#
# To Do:
#    
#

# Constants

# Cluster database constants
#LineDBName = "Databases/Atomic_Lines"
ClusterDBName = 'Databases/Cluster_Data'

# Line reference table
LineRefTableName = 'linerefs'
LineRefTableFormat = '\"{0}\" (RefNo INTEGER NOT NULL, FirstAuthor TEXT, Year TEXT, Ref TEXT, BibtexRef TEXT, Notes TEXT, PRIMARY KEY (RefNo))'.format(LineRefTableName)
LineRefGenFormat = '(?,?,?,?,?,?)'

# Line data table
LineTableName = 'lines'
LineTableFormat = '\"{0}\" (Wavelength REAL NOT NULL, Ion REAL NOT NULL, ExPot REAL, Loggf REAL, vdwDamping REAL, Refs TEXT, RefRes TEXT, notes TEXT, ForParams TEXT, Blended INT, HFS INT, PRIMARY KEY(Ion, Wavelength))'.format(LineTableName)
LineTableFieldNames = ['Wavelength','Ion','ExPot','Loggf','vdwDamping','Refs','RefRes','notes','ForParams', 'Blended', 'HFS']
LineGenFormat = '(?,?,?,?,?,?,?,?,?,?,?)'

# Hyper-fine split Lines Table
HFSTableName = 'hfs_lines'
HFSTableFormat = '\"{0}\" (Wavelength REAL NOT NULL, Ion REAL NOT NULL, BlendWL REAL NOT NULL, Loggf REAL)'.format(HFSTableName)
HFSGeneralFormat = '(?,?,?,?)'

# Cluster star table
StarTableName = 'cluster_stars'
StarTableFormat = '\"{0}\" (RA TEXT, DEC TEXT, Epoch INTEGER, Cluster TEXT, StarID TEXT, Teff REAL, LogG REAL, Vturb REAL, Simbad_VMag REAL, Simbad_BMag REAL, S2N REAL, PRIMARY KEY (Cluster, StarID))'.format(StarTableName, )
StarFieldNames = ['RA', 'DEC', 'Epoch', 'Cluster', 'StarID', 'Teff', 'LogG', 'Vturb', 'Simbad_VMag', 'Simbad_BMag', 'S2N']
StarGenFormat = '(?,?,?,?,?,?,?,?,?,?,?)'

# Spectra table
SpectraTableName = 'spectra'
SpectraTableFormat = '\"{0}\" (DateTime TEXT, Instrument TEXT, CCD TEXT, Filename TEXT, Resolution REAL, WLRange TEXT, Cluster TEXT, StarID TEXT, FOREIGN KEY(Cluster, StarID) REFERENCES {1}(Cluster, StarID), PRIMARY KEY (DateTime, Instrument, CCD)'.format(SpectraTableName, StarTableName)
SpectraFieldNames = ['DateTime', 'Instrument', 'CCD', 'Filename', 'Resolution', 'WLRange', 'Cluster', 'StarID']
SpectraGenFormat = '(?,?,?,?,?,?,?,?)'

# Measured lines table
MeasuredLinesTableName = 'measured_lines'
MeasuredLinesTableFormat = '\"{0}\" (Ion REAL, Wavelength REAL, Filename TEXT, MeasWL REAL, EQM REAL, FWHM REAL, FOREIGN KEY(Ion, Wavelength) REFERENCES {1}(Ion, Wavelength), FOREIGN KEY(Filename) REFERENCES {2}(Filename))'.format(MeasuredLinesTableName, LineTableName, SpectraTableName)
MeasuredLinesFieldNames = ['Ion', 'Wavelength', 'Filename', 'MeasWL', 'EQM', 'FWHM']
MeasuredLinesGenFormat = '(?,?,?,?,?,?)'

# Star/spectra atmospheric parameters
AtmParTableName = 'atmospheric_parms'
AtmParTableFormat = '\"{0}\" (Filename TEXT, Teff REAL, Loggf REAL, Vturb REAL, ParmSource TEXT, FOREIGN KEY(Filename) REFERENCES {1}(Filename))'.format(AtmParTableName, SpectraTableName)
AtmParFieldNames = ['Filename', 'Teff', 'Loggf', 'Vturb', 'ParmSource']
AtmParGenFormat = '(?,?,?,?,?)'

# Solar calibration lines table
CalibrationLinesTableName = 'calibration_lines'
CalibrationLinesTableFormat = '\"{0}\" (Wavelength REAL, Ion REAL, EQWList TEXT, DeltaEQW REAL, Star TEXT, FOREIGN KEY(Wavelength, Ion) REFERENCES {1}(Wavelength, Ion))'.format(CalibrationLinesTableName, LineTableName)
CalibrationLinesFieldNames = ['Wavelength', 'Ion', 'EQWList', 'DeltaEQW','Star']
CalibrationLinesGenFormat = '(?,?,?,?,?)'

# Pre-calculated reference spectra table
RefCorrTableName = 'ReferenceCorrections'
RefCorrTableFormat = '\"{0}\" (ReferenceStar TEXT, Ion REAL, Wavelength REAL, ExPot REAL, Loggf REAL, EQW REAL, LogRW REAL, Abundance REAL)'.format(RefCorrTableName)
RefCorrFieldNames = ['ReferenceStar', 'Ion', 'Wavelength', 'ExPot', 'Loggf', 'EQW','LogRW', 'Abundance']
RefCorrGenFormat = '(?,?,?,?,?,?,?,?)'


