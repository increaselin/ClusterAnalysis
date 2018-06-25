#! /usr/bin/env python
# 
# Program: constants.py
#
# Author: Mike Lum
#
# Usage: n/a
#
# Description: Contains constants used in the eqwm library
#
# Revision History:
#    Date          Vers.    Author        Description
#    2015-06-05    1.1F1    Lum       Added constants for MOOG log processing and the line database
#    2014-09-07    1.0f1    Lum        First checked in - consolidation of many 
#                                       different constant files into one for an 
#                                       overall project consolidation.
#
import sys

if (sys.version_info > (3, 0)):
    isPython3 = True
else:
    isPython3 = False
    
if isPython3:
    from constants.directories import *
else:
    from directories import *
    
# File and OS names (not sure that these are needed outside of command-line testing)
InterpreterName = 'python'
EQWMModuleName = 'eqwm.py'
SplotLogFilename = 'splot.log'
SkippedLineExtension = '.skipped.txt'
LineEvalExtension = '.eval.txt'

OutputFileHeader = '#Line WL.  Atom+Ion  Ex. Pot.  (Log)gf   C6 value          Meas. Equiv. Width\n'

# As of 9/10/16, this file name is used in STRPlots.py
TempAbOutFilename = 'Abundance.out'

tempFilename = 'temp.out'

# Statistical constants
#   Confidence intervals
sigmaFor68 = 1.000
sigmaFor90 = 1.645
sigmaFor95 = 1.960
sigmaFor98 = 2.330
sigmaFor99 = 2.575

# General constants:
NoCoefficientData = -99.99
EmptyList = []
MutipleFiles = 'Composite'
BadInputNameString = "XXX"
BadInputNumberString = "000"

InstDict = {'KH':'HIRES', 'HA':'HARPS', 'SO':'SOPHIE', 'EL':'ELODIE', 'SH':'HDS', 'UV':'UVES', 'NS':'NSO', 'NO':'NOAO'}

# A "cluster" name for calibration and testing stars - potentially
# also used for field stars. Note: May want to change to NGC-0000, in case 
# someone actually wants to do something with NGC-1 (spiral galaxy...)
CalibClusterName = 'NGC-0001'

# fits keywords for echelle spectra (does this work on non-Hires data?)
FitsKwdObjectName = 'OBJECT'
FitsKwdTelescope = 'TELESCOP'
FitsKwdInstrument = 'INSTRUME'
FitsKwdNumAxes = 'NAXIS'
FitsKwdLen1DAxis = 'NAXIS1'
FitsKwdNum2DAps = 'NAXIS2'
FitsKwdApertureSize = 'i_naxis1'
FitsKwdNumApertures = 'i_naxis2'
FitsKwdCtype = 'CTYPE1'
FitsKwdDisp = 'CDELT1'
FitsKwdRVal = 'CRVAL1'
FitsKwdUTTime = 'UT'
FitsKwdUTTime2 = 'UTC'
FitsKwdUTDate = 'DATE-OBS'
FitsKwdRA2000 = 'RA2000'
FitsKwdDEC2000 = 'DEC2000'
FitsKwdRA = 'RA'
FitsKwdDEC = 'DEC'
FitsKwdObsLat = 'LATITUDE'
FitsKwdObsLong = 'LONGITUD'
FitsKwdTargetName = 'TARGNAME'

# HIRES-specific keywords and data
FitsKwdKHOrderBase = 'CRVL1_'
FitsKwdKHOrderDisp = 'CDLT1_'
FitsKwdKHCCDLoc = 'CCDLOC'
HiresRedCCDLoc = '3'
HiresGreenCCDLoc = '2'
HiresBlueCCDName = '1'
HiresInstrument1 = 'HIRES:'
HiresInstrument2 = 'HIRES'
HiresInstrument3 = 'HIRES Spectrograph'
HiresCode = 'KH'
KeckTelescope = 'Keck I'
Keck1Latitude = 19.82658656
Keck1Longitude = -155.4722
Keck1Altitude = 4145. # Meters above mean sea level
KeckTimeZone = -10
HiresArchiveReduction1 = 'LAMBDA'
HiresArchiveReduction2 = 'LINEAR'
FitsKwdHIRESSN = 'SIG2NOIS'

# OHP- ELODIE-, and SOPHIE-specific keywords and data
FitsKwdSOBaseWL = FitsKwdRVal
FitsKwdSODisp = FitsKwdDisp
FitsKwdSODate = 'OHP OBS DATE START'
FitsKwdSORA = 'OHP TARG ALPHA'
FitsKwdSODEC = 'OHP TARG DELTA'
OHPLatitude = 43.93167
OHPLongitude = 5.7122
OHPAltitude = 650.
OHPTimeZone = 2
HDSInstrument = 'HDS'
ElodieInstrument = 'ELODIE'
ElodieCode = 'EL'
ElodieObjectName = 'Elodie Object'
SOPHIEInstrument = 'SOPHIE'
SophieCode = 'SO'
SOPHIEObjectName = 'SOPHIE Object'
OHPCCDName = '1' # Only 1 ccd for OHP

# HARPS-specific keywords and data
HarpsInstrument = 'HARPS'
HarpsObjectName = 'HARPS Object'
HarpsCCDName = '1'
HarpsSNPrefix = 'ESO DRS SPE EXT SN'

# UVES (VLT) keywords
UVESInstrument = 'UVES'
UVESObjectName = 'UVES Object'

# National Solar Observatory reference spectrum (Wallace 2011)
NSOInstrument = 'NSO'
NSOObjectName = 'Solar Reference Spectrum'

# National Optical Astronomy Observatory reference spectrum (Hinkle 2007)
NOAOInstrument = 'NOAO'
NOAOObjectName = 'Arcturus Reference Spectrum'

# The window  height and width doesn't matter too much - they're just for user interaction
HalfWindowWidth = 3.0        
SpecFluxMin = 0.10
SpecFluxMax = 2.00
# The following FWHM line measurement tolerances are in Angstroms
BadBaselineLimit = 2.500
BroadLineFWHMLimit = 0.400
# We can deal with some blending, but the line centers do need to be 
# separated by a certain amount (exact distance undetermined)
# So, we're going to require half of the broadest line width.
MinLineSeparation = BroadLineFWHMLimit/2.
# Must be at least three pixels wide, at a dispersion of 0.01A/pix
# Note: I have found most spectra are using a dispersion of around 0.02A/pix
NarrowLineFWHMLimit = 0.030
# Sometimes the Gaussian fit function decides on a ridiculous fit width...
WideFWHMLimit = 1.000
WLRegionTolerance = 0.500
BadWLValue = -1.0

# We set a minimum eqw measurement, based on a S/N of 40 (semi-arbitrary), a
# Angstrom/pixel dispersion of 0.03 (HiRes) and a multiplier (call it "sigma")
# of 2 to get: (1/20)*(0.0287)*2 = 2.9mA
MinSN = 20.
KeckDisp = 28.7 # mA/pix
EQWMeasSigma = (1./MinSN)*(KeckDisp)
MinEQWMeasurement = 2.*EQWMeasSigma
MinDetFakeStarName = 'Min-001'
MaxDetFakeStarName = 'Max-001'

# Reasons why a line fit would be rejected:
goodFit = 0
unknownReason = -1
unknownStr = 'Unknown'
noLines = 'No lines in region'
poorBaselineFit = -2    # Error code woule be the measured baseline (string)
lineTooBroad = -3
lineTooNarrow = -4
emissionLine = -5
centerWLOff = -6
badSpectralRegion = -7
badEQWMValue = -8
nanStr = 'nan'

# Constant for Doppler corrections: delta(lambda)/lambda = delta(v)/c -> delta(v) = c(km/s)*delta(lambda)/lambda
SpeedofLight = 299792
MinDopplerCorrection = 1.0
invcmperev = 1.239842*(10**(-4))

SolarTeff = 5777
SolarLogG = 4.44
SolarVturb = 1.52
SolarFeH = 7.50

# The logG limit to delineate dwarfs from giants:
# Note: The selection should be photometry-based...
dwarfGiantLimit = 3.5
giantLogGLimit = dwarfGiantLimit
giantTeffLimit = 5200.

# Range(s) for parameter determination, based on available ATLAS or MARCS models
TEffRange = [3800., 7500.]
LogGRange = [1.8, 5.0]
VturRange = [0., 5.]
TeffStep = 10.
LogGStep = 0.01
VTurbStep = 0.02
MetStep = 0.01
ParmStepArray = [(0.,0.,0.),(TeffStep,0.,0.), (-TeffStep,0.,0.),\
    (0.,LogGStep,0), (0.,-LogGStep,0.),\
    (0.,0.,VTurbStep), (0.,0.,-VTurbStep)]
# Statistically speaking, we need a minimum number of lines of a particular
# element/ion state to run our parameter analysis
MinNumStatLines = 20
# We accept a fewer number of lines for certain probability calculations
MinNumProbLines = 10

# Constants for baseline fitting
MinBaseline = 0.50
ContinuumBlanketFactor = 1.00
InvalidProbability = 7.E-77
BadProbability = 1.E-12

# Line variance array - The number of lines returned from an abundance measurement
# are limited by their variance from an average. Lines above each of the limits
# are eliminated until we have either a minimum number of lines, or all the measured
# lines vary from the mean by less than the minimum limit.
MinimumNumEQWLines = 5
# The code should exit when hitting a variance of 0.0
EQWAbVariances = [1.00, 0.75, 0.50, 0.33, 0.20]
# Phasing AbVariances out, since we are doing "Astrophysically determined" Loggf
# values. We instead assume that lines with (adjusted) abundances which vary by
# more than the values below, are blended, incorrect loggf values, overly
# affected by S/N (or Resolution) or NLTE effects, etc.
MaxBlendSTDev = 0.75
MinBlendSTDev = 0.50

# The Log(delta-Lambda/Lambda) limit for the high end of the 
# linear portion of the curve of growth. Noth: This actually
# includes a significant portion of the sqrt(ln) portion of the curve.
#LinearCOGLimit = -4.80
LinearCOGLimit = -4.50

# Measured line centers are never exact, so this designates a maximum 
# acceptable limit for how far a center can be from expected, and still 
# be accepted as a valid measurement
LambdaVarianceLimit = 0.050

# Variables for our Curve Of Growth (COG) plots.
cogPlotLimit = -4.70
evBinSize = 0.20
cogCorrFactor = 0.05    # How much to adjust the limit due to higher XP
BadAbundanceValue = -9.99

# Starting "guesses" for various atmospheric determinations
AtmGiantTGuess = 3800.0
AtmDwarfTGuess = 5800.0
AtmTGuessSigma = 200.0
AtmGiantGGuess = 2.0
AtmDwarfGGuess = 4.4
AtmGiantvGuess = 2.0
AtmDwarfvGuess = 1.0
AtmGiantMGuess = 0.0
AtmDwarfMGuess = -0.5

# MCMC parameters
NumWalkers = 12
#BSDims = 4  # Teff, LogG, Metal., Vel.
BSDims = 3  # Teff, LogG, Metal.

# Atmospheric parameter types
SpectroscopicParm = 'Spectroscopic'
PhotometricParm = 'Photometric'
StarfishParm = 'Starfish'
BayesianParm = 'Bayesian'

# Parsec Isochrone generator constants
MaxParsecMetal = 0.30
MaxParsecZ = 0.06



# Temporary MOOG file data for intermediate steps
MOOGScriptName = 'MOOGScript.temp'
MOOGTempModelName = 'MOOGTemp.m'
MOOGTempListName = 'MOOGTempList.lis'
MOOGTempLogName = 'MOOGTemp.log.tmp'
MOOGTempParName = 'MOOGTemp.par'
MOOGParFileHead = '''abfind
terminal       'xterm'
atmosphere     1
molecules      2
lines          1
freeform       0
flux/int       0
damping        1
plot           0
'''
MOOGBlendParHead = '''blends
terminal       'xterm'
atmosphere     1
molecules      2
lines          1
freeform       0
flux/int       0
damping        1
plot           0
'''
MOOGBlendParTail = '''blenlimits
    0.50  0.03    {0:2.1f}'''
    
MOOGParout1Head = 'standard_out   '
MOOGParout2Head = 'summary_out    '
MOOGParModelHead = 'model_in       '
MOOGParLinesHead = 'lines_in       '
MOOGLogFormat = '%8.3f %10.1f %10.5f %10.5f %10.1f %16.1f'
# MOOG format is: Wavelength; Ion; Ex. Pot.; LogGF; vdw factor; EQW; notes
MOOGListFormat = '{0:8.3f}  {1:9.1f}  {2:6.2f}  {3:7.2f}          {4:6.1f}          {5:8.1f}          {6}'


# Line list table output formats:
# Wavelength   Ion   X.P.   Loggf   # Meas.   Refs
LLTableFormat = '{0:8.3f}{1:6.1f}{2:7.3f}{3:6.2f}{4:3d} {5}'

# Generally, bad spectral regions (due to atmospheric lines, Strong (H/He) lines, UV noise, etc.)
BadSpectralRegions = [(2000., 4000.), (4099., 4102.5), (4858, 4864), (6559.4, 6566.5), (6274.5, 6314.), (6865., 6943.)]

# Due to a memory leak in pyspeckit, the maximum number of lines we can measure in one sitting is limited.
# We choose 75 lines here, which has worked with our low-memory (2GB) test system. If you have more
# working memory in your system, feel free to try larger line list files.
MaxLineListFileLength = 200
# We store our temporary line lists in a temporary directory, which we (plan to) delete when finished
TempListDirectory = '/TempLists'
TempListFileHead = 'templist_'
TempListFileExt = '.lis'

# Abundance correction weights and ranges if abs(correction) is > range, then
# weight by value.
SolarRefStarName = 'Sun'
DAOSolarRefStarName = 'Sund'
SolarCluster = 'NGC-0001'
SolarStarName = 'Sun-001'
SolarRefSpect = SolarCluster+' '+SolarStarName
DAOSolarStarName = 'Sun-002'
DAOSolarRefSpect = SolarCluster+' '+DAOSolarStarName
SolarRefParms = (5778.0, 4.44, 1.25, 0.00, 0.00)

ArcturusRefStarName = 'Arcturus'
ArcturusRefSpect = SolarCluster+' Arc-001'
# Parameters from Ramirez & Allende Prieto, 2011
ArcturusRefParams = (4286., 1.66, 1.74, -0.52, 1.08)
AldebaranRefStarName = 'Aldebaran'

# For a discussion of line weights and ranges, see Lum & Boesgaard 2018
AbWeightRanges = [(0.00,0.05), (0.05,0.10), (0.10,0.20), (0.20, 0.40), (0.40,1.00), (1.00,9.99)]
AbWeights = [10., 8., 6., 4., 2., 1.]
NoAbWeight = 0.
# How to determine how many lines to use:
# <=3 measured, use only the best quality
# 4-9, use all of rating 8.0 or better
# 10+, use all with Q>=6.0
abQRange = [3,10]

