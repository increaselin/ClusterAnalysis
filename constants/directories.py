#! /usr/bin/env python
# 
# Module: directories.py
#
# Author: Mike Lum
#
# Description: Set your local path and directory names in this file.
# EVERY USER WILL NEED TO SET THE APPROPRIATE DIRECTORIES FOR THEIR BUILD!
#
#
# Revision History:
#    Date        Vers.    Author        Description
#    07/12/16    1.0f0    Lum           First checked in
#
# To Do:
#    
#
mySourceHome = '/home/mikelum/Dropbox/CodeCloud/MyTools/ClusterAnalysis/'
useMARCSModels = False
# MARCS Model file path names
if useMARCSModels:
    DwarfModelPath = mySourceHome+'Models/MARCS_PlaneParallel'
    GiantModelPath = mySourceHome+'Models/MARCS_Spherical'
# ATLAS Models
else:
    DwarfModelPath = mySourceHome+'Models/ATLASModels'
    GiantModelPath = mySourceHome+'Models/ATLASModels'

# NLTE O database (from Takeda 2000)
NLTEODB = mySourceHome+'Databases/NLTEOxyCorrectionDB'

# Spectra data file directory tree
ToBeProcessed = mySourceHome+'SpectraData/Automation/ToDo'
ProcessedSpectra = mySourceHome+'SpectraData/Automation/Processed'
EQWMLogs = '/home/mikelum/Dropbox/DataCloud/EQWM_Logs'

# Output plots, in .ps form, are stored here
PlotDir = mySourceHome+'Plots/'
ParmPlotDir =  PlotDir+'SpectroscopicParms/'
XPAbPlotsDir = PlotDir+'XPvsAb/'
MCPlotsDir = PlotDir+'MCPlots/'
CoGPlotsDir = PlotDir+'CoGPlots/'
RelAbsDir = PlotDir+'RelElems/'
# Output Tables, in CDS (text) form here:
TableDir = mySourceHome+'Tables/'

# Line lists 
# (for import - possibly not needed, now that we have lines in the DB)
ListDirectory = '/home/mikelum/Dropbox/DataCloud/LineLists'

# Local output file for DB analysis - Not sure this is used anymore...
LineConflictFile = './LineConflicts.txt'

## MOOG Locations
#mspawnLoc = '/home/mikelum/Sources/mspawn72/mspawn72'
MOOGLoc = '/home/mikelum/Tools/Astro/moogjul2014/MOOGSILENT'

# DAOSpec Location
DAOSpecLoc = '/home/mikelum/Tools/Astro/DAOSpec'

