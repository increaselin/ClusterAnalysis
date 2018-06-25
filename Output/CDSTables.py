#! /usr/bin/env python
# 
# Module: CDSTables.py
#
# Author: Mike Lum
#
# Description: Outputs data in a CDS (Vizier) readable format
#
# Contents: List externally-callable functions here
#   Function xxx: (Description)
#
# Revision History:
#    Date        Vers.    Author        Description
#    2-21-18     1.0b0    Lum         First checked in
#

# Imports
import numpy as np

from Abundances import ReferenceCorrections as RC
from constants import constants as k
from constants import directories as d
from Databases import LineLookup as LL
from Databases import StarDB as SDB
from utilities import elements as el

# Constants
kEqualsDiv =\
'================================================================================'
kDashDiv =\
'--------------------------------------------------------------------------------'
kBbBTitle = 'Byte-by-byte Description of file: '
kBbBCols = '   Bytes Format  Units   Label    Explanations'
kNoLabel = '---'
kNoUnits = kNoLabel

# Functions
def BuildFormatString(DataFormat, DataLabels=[], DataUnits=[], DataExp=[],
                      MakeBbBDesc = False, FileName=''):
# Returns a Python format string, built from the passed CDS column format
# string. If "MakeBbBDesc" is True, an index block is also returned 
# (otherwise, the second returned variable is "None"
#
# DataFormat is ecpected to be a list of strings:
#       A##   : ASCII text string of length ##
#       I##   : A ## digit integer
#       F#.#  : Floating point with #.# digits.
    formatStr = ''
    if MakeBbBDesc:
        docBlock = kBbBTitle + FileName + '\n' +\
                   kDashDiv + '\n' + \
                   kBbBCols + '\n' + \
                   kDashDiv + '\n'
    else:
        docBlock = ''
        
    byteCounter = 0
    colCounter = 0
    for col in DataFormat:
        if len(DataLabels)> colCounter:
            label = DataLabels[colCounter]
        else:
            label = kNoLabel
            
        if len(DataUnits)> colCounter:
            units = DataUnits[colCounter]
        else:
            units = kNoUnits
            
        if len(DataExp)> colCounter:
            expStr = DataExp[colCounter]
        else:
            expStr = ''
            
        if col[0] == 'A':
            # We allow a zero-length string field to represent one with no
            # space between it and the next field (like a sign)
            fieldLen = int(col[1:])
            if fieldLen == 0:
                formatStr += '{{{0}:>1}}'.format(colCounter)
            else:
                formatStr += '{{{0}:>{1}}} '.format(colCounter, fieldLen)
        elif col[0] == 'I':
            fieldLen = int(col[1:])
            formatStr += '{{{0}:>{1}d}} '.format(colCounter, fieldLen)
        else: # Assume col[0] == 'F'
            intPart, decPart = col[1:].split('.')
            fieldLen = int(intPart)
            formatStr += '{{{0}:>{1:d}.{2}f}} '.format(\
                                            colCounter, int(intPart), decPart)
        colCounter += 1
        if fieldLen>1:
            docBlock += ' {0:>3}-{1:>3}  {2:<4}   {3:<6}  {4:<8} {5}\n'.format(
                     byteCounter+1, byteCounter+fieldLen, col, 
                     units, label, expStr)
        else:
            docBlock += '     {0:>3}  {1:<4}   {2:<6}  {3:<8} {4}\n'.format(
                     byteCounter+1, col[0]+'1', units, label, expStr)
        byteCounter += fieldLen+1
    docBlock += kEqualsDiv+'\n'
    return formatStr, docBlock
    
def MakeCDSTextFile(Data, DataFormat, DataLabels, DataUnits=[], DataExp=[], 
                    Notes=[], FileName=k.tempFilename, TableName='',\
                    BiblioInfo='', keywords=[], LongDesc='', ShortDesc=''):
    assert len(Data[0])==len(DataFormat)==len(DataLabels),\
        'Data columns must match formats and labels'
    
    fmt, dataBlock = BuildFormatString(DataFormat=Data, \
                      DataLabels=DataLabels, DataUnits=DataUnits, \
                      DataExp=DataExp, MakeBbBDesc = True, FileName=FileName)
    fp = open(FileName, 'w')
    for line in Data:
        fp.write(fmt.format(*line))
        fp.write('\n')
    fp.close()
    
# Should move these to Constants/constants.py
table1Format = ['I4','F4.1','F8.3','F5.1']   
table1Desc = ['Star identification in Platais (1991)',\
              'Element and Ionization state', 'Central wavelength of line',\
              'Measured Equivalent Width']
table1Labels=['Pla_ID','Ion','WL','EqW']
table1Units = ['---','---','0.1nm','0.1pm']
def MakeStarLineTable(starNames=[], clusterID='NGC-0752', ions=[], \
                      wls=[4500,8500], outFile=d.TableDir+'Table1.dat'):
    if len(starNames) == 0:
        starNames = SDB.GetStarsForCluster(clusterID)
    fp = open(outFile,'w')
    fmt, db = BuildFormatString(table1Format, DataLabels=table1Labels,
                                DataUnits=table1Units, DataExp=table1Desc, 
                                MakeBbBDesc = True)
    fp.write(db)
    for sn in starNames:
        lines = LL.getLinesForStar(clusterID, sn, elements=ions, 
                                   wavelengthRange=wls, dataFormat=None)
        starID = int(sn[4:])
        for l in lines:
            fp.write(fmt.format(starID,l[0],l[1],l[4])+'\n')
    fp.close()

table2Format = ['A2','A2','F8.3','F5.2','A0','F5.2','I2','A20']
table2Desc = ['Element','Ionization State','Central Wavelength',
              'Excitation Potential','Log(g_f) sign','Log(g_f)',
              'Quality indicator','References']
table2Labels =['El','Ion','WL','ExPot','','Logg_f','Qual','Refs.']
table2Units=['---','---','0.1nm','eV','---','---','---','---']
def MakeLineListTable(ions=[], wls=[], outFile=d.TableDir+'Table2.dat'):
    lines, refs = LL.getLines(elements=ions, wavelengthRange=wls,\
        dataFormat=None)
    unused1, unused2, weights = RC.GetSolarCorrections()

    fp = open(outFile,'w')
    fmt, db = BuildFormatString(table2Format, DataLabels=table2Labels,
                                DataUnits=table2Units, DataExp=table2Desc, 
                                MakeBbBDesc = True)
    fp.write(db)
    allIons = sorted(set(np.array(lines)[:,1]))
    lineData = []
    for i in allIons:
        ion = np.round(i,decimals=1)
        (elem, state, unused) = el.getIonState(el.getIonName(ion))
        ionLines = np.array(sorted([l for l in lines 
                            if np.round(l[1],decimals=1)==ion], 
                                key=lambda l:l[0]))
        for line in ionLines:
            wl = np.round(line[0],decimals=3)
            xp = line[2]
            gf = line[3]
            if gf<0:
                gfs = '-'
            else:
                gfs = ''
            gf = abs(gf)
            try:
                qu = int(weights[ion][wl])
            except KeyError:
                qu = 0
            rf = line[5]
            lineData.append([elem, state, wl, xp, gfs, gf, qu, rf])
    for line in lineData:
        fp.write(fmt.format(*line)+'\n')
    fp.close()

table3Format = []
table3Desc = []
table3Labels = []
table3Units = []
def MakeStarIDTable(clusterID='NGC-0752', outFile=d.TableDir+'Table3.dat'):
    stars = SDB.GetStarInfoForCluster(clusterName=clusterID)
    
    fp = open(outFile,'w')
    fmt, db = BuildFormatString(table3Format, DataLabels=table3Labels,
                                DataUnits=table3Units, DataExp=table3Desc, 
                                MakeBbBDesc = True)
    fp.write(db)

