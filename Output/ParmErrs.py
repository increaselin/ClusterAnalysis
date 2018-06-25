#! /usr/bin/env python
# 
# Module: ParmErrs.py
#
# Author: Mike Lum
#
# Description: Function to create a LaTeX table with parameter variation errors
#              for the passed cluster (stars)
#
# Contents: List externally-callable functions here
#   Function xxx: (Description)
#
# Revision History:
#    Date        Vers.    Author        Description
#    2-8-18      0.0a0    Lum        First checked in
#
# To Do:
#    
#

# Imports
import numpy as np

from Abundances import ReferenceCorrections as RC
from Output import STRPlots as STP
from utilities import elements as el

# Constants
RowHeaders = [r'$T_{\rm{eff}}+100K$', r'$T_{\rm{eff}}-100K$',\
              r'$log(g)+0.02$', r'$log(g)-0.02$',\
              r'$\xi+0.2km/s$', r'$\xi-0.2km/s$']
parmDeltas = [('',100,0,0,0,0,0),('',-100,0,0,0,0,0), ('',0,0.02,0,0,0,0),\
              ('',0,-0.02,0,0,0,0),('',0,0,0.2,0,0,0), ('',0,0,-0.2,0,0,0)]

latexHead = '''\\documentclass[preprint]{aastex62}
\\usepackage{rotating}
\\begin{document}

\\clearpage
\\begin{sidewaystable}
\\caption{NGC 752 Atmospheric Errors}
\\begin{tabular}{l'''

latexFoot='''\\hline
\\label{tab:atmErrs}
\\end{tabular}
\\end{sidewaystable}

\\end{document}'''


def DataOut(data, fp=None):
    if fp is None:
        print(data)
    else:
        fp.write(data)


def GetGiantAbs(parmList, abDict):
# Note: could use an XOR, and a isGiant flag to combine the Dwarf and Giant
# versions of this function
    accumulator = {}
    giantList = [s[0] for s in parmList if RC.isGiantStar(s[1:])]
    for gName in giantList:
        if gName in abDict.keys():
            for ion in abDict[gName].keys():
                if ion in accumulator.keys():
                    accumulator[ion].append(abDict[gName][ion][0])
                else:
                    accumulator[ion] = [abDict[gName][ion][0]]
    retDict = {}
    for ion in accumulator.keys():
        retDict[ion] = np.mean([ab for ab in accumulator[ion]])
        
    return retDict
    
def GetDwarfAbs(parmList, abDict):
# Note: could use an XOR, and a isGiant flag to combine the Dwarf and Giant
# versions of this function
    accumulator = {}
    dwarfList = [s[0] for s in parmList if not RC.isGiantStar(s[1:])]
    for dName in dwarfList:
        if dName in abDict.keys():
            for ion in abDict[dName].keys():
                if ion in accumulator.keys():
                    accumulator[ion].append(abDict[dName][ion][0])
                else:
                    accumulator[ion] = [abDict[dName][ion][0]]
    retDict = {}
    for ion in accumulator.keys():
        retDict[ion] = np.mean([ab for ab in accumulator[ion]])
        
    return retDict
    

def MakeParmErrorTable(clusterName='NGC-0752', starDataList=None, 
               ions=STP.kAllIonList, filterBlends=False, 
               referenceCorrect=False, outFilename=None):
    # We will vary around either around the passed parameter list, or if none
    # is passed, then the parms in the DB
    
    if starDataList==None:
        starDataList = STP.GetAllStarParms(clusterName=clusterName)
    
    baseParms = starDataList
    baseAbs = STP.GetAbTable(clusterName=clusterName,
                             starDataList=baseParms,
                             ions=ions, filterBlends=filterBlends,
                             referenceCorrect=referenceCorrect)
    baseGs = GetGiantAbs(baseParms, baseAbs)
    baseDs = GetDwarfAbs(baseParms, baseAbs)
    
    altAbs = {}
    altGs = {}
    altDs = {}
    for parmDelta, label in zip(parmDeltas, RowHeaders):
        temp = [e[0]+e[1] for b in baseParms for e in zip(parmDelta, b)]
        parms = [tuple(temp[i:i+7]) for i in range(0,len(temp),7)]
        altAbs[label] = STP.GetAbTable(clusterName=clusterName,
                             starDataList=parms,
                             ions=ions, filterBlends=filterBlends,
                             referenceCorrect=referenceCorrect)
        altGs[label] = GetGiantAbs(parms, altAbs[label])
        altDs[label] = GetDwarfAbs(parms, altAbs[label])
        
        
    MakeErrorTable(baseGs, baseDs, altGs, altDs, outfile=outFilename)
    return
    
    
def MakeErrorTable(baseGs, baseDs, altGs, altDs, outfile=None):
    if outfile is not None:
        fp = open(outfile, 'a')
    else:
        fp=None
        
    DataOut(latexHead, fp=fp)
    
    dElems = np.array(sorted(baseDs.keys()))
    gElems = np.array(sorted(baseGs.keys()))
    # 19 element/ion colums will fit on a sideways table page
    # Note: Assume the giants and dwarfs _should_ have the same number of ions
    halfway = int(round(len(dElems)/2+0.1,0))   # F*ing Python floating points!
    DataOut(r'r'*halfway+'} \n', fp=fp)
    
    for ionSet, sType, baseAbs, abSet in [(dElems[:halfway], 'Dwarf', baseDs, altDs), \
                                 (gElems[:halfway], 'Giant', baseGs, altGs), \
                                 (dElems[halfway:], 'Dwarf', baseDs, altDs), \
                                 (gElems[halfway:], 'Giant', baseGs, altGs)]:
        # Do the column headers once for the Dwarfs, only
        if sType == 'Dwarf':
            columnHead = r'\multicolumn{1}{l}{Parameter} '
            for ion in ionSet:
                columnHead = columnHead+r'& \multicolumn{{1}}{{c}}{{{0}}} '.format(el.getIonName(ion))
            DataOut(columnHead+r'\\{0}\hline{0}'.format('\n'), fp=fp)
        
        DataOut( r'\multicolumn{{{0:2d}}}{{|c|}}{{{1}s}}\\{2}\hline{2}'.format(len(ionSet)+1, sType,'\n'), fp=fp)
        
        for parmVar in abSet.keys():
            rowData = r'\multicolumn{{1}}{{|l|}}{{{0}}}'.format(parmVar)
            for ion in ionSet:
                rowData = rowData+r' & \multicolumn{{1}}{{r|}}{{{0:1.2f}}}'.format(abSet[parmVar][ion]-baseAbs[ion])
            DataOut(rowData+r'\\{0}'.format('\n'), fp=fp)
        
        DataOut(r'\hline {1}\multicolumn{{1}}{{|l|}}{{{0} Totals}}'.format(sType,'\n'), fp=fp)
        for ion in ionSet:
            subtotal = np.sqrt(sum([max((abSet[RowHeaders[2*i]][ion]-baseAbs[ion])**2, \
                                      (abSet[RowHeaders[2*i+1]][ion]-baseAbs[ion])**2 \
                                      ) for i in [0,1,2]]))
            DataOut(r' & \multicolumn{{1}}{{r|}}{{{0:1.2f}}}'.format(subtotal), fp=fp)
        DataOut(r'\\{0}\hline{0}'.format('\n'), fp=fp)
    DataOut(latexFoot, fp=fp)
    if fp: fp.close()
    return
