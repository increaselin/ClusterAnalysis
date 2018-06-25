#! /usr/bin/env python
# 
# Module: SolarComps.py
#
# Author: Mike Lum
#
# Description: Functions to create plots to assess line quality based on
#               Solar spectra measurements
#
# Contents: List externally-callable functions here
#   Function xxx: (Description)
#
# Revision History:
#    Date        Vers.    Author        Description
#    1-11-18     1.0a0    Lum        First checked in
#
# To Do:
#    
#

# Imports
import numpy as np

from Abundances import ReferenceCorrections as RC
from constants import constants as k
from constants import ptoe_abs as PA
from Output import STRPlots as STP
from utilities import elements as el
from utilities import utilities as u

def MakeSolarCompTable(outFilename=k.TempAbOutFilename, 
                       headerLeft='', headerRight=''):
# Function creates a LaTeX table with elemental abundances as calculated from
# our Solar reference spectrum, based on a "quality" measurement of each line.
# Each column of the table represents a measurement using lines of a particular
# quality or better. The final column is represents all of the measured lines.
    corrections, lines, weights = RC.GetSolarCorrections()
    
    # We're going to create one column for each "quality" setting from the
    # returned weights list. Basically, each weight is a quality level.
    # The columnDict contains:
    # {weight:{ion:[ab,std,#lines],...},...}
    columnDict = {}
    for wt in k.AbWeights:
        columnDict[wt] = {}
        
    allIonList = sorted(lines.keys())
    for ion in allIonList:
        ionLines = lines[ion]
        # A single line will not show up as a list of one item. Rather, the 
        # entire list will be that one line's parameters...grrr.
        if not u.is_list(ionLines[0]):
            ionLines = [ionLines]
        for wt in k.AbWeights:
            wtLines = np.array([line for line in ionLines \
                        if line[0] in weights[ion] and weights[ion][line[0]] >= wt])
            if len(wtLines) > 0:
                asplundCorr = PA.ptoe[int(ion)-1][PA.abIdx]
                columnDict[wt][ion] = [np.mean(wtLines[:,5])+asplundCorr, 
                                       np.std(wtLines[:,5]), len(wtLines)]
            
    # Create a nice column title for each quality type:
    nameList = [r'$\Delta\leq{0:1.2f}$'.format(r[1]) for r in k.AbWeightRanges[:-1]]
    nameList.append('All lines')

    outfile = open(outFilename, 'w')
    
    latexHeader = STP.kLaTexHeader1 + '\\rhead{\\textbf{' + headerRight\
                  + '}}\n\\lhead{\\textbf{' + headerLeft\
                  + '}}\n\\begin{document}\n'
    outfile.write(latexHeader)

    outfile.write('\\begin{landscape}\n'\
                  + '\\hspace*{-5cm}\n'\
                  + '\\begin{tabular}{|l|l|l|l|l|l|l|l|l|l|')
    outfile.write('}\n\\multicolumn{1}{l}{Ion} & \\multicolumn{1}{l}{Asplund (2009)}')

    for name in nameList:
        outfile.write(' & \\multicolumn{{1}}{{l}}{{{0}}}'.format(name))
    outfile.write('\\\\\n\\hline\n')
    
    for ion in allIonList:
        outfile.write('{0} & {1:1.2f}$\pm${2:1.2f}'.format(el.getIonName(ion),\
                                          PA.ptoe[int(ion)-1][PA.abIdx],\
                                          PA.ptoe[int(ion)-1][PA.solarErrIdx]))
        
        for wt in k.AbWeights:
            if ion in columnDict[wt].keys() and columnDict[wt][ion][2]>0:
                outfile.write(' & {0:1.2f}$\pm${1:1.2f} ({2:d})'.\
                  format(columnDict[wt][ion][0]-PA.ptoe[int(ion)-1][PA.abIdx], 
                         columnDict[wt][ion][1], 
                         columnDict[wt][ion][2]))
            else:
                outfile.write(' & \\multicolumn{1}{c|}{---} ')

        outfile.write('\\\\\n\\hline\n')

    outfile.write('\\label{{tab:SolarAbs}}\n'
                  + '\\end{tabular}\n\\end{landscape}\n'
                  + '\\clearpage\n')

    outfile.write('\\end{document}\n')
    outfile.close()


