#! /usr/bin/env python
# 
# Module: LineTables.py
#
# Author: Mike Lum
#
# Description: Functions to create line list tables for publication
#
# Contents: 
#   Function xxx: (Description)
#
# Revision History:
#    Date        Vers.    Author        Description
#    12-21-2017  1.0f1    Lum        First checked in
#
# To Do:
#    
#

# Imports
import numpy as np

from constants import constants as k
from Databases import StarDB as SDB
from Databases import LineLookup as LL

# Constants


# Functions
def PrintLineTableForCluster(cluster = 'NGC-0752'):
    starNames = SDB.GetStarsForCluster(cluster)
    
    allLineDict = {}
    
    for starName in starNames:
        starLines = LL.getLinesForStar(cluster,starName)
        for line in starLines:
            if (line[0], line[1]) not in allLineDict.keys():
                allLineDict[(line[0], line[1])] = 1
            else:
                allLineDict[(line[0], line[1])] += 1
    
    for line in sorted(allLineDict.keys()):
        lookup, refs = LL.getLines(elements=[line[0]], \
            wavelengthRange=[(line[1]-k.LambdaVarianceLimit, \
                              line[1]+k.LambdaVarianceLimit)],\
                              dataFormat=None)
        if len(lookup) == 0:
            continue

        luline = lookup[0]

        print(k.LLTableFormat.format(line[1], line[0], luline[2], luline[3], allLineDict[(line[0],line[1])], luline[5]))

PrintLineTableForCluster()

