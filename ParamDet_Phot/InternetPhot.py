#! /usr/bin/env python
# 
# Module: InternetPhot.py
#
# Author: Mike Lum
#
# Description: Function to fetch stellar photometric data
#           from various web sources. Currently configured
#           for: Simbad
#
#
# Revision History:
#    Date        Vers.    Author        Description
#    2015-10-15  0.1a1    Lum        First testing and functionality
#
# To Do:
#    
#

# Imports
import datetime
import numpy as np
import sqlite3
#Py2->3
#import urllib2
from urllib.request import urlopen

from constants import constants as k
from constants import DBFormats as DB
from Databases import LineLookup as LL
from Models import MakeAtmModel as mk
from MOOGInterface import MOOGInterface as MI
from ParamDet_Phot import CasagrandeT as CT
from ParamDet_Phot import RamMelT as RMT
from ParamDet_Spec import SpecParamDet as SPD
from utilities import utilities as u

# Functions
def MakeSimbadName(starName):
    simbadName = 'Cl*+'+starName.replace(' ','+')
    return simbadName.replace('-','+')
    
def ParseSimbadHTML(simbadBytes):
    fluxes = []
    simbadPage = simbadBytes.decode('UTF-8')
    if 'Fluxes' not in simbadPage:
        return fluxes
        
    # Get the table of fluxes
    fluxTable = simbadPage.split('Fluxes')[1].split('</TABLE>')[0]
    fluxStrs = fluxTable.split('<TR>')[1:]
    
    for fluxEntry in fluxStrs:
    # Each flux entry wants to have the form:
    # [filter, flux, [error], quality, flags, Ref]
        # First, some filters contain a space. Like "g (AB)" - so check
        # if the second field is a number, if not, assume it is a
        # continuation of a filter
        fe = fluxEntry.split('<TT>')[1].split('</TT>')[0].split()
        idx = 1
        if not u.is_number(fe[idx]):
            fil = ' '.join(fe[0:2])
            idx += 1
        else:
            fil = fe[0]
            
        mag = float(fe[idx])
        idx+=1
        try:
            err = float(fe[idx][1:-1])
        except ValueError:
            err = 0.
            
        idx+=1
        if fe[idx] == '<SPAN':
            grade = fe[idx+4]
            idx +=5
        else:
            grade = '-'
        flags = ' '.join(fe[idx:])
                      
        fluxes.append([fil, mag, err, grade, flags])
    return fluxes
    
def GetSimbadPhot(simbadName):
    simbadURLHead = 'http://simbad.u-strasbg.fr/simbad/sim-id?ident='
    simbadURL = simbadURLHead+simbadName

# Py2->3
#    urlReq = urllib2.Request(simbadURL)
#    urlRes = urllib2.urlopen(urlReq)
    urlRes = urlopen(simbadURL)
        
    starPage = urlRes.read()
    fluxes = ParseSimbadHTML(starPage)
    
    returnDict = {}
    for flux in fluxes:
#        print('Filter: {0} - {1:2.3f}'.format(flux[0], float(flux[1])))
        returnDict[flux[0]] = float(flux[1])
    return returnDict
'''
def pareToMagicLines(lineList):
    magicList = [4630.124,4635.853,4661.979,4728.549,4733.591,4735.848,4741.532,4802.887,4946.385,4962.576,5022.236,5044.211,5067.160,5079.227,5126.200,5127.367,5131.476,5141.747,5145.102,5285.120,5288.532,5294.530,5295.300,5298.770,5373.704,5379.574,5389.490,5412.788,5432.960,5466.400,5473.910,5491.835,5529.160,5538.510,5539.270,5543.940,5554.900,5560.214,5633.950,5638.270,5652.310,5679.025,5717.835,5731.770,5741.856,5763.000,5809.220,5814.810,5849.690,5852.225,5909.978,5916.250,5956.706,6003.020,6008.570,6012.230,6015.250,6024.070,6027.060,6056.010,6173.340,6180.208,6187.995,6213.440,6219.286,6229.232,6240.652,6270.231,6271.283,6335.339,6344.150,6355.032,6358.690,6419.960,6430.855,6518.373,6575.020,6581.214,6592.922,6593.876,6625.027,6677.999,6699.136,6710.320,6739.524,6750.155,6752.716,6810.267,6837.020,6945.200,6971.936,6978.860,7112.173,7401.689,7723.200,7742.720,7746.600,7751.113,7941.090,8075.158,8239.130,8327.061,8468.418,8699.461, 4385.387,4413.600,4416.829,4508.289,4515.339,4520.227,4576.340,4582.834,4620.521,4666.193,4731.453,5197.575,5264.807,5425.245,6247.560,6369.462,6432.680,6456.389,7224.487,7711.726] 
    
    ironLines = [x for x in lineList if (x[0] == 26.0 or x[0] == 26.1) and (x[1] in magicList)]
    editedList = [x for x in lineList if x[0] != 26.0 and x[0] != 26.1]
    editedList.extend(ironLines)
    
    return editedList
    
# Command line calls for testing only
if __name__ == '__main__':
    
    isGiant = False
    metal = 0.10
    extinct = 0.03

    # We'll be running multiple passes, but only want to load models once. Otherwise,
    # we could just call mk.MakeModelFile.
    # To implement giant atmosphere calculation, we would need to use the 
    # "spherical" MARCS models. In case it's not obvious: d... = dwarf model file, 
    # g... = giant model file.
    u.errorOut('Reading plane-parallel model files:\n')
    dModelPath = k.DwarfModelPath
    dModelFiles = mk.findDataFiles(dModelPath)
    dKuruczAtms, dPradks = mk.LoadModels(dModelFiles)
    
    u.errorOut('Reading spherical model files:\n')
    gModelPath = k.GiantModelPath
    gModelFiles = mk.findDataFiles(gModelPath)
    gKuruczAtms, gPradks = mk.LoadModels(gModelFiles)
       
    rawData = []
    with sqlite3.connect(DB.ClusterDBName) as DBConnection:
        DBcursor = DBConnection.cursor()
'''
#        DBcursor.execute('''SELECT Cluster, StarID FROM {0};'''.format(DB.StarTableName))
'''        rawData.extend(DBcursor.fetchall())
        
    starData = [x for x in set([y[0]+' '+y[1] for y in rawData])]
    
    parms = []
    with open("mylogfileA.txt", 'w') as myLogFile:
        myLogFile.write("# New data run starting:{0}\n\n".format(datetime.datetime.now().strftime("%Y-%m-%d  %H:%M:%S")))

        for starName in starData:
        #['701','828','1107','391','552','921','864','413','791','964','699','575','361','1365','475','1089','993','859','858','300','506']:

            if starName[12:] not in ['300', '429', '575', '786', '859', '1017']:
                continue
            fluxes = GetSimbadPhot(MakeSimbadName(starName))

            temps, sigmas = CT.CasagrandeT(fluxes, metal, extinct)
                
            tWeights = 0.
            weightedT = 0.
            tAvg = np.mean(temps.values())
            tStd = np.std(temps.values())
            
            # Throw out Ts which differ by 1.5 sigma or more
            for t in temps.keys():
                if np.absolute(temps[t] - tAvg) < 1.5*tStd:
                    weightedT += temps[t]*1./sigmas[t]
                    tWeights += 1./sigmas[t]
            TGuess = weightedT/tWeights
            # Should do a better check here than just T for giants.
            if TGuess < 5200:
                # If Cassagrande says we're too cool, use Ramirez/Martinez
                # as a giant.
                isGiant = True
                temps, sigmas = RMT.RamMelGiantT(fluxes, metal, extinct)
                kuruczAtms = gKuruczAtms
                pradks = gPradks
                mass = 1.50
                logGGuess = 2.8
                tAvg = np.mean(temps.values())
                tStd = np.std(temps.values())
                for t in temps.keys():
                    if np.absolute(temps[t] - tAvg) < 1.5*tStd:
                        weightedT += temps[t]*1./sigmas[t]
                        tWeights += 1./sigmas[t]
                TGuess = weightedT/tWeights
            else:
                kuruczAtms = dKuruczAtms
                pradks = dPradks
                mass = 0.0
                logGGuess = 4.4
               
            # We really can't expect better than 5K resolution...
            TGuess = np.around(TGuess/5.,0)*5.
            
            VGuess = u.EdvardssonVTurb(TGuess, logGGuess)
            myLogFile.write('\n**** {0} ****'.format(starName))
            myLogFile.write("Initial guess:\n    T = {0:5.1f}\n    LogG = {1:1.2f}\n    VTurb = {2:1.2f}\n".format(TGuess, logGGuess, VGuess))
            
            star = starName.split(' ')

            theLines = LL.getLinesForStar(star[0], star[1], elements=[22.0, 22.1, 26.0, 26.1], wavelengthRange=[4000, 8700])
            
            theLines = pareToMagicLines(theLines)
            
            # We require certain minimums for good statistics...
            numIonLines = LL.NumLinesOfIon(theLines, 26.0)

            if numIonLines < 5:
                myLogFile.write('    Insufficient FeI lines ({0:2d}) to perform spectroscopic parameter determination. Skipping: {1}.\n'.format(numIonLines, starName))
                continue

# We have *MAGIC* lines!
#            if  numIonLines < k.MinNumStatLines:
#                if numIonLines < k.MinNumProbLines:
#                    myLogFile.write('    Insufficient FeI lines ({0:2d}) to perform spectroscopic parameter determination. Skipping: {1}.\n'.format(numIonLines, starName))
#                    continue
#                else:
#                     myLogFile.write('    Low number of FeI lines ({0:2d}) for spectroscopic parameter determination.\n'.format(numIonLines))
            else:
                myLogFile.write('    Total number of FeI lines: {0:2d}\n'.format(numIonLines))
                     
            numIonLines = LL.NumLinesOfIon(theLines, 26.1)
            if numIonLines < 1:
                myLogFile.write('    Insufficient FeII lines ({0:2d}) to perform spectroscopic parameter determination. Skipping: {1}.\n'.format(numIonLines, starName))
                continue
            else:
                myLogFile.write('    Total number of FeII lines: {0:2d}\n'.format(numIonLines))

# We have *MAGIC* lines!
#            if  numIonLines < k.MinNumStatLines/2:
#                if numIonLines < k.MinNumProbLines/2:
#                    myLogFile.write('    Insufficient FeII lines ({0:2d}) to perform spectroscopic parameter determination. Skipping: {1}.\n'.format(numIonLines, starName))
#                    continue
#                else:
#                     myLogFile.write('    Low number of FeII lines ({0:2d}) for spectroscopic parameter determination.\n'.format(numIonLines))
#            else:
#                myLogFile.write('    Total number of FeII lines: {0:2d}\n'.format(numIonLines))

            if  numIonLines < k.MinNumStatLines:
#                if numIonLines < k.MinNumProbLines:
                if numIonLines < 1:
                    myLogFile.write('    Insufficient TiI lines ({0:2d}) to perform spectroscopic parameter determination. Skipping: {1}.\n'.format(numIonLines, starName))
                    continue
                else:
                     myLogFile.write('    Low number of TiI lines ({0:2d}) for spectroscopic parameter determination.\n'.format(numIonLines))
            else:
                myLogFile.write('    Total number of TiI lines: {0:2d}\n'.format(numIonLines))

            numIonLines = LL.NumLinesOfIon(theLines, 22.1)
            if  numIonLines < k.MinNumStatLines/2:
                if numIonLines < 1:
#                if numIonLines < k.MinNumProbLines/2:
                    myLogFile.write('    Insufficient TiII lines ({0:2d}) to perform spectroscopic parameter determination. Skipping: {1}.\n'.format(numIonLines, starName))
                    continue
                else:
                     myLogFile.write('    Low number of TiII lines ({0:2d}) for spectroscopic parameter determination.\n'.format(numIonLines))
            else:
                myLogFile.write('    Total number of TiII lines: {0:2d}\n'.format(numIonLines))
            lineLogFilename = MI.MakeMOOGEQWLogFile(theLines)
            myLogFile.flush()
                    
            bestT, bestG, bestV, bestM = SPD.determineParms(lineLogFilename, (TGuess, logGGuess, VGuess), kuruczAtms, pradks, mass, clustMetal = -0.11, clustAge = 1400)
    
            myLogFile.write("Final guess:\n    T = {0:5.1f}\n    LogG = {1:1.2f}\n    VTurb = {2:1.2f}\n    Met: {3:1.2f}\n".format(bestT, bestG, bestV, bestM))
            
            u.clearFiles([lineLogFilename])
'''  
