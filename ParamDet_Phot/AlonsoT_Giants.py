#! /usr/bin/env python
# 
# Program: AlonsoT_Giants.py
#
# Author: Mike Lum
#
# Usage: ./AlonsoT_Giants.py [-vh]
#
# Description: Using the relations from Alonso et al. 1999, this function 
#    returns the Teff values for all of the passed photometry.
#
# Revision History:
#    Date        Vers.    Author        Description
#    Today        0.0a0    Lum        First checked in
#
# To Do:
#    
#

from collections import deque
import os

thisScriptName = 'AlonsoT_Giants.py'
interpreterName = 'python'

# Constants
# These are the textual names of the potential filters
kFilterDict = {'U':'U (Cousins)','B':'B (Cousins)','V':'V (Cousins)','R':'R (Cousins)','J':'J (2MASS)','H':'H (2MASS)','K':'K (2MASS)', 'L':'L\' (2MASS)','I':'I (Johnson)','b':'b (Stromgren)','y':'y (Stromgren)','u':'u (Stromgren)'}

# Table 2 from Alonso 1996 contains the polynomial values for the 
# calculation. 
#    Key:
#        U,B,V,R,I     : Johnson-Cousins
#        J,H,K,L        : 2MASS
#        b,y,u        : Stromgren
kTable2Dict = {\
    'UmV_a':{'Mrange':(-3.0,0.2),'Crange':(0.35,1.35),'a0':0.6388,'a1':0.4065,'a2':-0.1117,'a3':-0.002308,'a4':-0.07783,'a5':-0.01200,'sigma':164},\
    'UmV_b':{'Mrange':(-2.5,0.2),'Crange':(1.35,3.50),'a0':0.8323,'a1':0.09374,'a2':0.01184,'a3':0.02351,'a4':-0.1392,'a5':-0.01944,'sigma':80},\
    'BmV_a':{'Mrange':(-3.0,0.2),'Crange':(0.20, 0.80),'a0':0.5716,'a1':0.5404,'a2':-0.06126,'a3':-0.04862,'a4':0.01777,'a5':-0.007969,'sigma':167},\
    'BmV_b':{'Mrange':(-3.0,0.2),'Crange':(0.70, 1.90),'a0':0.6177,'a1':0.4354,'a2':-0.004025,'a3':0.05204,'a4':-0.1127,'a5':-0.01385,'sigma':96},\
    'VmR':{'Mrange':(-3.0,0.2),'Crange':(0.15, 1.70),'a0':0.4972,'a1':0.8841,'a2':-0.1904,'a3':-0.01197,'a4':0.01025,'a5':-0.0055,'sigma':150},\
#    'VmI':{'Mrange':(-3.0, 0.2),'Crange':(0.20, 2.90),'a0':0.5379,'a1':0.3981,'a2':0.04432,'a3':-0.02693,'a4':0.0,'a5':0.0,'sigma':125},\
    'RmI':{'Mrange':(-3.0, 0.2),'Crange':(0.15, 1.40),'a0':0.4974,'a1':1.345,'a2':-0.5008,'a3':-0.08134,'a4':0.03705,'a5':-0.006184,'sigma':150},\
    'VmK':{'Mrange':(-3.0, 0.2),'Crange':(0.20, 2.50),'a0':0.5558,'a1':0.2105,'a2':0.001981,'a3':-0.009965,'a4':0.01325,'a5':-0.002726,'sigma':40},\
    'VmK':{'Mrange':(-3.0, 0.2),'Crange':(2.00, 4.90),'a0':0.3770,'a1':0.3660,'a2':-0.00317,'a3':-0.003074,'a4':-0.002765,'a5':-0.002973,'sigma':25},\
    'JmH':{'Mrange':(-3.0, 0.2),'Crange':(0.00, 0.90),'a0':0.5977,'a1':1.015,'a2':-0.1020,'a3':-0.01029,'a4':0.03006,'a5':0.01013,'sigma':170},\
    'JmK':{'Mrange':(-3.0, 0.2),'Crange':(0.00, 1.10),'a0':0.5816,'a1':0.9134,'a2':-0.1443,'a3':0.000,'a4':0.000,'a5':0.000,'sigma':125},\
#    'VmL':{'Mrange':(-0.5, 0.2),'Crange':(0.40, 5.00),'a0':0.5641,'a1':0.1882,'a2':0.0189,'a3':-0.004651,'a4':0.000,'a5':0.000,'sigma':65},\
    'ImK':{'Mrange':(-3.0, 0.2),'Crange':(0.00, 1.90),'a0':0.5859,'a1':0.4846,'a2':-0.02457,'a3':0.000,'a4':0.000,'a5':0.000,'sigma':130},\
    'bmy_a':{'Mrange':(-3.0, 0.2),'Crange':(0.00, 0.55),'a0':0.5815,'a1':0.7263,'a2':0.06856,'a3':-0.06832,'a4':-0.01062,'a5':-0.01079,'sigma':110},\
    'bmy_b':{'Mrange':(-3.0, 0.2),'Crange':(0.50, 1.00),'a0':0.4399,'a1':1.209,'a2':-0.3541,'a3':0.08443,'a4':-0.1063,'a5':-0.01686,'sigma':70},\
    'umb':{'Mrange':(-3.0, 0.2),'Crange':(1.60, 4.00),'a0':0.5883,'a1':0.2008,'a2':-0.005931,'a3':0.005319,'a4':-0.1000,'a5':0.01542,'sigma':110}
    }
kNoData = -99.99


# Globals
verboseMode = True
passedMags = {}
deltas = {}
Teffs = {}

# Functions
def metInTableRange(met, tableKey):
    return (met<=kTable2Dict[tableKey]['Mrange'][1] and met>=kTable2Dict[tableKey]['Mrange'][0])

def colorInTableRange(color, tableKey):
    return (color<=kTable2Dict[tableKey]['Crange'][1] and color>=kTable2Dict[tableKey]['Crange'][0])

def printHelpText():
    print('Program: CasagrandeT.py\n')
    print('This function calculates the effective surface temperature of a')
    print('star, using the relationship from Cassagrande, et al. 2010\n')
    print('Usage: ./CasagrandeT.py <photometry> -M x.xx [-E x.xx] [-hv]\n')
    print('Options:')
    print('General Options:')
    print('\t-h: Print this help text.')
    print('\t-v: Use verbose progress and error messages')
    print('Photometry (at least two required):')
    print('\t-U xx.xx: Johnson-Cousins U band magnitude')
    print('\t-B xx.xx: Johnson-Cousins B band magnitude')
    print('\t-V xx.xx: Johnson-Cousins V band magnitude')
    print('\t-R xx.xx: Johnson-Cousins R band magnitude')
    print('\t-I xx.xx: Johnson-Cousins I band magnitude')
    print('\t-J xx.xx: 2MASS J band magnitude')
    print('\t-H xx.xx: 2MASS H band magnitude')
    print('\t-K xx.xx: 2MASS K band magnitude')
    print('\t-L xx.xx: 2MASS L\' band magnitude')
    print('\t-u xx.xx: Stromgren u band magnitude')
    print('\t-b xx.xx: Stromgren b band magnitude')
    print('\t-y xx.xx: Stromgren y band magnitude')
    print('Stellar parameters:')
    print('\t-M x.xx: (required) Metallicity of the star')
    print('\t-E x.xx: Extinction for a given band')


def AlonsoT(magDict, metal, extinct):
    try:
        deltas['UmV'] = magDict['U'] - magDict['V'] - extinct
    except (KeyError, IndexError):
        deltas['UmV'] = kNoData
    try:
        deltas['BmV'] = magDict['B'] - magDict['V'] - extinct
    except (KeyError, IndexError):
        deltas['BmV'] = kNoData
    try:
        deltas['VmR'] = magDict['V'] - magDict['R'] - extinct
    except (KeyError, IndexError):
        deltas['VmR'] = kNoData
    try:
        deltas['VmI'] = magDict['V'] - magDict['I'] - extinct
    except (KeyError, IndexError):
        deltas['VmI'] = kNoData
    try:
        deltas['RmI'] = magDict['R'] - magDict['I'] - extinct
    except (KeyError, IndexError):
        deltas['RmI'] = kNoData
    try:
        deltas['VmK'] = magDict['V'] - magDict['K'] - extinct
    except (KeyError, IndexError):
        deltas['VmK'] = kNoData
    try:
        deltas['JmH'] = magDict['J'] - magDict['H'] - extinct
    except (KeyError, IndexError):
        deltas['JmH'] = kNoData
    try:
        deltas['JmK'] = magDict['J'] - magDict['K'] - extinct
    except (KeyError, IndexError):
        deltas['JmK'] = kNoData
    try:
        deltas['VmL'] = magDict['V'] - magDict['L'] - extinct
    except (KeyError, IndexError):
        deltas['VmL'] = kNoData
    try:
        deltas['ImK'] = magDict['I'] - magDict['K'] - extinct
    except (KeyError, IndexError):
        deltas['ImK'] = kNoData
    try:
        deltas['bmy'] = magDict['b'] - magDict['y'] - extinct
    except (KeyError, IndexError):
        deltas['bmy'] = kNoData
    try:
        deltas['umb'] = magDict['u'] - magDict['b'] - extinct
    except (KeyError, IndexError):
        deltas['umb'] = kNoData
    
    goodKeys = [x for x in deltas.keys() if deltas[x] is not kNoData]
    goodPhot = {k:v for k,v in deltas.iteritems() if k in goodKeys}
    for thisKey in goodPhot.keys() :
        if verboseMode: print("Have photometry for: {0} = {1:3.2f}:".format(thisKey, goodPhot[thisKey]))
        paramList = [kTable2Dict[parmkey] for parmkey in kTable2Dict.keys() if (parmkey[:len(thisKey)]==thisKey and metInTableRange(metal, parmkey) and colorInTableRange(goodPhot[thisKey], parmkey))]

        if len(paramList) == 0: continue
        
        params=paramList[0]
        
        Teffs[thisKey] = 5040/(params['a0']+params['a1']*goodPhot[thisKey]+params['a2']*(goodPhot[thisKey]**2)+params['a3']*goodPhot[thisKey]*metal+params['a4']*metal+params['a5']*(metal**2))
    return Teffs

if __name__ == '__main__':

    # Parse your command line call
    temp = os.sys.argv
    argv = deque(temp)
    
    Metallicity = kNoData
    Extinction = 0.0
    
    while len(argv) > 0:
        flag = argv.popleft()
        if flag == '-h':
            printHelpText()
            exit()
        elif flag == '-v':
            verboseMode = True
            showIRAFMessages = 1
            print('Verbose Mode enabled.')
        elif flag == '-M':
            try:
                Metallicity = float(argv.popleft())
            except ValueError:
                Metallicity = kNoData
        elif flag == '-E':
            try:
                Extinction = float(argv.popleft())
            except ValueError:
                Extinction = 0.0
        elif flag[1:] in kFilterDict.keys():
            temp = argv.popleft()
            try: 
                tempMag = float(temp)
            except ValueError:
                continue
            passedMags[flag[1:]] = tempMag
            if verboseMode: print('{0} Magnitude: {1:3.2f}'.format(kFilterDict[flag[1:]], tempMag))
        elif str.find(flag,thisScriptName) != -1 :
            if verboseMode: print('Executing command.')
        elif str.find(flag,interpreterName) != -1 :
            if verboseMode: print('...')
        #else: ...
        # Ignore everything else, or not...
    
    if (Metallicity == kNoData):
        print("No metallicity value passed. Use the -M flag.")
        exit()
        
    temps = AlonsoT(passedMags, Metallicity, Extinction)
    
    for t in temps.keys():
        print("Teff = {0:5.1f} ({1})".format(temps[t], t))
    
    exit()
