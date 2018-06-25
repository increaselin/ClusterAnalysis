#! /usr/bin/env python
# 
# Program: ProgramName.py
#
# Author: Mike Lum
#
# Usage: ./contAndStack.py [-vh]
#
# Description: This script searches the passed directory for reduced Keck Archive 
#               files (as designated by a "Flux-" prefix. It then flattens the 2-D
#               spectra to 1-D, using makee's combine function. The flattened spectra 
#               are then normalized to a continuum level, using IRAF's "continuum" function.
#               The normalized spectra for each star are then stacked into a single 
#               output file.
#
# Revision History:
#    Date        Vers.    Author        Description
#    Today        0.0a0    Lum        First checked in
#
# To Do:
#    
#

from collections import deque
import glob
from pyraf import iraf
import os
import subprocess as sub
from SpectraProcessing import dopplerCorrectLib as DCL
from SpectraProcessing import SpectraData as SD
from utilities import utilities as u

#thisScriptName = 'contAndStack.py'
#interpreterName = 'python'

# Constants
#kMakeeCombineCall = 'exec /home/mikelum/Tools/Astro/makee_5.2.4-sep08/bin/combine '
#kMakeeLinearCall = 'exec /home/mikelum/Tools/Astro/makee_5.2.4-sep08/bin/linear -1d '
kMakeeLinearCall = 'exec /home/mikelum/Tools/Astro/makee_5.2.4-sep08/bin/linear '
# Globals
gVerboseMode = False

# Functions
def test(directory='', fileTag=''):
    fluxFiles = sorted(glob.glob(directory+'/*_1_Flux.fits'))
    spectraCount = 0
    # Load the necessary IRAF packages
    iraf.onedspec()
    for thisFile in fluxFiles:
        # Resulting file names will use the self-reported object name
        # (if any) and an index based on the alphabetical order of
        # the original ("HI.xxxxx...") filename. Individual CCD files
        # (ie: HI.xxxx..._[1,2,3]_Flux.fits") will be left separate.
        objName = SD.GetObjectName(thisFile).strip().replace(' ','_')
        starName = '{0}{1}_{2:02d}'.format(fileTag, objName, spectraCount)
        spectraCount += 1
        
        # makee 'linearize' breaks a 2-d spectrum into it's component
        # orders, one order per file.
        makeeCall = kMakeeLinearCall + thisFile
        sub.call(makeeCall, shell=True)
        orderFiles = sorted(glob.glob(directory+'/*_Flux-??.fits'))
        for order in orderFiles:
            orderNo = int(order[-7:-5])
            if starName == '':
                fileRoot = order.split('/')[-1]
                orderFileName = directory+fileRoot
            else:
                orderFileName = directory+'/{0}_{1:02d}.fits'.format(starName, orderNo)
            if orderNo>6:
                u.clearFile(order)
            else:
                os.system('mv {0} {1}'.format(order, orderFileName))
            # Continuum the orders
#            if int(orderNo)<5:
#                iraf.continuum(input=order, output=orderFileName, functio='chebyshev', interac='no', order=5, low_rej=3.0, high_re=3.5, niterat=10)
#            u.clearFile(order)
#        contFileName = directory+'/aC_'+starName+'.fits'
        # ...and stack them back up.
#        iraf.onedspec.scopy(input=directory+'/temp_*', output=contFileName, format='multispec', renumber="yes")
#        filesToClean = glob.glob(directory+'/temp_*')
#        for fn in filesToClean:
#            u.clearFile(fn)
    return
    
    
    
def contAndStack(directory='', fileTag=''):
# Function goes through all the files in the passed directory, breaks them
# into their component orders, and continuum-normalizes each order. It 
# then re-combines the continuum-ed orders into a single 1-D spectrum.
# It then performs a user-assisted Doppler correction.
    fluxFiles = sorted(glob.glob(directory+'/*_Flux.fits'))
    spectraCount = 0
    # Load the necessary IRAF packages
    iraf.onedspec()
    
    # We're going to keep track of the files for each object separately.
    # Note: With the potential for multiple exposures, we may have multiple
    # images to weight and stack for a given object.
    fileDict = {}
    
    for thisFile in fluxFiles:
        # Reason #93 as to why IRAF sucks...80 character limit for
        # file paths and names...
        if len(thisFile) >80:
            tempName = directory+'/a_{0}.fits'.format(spectraCount)
            if len(tempName)>80:
                # make the user rename their path...
                print('File path/name too long: {0}'.format(thisFile))
                print('Skipping this file.')
                continue
            os.rename(thisFile, tempName)
            thisFile = tempName

        # Resulting file names will use the self-reported object name
        # (if any) and an index based on the alphabetical order of
        # the original ("HI.xxxxx...") filename. Individual CCD files
        # (ie: HI.xxxx..._[1,2,3]_Flux.fits") will be left separate.
        objName = SD.GetObjectName(thisFile).strip().replace(' ','_')
        if objName not in fileDict.keys():
            fileDict[objName] = []
        # Using a "spectraCount" ensures that if multiple exposures of
        # the same object exist, they are given unique filenames for the
        # combining process
        starName = '{0}{1}_{2:02d}'.format(fileTag, objName, spectraCount)
        spectraCount += 1
        
        # makee 'linearize' breaks a 2-d spectrum into it's component
        # orders, one order per file.
        makeeCall = kMakeeLinearCall + thisFile
        sub.call(makeeCall, shell=True)
        orderFiles = sorted(glob.glob(directory+'/*_Flux-??.fits'))
        for order in orderFiles:
            orderNo = order[-7:-5]
            if starName == '':
                fileRoot = order.split('/')[-1]
                orderFileName = directory+'/temp_'+fileRoot
            else:
                orderFileName = directory+'/temp_'+starName+'_'+orderNo+'.fits'
            # Continuum the orders
            try:
                iraf.continuum(input=order, output=orderFileName,
                    functio='chebyshev', interac='no', order=5, low_rej=3.0, 
                    high_re=3.5, niterat=10)
            except:
                print('Error in continuum {0}. Skipping.'.format(order))
                u.clearFile(order)
                continue
            u.clearFile(order)
            fileDict[objName].append(orderFileName)
    completedList=['H027771','H027859','H028099','H028205','H028344','H028462',\
    'H028992','Hy_vB59=HD28034','Hy_vB64=HD28099','Hy_vB65=HD28205',\
    'Hy_vB73=HD28344','Hy_vB92=HD28805','Hy_vB93=HD28878','Hy_vB97=HD28992',\
    'vA_294','vB_173','vB_180','vB_183','vB__46','vB__49','vB__64','vB__65',\
    'vB__73','vB__92','vB__93','vB__97','RHy_281','VA294','28344','HD_028344'\
    'HD28394','08-U20187']
    
    for objName in fileDict.keys():
        if len(fileDict[objName]) == 0:
            continue
        if objName in completedList:
            print('Skipping {0}'.format(objName))
            continue
        # Create a temporary input file list:
        fpInput = open('input.dat','w')
        fpWeight = open('weight.dat','w')
        print('---- {0} ----'.format(objName))
        for fn in fileDict[objName]:
            fpInput.write(fn+'\n')
            sn = SD.GetSpectraSN(fn)
            fpWeight.write('{0:6.1f}\n'.format(sn))
            print('    {0:40}: {1:5.1f}'.\
                format(fn.split('/')[-1], sn))
        fpInput.close()
        fpWeight.close()
        contFileName = directory+'/wc_'+objName+'.fits'
        # Stack all the images/orders for this object.
        iraf.onedspec.scombine(input='@input.dat', output=contFileName,\
            weight='@weight.dat', combine='median', \
            lthreshold=0.05, hthreshold=2.0, group="all")
#        iraf.stsdas.hst_calib.tomultispec(input=contFileName, output='ms_'+contFileName)
        u.clearFile('input.dat')
        u.clearFile('weight.dat')
        for fn in fileDict[objName]:
            u.clearFile(fn)

        # For now, the automated Doppler correction is good only for very small
        # corrections (ie: <5km/s), so use the manual version.
        dopCorVelocity, dopWLs = DCL.dopCor(contFileName, interactive=True)
#        print('Doppler correction = {0:2.3f}km/s'.format(dopCorVelocity))
        dopCorFilename = directory+'/D'+contFileName.split('/')[-1][1:]
        iraf.dopcor(input=contFileName, output=dopCorFilename, redshift=dopCorVelocity, isveloc='yes', verbose='yes')
        u.clearFile(contFileName)
        
    return
        
        
        
        
        
# Stuff for the command-line version (no longer used)
'''
def printHelpText():
    print('Program: contAndStack.py\n')
    print('Insert your \"help text\" here.\n')
    # Be sure to list your full list of flags and options here
    print('Usage: ./contAndStack.py [-i <file name> -s <file tag> -hv]\n')
    print('Options:')
    print('-i: Search directory (default ./)')
    print('-s: Optional tag for output files')
    print('-h: Print this help text.')
    print('-v: Use verbose progress and error messages')


if __name__ == '__main__':

    # Parse your command line call
    temp = os.sys.argv
    argv = deque(temp)
    
    while len(argv) > 0:
        flag = argv.popleft()
        if flag == '-h':
            printHelpText()
            exit()
        elif flag == '-v':
            gVerboseMode = True
            print('Verbose Mode enabled.')
        elif str.find(flag,thisScriptName) != -1 :
            if gVerboseMode: print('Executing command.')
        elif str.find(flag,interpreterName) != -1 :
            if gVerboseMode: print('...')
        elif flag == '-i':
            # Input directory name
            dirName = argv.popleft()
        elif flag == '-s':
            # Optional star name for output file
            starName = argv.popleft()
    
    # Call your program here
    contAndStack(directory=dirName, fileTag=starName)
    exit()'''
