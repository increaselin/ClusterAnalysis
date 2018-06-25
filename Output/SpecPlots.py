#! /usr/bin/env python
# 
# Module: SpecPlots.py
#
# Author: Mike Lum
#
# Description: Functions to create Postscript plots of spectra
#
# Contents: 
#   PlotSpec: Creates a plot of a single spectrum
#       Parms:  SpecFilename = Name (with relative or absolute path) 
#                   of file to plot
#               WLRange = Tuple with min and max wl to plot
#               OutFilename = Destination file (with relative or absolute path)
#               Title (optional) = Title to place on plot (default None)
#               
#   PlotMultiSpec: Creates a plot of multiple spectra, stacked one above 
#               another.
#       Parms:  SpecFilenames = List of names (with relative or absolute path)
#                   of file to plot
#               WLRange = Tuple with min and max wl to plot
#               OutFilename = Destination file (with relative or absolute path)
#               Title (optional) = Title to place on plot (default None)
#               Labels (optional) = List (same length as SpeFilenames list)
#                   of labels to place above each spectrum (default None)
#               
#
# Revision History:
#    Date        Vers.    Author        Description
#    3-20-2018   0.0a0    Lum        First checked in
#
# To Do:
#    
#

# Imports
from matplotlib import pyplot
from astropy.io import fits as pyfits
import numpy as np

from constants import constants as k
from utilities import utilities as u

def PlotSpec(SpecFilename, WLRange, OutFilename, Title=None):
    hdulist=pyfits.open(SpecFilename)
    hdr = hdulist[0].header
    data = hdulist[0].data
    WLs = np.arange(hdr[k.FitsKwdRVal], 
                    hdr[k.FitsKwdRVal]+(hdr[k.FitsKwdLen1DAxis])*hdr[k.FitsKwdDisp],
                    hdr[k.FitsKwdDisp])
    # np.where returns a tuple with a list of where the condition is true as
    # the first element. We want the first element of that list, thus the
    # double indexing
    WLStart = np.where(WLs==u.closest(WLRange[0], WLs))[0][0]
    WLEnd = np.where(WLs==u.closest(WLRange[1], WLs))[0][0]
    
    fig = pyplot.figure()
    ax = fig.gca()
    ax.plot(WLs[WLStart:WLEnd], data[WLStart:WLEnd])
    if Title != None:
        pyplot.title(Title)
    pyplot.savefig(OutFilename)
    pyplot.close()
    return


def PlotMultiSpec(SpecFilenames, WLRange, OutFilename, Title=None, Labels=None):

    numPlots = len(SpecFilenames)
    
    if Labels != None and len(Labels)==len(SpecFilenames):
        fnLabelPairs = zip(SpecFilenames, Labels)
    else:
        fnLabelPairs = zip(SpecFilenames, [None for i in range(len(SpecFilenames))])

    if numPlots > 5:
        linewidth = 0.25
        fontSize = 8
        XLabelCoord = 0.93
    else:
        linewidth = 0.5
        fontSize = 10
        XLabelCoord = 0.90
        
    plotHeight = min(10., 0.5*numPlots)
    plotWidth = 8.
    
    fig = pyplot.figure(figsize=(plotWidth, plotHeight))
    ax = fig.gca()

    heightCtr = 0.
    for fn, label in fnLabelPairs:
        hdulist=pyfits.open(fn)
        hdr = hdulist[0].header
        data = hdulist[0].data
        WLs = np.arange(hdr[k.FitsKwdRVal], 
                        hdr[k.FitsKwdRVal]+(hdr[k.FitsKwdLen1DAxis])*hdr[k.FitsKwdDisp],
                        hdr[k.FitsKwdDisp])
        # np.where returns a tuple with a list of where the condition is true as
        # the first element. We want the first element of that list, thus the
        # double indexing
        WLStart = np.where(WLs==u.closest(WLRange[0], WLs))[0][0]
        WLEnd = np.where(WLs==u.closest(WLRange[1], WLs))[0][0]
        if WLStart==WLEnd:
            continue
        ax.plot(WLs[WLStart:WLEnd], data[WLStart:WLEnd]+heightCtr, color='k',
                linewidth = linewidth)
        # The +0.2 is for a potential label
        heightCtr += max(data[WLStart:WLEnd])-min(data[WLStart:WLEnd])+0.2
        if label != None:
            ax.text(WLs[WLStart]+XLabelCoord*(WLs[WLEnd]-WLs[WLStart]), heightCtr+0.08, label, fontsize=fontSize)
    if Title != None:
        pyplot.title(Title)
    pyplot.savefig(OutFilename)
    pyplot.close()
    return

#PlotMultiSpec(['SpectraData/752Spectra/PLA_0300_temp.fits','SpectraData/752Spectra/PLA_0350_temp.fits','SpectraData/752Spectra/PLA_0356_temp.fits','SpectraData/752Spectra/PLA_0361_temp.fits','SpectraData/752Spectra/PLA_0475_temp.fits','SpectraData/752Spectra/PLA_0506_temp.fits','SpectraData/752Spectra/PLA_0520_temp.fits','SpectraData/752Spectra/PLA_0687_temp.fits','SpectraData/752Spectra/PLA_0699_temp.fits','SpectraData/752Spectra/PLA_0701_temp.fits'] , (6230, 6250), 'temp.ps', \
#Labels=['PLA-300','PLA-350','PLA-356','PLA-361',\
#        'PLA-475','PLA-506','PLA-520','PLA-687','PLA-699','PLA-701'])

