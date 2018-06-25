#! /usr/bin/env python
# 
# Module: ParmCompPlots.py
#
# Author: Mike Lum
#
# Description: Functions to compare atmospheric parameters from our database
#               to those passed to the function.
#
# Contents: 
#   Function xxx: (Description)
#
# Revision History:
#    Date        Vers.    Author        Description
#    2018-03-21  1.0f1    Lum           First checked in
#
# To Do:
#    
#

# Imports
from matplotlib import pyplot
from matplotlib import pylab
import numpy as np

from Databases import StarDB as SDB

# Constants
# For my sanity:
sestito = [ [5982,4.50],[6151,4.50],[6130,4.50],[5425,4.50],[5970,4.50],\
[5677,4.50],[5531,4.50],[6325,4.50],[5753,4.50],[6060,4.50],[6231,4.50],\
[6223,4.50],[6102,4.50],[5611,4.50],[6282,4.50],[5791,4.50],[5791,4.50],\
[5603,4.50]]
maderak = [ [5537,4.55],[5517,4.55],[6312,4.37],[5124,4.62],[6102,4.43],\
[5656,4.53],[5517,4.55],[6264,4.38],[6208,4.40],[5728,4.52],[6071,4.44],\
[6177,4.41],[6169,4.41],[5593,4.54],[5728,4.52],[5764,4.51],[5764,4.51],\
[5621,4.54],[5586,4.54]]
castro =  [ [5308,4.59],[6266,4.36],[5688,4.53],[5929,4.47],[5888,4.47],\
[5847,4.50],[5794,4.50],[5727,4.52],[5714,4.52],[5902,4.47],[5970,4.46],\
[6053,4.44],[5888,4.48],[5482,4.57],[6165,4.40],[5445,4.58],[5662,4.54]]
padovaDwarf = [[4726, 4.69], [4955, 4.63], [5181, 4.60], [5311, 4.59], [5395, 4.58], [5595, 4.55], [5640, 4.54], [5781, 4.51], [5917, 4.48], [5950, 4.47], [5994, 4.46], [6107, 4.43], [6179, 4.41], [6252, 4.39], [6253, 4.39], [6383, 4.35], [6390, 4.35], [6519, 4.31]]

reddy = [[4850,2.5],[5050,2.85],[4850,2.6]]
topcu=[[4839,2.42],[4966,2.73],[4832,2.51],[5039,2.88],[4874,2.68]]
castroGiant=[[4634,2.68],[4677,2.75],[4731,2.87]]
padovaGiant = [[3908, 1.152], [4856, 2.706], [4815, 2.647], [4804, 2.626], [4797, 2.615], [4797, 2.612], [4802, 2.617], [4807, 2.626], [4813, 2.633], [4814, 2.633], [4813, 2.631], [4812, 2.627], [4824, 2.647], [4853, 2.694], [4916, 2.797], [4931, 2.821], [4940, 2.834], [4941, 2.837], [4936, 2.830], [4953, 2.861], [4973, 2.894], [4974, 2.895], [4976, 2.888], [4980, 2.862], [4976, 2.826], [4972, 2.812], [4967, 2.794], [4952, 2.758], [4937, 2.735], [4922, 2.702], [4903, 2.679], [4866, 2.632], [4825, 2.577], [4785, 2.521], [4748, 2.468], [4708, 2.408], [4674, 2.356]]
dartmouthGiant = [[4330,1.72],[4490,1.98],[4493,1.99],[4502,2.00],[4523,2.04],[4558,2.09],[4599,2.16],[4639,2.23],[4669,2.28],[4685,2.31],[4694,2.32],[4699,2.33],[4738,2.39],[4811,2.52],[4955,2.76],[4987,2.81],[5002,2.83],[5006,2.84],[4982,2.80],[5013,2.81],[5022,2.79],[5025,2.77],[5024,2.73],[5016,2.70],[5008,2.68],[4989,2.64],[4976,2.62],[4960,2.59],[4958,2.58],[4907,2.52],[4863,2.46],[4819,2.40],[4779,2.34],[4745,2.29],[4714,2.25],[4688,2.21],[4676,2.19],[4655,2.16],[4642,2.13],[4643,2.14]]

markers=['P','d','^','o','v','s','*']
colors=['firebrick','navy','darkgreen','darkviolet','cyan','gold']

# Functions
def CompPlotTeffLogg(stars=None, cluster='NGC-0752', 
                     isoPoints=None, isoLabel=None, 
                     compPointSets=None, outfilename=None, plotTitle=None):
# Creates a Teff-Logg plot of the stars in the passed cluster with the points
# passed (if any).
# Data Formats:
#               stars = list of star ids which can be used as lookup into
#                       the database.
#               cluster = Name of cluster. If no stars ids are passed, then
#                       plot all stars from this cluster.
#               isoPoints = Isochrone points to overplot with a dashed line.
#               isoLabel = A text label for the isochrone. If None, then the
#                       isochrone points are labeled as "Isochrone"
#               compPointSets = A dictionary of reference points to also plot.
#                       format: {RefName:[[Teff, Logg],...],...}
#               outfilename: If not None, then write the result to this file.
#                       Otherwise, just plot onscreen.
#               plotTitle: Title printed on plot. Generally not used for
#                       publications.

    if stars==None:
        stars = SDB.GetStarsForCluster(cluster)
    
    dbTeffs = []
    dbLoggs = []
    for sName in stars:
        (Teff, Logg, unused) = \
                                SDB.GetStarParms(sName, clusterName=cluster)
        dbTeffs.append(Teff)
        dbLoggs.append(Logg)
    
    minX = min(dbTeffs)
    maxX = max(dbTeffs)
    minY = min(dbLoggs)
    maxY = max(dbLoggs)

    params = {'legend.fontsize': 'large',
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
    pylab.rcParams.update(params)

    fig = pyplot.figure()
    axes = pyplot.gca()
    pyplot.scatter(dbTeffs, dbLoggs, figure=fig,\
                    label='This work', marker=markers[0], color=colors[0])
    mCount = 1
    cCount=1
    if compPointSets is not None:
        for setName in compPointSets.keys():
            teffs = [i[0] for i in compPointSets[setName]]
            loggs = [i[1] for i in compPointSets[setName]]
            pyplot.scatter(teffs, loggs, figure=fig, label=setName, 
                            marker=markers[mCount], color=colors[cCount])
            mCount += 1
            if mCount == len(markers): mCount = 0
            cCount += 1
            if cCount == len(colors): cCount = 0
            
            minX = min(teffs+[minX])
            maxX = max(teffs+[maxX])
            minY = min(loggs+[minY])
            maxY = max(loggs+[maxY])

    if isoPoints is not None:
        teffs = [i[0] for i in isoPoints]
        loggs = [i[1] for i in isoPoints]
        if isoLabel is not None:
            pyplot.plot(teffs, loggs, 'b:', figure=fig, label=isoLabel)
        else:
            pyplot.plot(teffs, loggs, 'b:', figure=fig, label='Isochrone')

    # We want the plot to look like an HR diagram, so reverse both axes:
    axes.set_xlim(maxX+50., minX-50.)
    axes.set_ylim(maxY+0.05, minY-0.05)

    axes.set_xlabel(r'$T_{\rm{eff}}$')
    axes.set_ylabel(r'log g')
    axes.legend()

    if plotTitle is not None:
        axes.set_title(plotTitle)

    if outfilename is not None:
        pyplot.savefig(outfilename, figure=fig)
    else:
        pyplot.show()
    pyplot.close()
    
    return

def plotGiantTeffLogg(cluster='NGC-0752', 
                     isoPoints=None, isoLabel=None, 
                     compPointSets=None, outfilename=None, plotTitle=None):
    stars = [s[4] for s in SDB.GetGiantInfoForCluster(clusterName = cluster)]
    print(stars)
    CompPlotTeffLogg(stars=stars, cluster=cluster, 
                     isoPoints=isoPoints, isoLabel=isoLabel, 
                     compPointSets=compPointSets, outfilename=outfilename,
                     plotTitle=plotTitle)
                     
                     
def plotDwarfTeffLogg(cluster='NGC-0752', 
                     isoPoints=None, isoLabel=None, 
                     compPointSets=None, outfilename=None, plotTitle=None):
    stars = [s[4] for s in SDB.GetDwarfInfoForCluster(clusterName = cluster)]

    CompPlotTeffLogg(stars=stars, cluster=cluster, 
                     isoPoints=isoPoints, isoLabel=isoLabel, 
                     compPointSets=compPointSets, outfilename=outfilename,
                     plotTitle=plotTitle)
