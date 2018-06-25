#! /usr/bin/env python

# Astro- and 3rd-party Library Imports
from matplotlib import pyplot
import numpy as np
from scipy.stats import linregress

# Local Imports
from constants import directories as dirs
from constants import ptoe_abs as PA
from Output import STRPlots as STP
from utilities import elements as EL
from utilities import utilities as u

kGiants = ['PLA-0350', 'PLA-0356', 'PLA-0506', 'PLA-0687', 'PLA-1089', 'PLA-1172']

# Functions
def ColorLookup(starName):
    if starName in kGiants:
        return 'r'
    else:
        return 'b'


def PlotAllCombos():
    abTable=STP.GetAbTable()
    
    ionsToCombine=[6.0,7.0,8.0,11.0,12.0,13.0,14.0,20.0,22.0,26.0,63.1]
    allCombos = u.Cartesian(ionsToCombine,ionsToCombine)
#    allCombos = u.Cartesian(STP.kAllIonList, STP.kAllIonList)
    elemCombos = np.array([i for i in allCombos if i[0]<i[1]])

    SolarFe = PA.ptoe[25][PA.abIdx]

    for combo in elemCombos:
        Xs = []
        Ys = []
        XErrs = []
        YErrs = []
        colors = []
        for star in abTable.keys():
            starAbs = abTable[star]
            if combo[0] in starAbs.keys() and combo[1] in starAbs.keys():
                starFe = 0.75*starAbs[26.0][0] + 0.25*starAbs[26.1][0] - SolarFe
                
                SolarX = PA.ptoe[int(combo[0])-1][PA.abIdx]
                Xs.append(starAbs[combo[0]][0] - starFe - SolarX)
                XErrs.append(starAbs[combo[0]][1])
                
                SolarY = PA.ptoe[int(combo[1])-1][PA.abIdx]
                Ys.append(starAbs[combo[1]][0] - starFe - SolarY)
                YErrs.append(starAbs[combo[1]][1])

                colors.append(ColorLookup(star))
                
        slope, intercept, rVal, pVal, stderr = linregress(Xs, Ys)
        xElem = EL.getIonName(combo[0])
        yElem = EL.getIonName(combo[1])
        fig = pyplot.figure()
        ax = fig.gca()
        
        for pt in zip(Xs, Ys, XErrs, YErrs, colors):
            if pt[4] == 'm' or pt[4]=='c':
                shape='s'
            else:
                shape='o'
            eb = pyplot.errorbar(pt[0], pt[1], yerr=pt[3], xerr=pt[2],\
                linestyle='None', marker=shape, color=pt[4])
            eb[-1][0].set_linestyle('--')
            eb[-1][1].set_linestyle('--')

        minX = min(Xs)
        maxX = max(Xs)
        pyplot.plot([minX,maxX],[intercept+slope*minX, intercept+slope*maxX], 'k:')
        xMin = min(np.median(Xs)-0.5, minX)
        xMax = max(np.median(Xs)+0.5, maxX)
        ax.set_xlim(xMin, xMax)
        
        yMin = min(np.median(Ys)-0.5, min(Ys))
        yMax = max(np.median(Ys)+0.5, max(Ys))
        ax.set_ylim(yMin, yMax)
        
        pyplot.xlabel('['+xElem+'/Fe]')
        pyplot.ylabel('['+yElem+'/Fe]')
        ax.text(0.95, 0.01, 'Slope:{0:4.2f}, R-squared:{1:4.2f}, P-val:{2:4.2f}'.format(slope,rVal**2,pVal), fontsize=10, verticalalignment='bottom', horizontalalignment='right', transform=ax.transAxes,)
        pyplot.savefig(dirs.RelAbsDir+xElem+'_'+yElem+'.png')
        pyplot.close()

def PlotFeXx(clusterName='NGC-0752', starDataList=None, fileTag='', groupErrors=False,\
             ionList=STP.kAllIonList, filterBlends=False, referenceCorrect=False):


    # abTable = {starID:element/ion:[X/H, sigma, #, quality score]}
    abTable=STP.GetAbTable(clusterName=clusterName, starDataList=starDataList,\
             ions=ionList, filterBlends=filterBlends, referenceCorrect=referenceCorrect)

    SolarFe = PA.ptoe[25][PA.abIdx]

    for ion in ionList:
        Xs = []
        Ys = []
        XErrs = []
        YErrs = []
        colors = []
        allPts = []
        SolarY = float(PA.ptoe[int(ion)-1][PA.abIdx])
        for star in abTable.keys():
            starAbs = abTable[star]
            starFe = float(STP.GetStarFeH(starAbs))
            if ion in starAbs.keys():
                if not(u.is_number(starAbs[ion][0]) and \
                    u.is_number(starAbs[26.0][1]) and \
                    u.is_number(starAbs[ion][1])):
                    continue
                adj = SolarY
                if ion != 26.0:
                    adj+=starFe
                allPts.append([float(starFe), float(starAbs[ion][0])-adj, float(starAbs[26.0][1]), float(starAbs[ion][1]), ColorLookup(star)])
                
        yElem = EL.getIonName(ion)
        fig = pyplot.figure()
        axes = fig.gca()
        allPts = np.array(allPts)
#        allPts = np.array(zip(Xs, Ys, XErrs, YErrs, colors))

        gAbs = np.array([p[1] for p in allPts if p[4]=='r']).astype(np.float)
        gAvg = np.mean(gAbs)
        gStd = np.std(gAbs)
        gFeAbs = np.array([p[0] for p in allPts if p[4]=='r']).astype(np.float)
        gFeAvg = np.mean(gFeAbs)
        gFeStd = np.std(gFeAbs)
        dAbs = np.array([p[1] for p in allPts if p[4]=='b']).astype(np.float)
        dAvg = np.mean(dAbs)
        dStd = np.std(dAbs)
        dFeAbs = np.array([p[0] for p in allPts if p[4]=='b']).astype(np.float)
        dFeAvg = np.mean(dFeAbs)
        dFeStd = np.std(dFeAbs)
        
        axes.axhline(y=gAvg, ls='dashed', color='r', 
             label=r'Giant Avg. (${0:4.2f}\pm{1:4.2f}$)'.\
                     format(gAvg,gStd))
        axes.axhline(y=gAvg+gStd, ls='dotted', color='r')
        axes.axhline(y=gAvg-gStd, ls='dotted', color='r')
        
        axes.axhline(y=dAvg, ls='dashed', color='b', 
             label=r'Dwarf Avg. (${0:4.2f}\pm{1:4.2f}$)'.\
                     format(dAvg,dStd))
        axes.axhline(y=dAvg+dStd, ls='dotted', color='b')
        axes.axhline(y=dAvg-dStd, ls='dotted', color='b')

#        axes.axhline(y=SolarY, ls='dashed', color='g', 
#             label=r'Solar ({0:4.2f})'.format(SolarY))
        axes.axhline(y=0., ls='dashed', color='g', 
             label=r'Solar ({0:4.2f})'.format(SolarY))
        Xs = np.array(allPts[:,0]).astype(np.float)
        XErrs = np.array(allPts[:,2]).astype(np.float)
        xMin = max(np.median(Xs)-0.3, min(Xs)-0.1)
        xMax = min(np.median(Xs)+0.3, max(Xs)+0.1)
        axes.set_xlim(xMin, xMax)
        Ys = np.array(allPts[:,1]).astype(np.float)
        YErrs = np.array(allPts[:,3]).astype(np.float)
        yMin = min(np.median(Ys)-0.3, min(Ys))
        yMax = max(np.median(Ys)+0.3, max(Ys))
        axes.set_ylim(yMin, yMax)
        axes.set_xlabel('[Fe/H]', fontsize=18)
        if ion != 26.0:
            axes.set_ylabel('['+yElem+'/Fe]', fontsize=18)
        else:
            axes.set_ylabel('['+yElem+'/H]', fontsize=18)
        
        if groupErrors:
            axes.scatter(gFeAbs, gAbs, marker='o', color='r')
            eb = axes.errorbar(gFeAvg, gAvg, xerr=gFeStd,\
                linestyle='None', markersize=12.,  marker='H', color='r', ecolor='r', \
                label='Giant Mean', capsize=10., elinewidth=2, markerfacecolor='None')
            eb[-1][0].set_linestyle('--')
#            eb[-1][1].set_linestyle('--')
            eb[-1][0].set_linewidth(2.)
#            eb[-1][1].set_linewidth(2.)
            axes.scatter(dFeAbs, dAbs, marker='o', color='b')
            eb = axes.errorbar(dFeAvg, dAvg, xerr=dFeStd,\
                linestyle='None', markersize=12., marker='H', color='b', ecolor='b', \
                label='Dwarf Mean', capsize=10., elinewidth=2, markerfacecolor='None')
            eb[-1][0].set_linestyle('--')
#            eb[-1][1].set_linestyle('--')
            eb[-1][0].set_linewidth(2.)
#            eb[-1][1].set_linewidth(2.)
        else:
            eb = axes.errorbar(Xs, Ys, yerr=YErrs, xerr=XErrs,\
                linestyle='None', marker='o', ecolor=allPts[:,4])
            eb[-1][0].set_linestyle(':')
            eb[-1][1].set_linestyle(':')
            eb[-1][0].set_linewidth(0.5)
            eb[-1][1].set_linewidth(0.5)

        axes.legend()
        pyplot.savefig(dirs.RelAbsDir+'Fe_'+yElem+fileTag+'.png')
        pyplot.close()


