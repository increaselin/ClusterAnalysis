#! /usr/bin/env python

# Imports
import glob
import numpy as np
from matplotlib import pyplot

from constants import ptoe_abs as PA
from constants import constants as k
from Abundances import ReferenceCorrections as RC
from Databases import LineLookup as LL
from Output import STRPlots as STP
from utilities import elements as el

def MakeTAbPlots(clusterName='NGC-0752', starDataList=None, ionList=STP.kAllIonList, fileTag='', useSubfigures=False, referenceCorrect=False):
# Create a bunch of .ps plots with the star's effective temperature vs. [X/Fe]
# for all ions in the optional passed list. One ion per plot

    if starDataList == None:
        starDataList = STP.GetAllStarParms(clusterName=clusterName)
        
    abTable = STP.GetAbTable(clusterName=clusterName, 
                            starDataList=starDataList, ions=ionList,
                            referenceCorrect=referenceCorrect)

    solarFeH = PA.ptoe[25][PA.abIdx]

    clusterFeH = \
        0.75*np.mean([abTable[star][26.0][0] for star in abTable.keys()])+\
        0.25*np.mean([abTable[star][26.1][0] for star in abTable.keys() \
            if 26.1 in abTable[star].keys()])

    for ion in ionList:
    
        solarIonH = PA.ptoe[int(ion)-1][PA.abIdx]

        fig = pyplot.figure()
        ax = pyplot.gca()
        
        redGPoints = []
        magGPoints = []
        bluDPoints = []
        cyaDPoints = []
        for star in starDataList:
            starName = star[0]
            if starName not in abTable.keys():
                continue
            starIons = sorted(abTable[starName].keys())
            if ion not in starIons:
                continue
                
            starParms = star[1:]
            
            # Get the Fe/H for this star. Note: We don't do Fe/Fe comps.
            if int(ion) != 26:
                # If both ions are measured, weight 75/25 for FeI/FeII
                if 26.0 in starIons and 26.1 in starIons:
                    starFe = 0.75*abTable[starName][26.0][0]\
                             + 0.25*abTable[starName][26.1][0]\
                             - solarFeH
                else:
                    starFe = abTable[starName][26.0][0]-solarFeH
            else:
                starFe = 0.
                
            if RC.isGiantStar(starParms):
#                if starParms[5] < 50:
#                    magGPoints.append([starParms[0], 
#                                   abTable[starName][ion][0]-solarIonH-starFe, 
#                                   abTable[starName][ion][1]])
#                else:
                redGPoints.append([starParms[0], 
                               abTable[starName][ion][0]-solarIonH-starFe, 
                               abTable[starName][ion][1]])
            else:
#                if starParms[5] < 25:
#                    cyaDPoints.append([starParms[0], 
#                                   abTable[starName][ion][0]-solarIonH-starFe, 
#                                   abTable[starName][ion][1]])
#                else:
                bluDPoints.append([starParms[0], 
                               abTable[starName][ion][0]-solarIonH-starFe, 
                               abTable[starName][ion][1]])
                                   

        dMean = np.mean([l[1] for l in cyaDPoints+bluDPoints])
        dStd = np.std([l[1] for l in cyaDPoints+bluDPoints])
        ax.axhline(y=dMean, ls='dashed', color='b', 
                   label=r'Dwarf Avg. (${0:4.2f}\pm{1:4.2f}$)'.\
                   format(dMean,dStd))
        ax.axhline(y=dMean-dStd, ls='dotted', color='b')
        ax.axhline(y=dMean+dStd, ls='dotted', color='b')
        
        gMean = np.mean([l[1] for l in redGPoints+magGPoints])
        gStd = np.std([l[1] for l in redGPoints+magGPoints])
        ax.axhline(y=gMean, ls='dashed', color='r',
                   label=r'Giant Avg. (${0:4.2f}\pm{1:4.2f}$)'.\
                   format(gMean,gStd))
        ax.axhline(y=gMean-gStd, ls='dotted', color='r')
        ax.axhline(y=gMean+gStd, ls='dotted', color='r')

        ax.axhline(y=0., ls='dashed', color='g', linewidth=1.0,\
                   label='Solar (${0:4.2f}$)'.\
                   format(PA.ptoe[int(ion)-1][PA.abIdx]))
                           
        pyplot.rcParams.update({'legend.fontsize':10})
        ax.legend()

        if dMean > 0. and gMean > 0.:
            ax.set_ylim([min([dMean-0.5,gMean-0.5]), max([dMean+0.5,gMean+0.5])])
        elif dMean > 0.:
            ax.set_ylim([dMean-0.5, dMean+0.5])
        elif gMean > 0:
            ax.set_ylim([gMean-0.5, gMean+0.5])

        Xs = [l[0] for l in bluDPoints]
        Ys = [l[1] for l in bluDPoints]
        Es = [l[2] for l in bluDPoints]
        eb=ax.errorbar(Xs, Ys, yerr=Es, fmt='o', color='b')
        if len(eb)>0 and len(eb[-1])>0: eb[-1][0].set_linestyle('--')
        
#        Xs = [l[0] for l in cyaDPoints]
#        Ys = [l[1] for l in cyaDPoints]
#        Es = [l[2] for l in cyaDPoints]
#        eb=ax.errorbar(Xs, Ys, yerr=Es, fmt='s', color='c')
#        if len(eb)>0 and len(eb[-1])>0: eb[-1][0].set_linestyle('--')

        Xs = [l[0] for l in redGPoints]
        Ys = [l[1] for l in redGPoints]
        Es = [l[2] for l in redGPoints]
        eb=ax.errorbar(Xs, Ys, yerr=Es, fmt='o', color='r')
        if len(eb)>0 and len(eb[-1])>0: eb[-1][0].set_linestyle('--')

#        Xs = [l[0] for l in magGPoints]
#        Ys = [l[1] for l in magGPoints]
#        Es = [l[2] for l in magGPoints]
#        eb=ax.errorbar(Xs, Ys, yerr=Es, fmt='s', color='m')
#        if len(eb)>0 and len(eb[-1])>0: eb[-1][0].set_linestyle('--')
        
        ionStr = el.getIonName(ion)
        
        if int(ion) != 26:
            ax.set_ylim((-0.5,0.5))
        else:
            ax.set_ylim((-0.3,0.3))
        
        ax.set_xlabel(r'$T_{eff}$ (K)')
        if int(ion) != 26:
            ax.set_ylabel(r'$[{0}/Fe]$'.format(ionStr))
        else:
            ax.set_ylabel(r'$[{0}]$'.format(ionStr))
            
        pyplot.savefig('{0}{1}{2}.png'.format(k.PlotDir+'Elems/', ionStr, fileTag), figure=fig)
    #        pyplot.show()
        pyplot.close(fig)

        
        

def MakeTvsEPPlots(fileTag=''):
    feiFiles = glob.glob("/home/mikelum/Dropbox/CodeCloud/MyTools/ClusterAnalysis/Tables/FeI_{0}.tab".format(fileTag))
    tiiFile = glob.glob("/home/mikelum/Dropbox/CodeCloud/MyTools/ClusterAnalysis/Tables/TiI_{0}.tab".format(fileTag))

    ion = 26.0
    ionStr = 'FeI'
    if len(feiFiles) > 0:
        fp = open(feiFiles[0],'r')
        lines = fp.readlines()
        fp.close()
            
        starNames = [l.replace(" ","") for l in lines[0].split(',')[1:-1]]
        wls = [float(l.split(',')[0]) for l in lines[3:-1]]
        abundSTR = [l.split(',')[1:-1] for l in lines[3:-1]]
        abunds = np.array([[float(b) for b in a] for a in abundSTR])
        
        for ab,starName in zip(abunds.T,starNames):
            lines = [(w,x) for (w,x) in zip(wls, ab) if x>0]
            if len(lines) == 0:
                continue
            epAbs = []
            for wl,ab in lines:
                info, refs = LL.getLines(elements=[26.0], wavelengthRange=[[wl-0.100, wl+0.100]], dataFormat='')
                if len(info) > 0:
                    epAbs.append([info[0][2], ab])
#            logg = parmLUT[starName][1
#            teff = parmLUT[starName][0]

            fig = pyplot.figure()
            axes = pyplot.gca()
        
            Xs = [p[0] for p in epAbs]
            Ys = [p[1] for p in epAbs ]
            avg = np.mean(Ys)
            std = np.std(Ys)
            axes.scatter(Xs, Ys, marker='.', color='b')
            axes.axhline(y=avg, ls='dashed', color='b')
            axes.axhline(y=avg-std, ls='dotted', color='b')
            axes.axhline(y=avg+std, ls='dotted', color='b')
            axes.set_ylim([avg-1.1, avg+1.1])
            axes.set_xlabel(r'[FeI/H]')
            axes.set_ylabel('Ex.Pot. (eV)')
            axes.set_title(starName)
            pyplot.savefig('{0}{1}_{2}{3}.png'.format('/home/mikelum/Dropbox/CodeCloud/MyTools/ClusterAnalysis/Plots/TvsEP/', starName, ionStr, fileTag), figure=fig)
#            pyplot.show()
        pyplot.close(fig)


# This was the original plot-maker function. Since we want a more general 
# function, we re-wrote it, above.

'''def MakeSPLATPlots(fileTag=''):
    abFiles = glob.glob("/home/mikelum/Dropbox/CodeCloud/MyTools/ClusterAnalysis/Tables/*{0}.tab".format(fileTag))

    ionStr = ''
    ion = 0.0
    for fn in abFiles:
        dAvg = -1.
        dStd = -1.
        gAvg = -1.
        gStd = -1.
        dPoints = []
        gPoints = []
        d2Points = []
        g2Points = []
        ionStr = fn.split('/')[-1][0:-(len(fileTag)+4)]
        fp = open(fn,'r')
        lines = fp.readlines()
        fp.close()
        ion = el.ionStrToFloat(ionStr)
            
        starNames = [l.replace(" ","") for l in lines[0].split(',')[1:-1]]
        abundSTR = [l.split(',')[1:-1] for l in lines[3:-1]]
        abunds = np.array([[float(b) for b in a] for a in abundSTR])
        
        idx = 0
        
        for ab in abunds.T:
            goodAbs = [x for x in ab if x>0]
            if len(goodAbs) == 0:
                idx += 1
                continue
                
            if idx not in range(len(starNames)):
                break
                
            logg = parmLUT[starNames[idx]][1]
            teff = parmLUT[starNames[idx]][0]
            if logg < 3.5:
                gPoints.append([teff,np.mean(goodAbs),np.std(goodAbs)])
            else:
                dPoints.append([teff,np.mean(goodAbs),np.std(goodAbs)])
#                # REALLY BAD!
#                if ion == 6.0:
#                    goodAbs = [x-0.4 for x in ab if x>0]
#                    gPoints.append([teff,np.mean(goodAbs),np.std(goodAbs)])
#                else:
#                    if ion-int(ion) == 0.1:
#                        goodAbs = [x-0.4 for x in ab if x>0]
#                        g2Points.append([teff,np.mean(goodAbs),np.std(goodAbs)])
#                    else:
#                        goodAbs = [x-0.4 for x in ab if x>0]
#                        g1Points.append([teff,np.mean(goodAbs),np.std(goodAbs)])
#            else:
#                if ion == 22.0:
#                    d1Points.append([teff,np.mean(goodAbs),np.std(goodAbs)])
#                else:
#                    if ion-int(ion) != 0.0:
#                        d2Points.append([teff,np.mean(goodAbs),np.std(goodAbs)])
#                    else:
#                        d1Points.append([teff,np.mean(goodAbs),np.std(goodAbs)])
            idx += 1

        if False:#ion in elemFilters.keys():
            minAb = elemFilters[ion][0]
            maxAb = elemFilters[ion][1]
        else:
            minAb = 0.0
            maxAb = 12.0

        fig = pyplot.figure()
        axes = pyplot.gca()
        
        if len(dPoints)>0:
            blueXs = [p[0] for p in dPoints if p[1]>minAb and p[1]<maxAb]
            blueYs = [p[1] for p in dPoints if p[1]>minAb and p[1]<maxAb]
#            greenXs = [p[0] for p in d2Points if p[1]>minAb and p[1]<maxAb]
#            greenYs = [p[1] for p in d2Points if p[1]>minAb and p[1]<maxAb]
#            dAvg = np.mean(blueYs+greenYs)
#            dStd = np.std(blueYs+greenYs)
            dAvg = np.mean(blueYs)
            dStd = np.std(blueYs)
            blueYErrs = [p[2] for p in dPoints if p[1]>minAb and p[1]<maxAb]
            axes.errorbar(blueXs, blueYs, yerr=blueYErrs, fmt='o', color='b')
#            greenYErrs = [p[2] for p in d2Points if p[1]>minAb and p[1]<maxAb]
#            axes.errorbar(greenXs, greenYs, yerr=greenYErrs, fmt='o', color='g')
            #axes.scatter(blueXs, blueYs, marker='o', color='b')
            #axes.scatter(greenXs, greenYs, marker='o', color='g')
            axes.axhline(y=dAvg, ls='dashed', color='b')
            axes.axhline(y=dAvg-dStd, ls='dotted', color='b')
            axes.axhline(y=dAvg+dStd, ls='dotted', color='b')
        if len(gPoints)>0:
            redXs = [p[0] for p in gPoints if p[1]>minAb and p[1]<maxAb]
            redYs = [p[1] for p in gPoints if p[1]>minAb and p[1]<maxAb]
#            orangeXs = [p[0] for p in g2Points if p[1]>minAb and p[1]<maxAb]
#            orangeYs = [p[1] for p in g2Points if p[1]>minAb and p[1]<maxAb]
#            gAvg = np.mean(redYs+orangeYs)
#            gStd = np.std(redYs+orangeYs)
            gAvg = np.mean(redYs)
            gStd = np.std(redYs)
            redYErrs = [p[2] for p in gPoints if p[1]>minAb and p[1]<maxAb]
            axes.errorbar(redXs, redYs, yerr=redYErrs, fmt='o', color='r')
#            marYErrs = [p[2] for p in g2Points if p[1]>minAb and p[1]<maxAb]
#            axes.errorbar(orangeXs, orangeYs, yerr=marYErrs, fmt='o', color='m')
#            axes.scatter(redXs, redYs, marker='o', color='r')
#            axes.scatter(orangeXs, orangeYs, marker='o', color='m')
            axes.axhline(y=gAvg, ls='dashed', color='r')
            axes.axhline(y=gAvg-gStd, ls='dotted', color='r')
            axes.axhline(y=gAvg+gStd, ls='dotted', color='r')
        if dAvg > 0. and gAvg > 0.:
            axes.set_ylim([min([dAvg,gAvg])-1.1, max([dAvg,gAvg])+1.1])
        elif dAvg > 0.:
            axes.set_ylim([dAvg-1.1, dAvg+1.1])
        else:
            axes.set_ylim([gAvg-1.1, gAvg+1.1])
        axes.set_xlabel(r'$T_{eff}$ (K)')
        axes.set_ylabel('[{0}/H]'.format(ionStr))
        pyplot.savefig('{0}{1}{2}.png'.format('/home/mikelum/Dropbox/CodeCloud/MyTools/ClusterAnalysis/Plots/Elems/', ionStr, fileTag), figure=fig)
    #        pyplot.show()
        pyplot.close(fig)
'''

