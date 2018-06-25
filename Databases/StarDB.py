#! /usr/bin/env python
# 
# Module: StarDB.py
#
# Author: Mike Lum
#
# Description: Functions to look up cluster and individual star data in the 
#   Cluster DB.
#
# Contents: 
#   Function xxx: (Description)
#
# Revision History:
#    Date        Vers.    Author        Description
#    4/05/2018   1.0f1    Lum         Formatting and cleanup for github
#    7/14/2016   1.0b1    Lum         First checked in
#
# To Do:
#    
#

# Imports
import sqlite3

from constants import DBFormats as DB
from utilities import utilities as u

# Functions

def GetStarsForCluster(clusterName):
# Returns a list, containing the star IDs for all stars (in the DB) associated
# with the passed cluster name.
    with sqlite3.connect(DB.ClusterDBName) as DBConnection:
        DBcursor = DBConnection.cursor()
        
        DBcursor.execute(
        '''SELECT StarID FROM {0} WHERE Cluster is \"{1}\";'''.\
            format(DB.StarTableName, clusterName))
        returnedStars = DBcursor.fetchall()
    
    # sqlite returns tuples of unicode strings, so convert them into "normal"
    # strings for return values.
    return [str(x[0]) for x in returnedStars]
        
        
def GetStarInfo(starID, clusterName=None):
# Returns a tuple, containing all DB info
# (RA, DEC, Epoch, ClusterID, StarID, Teff, LogG, Vturb, Vmag, Bmag, s2n) of
# the passed star ID. in the case of multiple stars with the same ID
# (unlikely, but possible), the passed cluster name is used. If no cluster
# was passed, then the first star's parameters will be returned.
#

    with sqlite3.connect(DB.ClusterDBName) as DBConnection:
        DBcursor = DBConnection.cursor()
        
        if clusterName is None:
            DBcursor.execute(\
            '''SELECT * FROM {0} WHERE StarID is \"{1}\";'''.\
                format(DB.StarTableName, starID))
        else:
            DBcursor.execute(\
            '''SELECT * FROM {0} WHERE Cluster is \"{1}\" and StarID is \"{2}\";'''.\
                format(DB.StarTableName, clusterName, starID))
            
        starData = DBcursor.fetchall()

    return starData[0]
    
    
def GetStarParms(starID, clusterName=None):
# Returns a tuple, containing the (Teff, LogG, Vturb) of the passed
# star ID. in the case of multiple stars with the same ID (unlikely, but 
# possible), the passed cluster name is used. If no cluster was passed, then
# the first star's parameters will be returned.
#

    with sqlite3.connect(DB.ClusterDBName) as DBConnection:
        DBcursor = DBConnection.cursor()
        
        if clusterName is None:
            DBcursor.execute(\
            '''SELECT Teff, LogG, Vturb FROM {0} WHERE StarID is \"{1}\";'''.\
                format(DB.StarTableName, starID))
        else:
            DBcursor.execute(\
            '''SELECT Teff, LogG, Vturb FROM {0} WHERE Cluster is \"{1}\" and StarID is \"{2}\";'''.\
                format(DB.StarTableName, clusterName, starID))
            
        starData = DBcursor.fetchall()

    return starData[0]
    
    
def GetGiantInfoForCluster(clusterName=None):
# Looks up all stars in the passed cluster (or all stars in the DB for "None")
# and returns a list of giant stars with entries of the form:
# [RA, DEC, Epoch, clusterName, starID, Teff, LogG, vTurb, Vmag, Bmag, s2n]
    allStars = GetStarInfoForCluster(clusterName=clusterName)
    giantStars = []
    for star in allStars:
        if u.isGiant(star[5], star[6]):
            giantStars.append(star)
    return giantStars
    
    
def GetDwarfInfoForCluster(clusterName=None):
# As above, but for Dwarf stars.
    allStars = GetStarInfoForCluster(clusterName=clusterName)
    dwarfStars = []
    for star in allStars:
        if not u.isGiant(star[5], star[6]):
            dwarfStars.append(star)
    return dwarfStars
    
    
def GetStarInfoForCluster(clusterName=None):
# Looks up all stars in the passed cluster (or all stars in the DB for "None")
# and returns a list with entries of the form:
# [RA, DEC, Epoch, clusterName, starID, Teff, LogG, vTurb, Vmag, Bmag, s2n]
    with sqlite3.connect(DB.ClusterDBName) as DBConnection:
        DBcursor = DBConnection.cursor()
        
        if clusterName is None:
            DBcursor.execute(\
                '''SELECT * FROM {0} ORDER BY Cluster, StarID;'''.\
                    format(DB.StarTableName))
        else:
            DBcursor.execute(\
                '''SELECT * FROM {0} WHERE Cluster is \"{1}\" ORDER BY Cluster, StarID;'''.\
                    format(DB.StarTableName, clusterName))
            
        starData = DBcursor.fetchall()
    
    return starData


def SetStarAtmParms(starID, atmParms, clusterName=None):
# Set the parameters for the passed star in the passed cluster (optional).
# The parameter list order is:
# [Teff, LogG*, Vturb] *:See Note above
    
    # Quick validity check (May be able to do more, with some AI...)
    if len(atmParms) is not 3:
        return
        
    with sqlite3.connect(DB.ClusterDBName) as DBConnection:
        DBcursor = DBConnection.cursor()
        
        if clusterName is None:
            DBcursor.execute(\
            '''UPDATE {0} SET Teff={1}, LogG={2}, Vturb={3} WHERE StarID is \"{4}\";'''.\
            format(DB.StarTableName, atmParms[0], atmParms[1], \
                    atmParms[2], starID))
        else:
            DBcursor.execute(\
            '''UPDATE {0} SET Teff={1}, LogG={2}, Vturb={3} WHERE Cluster is \"{4}\" and StarID is \"{5}\";'''.\
            format(DB.StarTableName, atmParms[0], atmParms[1], \
                    atmParms[2], clusterName, starID))
        
        DBConnection.commit()
        
    return
