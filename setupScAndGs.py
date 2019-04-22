#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#Setup the spacecraft and ground station initial states
#SimpsonAerospace (c) 2019
from skyfield.api import Topos,  load
#------------------------------------------------------------------------------
def loadMostRecentTle(satgroup = 'noaa', satellite=33591):
    '''LOADMOSTRECENTTLE Day-old TLE and UA Ground Station State'''
    ##NOAA 15 and UAGroundStation state information
    rooturl = 'http://www.celestrak.com/NORAD/elements/'
    satgroup_tle = rooturl+satgroup+'.txt'

    #Download TLE and declare ground station
    satellites = load.tle(satgroup_tle)
    chosensat = satellites[satellite] 
    UAGroundStation = Topos('33.213271 N','87.544696 W')

    #Check age of TLE uploaded
    timetype = load.timescale()
    t = timetype.now()
    dt_tle = t - chosensat.epoch
    print('{:.2f} days away from the last epoch observed.'.format(dt_tle))
    #Update if older than a day
    if(abs(dt_tle)>1.00):
        satellites = load.tle(satgroup_tle,  reload=True)
        chosensat = satellites[satellite]
        dt_tle = t - chosensat.epoch
        print('Now {:.2f} days away from the last epoch observed.'.format(dt_tle))
    return chosensat,  UAGroundStation,  timetype
