#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#Generate an array of times for the ephemeris
#SimpsonAerospace (c) 2019
from datetime import datetime,  timedelta,  tzinfo
#------------------------------------------------------------------------------
def timeEphem(timetype,  days=1., dt=30.):
    '''TIMEEPHEM Points in time to predict state for ephemeris'''
    ##generate time for ephemeris
    ZERO = timedelta(0)

    # A UTC class.
    class UTC(tzinfo):
        """UTC"""
        def utcoffset(self, dt):
            return ZERO
        def tzname(self, dt):
            return "UTC"
        def dst(self, dt):
            return ZERO
    utc = UTC()

    t0 = datetime.utcnow().replace(tzinfo=utc)
    stepsize = timedelta(seconds=dt)
    stoploop = 86400/30 #30 sec timestep in terms of day
    trange = []
    for i in range(stoploop):
        trange.append(t0 + stepsize)
    tephem = timetype.utc(trange)
    return tephem,  t0
