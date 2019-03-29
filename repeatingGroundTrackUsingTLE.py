#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Repeating Ground Track Orbits in High-Fidelity Geopotential"""
import numpy as np
import numpy.matlib as npml
from astropy import constants as const
#local libraries
import readtle
from downloadfileurl import download_file
"""Python 3.7
   Simpson Aerospace, Copyright 2019
   Christopher R. Simpson
   simpsonchristo@gmail.com"""
#------------------------------------------------------------------------------

rooturl = 'http://www.celestrak.com/NORAD/elements/'
satnogs_tle = download_file(rooturl+'satnogs.txt')
tdrss_tle   = download_file(rooturl+'tdrss.txt')
noaa_tle    = download_file(rooturl+'noaa.txt')


satrequest = np.array([43013])
oe = readtle.readtle(noaa_tle)
