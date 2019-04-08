#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  7 22:21:48 2019

@author: willjpatton
"""

#Start of process to compare postion after a certain time to desire position
#continution of RK4 and will be added to it after RK4 is working
#from print(Xf) 

import numpy as np

#convert Xf to lat and long and h
def ecf2ecisimple(Xf):
    we = (2.0*np.pi/86164.0) #rad/sec, Earth avg rotational rate
    
    ag = 0.0
    
    x = Xf[0,0]*np.cos(ag)+Xf[1,0]*np.sin(ag)
    y = -(Xf[0,0]*np.sin(ag))+X[1,0]*np.cos(ag)
    z = Xf[2,0]
    
    Tecf2eci = np.array([[x],[y],[z]])
    
    return Tecf2eci

def ecf2spherical(Tecf2eci):
    #ECF2SPHERICAL Convert ITRF to geocentric coordinates
    #Convert position vector in earth fixed frame to latitude, longitude, and
    #altitude based coordinates or a spherical system. 
    # INPUTS:
    #   X0 - Given state in ITRF frame at time t. [meters, meters/seconds]
    # OUTPUTS:
    #  latlonalt - Geocentric latitude, longitude, and altitude. [rad,rad,meters]

    re = 6378137 #meters, spherical Earth radius
    rad = np.sqrt(Tecf2eci[0]^2 + Tecf2eci[1]^2 + Tecf2eci[2]^2)
    lat = np.arcsin(Tecf2eci[2]/rad)
    long = np.arctan2(Tecf2eci[1]/(rad*np.cos(lat)),Tecf2eci[0]/(rad*np.cos(lat)))
    h = rad - re
    
    latlonalt = np.array([[lat],[long],[h]])
    
    return latlonalt


Tecf2eci = ecf2ecisimple(Xf)
latlonalt = ecf2spherical(Tecf2eci)

latdiff = latlonalt[0]-intiallat
longdiff = latlonalt[1]-intiallong
hdiff = latlonalt[2]-intialh

#if differences are too great, alter intial elements and rerun
#how much difference should we care about
#how are we going to change 





