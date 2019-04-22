#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#Potential function for J2 spherical harmonic contribution
#SimpsonAerospace (c) 2019
import numpy as np
#------------------------------------------------------------------------------
def j2potential(t, Xskyfield):
    '''J2POTENTIAL Return derivative of motion for orbit'''
    #Considering only J2 perturbations and point mass force use the kinematics
    #and return the acceleration.
    #INPUTS:
    #   Xskyfield - Given state in J2000 frame at time t. [meters, meters/seconds]
    mu = 3.9860044e+14  #m^3/s^2, Earth gravitational parameter
    re = 6378137 #meters, spherical Earth radius
    J2 = 0.001082636 
    
    #pass subpoint values to X
    subpoint = Xskyfield.subpoint()
    x,  y,  z = Xskyfield.position.km
    r = np.linalg.norm(np.array([x*1000,  y*1000,  z*1000]))
    vx,  vy,  vz = Xskyfield.velocity.km_per_s
    X = np.array([[subpoint.latitude], [subpoint.longitude], [subpoint.elevation], [vx*1000], [vy*1000],  [vz*1000]])
    #eom
    fspherical = -(mu/(r**3))*X[0:2]
    fnsr = 3*mu*(re**2/r**4)*J2*((3*(np.sin(X[0])*np.sin(X[0]))-1)/2) #radial direction
    fnsphi = -mu*(re**2/r**4)*J2*3*np.sin(X[0])*np.cos(X[0]) #latitude
    fnslam = 0 #longitude
    
    d2Xdt2 = np.array([[0], [0], [0], [0], [0], [0]])
    d2Xdt2[0:2,0] = X[4:6]
    d2Xdt2[3:5,0] = fspherical + np.array([fnsr], [fnsphi], [fnslam])
    return d2Xdt2
