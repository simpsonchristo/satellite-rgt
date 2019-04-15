#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 22:11:12 2019

@author: willjpatton
"""

#eom
  fspherical = -(mu/(r**3))*X[0:1]
  Tecf2eci = ecf2ecisimple(t,ag0)
  x = Tecf2eci.transpose()*X[0:2] 
  
  latlonalt = ecf2spherical(x)
  Tsph2ecf = spherical2ecf(latlonalt)

  fnsr = 3.0*mu*(re**2/r**4)*J2*((3.0*(np.sin(latlonalt[0])*np.sin(latlonalt[0]))-1.0)/2.0) #radial direction
  fnsphi = -mu*(re**2/r**4)*J2*3.0*np.sin(latlonalt[0])*np.cos(latlonalt[0]) #latitude
  fnslam = 0.0 #longitude
  
  d2Xdt2[0:2] = X[3:5]
  d2Xdt2[3:5,0] = fspherical + Tecf2eci@Tsph2ecf@np.array([fnsr], [fnsphi], [fnslam])
  