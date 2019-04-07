#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 12:01:31 2019

@author: willjpatton
"""

import numpy as np
import numpy.matlib as npml
from numpy import linalg as lin

a = 7702187 #meters
inc = np.radians(66.081)
ecc = 0.000951
raan = np.radians(251.911)
aop = np.radians(302.034)
M = np.radians(58.33)


#Some initial oribital elements
oe = np.array([a, ecc, inc, raan, aop, M])


def oe2rv(oe, f):
  #TWOBODYEPHEMERIS Determine position and velocity at time t by solving
  #Kepler's Equation
  #Solve Kepler's Equation for eccentric anomaly and determine true anomaly
  #to update the initial state to some state at time t. This propagation can
  #only be used in the two-body case. E0 and f are solved using the orbital
  #elements at the initial time, t0. 
  #Earth orbit
  
  mu = 3.9860044e14 #m^3/s^2

  a     = oe[0]             #semimajor axis, meters
  e     = oe[1]             #eccentricity
  i     = oe[2] #inclination, degrees
  Si = np.sin(i)
  Ci = np.cos(i)
  raan  = oe[3] #ascending node, degrees
  Sraan = np.sin(raan)
  Craan = np.cos(raan)
  w     = oe[4] #aop, degrees
  Sw = np.sin(w)
  Cw = np.cos(w)
  M     = oe[5] #mean anomaly, degrees
  #Sm = np.sin(M)
  #Cm = np.cos(M)
  #Q = np.array([[((Cw*Craan)-(Sw*Sraan*Ci)), ((-Sw*Craan)-(Cw*Sraan*Ci)], \
               #[((Cw*Sraan)+(Sw*Craan*Ci)), ((-Sw*Sraan)+(Cw*Craan*Ci))],  \
               #[(Sw*Si), (Cw*Si)]])
  Q = np.array([[(Cw*Craan)-(Sw*Sraan*Ci), ((Sw*Craan)-(Cw*Sraan*Ci))],[((Cw*Sraan)+(Sw*Craan*Ci)), ((-Sw*Sraan)+(Cw*Craan*Ci))], [(Sw*Si), (Cw*Si)]])  

  p = a*(1-e**2);#semilatus rectum
  h = np.sqrt(mu*p);#mag of angular momentum
  r = p/(1+(e*np.cos(f))) #mag of position vector, meters

  Vr = (h*e/p)*np.sin(f) #radial velocity, meters/sec
  Vtheta = h/r         #angular velocity, meters/sec

  Xdotstar = (Vr*np.cos(f)) - (Vtheta*np.sin(f))
  Ydotstar = (Vr*np.sin(f)) + (Vtheta*np.cos(f))
  
  X = Q@np.array([[r*np.cos(f), Xdotstar], [r*np.sin(f), Ydotstar]])

  r = X[:,0].transpose()
  v = X[:,1].transpose()
  X = np.array([[r],[v]])
  
  return X,r


def ecf2ecisimple(r,t,step):
    we = (2.0*np.pi/86164.0) #rad/sec, Earth avg rotational rate
    ag0 = 0.0
    t0 = 0.0
    t = t + step
    
    ag = we*(t-t0) + ag0
    
    x = r[0,0]*np.cos(ag)+r[1,0]*np.sin(ag)
    y = -(X[0,0]*np.sin(ag))+X[1,0]*np.cos(ag)
    z = X[2,0]
    
    Tecf2eci = np.array([[x],[y],[z]])
    
    return Tecf2eci
  
    
def ecf2spherical(x):
    #ECF2SPHERICAL Convert ITRF to geocentric coordinates
    #Convert position vector in earth fixed frame to latitude, longitude, and
    #altitude based coordinates or a spherical system. 
    # INPUTS:
    #   X0 - Given state in ITRF frame at time t. [meters, meters/seconds]
    # OUTPUTS:
    #  latlonalt - Geocentric latitude, longitude, and altitude. [rad,rad,meters]

    re = 6378137 #meters, spherical Earth radius
    rad = np.sqrt(x[0]^2 + x[1]^2 + x[2]^2)
    lat = np.arcsin(x[2]/rad)
    long = np.arctan2(x[1]/(rad*np.cos(lat)),x[0]/(rad*np.cos(lat)))
    latlonalt = np.array([[lat],[long],[rad]])
    
    return latlonalt
    
    
def spherical2ecf(latlonalt):
    lat = latlonalt[0]
    long = latlonalt[1]
    rad = latlonalt[2]
    
    x = rad*np.cos(lat)*np.cos(long)
    y = rad*np.cos(lat)*np.sin(long)
    z = rad*np.sin(lat)
    
    Tsph2ecf = np.array([[x],[y],[z]])
    
    return Tsph2ecf


def j2potential(t,X,step):
  #J2POTENTIAL Return 2nd derivative of motion for orbit
  #Considering only J2 perturbations and point mass force use the kinematics and return the acceleration.
  #INPUTS:
  #   X - Given state in J2000 frame at time t. [meters, meters/seconds]
  #   t - State will be predicted for each point. [seconds]
  # OUTPUTS:
  #   d2Xdt2 - acceleration from point mass potential and J2 [meters/sec^2]
  mu = 3.9860044e+14 #m^3/s^2, Earth gravitational parameter
  re = 6378137  #meters, spherical Earth radius
  ag0 = 0.0 #Greenwich angle
  J2 = 0.001082636
  C2_0 = -0.0004841695  #normalized value of harmonic coefficient
  r = lin.norm(X[0:2])
  
  #eom
  fspherical = -(mu/(r**3))*X[0:2]
  Tecf2eci = ecf2ecisimple(X,t,step)
  x = Tecf2eci.transpose()*X[0:2] 
  latlonalt = ecf2spherical(x)
  Tsph2ecf = spherical2ecf(latlonalt)

  fnsr = 3.0*mu*(re**2/r**4)*J2*((3.0*(np.sin(latlonalt[0])*np.sin(latlonalt[0]))-1.0)/2.0) #radial direction
  fnsphi = -mu*(re^2/r^4)*J2*3.0*np.sin(latlonalt[0])*np.cos(latlonalt[0]) #latitude
  fnslam = 0.0 #longitude

  d2Xdt2[0:2,0] = X[3:5]
  d2Xdt2[3:5,0] = fspherical + Tecf2eci@Tsph2ecf@np.array([fnsr], [fnsphi], [fnslam])

  return d2Xdt2
    

def rk4(t,f,step):
  #RK4 Runge-Kutta Numerical integrator
  t0 = t
  tf = t0 + step
  
  k1 = npml.empty((7,1))
  k2 = npml.empty((7,1))
  k3 = npml.empty((7,1))
  k4 = npml.empty((7,1))
  
  #---Begin integration
  k1 = j2potential(t0,f,step)
  k2 = j2potential(t0 + step / 2,f + (step / 2) * k1,step)
  k3 = j2potential(t0 + step / 2,f + (step / 2) * k2,step)
  k4 = j2potential(tf,f + (step / 2) * k3,step)
  
  Xf = f + (step / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
  
  return Xf    

ag0 = 0.0
t = 0.0

#Kepler's Equation
error = 100
E = M
while (error < 1e-6):
  E_old = E
  E = E - (E-ecc*np.sin(E)-M)/(1.0-np.cos(E))
  error = abs(E_old-E)

f = 2*np.arctan((np.sqrt((1+ecc)/(1-ecc))*np.tan(E/2)))
X = oe2rv(oe, f)

step = 0.0
while (step < 500):
    
    i = 0
    if not(i==0):
        Xf[i,:] = np.append(rk4(t,Xf[i-1,:],5))
    elif (i==0):
        Xf[i,:] = rk4(t,X,5)
    i = i + 1        
    step = step + 0.01
       
print(Xf)

