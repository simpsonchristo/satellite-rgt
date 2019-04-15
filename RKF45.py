#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 12:01:31 2019

@author: willjpatton
"""

import numpy as np
import numpy.matlib as npml
from numpy import linalg as lin
import matplotlib.pyplot as plt

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
  i     = oe[2]             #inclination, degrees
  Si = np.sin(i)
  Ci = np.cos(i)
  raan  = oe[3]             #ascending node, degrees
  Sraan = np.sin(raan)
  Craan = np.cos(raan)
  w     = oe[4]             #aop, degrees
  Sw = np.sin(w)
  Cw = np.cos(w)
  M     = oe[5]             #mean anomaly, degrees
  #Sm = np.sin(M)
  #Cm = np.cos(M)
  #Q = np.array([[((Cw*Craan)-(Sw*Sraan*Ci)), ((-Sw*Craan)-(Cw*Sraan*Ci)], \
               #[((Cw*Sraan)+(Sw*Craan*Ci)), ((-Sw*Sraan)+(Cw*Craan*Ci))],  \
               #[(Sw*Si), (Cw*Si)]])
  Q = np.array([[(Cw*Craan)-(Sw*Sraan*Ci), ((Sw*Craan)-(Cw*Sraan*Ci))],[((Cw*Sraan)+(Sw*Craan*Ci)), ((-Sw*Sraan)+(Cw*Craan*Ci))], [(Sw*Si), (Cw*Si)]])  

  p = a*(1-e**2)          #semilatus rectum
  h = np.sqrt(mu*p)       #mag of angular momentum
  r = p/(1+(e*np.cos(f))) #mag of position vector, meters

  Vr = (h*e/p)*np.sin(f)  #radial velocity, meters/sec
  Vtheta = h/r            #angular velocity, meters/sec

  Xdotstar = (Vr*np.cos(f)) - (Vtheta*np.sin(f))
  Ydotstar = (Vr*np.sin(f)) + (Vtheta*np.cos(f))
  
  X = Q@np.array([[r*np.cos(f), Xdotstar], [r*np.sin(f), Ydotstar]])

  r = X[:,0].transpose()
  v = X[:,1].transpose()
  X = np.array([[r],[v]])
  #rj = X[1,0].transpose()
  #rk = X[2,0].transpose()
  #vi = X[0,1].transpose()
  #vj = X[1,1].transpose()
  #vk = X[2,1].transpose()
  #X = np.array([[ri,rj,rk],[vi,vj,vk]])
  
  return X

"""
def eci2ecfsimple(X,t,h):
    we = (2.0*np.pi/86164.0) #rad/sec, Earth avg rotational rate
    ag0 = 0.0
    t = t + h
    ag = we*t + ag0
    
    x = X[0,0]*np.cos(ag)+X[0,1]*np.sin(ag)
    y = -(X[0,0]*np.sin(ag))+X[0,1]*np.cos(ag)
    z = X[0,2]
    
    Tecf2eci = np.array([x,y,z])
    
    return Teci2ecf
"""

def ecf2ecisimple(t,ag0):
    #ECF2ECISIMPLE Convert from ecf to eci using Greenwich angle
    #Simplified conversion by rotating solely Z axis where Z and z axis are
    #coincident because of the angular velocity vector direction is constant.
    we = (2*np.pi/86164)   #rad/sec, Earth avg rotational rate
    ag = we*t + ag0
    
    Tecf2eci = np.array([[np.cos(ag),-np.sin(ag),0.0],[np.sin(ag),np.cos(ag),0.0],[0.0,0.0,1.0]])
    
    return Tecf2eci
    
def ecf2spherical(x):
    #ECF2SPHERICAL Convert ITRF to geocentric coordinates
    #Convert position vector in earth fixed frame to latitude, longitude, and
    #altitude based coordinates or a spherical system. 
    # INPUTS:
    #   X0 - Given state in ITRF frame at time t. [meters, meters/seconds]
    # OUTPUTS:
    #  latlonalt - Geocentric latitude, longitude, and altitude. [rad,rad,meters]

    rad = np.sqrt(x[0,0]**2 + x[0,1]**2 + x[0,2]**2)
    lat = np.arcsin(x[0,2]/rad)
    long = np.arctan2(x[0,1]/(rad*np.cos(lat)),x[0,0]/(rad*np.cos(lat)))
    latlonalt = np.array([[lat],[long],[rad]])
    
    return latlonalt
    
    
def spherical2ecf(latlonalt):
    lat = latlonalt[0]
    long = latlonalt[1]
    rad = latlonalt[2]
    
    x = rad*np.cos(lat)*np.cos(long)
    y = rad*np.cos(lat)*np.sin(long)
    z = rad*np.sin(lat)
    
    Tsph2ecf = np.array([x,y,z])
    
    return Tsph2ecf


def j2potential(t,f,h):
  #J2POTENTIAL Return 2nd derivative of motion for orbit
  #Considering only J2 perturbations and point mass force use the kinematics and return the acceleration.
  #INPUTS:
  #   X - Given state in J2000 frame at time t. [meters, meters/seconds]
  #   t - State will be predicted for each point. [seconds]
  # OUTPUTS:
  #   d2Xdt2 - acceleration from point mass potential and J2 [meters/sec^2]
  mu = 3.9860044e+14 #m^3/s^2, Earth gravitational parameter
  re = 6378137  #meters, spherical Earth radius
  J2 = 0.001082636
  C2_0 = -0.0004841695  #normalized value of harmonic coefficient
  r = lin.norm(X[0:2])
  
  #eom
  fspherical = -(mu/(r**3))*X[0:1]
  Tecf2eci = ecf2ecisimple(t,ag0)
  x = Tecf2eci.transpose()*X[0:2] 
  
  latlonalt = ecf2spherical(x)
  Tsph2ecf = spherical2ecf(latlonalt)

  fnsr = 3.0*mu*(re**2/r**4)*J2*((3.0*(np.sin(latlonalt[0])*np.sin(latlonalt[0]))-1.0)/2.0) #radial direction
  fnsphi = -mu*(re**2/r**4)*J2*3.0*np.sin(latlonalt[0])*np.cos(latlonalt[0]) #latitude
  fnslam = 0.0 #longitude
  
  #A = np.empty([1,3])
  #A[0,0] = fnsr[0]
  #A[0,1] = fnsphi[0]
  #A[0,2] = fnslam
  
  d2Xdt2[0:2] = X[3:5]
  d2Xdt2[3:5,0] = fspherical + Tecf2eci@Tsph2ecf@np.array([fnsr], [fnsphi], [fnslam])

  #d2Xdt2 = np.empty([2,3])  
  #d2Xdt2[0,0] = X[1,0]
  #d2Xdt2[0,1] = X[1,1]
  #d2Xdt2[0,2] = X[1,2]

  #d2Xdt2[[1,1,1],[0,1,2]] = fspherical + Tecf2eci@Tsph2ecf@A
  
  return d2Xdt2

def rkf45(t,X,h):
  #rkf45 Runge-Kutta-Fehlberg Numerical integrator  
  z = 0
  while (z < 10): #max step changes
      
      k1 = npml.empty((7,1))
      k2 = npml.empty((7,1))
      k3 = npml.empty((7,1))
      k4 = npml.empty((7,1))
      k5 = npml.empty((7,1))
      k6 = npml.empty((7,1))
  
      #---Begin integration
      k1 = j2potential(t,f,h)
      k2 = j2potential(t + h/5, f + (h/5)*k1, h)
      k3 = j2potential(t + (3/10)*h, f + (3/40)*h*k1 + (9/40)*h*k2, h)
      k4 = j2potential(t + (3/5)*h, f + (3/10)*h*k1 - (9/10)*h*k2 + (6/5)*h*k3, h)
      k5 = j2potential(t + h, f - (11/54)*h*k1 + (5/2)*h*k2 - (70/27)*h*k3 + (35/27)*h*k4, h)
      k6 = j2potential(t + (7/8)*h, f+(1631/55296)*h*k1+(175/512)*h*k2+(575/13824)*h*k3+(44275/110592)*h*k4+(253/4096)*h*k5, h)
  
      Xf4 = f + h*((37/378)*k1 + (250/621)*k3 + (125/594)*k4 + (512/1771)*k6)
      Xf5 = f + h*((2825/27648)*k1 + (18575/48384)*k3 + (13525/55296)*k4 + (277/14336)*k5 + (1/4)*k6) 

  
      r = (Xf5[0,0]-Xf4[0,0])/h
      check = 0.84*(1e-6/r)**(1/4)
      if r > 1e-6: 
          if check < 0.1:
              h = 0.1*h
          elif check > 4:
              h = 4*h 
          elif 1 <= check <= 4:
              h = check*h
          z = z + 1      
      elif r < 1e-6:
          z + 10
          break
          
  return Xf4

#Kepler's Equation
error = 100
E = M
while (error < 1e-6):
  E_old = E
  E = E - (E-ecc*np.sin(E)-M)/(1.0-np.cos(E))
  error = abs(E_old-E)

f = 2*np.arctan((np.sqrt((1+ecc)/(1-ecc))*np.tan(E/2)))

X = oe2rv(oe, f)

ag0 = 0.0 #Greenwich angle
t = 0.0 
h = 1.0     #initial step size
i = 0

while (t < 250):
    if not(i==0):
        Xf[i,:] = np.append(rkf45(t,Xf[i-1,:],h))             
    elif (i==0):
        Xf[i,:] = rkf45(t,X,h)
                  
    i = i + 1        
    t = t + h   
    
print(Xf)    
    
    