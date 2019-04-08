#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 09:18:30 2019

@author: willjpatton
"""
import numpy as np

#initial position
lat = np.pi
long = 2*np.pi
h = 450.0 #km

#Days of Orbit
Do = 3.0

ecc = 0.0

Re = 6378137.0 #meters Radius of Earth


h = h*1000.0
ao = h + Re

we = (2*np.pi*(365.256360+1))/31558149.504  #rad/s Angular Velocity of the Earth
mu = 3.9860044e+14 #m^3/s^2, Earth gravitational parameter

#Orbital Period
T = 2*np.pi*np.sqrt(ao**3/mu)

ND = (2*np.pi)/(T*we)
N = ND*Do
No = round(N)
ND = No/Do
T = (2*np.pi)/(ND*we)

a = ((T/(2*np.pi))**2*mu)**(1/3)


ws = 2*np.pi/31558149.504
J2 = 0.001082636

no = No/Do
#inc = np.arccos(-2/3*(a/Re)**2*(ws/(no*J2)))
kh = 10.10949
inc = np.arccos((-1/kh)*(a/Re)**(7/2))

incd = np.degrees(inc)

ALG = 0.0

RAAN = long - np.arcsin(np.tan(lat)/np.tan(inc))+((np.arcsin(np.sin(lat)/np.sin(inc))/(2*np.pi))*T*we)+ALG

RAANd = np.degrees(RAAN)

#if we have intial a,ecc,incd,RAANd can and lat,long,h can we convert all to oribital elements






