#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Predict Orbit Based On Tx Beacon at Known/Fixed Interval"""
import numpy as np
import astropy.units as u
u.imperial.enable()
from astropy.coordinates import EarthLocation
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from poliastro.core.angles import fp_angle
import matplotlib.pyplot as plt
#from poliastro.coordinates import Longitude
#custom
#from staged_launch import aim120ToSpace
from atmospherefunction import atmConditions

"""Python 3.7
   EH Group, Inc. (c) 2019
   Christopher Simpson: christopher.r.simpson@simpsonaerospace.com"""
#------------------------------------------------------------------------------
earthrotspeed = 464.5*u.m/u.s
#dv_required to reach orbit
orbit_alt = (np.linspace(250,600,num=7)*u.km).to(u.m) + Earth.R_mean
orbit_ecc = 0*u.one
#ToDo: 1500 m/s losses due to drag+gravity should be added here (SPAD)
orbit_inc = np.tile(np.deg2rad(np.linspace(0,90,num=91)), (len(orbit_alt),1))
Y_inc, X_alt = np.meshgrid(orbit_inc[0,:], orbit_alt)
orbit_speed = np.sqrt(Earth.k/X_alt) + 1500*(u.m/u.s)
#launch site location 
loc_equator = EarthLocation.from_geodetic(-87.545*u.deg,
                                     0*u.deg,
                                     height=(19.25*u.m + 13716*u.m),
                                     ellipsoid='WGS84')
launch_azI = np.arcsin(np.cos(orbit_inc)/np.cos(0))*u.rad
launch_sitevelocity = earthrotspeed*np.cos(0)
launch_correctaz = np.arctan2(launch_sitevelocity*np.cos(launch_azI),(-(earthrotspeed*np.cos(orbit_inc)) + orbit_speed))
launch_az = launch_correctaz + launch_azI
#aircraft
M = np.linspace(0.3,1.6)*u.one
ac_alt = 45000*u.imperial.ft
T, p, rho, R_air, g0 = atmConditions(ac_alt)
ac_speedofsound = ((np.sqrt(1.4*R_air*T)).value)*0.3048*u.m/u.s
ac_speed = ac_speedofsound*M
orbit_DV = np.ones((len(orbit_alt),orbit_inc.shape[1],len(ac_speed)))
for k in range(len(ac_speed)):
    for i in range(len(orbit_alt)):
        for j in range(orbit_inc.shape[1]):
            orbit_target = Orbit.circular(Earth, (orbit_alt[i] - Earth.R_mean), (orbit_inc[i,j]*u.rad).to(u.deg))
            orbit_target = orbit_target.propagate(47.30*u.s)
            orbit_fp = fp_angle(orbit_target.nu, orbit_target.ecc)
            holder = [[(-(orbit_speed[i,j]-ac_speed[k])*np.cos(launch_az[i,j])*np.cos(orbit_fp)).value,
                      ((orbit_speed[i,j]-ac_speed[k])*np.sin(launch_az[i,j])*np.cos(orbit_fp) - launch_sitevelocity).value,
                      ((orbit_speed[i,j]-ac_speed[k])*np.sin(orbit_fp)).value]]
            holder = np.linalg.norm(holder)
            orbit_DV[i,j,k] = holder

#missile_DV = aim120ToSpace()

#plotting
inc_Y, alt_X = np.meshgrid(np.rad2deg(orbit_inc[0,:]), orbit_alt)
fig, ax = plt.subplots()
dv = ax.contourf(alt_X, inc_Y, orbit_DV[:,:,0])
ax.clabel(dv, inline=1, fontsize=8)
ax.set_title('$\Delta$V Required to Achieve Orbit')

plt.show()