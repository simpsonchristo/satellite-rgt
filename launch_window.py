#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Predict Orbit Based On Tx Beacon at Known/Fixed Interval"""
import numpy as np
import astropy.units as u
u.imperial.enable()
from astropy.coordinates import EarthLocation, SkyCoord, AltAz, GCRS, ICRS
from astropy.time import Time, TimeDelta
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from poliastro.twobody.propagation import mean_motion
#custom
import readtle
from downloadfileurl import download_file
from atmospherefunction import atmConditions

"""Python 3.7
   EH Group, Inc. (c) 2019
   Christopher Simpson: christopher.r.simpson@simpsonaerospace.com"""
#------------------------------------------------------------------------------
#COE of ISS
rooturl = 'http://www.celestrak.com/NORAD/elements/supplemental/'
iss_tles = download_file(rooturl+'iss.txt')
coe = readtle.readtle(iss_tles)
#a,    semi-major axis [m]
#e,    eccentricity []
#i,    inclination [deg]
#RAAN, right angular ascending node [deg]
#AOP,  argument of perigee [deg]
#M,    mean anomaly [deg]

#Launched over Hardaway Hall
loc_hh = EarthLocation.from_geodetic(-87.545*u.deg,
                                     33.213*u.deg,
                                     height=19.25*u.m,
                                     ellipsoid='WGS84')

#time
t_ob = Time(Time.now(),scale='tai')
delta_t = TimeDelta(np.linspace(0,86400,num=86401)*u.second,scale='tai')
t_del_nounit = np.linspace(0,86400,num=86401)*u.second

#topocentric frame at Hardaway Hall
frame_now_hh = AltAz(obstime=t_ob+delta_t,location=loc_hh)

#Velocity of Earth at equator
Veq = 464.5 *(u.m/u.s)
#Earth inertial velocity addition to launch
Vgs_lat = Veq*np.cos(np.deg2rad(33.213))*(u.m/u.s)
#aircraft speed
V_ac = np.linspace(0.3,1.6,num=14) * u.one #mach number
h_ac = 45000*u.imperial.ft
Temp, press, density, R_air = atmConditions(h_ac)
speed_of_sound = np.sqrt(1.4*R_air*Temp)
V_ac = speed_of_sound*V_ac
#speed of satellite immediately after launch
V_sat = V_ac + (7800*u.m/u.s)
#inertial launch azimuth
beta_i = np.arcsin(np.cos(np.deg2rad(coe[2,0]))/np.cos(np.deg2rad(33.213)))
#correction to account for Earth rotation velocity addition
gam_correction = np.arctan2(((Vgs_lat*np.cos(beta_i))+V_ac)/(V_sat - (Veq*np.cos(np.deg2rad(coe[2,0])))))
#launch azimuth (ascending node)
beta = beta_i + gam_correction

