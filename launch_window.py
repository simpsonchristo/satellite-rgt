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
from poliastro.core.angles import D_to_nu
from poliastro.coordinates import Longitude
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
Vgs_lat = Veq*np.cos(np.deg2rad(33.213))
#aircraft speed
V_ac = (np.array([np.linspace(0.3,1.6,num=14)])).transpose() * u.one #mach number

h_ac = 45000*u.imperial.ft
Temp, press, density, R_air = atmConditions(h_ac)
speed_of_sound = np.sqrt(1.4*R_air*Temp)
V_ac = speed_of_sound*V_ac
#speed of satellite immediately after launch
V_sat = V_ac + (7800*u.m/u.s)
#inertial launch azimuth
beta_i = np.arcsin(np.cos(np.deg2rad(coe[2,0]))/np.cos(np.deg2rad(33.213)))*u.rad
#correction to account for Earth rotation velocity addition
gam_correction = np.arctan2(((Vgs_lat*np.cos(beta_i))+V_ac),(V_sat - (Veq*np.cos(np.deg2rad(coe[2,0])))))
#launch azimuth (ascending node)
beta = beta_i + gam_correction
#local sidereal time
theta = Longitude(coe[3,0]*u.deg + t_ob.sidereal_time('apparent',longitude=360-87.545))

nu = D_to_nu(coe[5,0]*u.deg)
iss_target = Orbit.from_classical(Earth, coe[0,0]*u.m,
                                  coe[1,0]*u.one,
                                  coe[2,0]*u.deg,
                                  coe[3,0]*u.deg,
                                  coe[4,0]*u.deg,
                                  nu*u.deg,
                                  epoch = t_ob)
X = mean_motion(Earth.k, iss_target.r, iss_target.v, t_del_nounit)
#iss_SkyCoord = SkyCoord(X[0][:,0], X[0][:,1], X[0][:,2], frame = 'gcrs', representation_type='cartesian')
#ssaltaz_hh = iss_SkyCoord.transform_to(frame_now_hh)
#iss_target_SEZ = transform(iss_target, 'gcrs', frame_now_hh)
V_req = np.array([[-np.cos(beta)*X[1][:,0]],[(X[1][:,1]*np.sin(beta))-(Vgs_lat+V_ac)],[0*X[1][:,2]]])
