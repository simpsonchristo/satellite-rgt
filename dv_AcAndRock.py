#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Predict Orbit Based On Tx Beacon at Known/Fixed Interval"""
import numpy as np
import astropy.units as u
u.imperial.enable()
#from astropy.coordinates import EarthLocation, SkyCoord, AltAz, GCRS, ICRS
#from astropy.time import Time, TimeDelta
#from poliastro.bodies import Earth
#from poliastro.twobody import Orbit
#from poliastro.twobody.propagation import mean_motion
#from poliastro.core.angles import D_to_nu
#from poliastro.coordinates import Longitude
#custom
from atmospherefunction import atmConditions

"""Python 3.7
   EH Group, Inc. (c) 2019
   Christopher Simpson: christopher.r.simpson@simpsonaerospace.com"""
#------------------------------------------------------------------------------
"""AMRAAM (AIM-120)"""
L = (12 * u.imperial.ft).to(u.m)
D = (7.0 * u.imperial.inch).to(u.m)
Vol =  (L/2)*(np.pi*(D/2)**2) #motor section is about half of AMRAAM (6.11 ft)
W = 335 * u.imperial.lbf
Temp, press, air_density, air_R, g0 = atmConditions(45000*u.imperial.ft)
a = np.sqrt(1.4*air_R*Temp)
V_ac = a*np.linspace(0.3,1.6)
V_amraam = np.linspace(V_ac,(a*4.0))
mpay_sparrow = (30*u.kg).to(u.imperial.slug)
#Prop is HTPB --> Isp = 225.7 s <https://pdfs.semanticscholar.org/a15a/a1b29fffeebfcdb70965e4241af07eb5f00f.pdf>
Isp_amraam = 225.7 * u.s
DV_amraam = -Isp_amraam*g0*np.log(W/(W*0.9 - (mpay_sparrow.value*(-g0.value)*u.imperial.lbf)))

htpb_density = 930*u.kg/(u.m*u.m*u.m)
htpb_mass = htpb_density*Vol

check_DV = DV_amraam.value/a.value



#
#"""Sparrow (AIM-7M)"""
#L = 12 * u.imperial.ft
#D = 8.0 * u.imperial.inch
#W = 508 * u.imperial.lbf
#
#Tmax_sparrow = (34.30e+3 * u.N).to(u.imperial.lbf)
#Mmotor_sparrow = 6.83781*u.imperial.slug
#tburn_sparrow = 1.80*u.s
#mdot_sparrow = Mmotor_sparrow*0.9/tburn_sparrow
#Isp_sparrow = (Tmax_sparrow.value/(-g0.value*mdot_sparrow.value)) * u.s
#
#DV_sparrow = Isp_sparrow*g0*np.log(W/(W*0.9 - (mpay_sparrow.value*(-g0.value)*u.imperial.lbf)))