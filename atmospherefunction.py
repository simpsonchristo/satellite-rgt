#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Predict Orbit Based On Tx Beacon at Known/Fixed Interval"""
import numpy as np
import astropy.units as u
#import scipy.constants as con
u.imperial.enable()
from poliastro.bodies import Earth

"""Python 3.7
   EH Group, Inc. (c) 2019
   Christopher Simpson: christopher.r.simpson@simpsonaerospace.com"""
#------------------------------------------------------------------------------
def atmConditions(h):
    """ATMCONDITIONS Gradient method (lapse rate) to determine P, rho, T @ h
       MIL-STD-210A Assumptions:
       1. Pressure at SSL is same for cold and hot atmospheres (29.92 in Hg)
       2. Constant value of gravity to 100,000 ft (32.174 ft/sec**2)
       (Constant gravity will not be used in this case.)
       3. Constant composition of atmosphere throughout altitude range.
       4. Relative humidity not considered (dry air!)."""
#    g0 = -32.174*(u.imperial.ft/(u.s*u.s)) # ft/sec**2
    g = Earth.k/((Earth.R_mean+h)**2)
    deg_R = u.def_unit('deg_R', represents=u.K, format={'latex': r'\degree R'})
#    Rankine = u.def_unit('Rankine',u.deg_F + 459.67)
    T_ssl = (15 + 273.15)*deg_R 
    p_ssl = 14.695332*u.imperial.psi #psi
    rho_ssl = 0.00237*(u.imperial.slug/u.imperial.ft*u.imperial.ft*u.imperial.ft)
    R_air = 1716.5*(u.imperial.ft*u.imperial.lbf/(deg_R*u.imperial.slug))
    
    a_36089 = -0.003566 * (deg_R/u.imperial.ft) #Rankine per ft
    T_36089 = T_ssl + a_36089*((36089)*u.imperial.ft) #Rankine
    p_36089 = p_ssl * ((T_36089/T_ssl)**(-g/(a_36089*R_air))) #psi
    rho_36089 = rho_ssl * (T_36089/T_ssl)**(-((g/(a_36089*R_air))+1))
    
    a_65000 = 0* (deg_R/u.imperial.ft)
    h_65000 = 65000*u.imperial.ft
    T_65000 = T_36089 + a_65000*h_65000
    p_65000 = p_36089*np.exp((-g/(R_air*T_65000))*h_65000)
    rho_65000 = rho_36089 * np.exp((-g/(R_air*T_65000))*h_65000)
    
    h_104987 = 104987*u.imperial.ft
    a_104987 = 0.000549 * (deg_R/u.imperial.ft) #Rankine per ft
    T_104987 = T_65000 + a_104987*((h_104987 - h_65000)) #Rankine
    p_104987 = p_65000 * ((T_104987/T_65000)**(-g/(a_104987*R_air))) #psi
    rho_104987 = rho_65000 * (T_104987/T_65000)**(-((g/(a_104987*R_air))+1))
        
    if(h==0*u.imperial.ft):
        T = T_ssl
        p = p_ssl
        rho = rho_ssl
    
    elif(h <= 36089*u.imperial.ft):
        a = -0.003566 * (deg_R/u.imperial.ft) #Rankine per ft
        T = T_ssl + a*(h-0) #Rankine
        p = p_ssl * pow((T/T_ssl),(-g/a*R_air)) #psi
        rho = rho_ssl * pow((T/T_ssl),-((g/(a*R_air))+1))
    
    elif(h <=65000*u.imperial.ft):   
        a = 0* (deg_R/u.imperial.ft)
        dh = (h - (36089*u.imperial.ft))
        T = T_36089 + a*dh
        p = p_36089*np.exp((-g/(R_air*T))*dh)
        rho = rho_36089 * np.exp((-g/(R_air*T))*dh)
    
    elif(h <=104987*u.imperial.ft):        
        a = 0.000549 * (deg_R/u.imperial.ft) #Rankine per ft
        T = T_65000 + a*((h - h_65000)) #Rankine
        p = p_65000 * ((T/T_65000)**(-g/(a*R_air))) #psi
        rho = rho_65000 * (T/T_65000)**(-((g/(a*R_air))+1))
    
    elif(h <=160000*u.imperial.ft):        
        a = 0.001536 * (deg_R/u.imperial.ft) #Rankine per ft
        T = T_104987 + a*((h - h_104987)) #Rankine
        p = p_104987 * ((T/T_104987)**(-g/(a*R_air))) #psi
        rho = rho_104987 * (T/T_104987)**(-((g/(a*R_air))+1))
        
    return T,p,rho,R_air,g
        
        
       