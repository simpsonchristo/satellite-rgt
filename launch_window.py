#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Predict Orbit Based On Tx Beacon at Known/Fixed Interval"""
import numpy as np
import astropy.units as u
u.imperial.enable()
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from poliastro.core.angles import fp_angle
import matplotlib.pyplot as plt
from scipy.integrate import odeint
#custom
from staged_launch import rocketDV
from atmospherefunction import atmConditions

"""Python 3.7
   EH Group, Inc. (c) 2019
   Christopher Simpson: christopher.r.simpson@simpsonaerospace.com"""
#------------------------------------------------------------------------------
aim120_DV, aim120_pay, aim120_s1, aim120_s2 = rocketDV()
aim7m_DV, aim7m_pay, aim7m_s1, aim7m_s2 = rocketDV(Dia=0.2032*u.m)
agm84HK_DV, agm84HK_pay, agm84HK_s1, agm84HK_s2 = rocketDV(totalLen=4.38912*u.m, Dia=0.3429*u.m)
tb_s1 = 62.8*u.s
T_s1 = 2855.08*u.N
mdot_s1 = (np.array([[aim120_s1.monoprop().value, aim7m_s1.monoprop().value, agm84HK_s1.monoprop().value]])/tb_s1.value)*(u.g/u.s)
#cross-sectional area and Cd of missile
aim120_A  = np.pi*(aim120_s1.Dia/2)**2
aim7m_A   = np.pi*(aim7m_s1.Dia/2)**2
agm84HK_A = np.pi*(agm84HK_s1.Dia/2)**2
Cd = 0.1

#aircraft
#M = np.linspace(0.3,1.6, num=(16 - 3))*u.one
M = 0.6*u.one
ac_alt = 45000*u.imperial.ft
T, p, rho, R_air, g = atmConditions(ac_alt)
ac_speedofsound = ((np.sqrt(1.4*R_air*T)).value)*0.3048*u.m/u.s
ac_speed = ac_speedofsound*M

#'ToDo: h0, v0 from a/c (launch platform) + earth contribution'
#launch site (lon, lat, alt)
loc_equator = np.array([[np.deg2rad(-87.545), 0, (ac_alt.value*12*2.54/100)]])
earthrotspeed = 464.5*u.m/u.s
v_launchsite = earthrotspeed*np.cos(loc_equator[0,1])
#orbital parameters
orbit_alt = (np.linspace(250,650,num=7)*u.km).to(u.m) + Earth.R_mean
orbit_ecc = 0*u.one
orbit_inc = np.tile(np.deg2rad(np.linspace(0,90,num=91)), (len(orbit_alt),1))
Y_inc, X_alt = np.meshgrid(orbit_inc[0,:], orbit_alt - Earth.R_mean)
#ToDo: 1500 m/s losses due to drag+gravity should be added if phasing used (SPAD)
#orbit_speed = np.sqrt(Earth.k/X_alt)

'ToDo: Determine gamma adjust value to find turnover rate of flight path angle'

#launch_azI = np.arcsin(np.cos(orbit_inc)/np.cos(loc_equator[0,1]))*u.rad
#launch_correctaz = np.arctan2(v_launchsite*np.cos(launch_azI),(-(earthrotspeed*np.cos(orbit_inc)) + orbit_speed))
#launch_az = launch_correctaz + launch_azI
launch_az = np.deg2rad(90)


#orbit_DV = np.ones((len(orbit_alt),orbit_inc.shape[1],len(ac_speed)))
v0 = np.array([[0], [ac_speed.value+v_launchsite.value]])
x0 = np.array([[0], [(ac_alt.value*12*2.54/100)]])

theta0 = np.arccos(v0[0,0]/np.linalg.norm(v0))
aim120_m0  = aim120_s1.monoprop() + aim120_s1.m_i() + aim120_s2.mprop() + aim120_s2.m_i()
aim7m_m0   = aim7m_s1.monoprop() + aim7m_s1.m_i() + aim7m_s2.mprop() + aim7m_s2.m_i()
agm84HK_m0 = agm84HK_s1.monoprop() + agm84HK_s1.m_i() + agm84HK_s2.mprop() + agm84HK_s2.m_i()

X_aim120  = np.array([[x0],[v0], [aim120_m0.value,  aim120_s1.monoprop().value/tb_s1.value],  [aim120_A,   theta0]])
X_aim7m   = np.array([[x0],[v0], [aim7m_m0.value,   aim7m_s1.monoprop().value/tb_s1.value],   [aim7m_A,    theta0]])
X_agm84HK = np.array([[x0],[v0], [agm84HK_m0.value, agm84HK_s1.monoprop().value/tb_s1.value], [agm84HK_A,  theta0]])
t = np.linspace(0,tb_s1.value + 100, num=1628)
'''Integrate over first stage burn time'''
'''ToDo: (Integration) 
        Calculate drag force at each step. (Requires calling atmosphere func)
        Calculate velocity using a_g, a_thrust, and a_D. Remember to use flight-path angle.
        Mass loss will be mdot*dt.
        Altitude will be h(i) = 0.5*V*dt*sin(fp_angle)+h(i-1)'''
def fprime(X, t):
    #constant
    T_s1 = 2855.08*u.N
    Cd = 0.1*u.one
    #pull out state info
    h = X[0,1]*100/(2.54*12)*u.imperial.ft
    v = X[1,:]
    mass = X[2,0]*u.g
    mdot = X[2,1]*u.g/u.s
    A = X[3,0]*u.m*u.m
    theta = X[3,1]*u.rad
    #use state info for operations
    theta = np.arccos(v[0,0]/np.linalg.norm(v)) #angle of path with horizontal
    Temp, press, rho, R_air, g = atmConditions(h)
    psi = np.deg2rad(90)*u.rad #thrust direction
    D = 0.5 * Cd * (rho.value*515.379*u.kg/(u.m*u.m*u.m)) * A * (np.linalg.norm(v))**2 #drag
    ay = (T_s1/mass)*np.cos(psi-theta) - (D/mass)*np.sin(theta) - g*np.sin(theta)
    ax = (T_s1/mass)*np.sin(psi-theta) - (D/mass)*np.cos(theta) - g*np.cos(theta)
    a = np.array([[ax.value, ay.value]])*ax.unit
    m_new = mass - mdot*t
    
    Xdot = np.array([[v], [a], [m_new, mdot], [A, theta]])
    return Xdot

X_final = odeint(fprime, X_aim120, t)
    
    
#for k in range(len(ac_speed)):
#    for i in range(len(orbit_alt)):
#        for j in range(orbit_inc.shape[1]):
#            orbit_target = Orbit.circular(Earth, (orbit_alt[i] - Earth.R_mean), (orbit_inc[i,j]*u.rad).to(u.deg))
#            orbit_target = orbit_target.propagate(47.30*u.s)
#            orbit_fp = fp_angle(orbit_target.nu, orbit_target.ecc)
#            holder = [[(-(orbit_speed[i,j]-ac_speed[k])*np.cos(launch_az[i,j])*np.cos(orbit_fp)).value,
#                      ((orbit_speed[i,j]-ac_speed[k])*np.sin(launch_az[i,j])*np.cos(orbit_fp) - v_launchsite).value,
#                      ((orbit_speed[i,j]-ac_speed[k])*np.sin(orbit_fp)).value]]
#            holder = np.linalg.norm(holder)
#            orbit_DV[i,j,k] = holder
#orbit_target = Orbit.circular(Earth, (650*u.km), (0*u.rad).to(u.deg))
#orbit_target = orbit_target.propagate(tb_s1)
#orbit_fp = fp_angle(orbit_target.nu, orbit_target.ecc)
#holder = [[((-(np.linalg.norm(orbit_target.v))*np.cos(launch_az)*np.cos(orbit_fp))),
#          ((np.linalg.norm(orbit_target.v))*np.sin(launch_az)*np.cos(orbit_fp) - 0.4645),
#          ((np.linalg.norm(orbit_target.v))*np.sin(orbit_fp))]]
#holder = np.linalg.norm(holder)
#orbit_DV[i,j,k] = holder
            
#    #plotting
#    fig, ax = plt.subplots()
##    dv = ax.contour(X_alt/1000, np.rad2deg(Y_inc), orbit_DV[:,:,k], levels=missile_DV[0,:], colors=('k',), linestyles=('-',), linewidths=(2,))
#    dvfill = ax.contourf(X_alt, np.rad2deg(Y_inc), orbit_DV[:,:,k])
##    ax.clabel(dv, fmt = '%2.1d 5 kg', colors = 'k', inline=1, fontsize=8)
#    plt.colorbar(dvfill, shrink=0.8, extend='both')
#    ax.set_title('$\Delta$V Required to Achieve Orbit')
#    
#    plt.show()