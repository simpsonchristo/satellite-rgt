#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Predict Orbit Achieved From Different Second Stages on Launch Vehicle"""
import numpy as np
import astropy.units as u
u.imperial.enable()
from poliastro.bodies import Earth
import matplotlib.pyplot as plt
from scipy.optimize import newton
#custom
from staged_launch import rocketDV
from atmospherefunction import atmConditions
from gamma_adjust import gammaAdjust

"""Python 3.7
   EH Group, Inc. (c) 2019
   Christopher Simpson: christopher.r.simpson@simpsonaerospace.com"""
#------------------------------------------------------------------------------
'''Rocket'''
'Rocket stages, payload, and DV for AIM-120, AIM-7M, and AGM-84HK sizing.'
DV_s2_aim120, aim120_DV, aim120_pay, aim120_s1, aim120_s2 = rocketDV(totalLen = 4.389*u.m, Dia=0.3505*u.m,)
DV_s2_aim7m, aim7m_DV, aim7m_pay, aim7m_s1, aim7m_s2 = rocketDV(Dia=0.2032*u.m)
DV_s2_agm84HK, agm84HK_DV, agm84HK_pay, agm84HK_s1, agm84HK_s2 = rocketDV(totalLen=4.38912*u.m, Dia=0.3429*u.m)
DV_s2_10kg = np.array([DV_s2_aim120[DV_s2_aim120.shape[0] - 1,:], DV_s2_aim7m[DV_s2_aim7m.shape[0] - 1,:], DV_s2_agm84HK[DV_s2_agm84HK.shape[0] - 1,:]]).flatten()*(u.m/u.s)
tb_s1 = 62.8*u.s
tb_s2 = 60.0*u.s
T_s1 = np.repeat(2855.08*u.N, 12)
m_s2 = (np.array([(aim120_s2.mprop() + aim120_s2.m_i()), (aim7m_s2.mprop() + aim7m_s2.m_i()), (agm84HK_s2.mprop() + agm84HK_s2.m_i())]).flatten()*u.g).to(u.kg)
T_s2 = (m_s2*DV_s2_10kg/tb_s2).to(u.N)
tb_s2 = np.array([20, 40, 20, 40,
                  20, 40, 20, 40,
                  10, 30, 10, 30])*u.s

mdot_s1 = (np.array([[aim120_s1.monoprop().to(u.kg).value, aim7m_s1.monoprop().to(u.kg).value, agm84HK_s1.monoprop().to(u.kg).value]])/tb_s1.value)
mdot_s1 = np.array([mdot_s1[0,0],mdot_s1[0,0],mdot_s1[0,0],mdot_s1[0,0],
                    mdot_s1[0,1],mdot_s1[0,1],mdot_s1[0,1],mdot_s1[0,1],
                    mdot_s1[0,2],mdot_s1[0,2],mdot_s1[0,2],mdot_s1[0,2]])
mdot_s2 = (np.array([aim120_s2.mprop().to(u.kg).value, aim7m_s2.mprop().to(u.kg).value, agm84HK_s2.mprop().to(u.kg).value]).flatten()/tb_s2.value)
#mdot_s2 = mdot_s2.flatten()
'Cross-sectional area/Cd of missile.'
aim120_A  = np.repeat(np.pi*(aim120_s1.Dia/2)**2, 4)
aim7m_A   = np.repeat(np.pi*(aim7m_s1.Dia/2)**2, 4)
agm84HK_A = np.repeat(np.pi*(agm84HK_s1.Dia/2)**2, 4)
Area = np.array([aim120_A, aim7m_A, agm84HK_A]).flatten()
Cd = 0.1

'''Aircraft'''
M = 0.6*u.one
ac_alt = 45000*u.imperial.ft
T, p, rho, R_air = atmConditions(ac_alt)
ac_speedofsound = ((np.sqrt(1.4*R_air*T)).value)*0.3048*u.m/u.s
ac_speed = ac_speedofsound*M
'Launch site (lon, lat, alt)'
loc_equator = np.array([[np.deg2rad(-87.545), 0, (ac_alt.value*12*2.54/100)]])
earthrotspeed = 464.5*u.m/u.s
v_launchsite = earthrotspeed*np.cos(loc_equator[0,1])
#launch_azI = np.arcsin(np.cos(orbit_inc)/np.cos(loc_equator[0,1]))*u.rad
#launch_correctaz = np.arctan2(v_launchsite*np.cos(launch_azI),(-(earthrotspeed*np.cos(orbit_inc)) + orbit_speed))
#launch_az = launch_correctaz + launch_azI
launch_az = np.deg2rad(89.5)

h0 = loc_equator[0,2]
v0 = v_launchsite + ac_speed

'''Target Orbits'''
orbit_alt = (np.linspace(250,650,num=7)*u.km).to(u.m) + Earth.R_mean
orbit_ecc = 0*u.one
orbit_inc = np.tile(np.deg2rad(np.linspace(0,90,num=91)), (len(orbit_alt),1))
Y_inc, X_alt = np.meshgrid(orbit_inc[0,:], orbit_alt - Earth.R_mean)

'''Initial state'''
X = np.empty((int(tb_s1.value*2),12,4))
X[0,:,0] = 0*u.m #range
X[0,:,1] = h0    #altitude
X[0,:,2] = 0     #vx
X[0,:,3] = v0    #vy
'flight-path angle'
gamma = np.array([np.deg2rad(90)])
dgamma = np.deg2rad(np.linspace(90,9,num=int(tb_s1.value*2)))
'initial dv'
#DV_1 = np.empty((int(tb_s1.value),12))
#DV_2 = np.empty((int(tb_s1.value*2) - int(tb_s1.value + 3),12))
DV = np.zeros((int(tb_s1.value*2), mdot_s2.shape[0]))
'initial mass'
mass = np.empty((int(tb_s1.value*2),12))
mass[0,:] = np.array([(aim120_s2.mprop() + aim120_s2.m_i() + aim120_s1.monoprop() + aim120_s1.m_i() + aim120_pay.m[len(aim120_pay.m)-1]),
                      (aim7m_s2.mprop() + aim7m_s2.m_i() + aim7m_s1.monoprop() + aim7m_s1.m_i() + aim7m_pay.m[len(aim7m_pay.m)-1]),
                      (agm84HK_s2.mprop() + agm84HK_s2.m_i() + agm84HK_s1.monoprop() + agm84HK_s1.m_i() + agm84HK_pay.m[len(agm84HK_pay.m)-1])]).flatten()
mass[int(tb_s1.value),:] = np.array([(aim120_s2.mprop() + aim120_s2.m_i() + aim120_pay.m[len(aim120_pay.m)-1]),
                                     (aim7m_s2.mprop() + aim7m_s2.m_i() + aim7m_pay.m[len(aim7m_pay.m)-1]),
                                     (agm84HK_s2.mprop() + agm84HK_s2.m_i() + agm84HK_pay.m[len(agm84HK_pay.m)-1])]).flatten()
mass[int(tb_s1.value + 4),:] = mass[int(tb_s1.value),:]
mass = (mass*u.g).to(u.kg).value

masshold_s1 = mass[0,:] - np.array([np.repeat(aim120_s1.monoprop().to(u.kg).value, 4), np.repeat(aim7m_s1.monoprop().to(u.kg).value, 4), np.repeat(agm84HK_s1.monoprop().to(u.kg).value, 4)]).flatten()
masshold_s2 = mass[int(tb_s1.value),:] - np.array([(aim120_s2.mprop().to(u.kg).value), (aim7m_s2.mprop().to(u.kg).value), (agm84HK_s2.mprop().to(u.kg).value)]).flatten()

rho = np.empty((1,12))
g = np.empty((1,12))
'''For loops'''
for t in range(0,int(tb_s1.value*2)-1):
    dt = abs(1)
    'atmospheric conditions at current altitude'
    alt = X[t,:,1]*(100/(2.54*12))*u.imperial.ft
    for i in range(mass.shape[1]):
        Temp, press, p, R_air = atmConditions(alt[i])
        rho[0,i] = p.value/515.379 #* u.kg/(u.m*u.m*u.m)
        g[0,i] = Earth.k.value/((Earth.R_mean.value+X[t,i,1])**2)
    
    'calculate drag losses'
    Vsqrd = (X[t,:,2]**2 + X[t,:,3]**2) 
    Dovermass = 0.5*rho*Cd*Area*Vsqrd/mass[t,:]
    
    
    'flight path angle'
    if(t==3):
      gamma = np.array([np.deg2rad(89.5)])  
    elif(t>3):
#        fprime = lambda y, a: -2*(Earth.k.value/(Earth.R_mean.value + X[a[0],:,1])**2) * (np.cos(a[3])/(a[1]*a[2]))*1*y
        gamma_adjust = 10.2
        Vsqrdminus1 = (X[t - 1,:,2]**2 + X[t - 1,:,3]**2)
        dgamma_dt = -2*(Earth.k.value/(Earth.R_mean.value + X[t,:,1])**2) * (np.cos(gamma)/(np.sqrt(Vsqrd)*np.sqrt(Vsqrdminus1)))*dt*gamma_adjust
        gamma = gamma + dgamma_dt
        if(any(gamma)<0):
            flag = 1
#        gamma = dgamma[t]
    
    'calculate thrust and steering contribution'
    if (t<int(tb_s1.value)):
        thrust = -(T_s1/mass[t,:])*np.cos(gamma) + (T_s1/mass[t,:])
#        DV_1[t,:] = (T_s1/mass[t,:])
        DV[t + 1,:] = (T_s1/mass[t,:])
    elif(t==int(tb_s1.value)):
        for k in range(1, 4):
            DV[t + k,:] = (0/mass[t,:])
    elif(t>int(tb_s1.value + 3)):
        for k in range(mass.shape[1]):
            if (t>int(tb_s2[k].value + tb_s1.value + 3)):
                T_s2[k] = 0
                mass[t, k] = masshold_s2[k]
                
        thrust = -(T_s2/mass[t,:])*np.cos(gamma) + (T_s2/mass[t,:])
        DV[t + 1,:] = (T_s2/mass[t,:])
   
    'calculate gravity loss'
    g_acc = g*np.sin(gamma)
    
    'determine velocity at this point'
    X[t + 1, :, 3] = X[t, :, 3] + ((thrust.value - Dovermass - g_acc)*(dt))*np.sin(gamma)
    X[t + 1, :, 2] = X[t, :, 2] + ((thrust.value - Dovermass - g_acc)*(dt))*np.cos(gamma)
    
    'determine new height and range'
    X[t + 1, :, 1] = X[t, :, 1] + X[t,:,3]*dt + 0.5*((thrust.value - Dovermass - g_acc)*(dt**2))*np.sin(gamma)
    X[t + 1, :, 0] = X[t, :, 0] + X[t,:,2]*dt + 0.5*((thrust.value - Dovermass - g_acc)*(dt**2))*np.cos(gamma)
    
    'update mass'
    if (t<int(tb_s1.value + 1)):
        mass[t + 1, :] = mass[t,:] - mdot_s1*dt
        for k in range(0, mass.shape[1]):
            if (mass[t + 1, k] < masshold_s1[k]):
                mass[t + 1, k] = masshold_s1[k]
    elif(t==int(tb_s1.value + 1)):
        for k in range(0, 4):
            mass[t + k, :] = m_s2.value
    elif(t>int(tb_s1.value + 3)):
        mass[t + 1, :] = mass[t,:] - mdot_s2*dt
        for k in range(0, mass.shape[1]):
            if (mass[t + 1, k] < masshold_s2[k]):
                mass[t + 1, k] = masshold_s2[k]
            if (t>int(tb_s2[k].value + tb_s1.value + 3)):
                mass[t, k] = masshold_s2[k]
    
    
#DV_1 = DV_1.sum(axis=0)
#DV_2 = DV_2.sum(axis=0)# + DV_1        
'''Plotting'''
'HTPB/LOX'
fig, ax = plt.subplots()
plt.title('$2^{nd}$ Stage: HTPB/LOX')
rangeX = np.array([X[:,1,0], X[:,3,0], X[:,5,0], X[:,7,0], X[:,9,0], X[:,11,0]]).transpose()
altitY = np.array([X[:,1,1], X[:,3,1], X[:,5,1], X[:,7,1], X[:,9,1], X[:,11,1]]).transpose()
ax.plot(rangeX, altitY)
plt.ylabel('Altitude [m]')
plt.xlabel('Range [m]')
for i in range(0, rangeX.shape[1]):
    plt.annotate('MBO', (rangeX[int(tb_s1.value),i], altitY[int(tb_s1.value),i]), xytext=(3, 1.5), 
                 arrowprops=dict(facecolor='black', shrink=0.0005))
    plt.annotate('SECO', (rangeX[int(tb_s1.value + tb_s2[i].value + 5),i], altitY[int(tb_s1.value + tb_s2[i].value + 5),i]), xytext=(2000, 8000), 
                     arrowprops=dict(facecolor='blue', shrink=0.0005))
ax.axis('tight')
#dv = ax.contour(X_alt/1000, np.rad2deg(Y_inc), X[:,:,], levels=missile_DV[0,:], colors=('k',), linestyles=('-',), linewidths=(2,))
#dvfill = ax.contourf(X_alt, np.rad2deg(Y_inc), orbit_DV[:,:,k])
#ax.clabel(dv, fmt = '%2.1d 5 kg', colors = 'k', inline=1, fontsize=8)
#plt.colorbar(dvfill, shrink=0.8, extend='both')
#ax.set_title('$\Delta$V Required to Achieve Orbit')

plt.show()