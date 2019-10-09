#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Predict Total DV Generated By Different Staged Motors"""
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
#custom

"""Python 3.7
   EH Group, Inc. (c) 2019
   Christopher Simpson: christopher.r.simpson@simpsonaerospace.com"""
#------------------------------------------------------------------------------
def aim120ToSpace(stage_n = 2):
    #staging
    '''Motors are HTPB+N2O, HTPB+LOX, Paraffin+N2O, Paraffin+LOX'''
    stage_motor_Isp = [[225.7*u.s, 368*u.s, 226*u.s, 372*u.s]]
    stage_motor_OF = [[(0.300/0.077), 2.5, (0.300/0.077), 2.7]]
    #density of HTPB, SP-1a (paraffin), LOX, N2O
    stage_motor_rho = [[915.0e+3 * (u.g/(u.m*u.m*u.m)), 920.0e+3 * (u.g/(u.m*u.m*u.m)), 1141.0e+3 * (u.g/(u.m*u.m*u.m)), 772.3e+3 * (u.g/(u.m*u.m*u.m))]]
    #missile sizing
    stage_all_length = 3.9624*u.m #13 ft
    stage_all_dia = 0.254*u.m #10 in
    stage_all_vol = stage_all_length*np.pi*(stage_all_dia/2)**2
    #stage 1
    stage_1_length = 1.8288*u.m #6 ft
    stage_1_dia = 0.254*u.m #10 in
    stage_1_vol = stage_1_length*np.pi*(stage_1_dia/2)**2
    intermediate = 1 + (stage_motor_rho[0][3]/(stage_motor_OF[0][0]*stage_motor_rho[0][0]))
    stage_1_mHTPB = (stage_motor_rho[0][3]/stage_motor_OF[0][0])*(1/intermediate)*stage_1_vol
    stage_1_mN2O = stage_motor_rho[0][3]*(stage_1_vol - (stage_1_mHTPB/stage_motor_rho[0][0]))
    stage_1_mprop = stage_1_mHTPB + stage_1_mN2O
    
    #payload section
    pay_length = 0.36*u.m
    pay_height = 0.12*u.m
    pay_width = 0.24*u.m
    pay_vol = pay_length*pay_width*pay_height
    stage_pay_vol = pay_length*np.pi*(stage_all_dia/2)**2
    stage_pay_m = np.linspace(5.0e+3,10.0e+3)*u.g
    #stage_pay_m = 10.0e+3 * u.g
    
    #stage 2
    stage_2_length = stage_all_length - stage_1_length - pay_length
    stage_2_vol = stage_all_vol - stage_1_vol - stage_pay_vol
    #motor 1 --> HTPB+N2O
    intermediate = 1 + (stage_motor_rho[0][3]/(stage_motor_OF[0][0]*stage_motor_rho[0][0]))
    stage_2_motor1_mHTPB = (stage_motor_rho[0][3]/stage_motor_OF[0][0])*(1/intermediate)*stage_2_vol
    stage_2_motor1_mN2O = stage_motor_rho[0][3]*(stage_2_vol - (stage_2_motor1_mHTPB/stage_motor_rho[0][0]))
    stage_2_motor1_mprop = stage_2_motor1_mHTPB + stage_2_motor1_mN2O
    #motor 2 --> HTPB+LOX
    intermediate = 1 + (stage_motor_rho[0][2]/(stage_motor_OF[0][1]*stage_motor_rho[0][0]))
    stage_2_motor2_mHTPB = (stage_motor_rho[0][2]/stage_motor_OF[0][1])*(1/intermediate)*stage_2_vol
    stage_2_motor2_mLOX = stage_motor_rho[0][2]*(stage_2_vol - (stage_2_motor2_mHTPB/stage_motor_rho[0][0]))
    stage_2_motor2_mprop = stage_2_motor2_mHTPB + stage_2_motor2_mLOX
    #motor 3 --> SP1-a+N20
    intermediate = 1 + (stage_motor_rho[0][3]/(stage_motor_OF[0][2]*stage_motor_rho[0][1]))
    stage_2_motor3_mSP1a = (stage_motor_rho[0][3]/stage_motor_OF[0][2])*(1/intermediate)*stage_2_vol
    stage_2_motor3_mN2O = stage_motor_rho[0][3]*(stage_2_vol - (stage_2_motor3_mSP1a/stage_motor_rho[0][1]))
    stage_2_motor3_mprop = stage_2_motor3_mSP1a + stage_2_motor3_mN2O
    #motor 4 --> SP1-a+LOX
    intermediate = 1 + (stage_motor_rho[0][2]/(stage_motor_OF[0][3]*stage_motor_rho[0][1]))
    stage_2_motor4_mSP1a = (stage_motor_rho[0][2]/stage_motor_OF[0][3])*(1/intermediate)*stage_2_vol
    stage_2_motor4_mLOX = stage_motor_rho[0][2]*(stage_2_vol - (stage_2_motor4_mSP1a/stage_motor_rho[0][1]))
    stage_2_motor4_mprop = stage_2_motor4_mSP1a + stage_2_motor4_mLOX
    
    stage_finert = [[0.100, 0.08]]  #inert mass fraction (solid first+second stage) (SPAD)
    missile_DV = np.ones((1,4))
    for i in range(len(stage_pay_m)):
        stage_2_motor1_minert = stage_finert[0][1]*(stage_2_motor1_mprop + stage_pay_m[i]) + stage_2_motor1_mprop + stage_pay_m[i]
        stage_2_motor2_minert = stage_finert[0][1]*(stage_2_motor2_mprop + stage_pay_m[i]) + stage_2_motor2_mprop + stage_pay_m[i]
        stage_2_motor3_minert = stage_finert[0][1]*(stage_2_motor3_mprop + stage_pay_m[i]) + stage_2_motor3_mprop + stage_pay_m[i]
        stage_2_motor4_minert = stage_finert[0][1]*(stage_2_motor4_mprop + stage_pay_m[i]) + stage_2_motor4_mprop + stage_pay_m[i]
        
        stage_2_DV = np.atleast_2d(np.empty((1,4)))
        mi_mf = (stage_pay_m[i] + stage_2_motor1_minert + stage_2_motor1_mprop)/(stage_pay_m[i])# + stage_2_motor1_minert)
        stage_2_DV[0][0] = (stage_motor_Isp[0][0]*(9.81*u.m/(u.s*u.s))*np.log(mi_mf)).value
        
        mi_mf = (stage_pay_m[i] + stage_2_motor2_minert + stage_2_motor2_mprop)/(stage_pay_m[i])# + stage_2_motor1_minert)
        stage_2_DV[0][1] = (stage_motor_Isp[0][1]*(9.81*u.m/(u.s*u.s))*np.log(mi_mf)).value
        
        mi_mf = (stage_pay_m[i] + stage_2_motor3_minert + stage_2_motor3_mprop)/(stage_pay_m[i])# + stage_2_motor1_minert)
        stage_2_DV[0][2] = (stage_motor_Isp[0][2]*(9.81*u.m/(u.s*u.s))*np.log(mi_mf)).value
        
        mi_mf = (stage_pay_m[i] + stage_2_motor4_minert + stage_2_motor4_mprop)/(stage_pay_m[i])# + stage_2_motor1_minert)
        stage_2_DV[0][3] = (stage_motor_Isp[0][3]*(9.81*u.m/(u.s*u.s))*np.log(mi_mf)).value
        
        stage_1_minert1 = stage_finert[0][0]*(stage_1_mprop + stage_pay_m[i] + stage_2_motor1_mprop + stage_2_motor1_minert) + stage_1_mprop + stage_pay_m[i] + stage_2_motor1_mprop + stage_2_motor1_minert
        stage_1_minert2 = stage_finert[0][0]*(stage_1_mprop + stage_pay_m[i] + stage_2_motor2_mprop + stage_2_motor2_minert) + stage_1_mprop + stage_pay_m[i] + stage_2_motor2_mprop + stage_2_motor2_minert
        stage_1_minert3 = stage_finert[0][0]*(stage_1_mprop + stage_pay_m[i] + stage_2_motor3_mprop + stage_2_motor3_minert) + stage_1_mprop + stage_pay_m[i] + stage_2_motor3_mprop + stage_2_motor3_minert
        stage_1_minert4 = stage_finert[0][0]*(stage_1_mprop + stage_pay_m[i] + stage_2_motor4_mprop + stage_2_motor4_minert) + stage_1_mprop + stage_pay_m[i] + stage_2_motor4_mprop + stage_2_motor4_minert
        
        stage_1_DV = np.atleast_2d(np.empty((1,4)))
        mi_mf = (stage_1_mprop + stage_pay_m[i] + stage_2_motor1_mprop + stage_2_motor1_minert + stage_1_minert1)/(stage_pay_m[i] + stage_2_motor1_minert + stage_2_motor1_mprop)
        stage_1_DV[0][0] = (stage_motor_Isp[0][0]*(9.81*u.m/(u.s*u.s))*np.log(mi_mf)).value
        
        mi_mf = (stage_1_mprop + stage_pay_m[i] + stage_2_motor2_mprop + stage_2_motor2_minert + stage_1_minert2)/(stage_pay_m[i] + stage_2_motor2_minert + stage_2_motor2_mprop)
        stage_1_DV[0][1] = (stage_motor_Isp[0][0]*(9.81*u.m/(u.s*u.s))*np.log(mi_mf)).value
        
        mi_mf = (stage_1_mprop + stage_pay_m[i] + stage_2_motor3_mprop + stage_2_motor3_minert + stage_1_minert3)/(stage_pay_m[i] + stage_2_motor3_minert + stage_2_motor3_mprop)
        stage_1_DV[0][2] = (stage_motor_Isp[0][0]*(9.81*u.m/(u.s*u.s))*np.log(mi_mf)).value
        
        mi_mf = (stage_1_mprop + stage_pay_m[i] + stage_2_motor4_mprop + stage_2_motor4_minert + stage_1_minert4)/(stage_pay_m[i] + stage_2_motor4_minert + stage_2_motor4_mprop)
        stage_1_DV[0][3] = (stage_motor_Isp[0][0]*(9.81*u.m/(u.s*u.s))*np.log(mi_mf)).value
        
        holder = (stage_1_DV + stage_2_DV)
        missile_DV = np.concatenate((missile_DV, holder), axis=0)
    
    missile_DV = np.delete(missile_DV, 0, 0)
    fig, ax = plt.subplots()
    motor11 = ax.plot(stage_pay_m, missile_DV[:,0], 'b--', label='F: HTPB O: N$_2$O')
    motor12 = ax.plot(stage_pay_m, missile_DV[:,1], 'ro-', label='F: HTPB O: LOX')
    motor13 = ax.plot(stage_pay_m, missile_DV[:,2], 'k*-', label='F: SP1$_a$ O:N$_2$O')
    motor14 = ax.plot(stage_pay_m, missile_DV[:,3], 'g+-', label='F: SP1$_a$ O:LOX')
    
    plt.xlabel('$m$ [g]')
    plt.ylabel('$\Delta$V [m/s]')
    plt.title('$\Delta$V of TACAIR for different 2$^{nd}$ Stages over payload mass')
    ax.legend()
    plt.show()
    return missile_DV, stage_pay_m

