#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Predict Orbit Achieved From Different Second Stages on Launch Vehicle"""
import numpy as np
from poliastro.bodies import Earth
from scipy.optimize import newton
#custom

"""Python 3.7
   EH Group, Inc. (c) 2019
   Christopher Simpson: christopher.r.simpson@simpsonaerospace.com"""
#------------------------------------------------------------------------------
def gammaAdjust(X, gamma, t):
    def fder(yadjust, X, gamma, t):
       Vi1 = np.sqrt(X[t - 1,2]**2 + X[t - 1,3]**2)
       Vi = np.sqrt(X[t,2]**2 + X[t,3]**2)
       return -2*(Earth.k.value/(Earth.R_mean.value + X[t,1])**2) * (np.cos(gamma)/(Vi1*Vi))*yadjust*t
    def func(yadjust, X, gamma, t):
        Vi1 = np.sqrt(X[t - 1,2]**2 + X[t - 1,3]**2)
        Vi = np.sqrt(X[t,2]**2 + X[t,3]**2)
        return ((-2*(Earth.k.value/(Earth.R_mean.value + X[t,1])**2) * (np.cos(gamma)/(Vi1*Vi))*t*yadjust) + gamma)
    
    yadjust = newton(func, 10.2, fprime=fder, args=(X, gamma, t))
    
    return yadjust


