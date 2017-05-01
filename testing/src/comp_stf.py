#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 10:11:13 2017

@author: william
"""

import numpy as np
import matplotlib.pyplot as plt

import compare_stf_kin_dyn_functions_v2_test as f

Ts = 1.1
Tr = 2.0
t0 = 1.0

t = np.arange(0,10.0,0.001)
stf_c1 = np.loadtxt('case1.txt')
stf_c2 = np.loadtxt('case2.txt')
#stf_an = f.tinti(t, Ts, Tr, t0)
Ts = 0.1
stf_nu1 = f.regularized_yoffe_triangle(t, Ts, Tr, t0)
Ts = 1.1
stf_nu2 = f.regularized_yoffe_triangle(t, Ts, Tr, t0)
plt.figure()
plt.plot(t, stf_c1, label='(C) : Case 1')
plt.plot(t, stf_nu1,'r--', label='(Py) : Case 1')
plt.plot(t, stf_c2, label='(C) : Case 2')
plt.plot(t, stf_nu2,'r--', label='(Py) : Case 2')
plt.xlabel('Time [s]')
plt.ylabel('Slip Velocity [m/s]')
plt.legend()
plt.show()
