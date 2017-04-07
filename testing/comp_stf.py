#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 10:11:13 2017

@author: william
"""

import numpy as np
import matplotlib.pyplot as plt

import compare_stf_kin_dyn_functions_v2 as f

Ts = 0.1
Tr = 2.0
t0 = 1.0

t = np.arange(0,10.0,0.001)
stf_c = np.loadtxt('stf_test9.txt')
stf_an = f.tinti(t, Ts, Tr, t0)
plt.plot(t, stf_c, label='From C')
plt.plot(t, stf_an,'r--', label='From Python')
plt.xlabel('Time [s]')
plt.ylabel('Slip Velocity [m/s]')
plt.legend()
plt.show()
