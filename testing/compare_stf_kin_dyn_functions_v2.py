#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 25 13:52:54 2017
 
@author: william
"""
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

def yoffe_singularity_integral(Tr, dt):
    # integral evaluated numerically
    nu=2*(np.sqrt(dt*(-dt+Tr))+Tr*np.arcsin(np.sqrt(dt/Tr)))
    de=dt*np.pi*Tr 
    return nu/de

# yoffe function
def yoffe(t, Tp, Tr, t0):
    dt = np.real(t[1])-np.real(t[0])
    num = np.sqrt(t0 + Tr - Tp - t)
    deno = np.sqrt(t-t0)
    sft = np.real(num/deno).astype(np.float64)
    # how to deal with this infinity?
    val = yoffe_singularity_integral(Tr, dt)
    sft[np.where(sft == np.inf)] = val
    return sft
 
def heaviside(t):
    return 0.5 * (np.sign(t) + 1)
 
def half_sin(t, Tp):
    return np.sin(np.pi*t/Tp)
 
# regularized yoffe function
def regularized_yoffe(t, Tp, Tr, t0):
    dt = t[1] - t[0]
    tc = t.astype(np.complex64)
     
    # compute individual components
    srf = yoffe(tc, Tp, Tr, t0)
    hs = half_sin(t, Tp)*heaviside(Tp-t)
     
    # remove zero-padding from half-sine function
    n = int(np.ceil(Tp/dt))
    hs = hs[:n]
    stf = sp.signal.fftconvolve(srf, hs, mode='full')
    stf = stf[:len(srf)]
     
    plt.figure()
    plt.plot(stf)
    plt.show()
    # remove shift in convolution
#    stf = stf[:-n]
     
    # normalize so that integral = 1
    stf = stf/np.trapz(stf, dx=dt)
 
    # guarantee that t0 is the same in both cases
    t0_est = dt*(np.where(stf > 0.001)[0]).min()
    if t0 != t0_est:
        d = int((t0 - t0_est)/dt)
        stf = np.roll(stf, d)
        t0_est = dt*(np.where(stf > 0.001)[0]).min() 
    return stf
 
def i1(t,Tr):
    p1 = np.pi*(4*t - Tr)*Tr + 8*t**1.5*np.sqrt(-t + Tr)
    p2 = 4*Tr*np.sqrt(t*(-t + Tr)) 
    p3 = - 2*(4*t - Tr)*Tr*np.arctan((-2*t + Tr)/(2.*np.sqrt(t*(-t + Tr))))
    return (p1+p2+p3)/16.0

def i2(t,Tr,Ts):
    p1 = (Tr*np.sqrt(t*(-t + Tr)) + 2*np.sqrt(t**3*(-t + Tr)) - t*Tr*np.sqrt((-t + Tr + Ts)/(t - Ts)) + 
          Tr*Ts*np.sqrt((-t + Tr + Ts)/(t - Ts)) - 2*t*np.sqrt((t - Ts)*(-t + Tr + Ts)) - 
          2*Ts*np.sqrt((t - Ts)*(-t + Tr + Ts)) + (4*t - Tr)*Tr*np.arcsin(np.sqrt(t/Tr)) + 
          Tr*(-4*t + Tr)*np.arcsin(np.sqrt((t - Ts)/Tr)))/4. 
          
    p2 = (-((2*t + Tr - 6*Ts)*np.sqrt((t - Ts)*(-t + Tr + Ts))) + 
          Tr*(-4*t + Tr + 8*Ts)*np.arcsin(np.sqrt((t - Ts)/Tr)))/4.
    return p1+p2

def i3(t,Tr,Ts):
    p1 = (Tr*np.sqrt(t*(-t + Tr)) + 2*np.sqrt(t**3*(-t + Tr)) - t*Tr*np.sqrt((-t + Tr + Ts)/(t - Ts)) + 
          Tr*Ts*np.sqrt((-t + Tr + Ts)/(t - Ts)) - 2*t*np.sqrt((t - Ts)*(-t + Tr + Ts)) - 
          2*Ts*np.sqrt((t - Ts)*(-t + Tr + Ts)) + (4*t - Tr)*Tr*np.arcsin(np.sqrt(t/Tr)) + 
          Tr*(-4*t + Tr)*np.arcsin(np.sqrt((t - Ts)/Tr)))/4. 
          
    p2 = (2*(-2*t*np.sqrt((t - Ts)*(-t + Tr + Ts)) - Tr*np.sqrt((t - Ts)*(-t + Tr + Ts)) + 
        6*Ts*np.sqrt((t - Ts)*(-t + Tr + Ts)) + 2*t*np.sqrt((t - 2*Ts)*(-t + Tr + 2*Ts)) + 
        Tr*np.sqrt((t - 2*Ts)*(-t + Tr + 2*Ts)) - 4*Ts*np.sqrt((t - 2*Ts)*(-t + Tr + 2*Ts))) - 
        Tr*(-4*t + Tr + 8*Ts)*np.arctan((-2*t + Tr + 2*Ts)/(2.*np.sqrt((t - Ts)*(-t + Tr + Ts)))) + 
        Tr*(-4*t + Tr + 8*Ts)*np.arctan((-2*t + Tr + 4*Ts)/(2.*np.sqrt((t - 2*Ts)*(-t + Tr + 2*Ts)))))/8.
    return p1+p2
          

def i4(t,Tr,Ts):
    p1 = t*Tr*np.arccos(np.sqrt((t - Ts)/Tr)) 
    
    p2 = (-(np.sqrt((t - Ts)*(-t + Tr + Ts))*(2*t + Tr + 2*Ts)) - 
        Tr**2*np.arccos(1.0/(np.sqrt(Tr/(t - Ts)))))/4. 
                        
    p3 = (2*(-2*t*np.sqrt((t - Ts)*(-t + Tr + Ts)) - Tr*np.sqrt((t - Ts)*(-t + Tr + Ts)) + 
        6*Ts*np.sqrt((t - Ts)*(-t + Tr + Ts)) + 2*t*np.sqrt((t - 2*Ts)*(-t + Tr + 2*Ts)) + 
        Tr*np.sqrt((t - 2*Ts)*(-t + Tr + 2*Ts)) - 4*Ts*np.sqrt((t - 2*Ts)*(-t + Tr + 2*Ts))) - 
        Tr*(-4*t + Tr + 8*Ts)*np.arctan((-2*t + Tr + 2*Ts)/(2.*np.sqrt((t - Ts)*(-t + Tr + Ts)))) + 
        Tr*(-4*t + Tr + 8*Ts)*np.arctan((-2*t + Tr + 4*Ts)/(2.*np.sqrt((t - 2*Ts)*(-t + Tr + 2*Ts)))))/8.
    return p1+p2+p3

def i5(t,Tr,Ts):
    p1 = (4*(2*t + Tr - 4*Ts)*np.sqrt(-((t - 2*Ts)*(t - Tr - 2*Ts))) + np.pi*Tr*(-4*t + Tr + 8*Ts) + 
    2*Tr*(-4*t + Tr + 8*Ts)*np.arctan((-2*t + Tr + 4*Ts)/(2.*np.sqrt((t - 2*Ts)*(-t + Tr + 2*Ts)))))/16.
    return p1

def tinti(t, Ts, Tr, t0):
    n = len(t)
    dt = t[1]-t[0]
    
    # from the derivation in mathematica
    k = 2 / (np.pi * Tr * Ts**2)
    
    stf = np.zeros(n)
    
    if Tr > 2*Ts:
        stf[(t >= 0) & (t <= Ts)] = i1(t[(t >= 0) & (t <= Ts)], Tr)
        stf[(t > Ts) & (t <= 2*Ts)] = i2(t[(t > Ts) & (t <= 2*Ts)], Tr, Ts)
        stf[(t >= 2*Ts) & (t < Tr)] = i3(t[(t >= 2*Ts) & (t < Tr)], Tr, Ts)
        stf[(t >= Tr) & (t < Tr+Ts)] = i4(t[(t >= Tr) & (t < Tr+Ts)], Tr, Ts)
        stf[(t >= Tr+Ts) & (t < Tr + 2*Ts)] = i5(t[(t >= Tr+Ts) & (t < Tr + 2*Ts)], Tr, Ts)
    
    d = int(np.ceil((t0)/dt))
    stf = np.roll(stf, d)
    return k*stf    


def triangle(t,Ts):
    """ assuming equal spaced time-series."""
    dt = t[1]-t[0]
    stf = np.zeros(len(t))
    stf[t < Ts] = t[t < Ts]
    stf[((t>=Ts)&(t<=2*Ts))] = 2*Ts - t[((t>=Ts)&(t<=2*Ts))]
#    stf = stf[:int(2*Ts/dt)+10]/np.trapz(stf[:int(2*Ts/dt)+10],dx=dt)
    stf = stf[:int(2*Ts/dt)]/np.trapz(stf[:int(2*Ts/dt)],dx=dt)
    return stf
 
# regularized yoffe function
def regularized_yoffe_triangle(t, Tp, Tr, t0):
    dt = t[1] - t[0]
    tc = t.astype(np.complex64)
     
    # compute individual components
    srf = yoffe(tc, Tp, Tr, t0)
    hs = triangle(t, Tp)
    
    stf = sp.signal.fftconvolve(srf, hs, mode='full')
    stf = stf[:len(srf)]
     
    # normalize so that integral = 1
    stf = stf/np.trapz(stf, dx=dt)
 
    # guarantee that t0 is the same in both cases
    t0_est = dt*(np.where(stf > 0.001)[0]).min()
    while t0 != t0_est:
        d = int((t0 - t0_est)/dt)
        stf = np.roll(stf, d)
        t0_est = dt*(np.where(stf > 0.001)[0]).min() 
    return stf
 
def compute_tarr( svm, dt ):
    """
    t1: first occurance where svm < 0.001 m/s
    t2: time after t0 where 95% of slip is constrained
    t3: entire rupture time (force retaining dynamic slip)
    t4: used for fitting tr during optimization.
    """
     
    t0_ind = int((np.where(svm > 0.001)[0]).min())
    t0 = dt * t0_ind
     
    # search from t0 + 1, return first occurence where svm is less than 0.001
    te1 = dt * (int((np.where(svm[t0_ind+1:] < 0.001)[0]).min()) + t0_ind)
     
    # from mai et al., 95% of slip is constrained
    svm_norm = svm / np.trapz(svm, dx=dt)
    cum_slip = sp.integrate.cumtrapz(svm_norm, dx=dt)
    te2 = dt * np.where(cum_slip >= 0.95)[0].min()
     
    # entire time
    te3 = dt * np.where(svm > 0.001)[0].max()
     
    # first occurance where psv reaches 2.5% of peak value, using for sliding
    psv = svm.max()
    psv_ind = np.argmax(svm)
    te4 = dt * (np.where(svm[psv_ind+1:] < 0.025*psv)[0].min() + psv_ind)
     
    return (te1-t0,te2-t0,te3-t0,te4-t0)
     
def exponential_stf( time, slip, psv, trup ):
    """accepts a time index and returns the slip-rate for specified sub-fault parameters."""
    tpeak = slip / (np.exp( 1 ) * psv)
    stf = slip / tpeak * ((time - trup) / tpeak) * np.exp( -(time-trup)/tpeak )
    stf[np.where(time<trup)] = 0.0
    return stf
 
def objective_function_slip_psv(x, *params):
    t, svm = params
    tp, tr, slip = x
    dt = t[1]-t[0]
     
    # compute t0 first, as that parameter does not need to be fitted.
    t0 = dt*(np.where(svm > 0.001)[0]).min()
     
    t0_ind = int((np.where(svm > 0.001)[0]).min())
    tend_ind = int((np.where(svm > 0.001)[0]).max())
     
    test_function = regularized_yoffe(t, tp, tr, t0)*slip
     
    psv = svm[t0_ind:tend_ind].max()
    slip = np.trapz(svm[t0_ind:tend_ind], dx=dt)
    psv_t = test_function[t0_ind:tend_ind].max()
    slip_t = np.trapz(test_function[t0_ind:tend_ind], dx=dt)
     
    err = 0.5*(psv-psv_t)**2+0.5*(slip-slip_t)**2
     
#    res = np.sqrt(np.mean((test_function[t0_ind:tend_ind] - svm[t0_ind:tend_ind])**2))
#    res = np.sum(np.square(test_function[t0_ind:tend_ind] - svm[t0_ind:tend_ind]))
 
    return err
 
def objective_function(x, *params):
    """
    Parameters
    ==========
    x[0] = Tp
    x[1] = Tr
    x[2] = Slip
    """
    t, svm = params
    tp, tr, slip = x
    dt = t[1]-t[0]
    # compute t0 first, as that parameter does not need to be fitted.
    t0 = dt*(np.where(svm > 0.001)[0]).min()
     
    t0_ind = int((np.where(svm > 0.001)[0]).min())
    tend_ind = int((np.where(svm > 0.001)[0]).max())
     
    test_function = regularized_yoffe(t, tp, tr, t0)*slip
     
#    res = np.sqrt(np.mean((test_function[t0_ind:tend_ind] - svm[t0_ind:tend_ind])**2))
    res = np.sum(np.square(test_function[t0_ind:tend_ind] - svm[t0_ind:tend_ind]))
 
    return res
 
def objective_function_windowed_fit_tp(x, *params):
    """
    Parameters
    ==========
    x[0] = Tp
    x[1] = Tr
    x[2] = Slip
    """
#    t, svm, tr = params
#    tp, slip = x
    t, svm, tr, slip = params
    tp = x[0]
    dt = t[1]-t[0]
     
    # compute t0 first
    t0_ind = int((np.where(svm > 0.001)[0]).min())
    t0 = dt * t0_ind
     
    # compute index of ending point
    tend_ind = int((t0+tr)/dt)
     
    # compute test function
    test_function = regularized_yoffe(t, tp, tr, t0)*slip
     
    # return the square of the errors
    sqerr = np.sum(np.square(test_function[t0_ind:tend_ind] - svm[t0_ind:tend_ind]))
 
    return sqerr
 
def objective_function_windowed_fit_tp_and_tr(x, *params):
    """
    Parameters
    ==========
    x[0] = Tp
    x[1] = Tr
    x[2] = Slip
    """
#    t, svm, tr = params
#    tp, slip = x
    t, svm, slip = params
    tp, tr = x
    dt = t[1]-t[0]
     
    # compute t0 first
    t0_ind = int((np.where(svm > 0.001)[0]).min())
    t0 = dt * t0_ind
     
    # compute index of ending point
    tend_ind = int((t0+tr)/dt)
     
    # compute test function
    test_function = regularized_yoffe(t, tp, tr, t0)*slip
     
    # return the square of the errors
    sqerr = np.sum(np.square(test_function[t0_ind:tend_ind] - svm[t0_ind:tend_ind]))
 
    return sqerr
 
def objective_function_windowed_fit_tp_and_tr_and_slip(x, *params):
    """
    Parameters
    ==========
    x[0] = Tp
    x[1] = Tr
    x[2] = Slip
    """
#    t, svm, tr = params
#    tp, slip = x
    t, svm = params
    tp, tr, slip = x
    dt = t[1]-t[0]
     
    # compute t0 first
    t0_ind = int((np.where(svm > 0.001)[0]).min())
    t0 = dt * t0_ind
     
    # compute index of ending point
    tend_ind = int((t0+tr)/dt)
     
    # compute test function
    test_function = regularized_yoffe(t, tp, tr, t0)*slip
     
    # return the square of the errors
    sqerr = np.sum(np.square(test_function[t0_ind:tend_ind] - svm[t0_ind:tend_ind]))
 
    return sqerr
 
def objective_function_windowed_fit_tp_and_slip(x, *params):
    """
    Parameters
    ==========
    x[0] = Tp
    x[1] = Tr
    x[2] = Slip
    """
#    t, svm, tr = params
#    tp, slip = x
    t, svm, tr = params
    tp, slip = x
    dt = t[1]-t[0]
     
    # compute t0 first
    t0_ind = int((np.where(svm > 0.001)[0]).min())
    t0 = dt * t0_ind
     
    # compute index of ending point
    tend_ind = int((t0+tr)/dt)
     
    # compute test function
    test_function = regularized_yoffe(t, tp, tr, t0)*slip
     
    # return the square of the errors
    sqerr = np.sum(np.square(test_function[t0_ind:tend_ind] - svm[t0_ind:tend_ind]))
 
    return sqerr#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 16:15:48 2017

@author: william
"""

