#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May  2 08:53:21 2017

@author: william
"""
from numpy import cos, pi, arange, column_stack, matrix
from matplotlib.pyplot import figure, plot, show, imshow

def rastr(xx):
    d = len(xx)
    sum = 0;
    for ii in range(d):
        xi = xx[ii]
        sum = sum + (xi**2 - 10*cos(2*pi*xi))
        
    y = 10*d + sum
    return y


x1 = arange(-5.12,5.12,0.01)
x2 = arange(-5.12,5.12,0.01)
xx = column_stack([x1,x2])

y = rastr(xx.T)
print y.shape
yy = matrix(y).T*matrix(y)
print yy.shape

figure()
imshow(yy)
show()