#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May  2 09:55:57 2017

@author: william
"""

from pylab import *

n = 1025
r = fromfile('rastr.bin', '<f4').reshape([n,n])
print r.min()
print unravel_index(argmin(r), [n,n])

figure()
imshow(r)
show()