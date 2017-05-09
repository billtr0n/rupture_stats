# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from pylab import *

n = 1025
# computed from optimization algorithm
fd = fromfile('test_fdmin_output.bin', '<f4').reshape([n,n])

# analytical solution
r = fromfile('rastr.bin', '<f4').reshape([n,n])


figure()
imshow(fd)
colorbar()
savefig('rastr.pdf')

figure()
imshow(fd/r)
colorbar()
savefig('fdmin_ratio.pdf')

show()