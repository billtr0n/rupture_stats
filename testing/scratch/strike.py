from pylab import *

# fault normal vectors
nhat1 = fromfile("nhat1", ">f4").reshape([801,2601])
nhat2 = fromfile("nhat2", ">f4").reshape([801,2601])
nhat3 = fromfile("nhat3", ">f4").reshape([801,2601])

# slip vectors
su1 = fromfile("su1", ">f4").reshape([801,2601])
su2 = fromfile("su2", ">f4").reshape([801,2601])
su3 = fromfile("su3", ">f4").reshape([801,2601])

norm = lambda v: sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])

c = 0
strike = zeros([801,2601])
for i in range(200,300):
    print i
    for j in range(1000,1500):
        nproj = (nhat1[i,j], 0, nhat3[i,j])
        scaling = 1.0 / norm(nproj)
        strike[i,j] = rad2deg(scaling * arccos( nproj[0] )) - 90






