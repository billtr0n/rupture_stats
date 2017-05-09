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
rake2 = zeros([801,2601])
for i in range(200,300):
    print i
    for j in range(1000,1500):
        n = matrix([nhat1[i,j], nhat2[i,j], nhat3[i,j]]).T
        u = matrix([su1[i,j], su2[i,j], su3[i,j]]).T

        if isclose(norm(u), 0): 
            rake[i,j] = nan

        else:
            # vector along strike
            m = (identity(3) - n*n.T)[:,0]

            # compute angle of rake
            argc = (u.T*m/(norm(u)*norm(m))).A1
            args = norm(cross(u.A1, m.A1))/(norm(u)*norm(m))

            rake2[i,j] = rad2deg(arctan2(args, argc))
            
        




