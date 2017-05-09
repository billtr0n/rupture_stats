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

def get_rake(su, nhat):
    n = matrix([nhat[0], nhat[1], nhat[2]]).T
    u = matrix([su[0], su[1], su[2]]).T

    if isclose(norm(u), 0): 
        rake[i,j] = nan

    else:
        # vector along strike
        m = (identity(3) - n*n.T)[:,0]

        # compute angle of rake
        arg = u.T*m/(norm(u)*norm(m)).A1

        if isclose(arg, 1):
            arg = 1
        rake = rad2deg(arccos(arg))
    return rake


print get_rake([1,0,0],[0,0,1])




