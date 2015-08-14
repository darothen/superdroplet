#cython: cdivision=False
#cython: nonecheck=False
#cython: boundscheck=False
#cython: wraparound=False

cimport cython
import numpy as np
cimport numpy as np

from numpy import double as np_double
from numpy cimport double_t as np_double_t

cdef inline double dmin(double a, double b): return a if a <= b else b

effic_array = np.array(
    [0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,
     0.001,0.001,0.001,0.001,0.001,0.003,0.003,0.003,0.004,0.005,
     0.005,0.005,0.010,0.100,0.050,0.200,0.500,0.770,0.870,0.970,
     0.007,0.007,0.007,0.008,0.009,0.010,0.010,0.070,0.400,0.430,
     0.580,0.790,0.930,0.960,1.000,0.009,0.009,0.009,0.012,0.015,
     0.010,0.020,0.280,0.600,0.640,0.750,0.910,0.970,0.980,1.000,
     0.014,0.014,0.014,0.015,0.016,0.030,0.060,0.500,0.700,0.770,
     0.840,0.950,0.970,1.000,1.000,0.017,0.017,0.017,0.020,0.022,
     0.060,0.100,0.620,0.780,0.840,0.880,0.950,1.000,1.000,1.000,
     0.030,0.030,0.024,0.022,0.032,0.062,0.200,0.680,0.830,0.870,
     0.900,0.950,1.000,1.000,1.000,0.025,0.025,0.025,0.036,0.043,
     0.130,0.270,0.740,0.860,0.890,0.920,1.000,1.000,1.000,1.000,
     0.027,0.027,0.027,0.040,0.052,0.200,0.400,0.780,0.880,0.900,
     0.940,1.000,1.000,1.000,1.000,0.030,0.030,0.030,0.047,0.064,
     0.250,0.500,0.800,0.900,0.910,0.950,1.000,1.000,1.000,1.000,
     0.040,0.040,0.033,0.037,0.068,0.240,0.550,0.800,0.900,0.910,
     0.950,1.000,1.000,1.000,1.000,0.035,0.035,0.035,0.055,0.079,
     0.290,0.580,0.800,0.900,0.910,0.950,1.000,1.000,1.000,1.000,
     0.037,0.037,0.037,0.062,0.082,0.290,0.590,0.780,0.900,0.910,
     0.950,1.000,1.000,1.000,1.000,0.037,0.037,0.037,0.060,0.080,
     0.290,0.580,0.770,0.890,0.910,0.950,1.000,1.000,1.000,1.000,
     0.037,0.037,0.037,0.041,0.075,0.250,0.540,0.760,0.880,0.920,
     0.950,1.000,1.000,1.000,1.000,0.037,0.037,0.037,0.052,0.067,
     0.250,0.510,0.770,0.880,0.930,0.970,1.000,1.000,1.000,1.000,
     0.037,0.037,0.037,0.047,0.057,0.250,0.490,0.770,0.890,0.950,
     1.000,1.000,1.000,1.000,1.000,0.036,0.036,0.036,0.042,0.048,
     0.230,0.470,0.780,0.920,1.000,1.020,1.020,1.020,1.020,1.020,
     0.040,0.040,0.035,0.033,0.040,0.112,0.450,0.790,1.010,1.030,
     1.040,1.040,1.040,1.040,1.040,0.033,0.033,0.033,0.033,0.033,
     0.119,0.470,0.950,1.300,1.700,2.300,2.300,2.300,2.300,2.300,
     0.027,0.027,0.027,0.027,0.027,0.125,0.520,1.400,2.300,3.000,
     4.000,4.000,4.000,4.000,4.000], dtype=np_double)
r0_arr = np.array([6., 8., 10., 15., 20., 25., 30., 40., 50.,
               60., 70., 100., 150., 200., 300.], dtype=np_double)
rat_arr = np.array([0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,
                0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0], 
               dtype=np_double)

cdef int n_r0 = 15
cdef int n_rat = 21

assert len(r0_arr) == n_r0
assert len(rat_arr) == n_rat

effic_array = effic_array.reshape((n_r0, n_rat), order='F')

# Memory view into efficiency lookup table
cdef double[:,:] effic = effic_array
cdef double[:] rat = rat_arr
cdef double[:] r0 = r0_arr

cdef double hall_effic(double r_coll, double r_small):

    cdef:
        unsigned int k, ir, kk, irat
        double ratio
        double p, q, ec, ek
    
    ## Seek for row in effic table
    for k in xrange(1, n_r0):
        if ((r_coll <= r0[n_r0-1]) and (r_coll >= r0[k-1])):
            ir = k
        elif (r_coll > r0[n_r0-1]):
            ir = n_r0
        elif (r_coll < r0[0]):
            ir = 0
    
    ## Seek column; simpler because we have an invariant that 
    ## 0 < ratio < 1
    ratio = r_small / r_coll
    for kk in xrange(1, n_rat):
        if ((ratio >= rat[kk-1]) and (ratio < rat[kk])):
            irat = kk
            break
            
    # print ir, irat
    
    ## Interpolation logic - 4 point, linear weighting
    if ir < n_r0:
        if ir >= 1:
            p = (r_coll - r0[ir-1]) / (r0[ir] - r0[ir-1])
            q = (ratio - rat[irat-1]) / (rat[irat] - rat[irat-1])
            ec = ( (1. - p)*(1. - q)*effic[ir-1, irat-1] + 
                          p*(1. - q)*effic[ir, irat-1]   + 
                          q*(1. - p)*effic[ir-1, irat]   + 
                                  p*q*effic[ir, irat]      )
        else:
            q = (ratio - rat[irat-1]) / (rat[irat] - rat[irat-1])
            ec = (1. - q)*effic[0, irat-1] + q*effic[0, irat]
    else:
        q = (ratio - rat[irat-1]) / (rat[irat] - rat[irat-1])
        ek = (1. - q)*effic[n_r0-1, irat-1] + q*effic[n_r0-1, irat]
        ec = dmin(ek, 1.0)
        
    return ec