""" Common type definitions and functions for including
in model modules.

"""

## Imports and type definitions
from numpy import double as np_double
from numpy cimport double_t as np_double_t
# ctypedef double real_t

## Constants
DEF PI = 3.1415926535897932384626433832
DEF RHO_WATER = 1e3 # kg/m3
DEF RHO_AIR = 1.0 
DEF THIRD = 1./3.
DEF MULTI_THRESH = 1e4

## Inline alias functions
cdef inline double dmax(double a, double b) nogil: return a if a >= b else b
cdef inline double dmin(double a, double b) nogil: return a if a <= b else b
cdef inline long lmax(long a, long b) nogil: return a if a >= b else b
