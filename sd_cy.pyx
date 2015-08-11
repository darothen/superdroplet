#cython: cdivision=True
#cython: nonecheck=False
#cython: boundscheck=False
#cython: wraparound=False

import numpy as np
cimport numpy as cnp

from libc.stdlib cimport rand, RAND_MAX
from libc.math cimport abs, sqrt

from cpython.object cimport Py_LT, Py_LE, Py_EQ, Py_GE, Py_GT, Py_NE

cdef extern from "math.h":
    double floor(double x)

cdef extern from "stdlib.h":
    void qsort(void *array, size_t count, size_t size,
                int (*compare)(const void *, const void *))
    void *malloc(size_t size) 
    void free(void *ptr)

cdef inline int ifloor(double x): return int(floor(x))
cdef inline double dmin(double a, double b): return a if a <= b else b

DEF VERBOSITY = 1
DEF RHO_WATER = 1e3
DEF RHO_AIR = 1.0
DEF PI = 3.1415926535897932384626433832
DEF MULTI_THRESH = 1e4
cdef double RAND_FACT = float(RAND_MAX)

cdef int superdroplet_count

ctypedef Superdroplet Superdroplet_t
cdef class Superdroplet:

    cdef public int multi
    cdef public double rcubed, solute
    cdef public double density
    cdef public int id

    def __init__(self, int multi, double rcubed, double solute):
        global superdroplet_count
        superdroplet_count += 1

        self.multi  = multi
        self.rcubed = rcubed
        self.solute = solute

        self.density = RHO_WATER
        self.id = superdroplet_count

    def get_terminal_v(self):
        return self._get_terminal_v()
    cdef double _get_terminal_v(self):        
        cdef double diameter, g, C_D, t_v

        diameter = 2*self.rcubed**(1./3.)
        g   = 9.8 # gravitational acceleration, m/s
        C_D = 0.6 # drag coefficient, unitless
        tv  = (4./3.)*(g*diameter/C_D)*(self.density/RHO_AIR)
        tv  = sqrt(tv)
        return tv

    def get_volume(self):
        return self._get_volume()
    cdef double _get_volume(self):
        return self.rcubed*4.*PI/3.

    def get_mass(self):
        return self._get_mass()
    cdef double _get_mass(self):
        return self.density*self._get_volume()

    def attributes(self):
        return { 
            'multi': self.multi, 
            'rcubed': self.rcubed,
            'solute': self.solute,
        }

    def copy(self):
        return Superdroplet(**self.attributes())

    def __repr__(self):
        return "%d" % self.id

cdef int sd_compare(Superdroplet_t a, Superdroplet_t b):
    return a.multi - b.multi

## Collision Kernel
cdef double golovin(Superdroplet sd_j, Superdroplet sd_k):
    cdef double b = 1.5e3
    cdef double rj3 = sd_j.rcubed, rk3 = sd_k.rcubed
    return b*(rj3 + rk3)*4.*PI/3.

cdef double hydro(Superdroplet sd_j, Superdroplet sd_k):
    cdef double E = 1.0
    cdef double p, r_j, r_k, tv_j, tv_k
    p = 1./3.
    r_j = sd_j.rcubed**p
    r_k = sd_k.rcubed**p
    tv_j = sd_j.get_terminal_v()
    tv_k = sd_k.get_terminal_v()
    return E*PI*((r_j + r_k)**2)*abs(tv_j - tv_k)

cpdef double detect_collision(Superdroplet_t sd_a, Superdroplet_t sd_b, 
                              double scaling, double t_c, double delta_V):
    cdef double phi, p_alpha, gamma

    #phi = c_libc_rand()
    #phi = np.random.uniform()
    phi = rand() / RAND_FACT
    p_alpha = prob_collision(sd_a, sd_b, 
                             scaling=scaling,t_c=t_c,delta_V=delta_V)
    if phi < p_alpha - floor(p_alpha):
        gamma = floor(p_alpha) + 1
    else:
        gamma = floor(p_alpha)

    return gamma

cdef double prob_collision(Superdroplet sd_j, Superdroplet sd_k, 
                           double scaling, double t_c, double delta_V):
    
    cdef int xi_j, xi_k, max_xi
    cdef double K_ij

    xi_j = sd_j.multi
    xi_k = sd_k.multi

    K_ij = hydro(sd_j, sd_k)

    #if kernel == 'golovin':
    #    K_ij = golovin(sd_j, sd_k)
    #elif kernel == 'hydro':
    #    K_ij = hydro(sd_j, sd_k)
    #else:
    #    raise ValueError("Undefined collision kernel (%s)" % kernel)
    if xi_j > xi_k:
        max_xi = xi_j
    else:
        max_xi = xi_k

    return scaling*max_xi*(t_c/delta_V)*K_ij

cdef list multi_coalesce(Superdroplet sd_j, Superdroplet sd_k, double gamma):
    """
    Coalesce two superdroplets with one another.
    """

    cdef Superdroplet sd_temp, sd_recycle
    cdef double gamma_tilde
    cdef int multi_j_p, multi_k_p, excess
    cdef double solute_j_p, solute_k_p, rcubed_j_p, rcubed_k_p

    if sd_j.multi < sd_k.multi:
        sd_temp = sd_j.copy()
        sd_j = sd_k.copy()
        sd_k = sd_temp

    # gamma_tilde = dmin(gamma, floor(sd_j.multi/sd_k.multi))
    gamma_tilde = 1.
    excess = sd_j.multi - int(gamma_tilde*sd_k.multi)

    # print "   ", gamma, gamma_tilde, excess

    if excess > 0:

        multi_j_p = excess # = sd_j.multi - int(gamma_tilde*sd_k.multi)
        multi_k_p = int(gamma_tilde*sd_k.multi)

        rcubed_j_p = sd_j.rcubed
        rcubed_k_p = gamma_tilde*sd_j.rcubed + sd_k.rcubed

        solute_j_p = sd_j.solute
        solute_k_p = gamma_tilde*sd_j.solute + sd_k.solute

        return [
            Superdroplet(multi_j_p, rcubed_j_p, solute_j_p),
            Superdroplet(multi_k_p, rcubed_k_p, solute_k_p)
        ]

    else:

        multi_j_p = ifloor(sd_k.multi/2)
        multi_k_p = sd_k.multi - multi_j_p

        sd_temp = Superdroplet(multi_k_p, 
                               gamma_tilde*sd_j.rcubed + sd_k.rcubed,
                               gamma_tilde*sd_j.solute + sd_k.solute)
        sd_recycle = Superdroplet(0, 
                                  1., 1.)

        return [ sd_temp, sd_recycle ]

def recycle(list sds):
    """ For a list of superdroplets, identify which ones have 0 
    multiplicities; for each *i* of these, pick the superdroplet 
    with the *i*th-most multiplicity, split it in half, and copy 
    the properties to the remaining superdroplet. """

    print "RECYCLE", 

    sds.sort(cmp=sd_compare)

    for i, sd in enumerate(sds):
        # Short circuits:
        # 1) Did we already encounter non-zero superdroplets?
        if sd.multi > 0: break

        sd_donor = sds[-i]
        # 2) Does the donor superdroplet have data to spare?
        if sd_donor.multi <= 0: break

        sd.multi = ifloor(sd_donor.multi/2)
        sd_donor.multi -= ifloor(sd_donor.multi/2)

        sd.rcubed = sd_donor.rcubed
        sd.solute = sd_donor.solute
    print " %d superdroplets" % i

    return sds

def step(list sd_list, double t_c, double delta_V):
    cdef list diag_msgs = []

    print "PREP STEPS"
    # 1) Make the random permutation of the super-droplet list
    print "   SHUFFLE LIST"
    np.random.shuffle(sd_list)

    # 2) Make the candidate pairs
    print "   GEN PAIRS"
    cdef int n_part = len(sd_list)

    # 3) Generate the uniform random numbers
    print "PROBABILITY LOOP"
    cdef double scaling = (n_part*(n_part - 1)/2.)/ifloor(n_part/2)    


    print "PROB / COLLISION LOOP"
    cdef list output = []
    cdef bint collisions = False
    cdef int counter = 0
    cdef:
        unsigned int i
        Superdroplet_t sd_j, sd_k
        double gamma
        list new_pair
        double phi, p_alpha
        int xi_j, xi_k, max_xi
        double K_ij

    for i in xrange(n_part/2):
        sd_j, sd_k = sd_list[i], sd_list[i + n_part/2]

        phi = rand() / RAND_FACT
        xi_j = sd_j.multi
        xi_k = sd_k.multi

        K_ij = hydro(sd_j, sd_k)
        if xi_j > xi_k:
            max_xi = xi_j
        else:
            max_xi = xi_k
        prob = scaling*max_xi*(t_c/delta_V)*K_ij

        if prob < phi: # no collision!
            output.extend([sd_j, sd_k])

        else:
            # print sd_j, sd_k, prob
    
            new_pair = multi_coalesce(sd_j, sd_k, gamma)
            output.extend(new_pair)

            if not collisions:
                collisions = True
            counter += 1

    if collisions:
        diag_msgs.append("%5d collisions simulated" % counter)

    return output, diag_msgs

