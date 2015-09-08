#cython: cdivision=True
#cython: nonecheck=False
#cython: boundscheck=False
#cython: wraparound=False
#cython: profile=True
#cython: linetrace=True

STUFF = "Hi"

cimport cython

import numpy as np
cimport numpy as cnp

from libc.stdlib cimport RAND_MAX

from hall cimport hall_effic

cdef extern from "math.h" nogil:
    double floor(double x)
    double sqrt(double x)
    double fabs(double x)

cdef extern from "stdlib.h" nogil:
    double rand()

## Inline alias functions
cdef inline double dmax(double a, double b): return a if a >= b else b
cdef inline double dmin(double a, double b): return a if a <= b else b
cdef inline long lmax(long a, long b): return a if a >= b else b

ctypedef kernel_id kernel_id_t
cpdef enum kernel_id:
    GOLOVIN = 1 
    LONG    = 2 
    HYDRO   = 3 
    HALL    = 4    

# DEF KERNEL_ID = 3 # 1 = Golovin, 
#                   # 2 = Hydro w/ E_coll=1
#                   # 3 = Hydro w/ Long collection
#                   # 4 = Hydro w/ Hall efficiencies (from Bott)
DEF VERBOSITY = 1
#: Density of water in kg/m3
DEF RHO_WATER = 1e3 
DEF RHO_AIR = 1.0 
DEF THIRD = 1./3.
DEF PI = 3.1415926535897932384626433832
DEF MULTI_THRESH = 1e4
cdef double RAND_FACT = float(RAND_MAX)

cdef int superdroplet_count

cdef class Superdroplet:

    cdef readonly long multi
    cdef readonly double rcubed, solute, density
    cdef readonly int id

    def __init__(self, long multi, double rcubed, double solute):
        global superdroplet_count
        superdroplet_count += 1

        self.multi  = multi
        self.rcubed = rcubed
        self.solute = solute
        self.density = RHO_WATER # kg/m3
        self.id = superdroplet_count

    cdef void set_properties(self, multi, rcubed, solute):
        self.multi = multi
        self.rcubed = rcubed
        self.solute = solute

    property terminal_v:
        def __get__(self):      
            return self.calc_terminal_v()
    cdef double calc_terminal_v(self):
        ## Physical formula - DOESN'T WORK WELL!
        # cdef double diameter, g, C_D, t_v
        # diameter = 2*self.rcubed**(1./3.)
        # g   = 9.8 # gravitational acceleration, m/s
        # C_D = 0.6 # drag coefficient, unitless
        # tv  = (4./3.)*(g*diameter/C_D)*(self.density/RHO_AIR)
        # tv  = sqrt(tv)
        # return tv

        ## BEARD, 1976
        cdef double alpha, r, d, x, x_to_beta

        r = self.rcubed**(1./3.)
        d = 2.*r*1e6 # diameter, m -> micron
        x = self.calc_mass() * 1e3 # convert kg -> g 

        if d <= 134.43:
            alpha = 4.5795e5
            x_to_beta = x**(2./3.)
        elif 134.43 < d <= 1511.64:
            alpha = 4962.0
            x_to_beta = x**(1./3.)
        elif 1511.64 < d <= 3477.84:
            alpha = 1732.0
            x_to_beta = x**(1./6.)
        else:
            alpha = 917.0
            x_to_beta = 1.0

        return 1e-2 * alpha * x_to_beta # from cm/s -> m/s

    property volume:
        def __get__(self):
            return self.calc_volume()
    cdef double calc_volume(self) nogil:
        return self.rcubed*4.*PI/3.

    property mass:
        def __get__(self):
            return self.calc_mass()
    cdef double calc_mass(self) nogil:
        return self.density*self.calc_volume()

    def attributes(self):
        return { 
            'multi': self.multi, 
            'rcubed': self.rcubed,
            'solute': self.solute,
        }

    cpdef copy(self):
        return Superdroplet(self.multi, self.rcubed, self.solute)

    def __repr__(self):
        return "%d" % self.id

ctypedef Superdroplet Superdroplet_t

cdef double tv(double r, double mass):
    ## BEARD, 1976
    cdef double alpha, d, x, x_to_beta

    d = 2.*r*1e6 # diameter, m -> micron
    x = mass * 1e3 # convert kg -> g 

    if d <= 134.43:
        alpha = 4.5795e5
        x_to_beta = x**(2./3.)
    elif 134.43 < d <= 1511.64:
        alpha = 4962.0
        x_to_beta = x**(1./3.)
    elif 1511.64 < d <= 3477.84:
        alpha = 1732.0
        x_to_beta = x**(1./6.)
    else:
        alpha = 917.0
        x_to_beta = 1.0

    return 1e-2 * alpha * x_to_beta # from cm/s -> m/s

cdef int sd_compare(Superdroplet_t a, Superdroplet_t b):
    return a.multi - b.multi

cdef double kernel(Superdroplet_t sd_j, Superdroplet_t sd_k,
                    kernel_id_t kern):
    cdef double b = 1.5e3 # constant for Golovin
    cdef double E_coll, E_coal
    cdef double p, r_j, r_k, x_j, x_k
    cdef double tv_j, tv_k
    cdef double r_small, r_large
    cdef double tv_diff, r_sum

    if kern == GOLOVIN:
        return b*(sd_j.rcubed + sd_k.rcubed)*4.*PI/3.

    p = 1./3.
    r_j = sd_j.rcubed**p
    r_k = sd_k.rcubed**p
    x_j = sd_j.calc_mass()
    x_k = sd_k.calc_mass()

    tv_j = sd_j.calc_terminal_v() 
    tv_k = sd_k.calc_terminal_v()
    
    tv_diff = tv_j - tv_k
    r_sum = r_j + r_k

    # Basic setting for coalescence efficiency
    E_coal = 1.0

    if kern == HYDRO:
        E_coll = 1.0

    elif kern == LONG:
        ## Long (1974) collision kernel
        r_small = dmin(r_j, r_k)*1e6 # convert to micron
        r_large = dmax(r_j, r_k)*1e6

        if r_large >= 50.0: # microns
            E_coll = 1.0
        else:
            ## Simmel et al, 2002 # REMOVE MICRON CONVERSION! 
            # E_coll = dmax(4.5e4 * (r_large*r_large*1e4) * \
            #               (1. - 3e-4/(r_small*1e2)),
            #               1e-3 )
            ## Bott code
            E_coll = 4.5e-4*(r_large*r_large)* \
                     (1.0 - 3.0/(dmax(3., r_small) + 1e-2))

    elif kern == HALL:
        # Hall (1980) collection kernel as collatd by Bott (1998)
        r_small = dmin(r_j, r_k)*1e6 # convert to micron
        r_large = dmax(r_j, r_k)*1e6

        E_coll = hall_effic(r_large, r_small)
        # print r_large, r_small / r_large, E_coll

    else: 
        E_coll = 1.0

    # Limit collection efficiency to 0 <= E_coll <= 1.0
    E_coll = dmin(E_coll, 1.0)
    E_coll = dmax(0.0, E_coll)

    return (E_coll*E_coal)*PI*(r_sum*r_sum)*fabs(tv_diff)

cdef void multi_coalesce(Superdroplet_t sd_j, 
                         Superdroplet_t sd_k, 
                         double gamma):
    """
    Coalesce two superdroplets with one another. Assume
    sd_j.multi > sd_k.multi

    """

    cdef Superdroplet_t sd_temp, sd_recycle
    cdef double gamma_tilde
    cdef long multi_j_p, multi_k_p, excess
    cdef double solute_j_p, solute_k_p, rcubed_j_p, rcubed_k_p

    gamma_tilde = dmin(gamma, <long> floor(sd_j.multi/sd_k.multi))
    excess = sd_j.multi - <long> floor(gamma_tilde*sd_k.multi)

    if excess > 0:

        multi_j_p = excess 
        multi_k_p = sd_k.multi

        rcubed_j_p = sd_j.rcubed
        rcubed_k_p = gamma_tilde*sd_j.rcubed + sd_k.rcubed

        solute_j_p = sd_j.solute
        solute_k_p = gamma_tilde*sd_j.solute + sd_k.solute

        sd_j.set_properties(multi_j_p, rcubed_j_p, solute_j_p)
        sd_k.set_properties(multi_k_p, rcubed_k_p, solute_k_p)

    else: # implies excess == 0

        multi_j_p = <long> floor(sd_k.multi / 2)
        multi_k_p = sd_k.multi - multi_j_p

        rcubed_j_p = gamma_tilde*sd_j.rcubed + sd_k.rcubed
        rcubed_k_p = gamma_tilde*sd_j.rcubed + sd_k.rcubed

        solute_j_p = gamma_tilde*sd_j.solute + sd_k.solute
        solute_k_p = gamma_tilde*sd_j.solute + sd_k.solute

        sd_j.set_properties(multi_k_p, rcubed_j_p, solute_j_p)
        sd_k.set_properties(multi_j_p, rcubed_k_p, solute_k_p)

def recycle(list sds):
    """ For a list of superdroplets, identify which ones have 0 
    multiplicities; for each *i* of these, pick the superdroplet 
    with the *i*th-most multiplicity, split it in half, and copy 
    the properties to the remaining superdroplet. """

    print "RECYCLE"

    cdef:
        int i, n_parts
        Superdroplet sd, sd_donor

    sds.sort(cmp=sd_compare)
    n_parts = len(sds)

    for i in range(n_parts / 2):
        sd = sds[i]

        # Short circuits:
        # 1) Did we already encounter non-zero superdroplets?
        if sd.multi > 0: break

        sd_donor = sds[n_parts-i-1]

        # 2) Does the donor superdroplet have data to spare?
        if sd_donor.multi <= 0: break

        sd.multi = <long> floor(sd_donor.multi/2)
        sd_donor.multi -= sd.multi

        sd.rcubed = sd_donor.rcubed
        sd.solute = sd_donor.solute

    print "   %d superdroplets recycled" % i

    return sds

def step(list sd_list, 
         double t_c, double delta_V, kernel_id kern):

    print "PREP STEPS"
    # 1) Make the random permutation of the super-droplet list
    print "   SHUFFLE LIST"
    np.random.shuffle(sd_list)

    # 2) Make the candidate pairs
    print "   GEN PAIRS"
    cdef long n_part = <long> len(sd_list)

    # 3) Generate the uniform random numbers
    print "PROBABILITY LOOP"
    cdef double scaling = (n_part*(n_part - 1)/2.)/(<int> floor(n_part/2))

    print "PROB / COLLISION LOOP"
    cdef bint collisions = False
    cdef int counter = 0
    cdef:
        unsigned int i
        Superdroplet_t sd_j, sd_k
        double gamma, phi, prob
        long xi_j, xi_k, max_xi
        double K_ij

        cdef int big_probs = 0
        cdef double max_prob = 0.0, min_prob = 1.0

    for i in xrange(n_part/2):
        sd_j = (<Superdroplet_t> sd_list[i])
        sd_k = (<Superdroplet_t> sd_list[i + n_part/2])

        phi = rand() / RAND_FACT
        xi_j = sd_j.multi
        xi_k = sd_k.multi

        K_ij = kernel(sd_j, sd_k, kern)
        max_xi = lmax(xi_j, xi_k)
        # print xi_j, xi_k, max_xi
        prob = scaling*max_xi*(t_c/delta_V)*K_ij

        if prob > max_prob: max_prob = prob
        if prob < min_prob: min_prob = prob
        if prob > 1: big_probs += 1

        # Check for collision and coalesce if necessary
        if ( prob - floor(prob) ) >= phi:
            gamma = floor(prob) + 1

            if xi_j > xi_k:
                multi_coalesce(sd_j, sd_k, gamma)
            else:
                multi_coalesce(sd_k, sd_j, gamma)

            if not collisions:
                collisions = True
            counter += 1

    print "%5d collisions simulated" % counter
    print " Max/min probabilities (count): ", min_prob, max_prob, big_probs

    return sd_list
