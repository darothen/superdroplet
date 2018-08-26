#cython: cdivision=True

include "common.pxi"

cdef int superdroplet_count

cdef class Superdroplet:

    def __init__(self, long multi, double rcubed, double solute):
        global superdroplet_count
        superdroplet_count += 1

        self.multi  = multi
        self.rcubed = rcubed
        self.solute = solute
        self.density = RHO_WATER # kg/m3
        self.id = superdroplet_count

    cdef void set_properties(Superdroplet self,
                             long multi,
                             double rcubed,
                             double solute) nogil:
        self.multi = multi
        self.rcubed = rcubed
        self.solute = solute

    property terminal_v:
        def __get__(self):
            return self.calc_terminal_v()
    cdef double calc_terminal_v(self) nogil:
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

    def __repr__(self):
        return "%d" % self.id
