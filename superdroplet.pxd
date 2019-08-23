#cython: cdivision=True
#cython: nonecheck=False

cdef class Superdroplet:

    cdef readonly long multi
    cdef readonly double rcubed, solute, density
    cdef readonly int id

    cdef void set_properties(Superdroplet self, 
                             long multi, 
                             double rcubed, 
                             double solute) nogil
    cdef double calc_terminal_v(Superdroplet self) nogil
    cdef double calc_volume(Superdroplet self) nogil
    cdef double calc_mass(Superdroplet self) nogil
        
ctypedef Superdroplet Superdroplet_t