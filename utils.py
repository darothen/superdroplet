import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

RHO_WATER = 1e3 # kg / m3

def ifloor(x):
    return int(np.floor(x))

def exp_dist_moments(x, n0, x0, l=0.):
    """ Exponential distribution PDF multiplied by one
    of its ordinate moments. """
    return (x**l)*(n0/x0)*np.exp(-x/x0)

def estimate_n0(R_0, L):
    """ Given a liquid water content (L) and a mean radius, estimate
    the number concentration of a distribution. """

    X_0 = (4.*np.pi/3.)*(R_0**3)
    M_0 = X_0*RHO_WATER

    r_lo, r_hi = 1e-6, 1e-4
    m_lo, m_hi = (4.*np.pi/3.)*(np.array([r_lo, r_hi])**3)*RHO_WATER

    # Objective function minimize - total liquid water content given n_0
    f = lambda n_0: quad(exp_dist_moments, m_lo, m_hi, 
                         args=(n_0, M_0, 1.))[0]*1e3 # kg/m3 -> g/m3
    g = lambda n_0: f(n_0) - L
    n_opt = brentq(g, 2**20, 5e9, xtol=1e-1, maxiter=1000)

    return n_opt
