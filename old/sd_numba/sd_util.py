from enum import Enum

from numba import jit

import numpy as np

from hall import hall_effic
from utils import *

class Kernel(Enum):
    GOLOVIN = 1
    LONG    = 2
    HYDRO   = 3
    HALL    = 4

def sd_compare(a, b):
    return a.multi - b.multi

# @jit(nopython=True)
# @profile
def _kernel_golovin(rcubed_j, rcubed_k):
    b = 1.5e3
    return b*(rcubed_j + rcubed_k)*4.*PI/3.

@jit(nopython=True)
def _kernel(rcubed_j, rcubed_k, tv_diff, kern):

    p = 1./3.
    r_j = rcubed_j**p
    r_k = rcubed_k**p

    tv_diff = np.abs(tv_diff)
    r_sum = r_j + r_k

    # Basic setting for coalescence efficiency
    E_coal = 1.0

    if kern == Kernel.HYDRO:
        E_coll = 1.0

    elif kern == Kernel.LONG:
        ## Long (1974) collision kernel
        r_small = dmin(r_j, r_k)*1e6 # convert to micron
        r_large = dmax(r_j, r_k)*1e6

        if r_large >= 50.0: # microns
            E_coll = 1.0
        else:
            ## Bott code
            E_coll = 4.5e-4*(r_large*r_large)* \
                     (1.0 - 3.0/(dmax(3., r_small) + 1e-2))

    elif kern == Kernel.HALL:
        # Hall (1980) collection kernel as collatd by Bott (1998)
        r_small = dmin(r_j, r_k)*1e6 # convert to micron
        r_large = dmax(r_j, r_k)*1e6

        E_coll = hall_effic(r_large, r_small)

    else:
        E_coll = 1.0

    # Limit collection efficiency to 0 <= E_coll <= 1.0
    E_coll = dmin(E_coll, 1.0)
    E_coll = dmax(0.0, E_coll)

    return (E_coll*E_coal)*PI*(r_sum*r_sum)*tv_diff


# @jit
# @profile
def kernel(sd_j, sd_k, kern):

    if kern == Kernel.GOLOVIN:
        b = 1.5e3
        return b*(sd_j.rcubed + sd_k.rcubed)*4.*PI/3.

    else: # Hydrodynamic base kernel

        tv_diff = sd_j.terminal_v - sd_k.terminal_v

        return _kernel(sd_j.rcubed, sd_k.rcubed, tv_diff, kern)

# @profile
def multi_coalesce(sd_j, sd_k, gamma):
    """
    Coalesce two superdroplets with one another. Assume
    sd_j.multi > sd_k.multi

    """

    gamma_tilde = dmin(gamma, np.floor(sd_j.multi/sd_k.multi))
    excess = sd_j.multi - np.floor(gamma_tilde*sd_k.multi)

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

        multi_j_p = np.floor(sd_k.multi / 2)
        multi_k_p = sd_k.multi - multi_j_p

        rcubed_j_p = gamma_tilde*sd_j.rcubed + sd_k.rcubed
        rcubed_k_p = gamma_tilde*sd_j.rcubed + sd_k.rcubed

        solute_j_p = gamma_tilde*sd_j.solute + sd_k.solute
        solute_k_p = gamma_tilde*sd_j.solute + sd_k.solute

        sd_j.set_properties(multi_k_p, rcubed_j_p, solute_j_p)
        sd_k.set_properties(multi_j_p, rcubed_k_p, solute_k_p)

    return sd_j, sd_k

def recycle(sds):
    """ For a list of superdroplets, identify which ones have 0
    multiplicities; for each *i* of these, pick the superdroplet
    with the *i*th-most multiplicity, split it in half, and copy
    the properties to the remaining superdroplet. """

    print("RECYCLE")

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

        sd.multi = np.floor(sd_donor.multi/2)
        sd_donor.multi -= sd.multi

        sd.rcubed = sd_donor.rcubed
        sd.solute = sd_donor.solute

    print("   %d superdroplets recycled" % i)

    return sds

# @profile
def step(sd_list, t_c, delta_V, kern):

    print("PREP STEPS")
    # 1) Make the random permutation of the super-droplet list
    print("   SHUFFLE LIST")
    np.random.shuffle(sd_list)

    # 2) Make the candidate pairs
    print("   GEN PAIRS")
    n_part = len(sd_list)

    # 3) Generate the uniform random numbers
    print("PROBABILITY LOOP")
    scaling = (n_part*(n_part - 1)/2.)/(np.floor(n_part/2))

    print("PROB / COLLISION LOOP")
    collisions = False
    counter = 0
    big_probs = 0
    max_prob = 0.0
    min_prob = 1.0

    for i in range(n_part//2):
        sd_j = sd_list[i]
        sd_k = sd_list[i + n_part//2]

        phi = np.random.rand()
        xi_j = sd_j.multi
        xi_k = sd_k.multi

        K_ij = kernel(sd_j, sd_k, kern)
        max_xi = lmax(xi_j, xi_k)
        prob = scaling*max_xi*(t_c/delta_V)*K_ij

        if prob > max_prob: max_prob = prob
        if prob < min_prob: min_prob = prob
        if prob > 1: big_probs += 1

        # Check for collision and coalesce if necessary
        if (phi < (prob - np.floor(prob))):
            gamma = np.floor(prob) + 1
        else:
            gamma = np.floor(prob)

        if gamma <= 0:
            continue

        # if ( prob - np.floor(prob) ) >= phi:
            # gamma = np.floor(prob) + 1
        if xi_j > xi_k:
            sd_j, sd_k = multi_coalesce(sd_j, sd_k, gamma)
        else:
            sd_j, sd_k = multi_coalesce(sd_k, sd_j, gamma)
        sd_list[i] = sd_j
        sd_list[i + n_part//2] = sd_k

        if not collisions:
            collisions = True
        counter += 1

    print("%5d collisions simulated" % counter)
    print(" Max/min probabilities (count): ", min_prob, max_prob, big_probs)

    return sd_list
