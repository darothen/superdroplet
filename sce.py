import numpy as np
from numpy import random as npr
from scipy.integrate import quad
from scipy.stats import expon as sse

import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style='ticks')

#from superdroplet import Superdroplet
import pyximport; pyximport.install(setup_args={'include_dirs': np.get_include()})
from sd_cy import *
from sd_cy import step as cython_step

from numba import jit
from functools import partial

import pstats, cProfile

RHO_WATER = 1e3 # kg/m^3

DIAG_PLOTS = True

def ifloor(x):
    return int(np.floor(x))

## Cell/experiment setup
delta_V = 1.0  # Cell volume, m^3
t_c     = 1.0  # timestep, seconds    
t_end   = 1801
plot_dt = 60

# Initial size distribution
n_0     = 2.**23. # initial number density of droplets, m^-3
#n_0     = 3e8 #
R_0     = 30.531e-6 # base drop radius, meters
#R_0     = 9.3e-6
X_0     = (4.*np.pi/3.)*(R_0**3.) # base drop volume, m^3
M_0     = X_0*RHO_WATER

m_tot_ana = 1.0 # g m^-3

def mom_dist(x, l=0., x0=X_0):
    return (x**l)*(n_0/x0)*np.exp(-x/x0)

dist = sse(scale=X_0)
size_dist = lambda x: dist.pdf(x) # m^-3
number_dens = lambda x: (n_0/RHO_WATER) * size_dist(x) # m^-6
mass_dens = lambda x: (x*RHO_WATER)*(n_0/RHO_WATER)*size_dist(x) # is volume

n_part = 2**13
x_grid = np.sort( dist.rvs(size=n_part) ) # m^3
r_grid = np.power( x_grid*3./np.pi/4., 1./3. ) # m
m_grid = x_grid*RHO_WATER # kg

print "GRID SETUP"
print "   radii: %1.3e - %1.3e m" % (r_grid[0], r_grid[-1]) 
print "  volume: %1.3e - %1.3e m^3" % (x_grid[0], x_grid[-1])
print "    mass: %1.3e - %1.3e kg" % (m_grid[0], m_grid[-1])
print

m_tot_cal = quad(mom_dist, m_grid[0], m_grid[-1], args=(1., M_0))[0]*1e3

print "TOTAL MASS"
print "   estimate:", m_tot_ana
print "    numeric:", m_tot_cal
print

total_droplets = delta_V*n_0 # unitless 

## MULTIPLICITY BASED ON TOTAL NUMBER OF PARTICLES
xi_i = ifloor(total_droplets / float(n_part))

## MULTIPLICITY BASED ON SIZE DISTRIBUTION

## MULTIPLICITY BASED ON MASS
# xi_i = ifloor(m_tot_cal*delta_V/np.sum((4.*np.pi/3.)*(RHO_WATER*1e3)*x_grid))
# sd_mass = np.sum(xi_i*(4.*np.pi/3.)*(RHO_WATER*1e3)*x_grid)/delta_V
# print "sd_mass:", sd_mass, m_tot_cal

print "SD SETUP"
print "   N_s:", n_part
print "  xi_i:", xi_i
total_xi = xi_i*n_part # unitless
N_per_SD = total_droplets / total_xi
print " N per SD_xi: ", N_per_SD

Rs = np.logspace(0., np.log10(5e3), 50)

@jit("f8[:](f8[:],f8[:],f8[:],f8)")
def gtilde_jit(R, r_grid, xi, sigma):
    nr = len(R)
    nrg = len(r_grid)
    assert len(r_grid) == len(xi)

    ln_R = np.log(R)
    ln_r_grid = np.log(r_grid)

    A = (1./np.sqrt(2.*np.pi)/sigma)

    result = np.zeros_like(Rs)
    for i in xrange(nr):
        s = 0.
        for j in xrange(nrg):
            Y = ln_R[i] - ln_r_grid[j]
            expon = (Y**2.)/2./(sigma**2.)
            W = A*np.exp(-expon)

            s += xi[j]*(r_grid[j]**3.)*W
        s *= RHO_WATER*4.*np.pi/3./delta_V
        result[i] = s

    return result

def kde_plot(sds_list=None, r_grid=None, xi=None):

    if r_grid is None:
        r_grid = np.array([sd.rcubed for sd in sds_list])
        r_grid = np.power( r_grid, 1./3. )

        xi = np.array([sd.multi for sd in sds_list])

    else:
        assert xi is not None
        if type(xi) is int:
            print "converting xi -> array"
            xi = np.ones_like(r_grid, dtype=float)*xi
        else:
            print "xi_i type issue", type(xi), xi

    sigma_0 = 0.62
    sigma = sigma_0*(len(r_grid)**(-1./5.))

    # @np.vectorize
    # def W_sigma(Y):
    #     A = (1./np.sqrt(2.*np.pi)/sigma)
    #     expon = (Y**2.)/2./(sigma**2.)
    #     return A*np.exp(-expon)

    # @np.vectorize
    # def gtilde(R): # kg m^-3 
    #     return (1./delta_V) *                            \
    #             np.sum(                                  \
    #                xi                                    \
    #              * (4.*np.pi/3.)*(R**3.)*RHO_WATER       \
    #              * W_sigma( np.log(R) - np.log(r_grid) ) \
    #            )

    xx = Rs
    print len(Rs), len(r_grid), len(xi)
    yy = gtilde_jit(Rs*1e-6, r_grid, xi, sigma) * 1e3 # g m^-3

    return xx, yy

print
print "SMOOTHING PARTICLES"
xx, yy = kde_plot(r_grid=r_grid, xi=xi_i)

if DIAG_PLOTS:
    fig = plt.figure(1)
    plt.clf()
    plt.plot(xx, yy, '-k', lw=2, label='original')
    plt.semilogx()
    plt.xlim(10, 5000)
    plt.ylim(0, 1.8)

# ## Bin the particles based on Rs
# print
# print "BINNING PARTICLES"
# particles = list(np.copy(r_grid*1e6))
# r_bins = np.zeros_like(Rs)
# for k, particle in enumerate(particles):
#     binned = False
#     for i in xrange(len(r_bins) - 1):
#         if Rs[i] <= particle <= Rs[i+1]:
#             r_bins[i] += 1
#             binned = True
#             break
#     if not binned: print "  could not bin particle %d (r = %1.3e)" % (k, particle)

# bar_left = Rs[:-1]
# bar_right = Rs[1:]
# bar_geomid = np.sqrt(bar_left*bar_right)

# x = (4.*np.pi/3.)*((bar_geomid*1e-6)**3)*RHO_WATER
# Ns = xi_i*r_bins[:-1]*delta_V
# bar_height = 3*(x**2)*Ns

# plt.bar(bar_left, bar_height, bar_right - bar_left, 
#         color='g', alpha=0.5)
# plt.xlabel("Radius ($\mu$m)")
# plt.semilogx()
# plt.plot(r_grid*1e6, 
#          3*((m_grid)**2.)*(number_dens(x_grid)/RHO_WATER), 
#          '-k')
# plt.ylabel("g (g m$^{-3}$)")
# plt.xlim(10, 5000)

# ## Collision Kernel
# @np.vectorize
# def golovin(sd_j, sd_k, b=1.5e3):
#     R3_j = sd_j.rcubed
#     R3_k = sd_k.rcubed
#     return b*(R3_j + R3_k)*4.*np.pi/3.

# @np.vectorize
# def hydro(sd_j, sd_k, E=1.):
#     p = 1./3.
#     r_j = sd_j.rcubed**p
#     r_k = sd_k.rcubed**p
#     tv_j = sd_j.get_terminal_v(RHO_WATER, 1.0)
#     tv_k = sd_k.get_terminal_v(RHO_WATER, 1.0)
#     return E*np.pi*((r_j + r_k)**2)*np.abs(tv_j - tv_k)

# def prob_collision(sd_j, sd_k, t_c=t_c, delta_V=delta_V, 
#                    kernel='golovin'):
#     xi_j = sd_j.multi
#     xi_k = sd_k.multi

#     if kernel == 'golovin':
#         K_ij = golovin(sd_j, sd_k)
#     elif kernel == 'hydro':
#         K_ij = hydro(sd_j, sd_k)
#     else:
#         raise ValueError("Undefined collision kernel (%s)" % kernel)

#     return np.max([xi_j, xi_k])*(t_c/delta_V)*K_ij#*100

#@profile
def step(sd_list):
    diag_msgs = []

    print "PREP STEPS"
    # 1) Make the random permutation of the super-droplet list
    print "   SHUFFLE LIST"
    np.random.shuffle(sd_list)

    # 2) Make the candidate pairs
    print "   GEN PAIRS"
    n_part = len(sd_list)
    #pairs = zip(sd_list[::2], sd_list[1::2])
    #droplets = np.array(pairs).flatten()

    # 3) Generate the uniform random numbers
    print "PROBABILITY LOOP"
    # gammas = detect_collisions(droplets, t_c, delta_V)
    #gammas = []
    scaling = (n_part*(n_part - 1)/2.)/ifloor(n_part/2)    
    # for pair in pairs:
    # for i in xrange(n_part/2):
    #     sd_a, sd_b = sd_list[i], sd_list[i + n_part/2]
    #     phi = npr.uniform()
    #     p_alpha = prob_collision(sd_a, sd_b,
    #                              scaling=scaling,t_c=t_c,delta_V=delta_V)
    #     if phi < p_alpha - np.floor(p_alpha):
    #         gammas.append(int(np.floor(p_alpha) + 1))
    #     else:
    #         gammas.append(int(np.floor(p_alpha)))

    print "PROB / COLLISION LOOP"
    output = []
    collisions = False
    counter = 0
    # for i, (gamma, pair) in enumerate(zip(gammas, pairs)):
    for i in xrange(n_part/2):
        sd_a, sd_b = sd_list[i], sd_list[i + n_part/2]
        gamma = detect_collision(sd_a, sd_b, scaling, t_c, delta_V)

        # phi = npr.uniform()
        # p_alpha = prob_collision(sd_a, sd_b,
        #                          scaling=scaling,t_c=t_c,delta_V=delta_V)
        # gamma = ifloor(p_alpha) + 1 if phi < p_alpha - ifloor(p_alpha) else \
        #         ifloor(p_alpha)

        #print "   ", i, gamma
        if np.abs(gamma) > 0:
            new_pair = multi_coalesce(sd_a, sd_b, gamma)
            output.extend(new_pair)

            # sd_j, sd_k = pair[0].copy(), pair[1].copy()
            # if sd_j.multi < sd_k.multi:
            #     sd_temp = sd_j.copy()
            #     sd_j = sd_k.copy()
            #     sd_k = sd_temp

            # gamma_tilde = np.min([gamma, np.floor(sd_j.multi/sd_k.multi)])
            # #print gamma, gamma_tilde,
            # excess = sd_j.multi - gamma_tilde*sd_k.multi
            # #print excess
            # if excess < 0:
            #     raise ValueError("Excess multiplicity calc between %r and %r yielded \
            #                       excess=%d" % (sd_j, sd_k, excess))

            # if excess > 0:

            #     multi_j_p = sd_j.multi - gamma_tilde*sd_k.multi
            #     multi_k_p = sd_k.multi

            #     rcubed_j_p = sd_j.rcubed
            #     rcubed_k_p = gamma_tilde*sd_j.rcubed + sd_k.rcubed

            #     solute_j_p = sd_j.solute
            #     solute_k_p = gamma_tilde*sd_j.solute + sd_k.solute

            # else:

            #     multi_j_p = ifloor(sd_k.multi / 2.)
            #     multi_k_p = sd_k.multi - multi_j_p

            #     rcubed_j_p = rcubed_k_p = \
            #         gamma_tilde*sd_j.rcubed + sd_k.rcubed

            #     solute_j_p = solute_k_p = \
            #         gamma_tilde*sd_j.solute + sd_k.solute

            # ## Book-keeping

            # if multi_j_p > 0:
            #     output.extend([
            #         Superdroplet(multi_k_p, rcubed_k_p, solute_k_p),
            #         Superdroplet(multi_j_p, rcubed_j_p, solute_j_p)
            #     ])
            # else:
            #     output.append(Superdroplet(multi_k_p, rcubed_k_p, solute_k_p))

            if not collisions:
                collisions = True
            counter += 1

        else:
            output.extend([sd_a, sd_b])

    if collisions:
        diag_msgs.append("%5d collisions simulated" % counter)

    return output, diag_msgs

def main():

    print
    print "BEGINNING MAIN ROUTINE"

    # Generate list of Superdroplets
    sds = [Superdroplet(xi_i, r**3., 0.) for r in r_grid]
    wm0 = np.sum([s.get_mass() for s in sds])/1e3
    print "Initial water mass = ", wm0

    #raw_input("Continue? ")

    # add a few big drops to seed collection
    # big_sd1 = Superdroplet(xi_i/1e8, (2500e-6)**3., 0.)
    # big_sd2 = Superdroplet(xi_i/1e8, (5000e-6)**3., 0.)
    # sds.extend([big_sd1, big_sd2])

    results = []

    # xx, yy = kde_plot(sds)
    # plt.figure(1)
    # plt.plot(xx, yy, '--r', label='timestep 0')
    # plt.draw()
    # results.append(yy)

    # step(sds)
    # cProfile.runctx("step(sds)", globals(), locals(), "step.prof")
    # s = pstats.Stats("step.prof")
    # s.strip_dirs().sort_stats("time").print_stats()

    n_steps = int(t_end / t_c)

    c_step = partial(cython_step, t_c=t_c, delta_V=delta_V)

    t, ti = 0., 0
    n_drops = len(sds)
    n_init = n_drops*1
    wms = [np.sum([s.get_mass() for s in sds])/1e3, ]
    xi_s = [np.sum([s.multi for s in sds]), ]
    while t < t_end:
        t += t_c
        ti += 1

        print "STEP %d (%5.1f s)" % (ti, t)
        sds, msgs = c_step(sds)
        print len(sds)
        if len(sds) < n_drops:
            print "   Lost %d superdrops" % (n_drops - len(sds), )
            n_drops = len(sds)
        else:
            print "   No change in superdrop number"

        print
        for msg in msgs:
            print "   %s" % msg

        if t % plot_dt == 0:
            if DIAG_PLOTS:
                xx, yy = kde_plot(sds)
                plt.figure(1)
                plt.plot(xx, yy, label='t = %4d' % t)
                plt.draw()
                results.append(yy)

        if n_drops < n_init/(2**3):            
            if DIAG_PLOTS:
                xx, yy = kde_plot(sds)
                plt.figure(1)
                plt.plot(xx, yy, label='crash')
                plt.draw()


            print "Superdroplet number dropped too low"
            break

        wms.append(np.sum([s.get_mass() for s in sds])/1e3)
        xi_s.append(np.sum([s.multi for s in sds]))
        print "--"*40

    if DIAG_PLOTS:
        xx, yy = kde_plot(sds)
        plt.figure(1)
        plt.plot(xx, yy, '-b', label='final')
        plt.legend()
        plt.draw()

        plt.show()

    return wms, xi_s

if __name__ == "__main__":

    DEBUG = False

    if DEBUG:
        from math import floor
        ifloor = lambda x: int(floor(x))
        n_part = len(r_grid)
        scaling = (n_part*(n_part - 1)/2.)/ifloor(n_part/2)

        # Generate list of Superdroplets
        sds = [Superdroplet(xi_i, r**3., 0.) for r in r_grid]
        wm0 = np.sum([s.get_mass() for s in sds])/1e3
        print "Initial water mass = ", wm0

        results = []
        n_steps = int(t_end / t_c)
        c_step = partial(cython_step, t_c=t_c, delta_V=delta_V)

    ## DIAGNOSTICS
    out = main()
    # cProfile.run("main()", sort='time')


