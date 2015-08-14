import sys

import pandas as pd
import numpy as np
from numpy import random as npr
from scipy.integrate import quad
from scipy.stats import expon as sse

import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style='ticks')

PYXIMPORT = False
if PYXIMPORT:
    import pyximport
    pyximport.install(
        setup_args = { 'include_dirs': np.get_include() },
    )

# from sd_cy import *
from sd_cy import step as cython_step, recycle, Superdroplet

from cases import get_case

from numba import jit
from functools import partial
from operator import attrgetter

import pstats, cProfile

RHO_WATER = 1e3 # kg/m^3

DIAG_PLOTS = True    
DEBUG      = False
PROFILE    = False

# Don't do diag plots and profile simultaneously 
if PROFILE and DIAG_PLOTS: DIAG_PLOTS = False

def ifloor(x):
    return int(np.floor(x))

## Cell/experiment setup
delta_V = 1e6  # Cell volume, m^3
t_c     = 1.0  # timestep, seconds 
n_part  = 2**13 # number of superdroplets to use in simulation
casename = "shima_golo"   

settings = get_case(casename)
t_end, plot_dt, n_0, R_0, X_0, M_0, m_tot_ana = settings
out_dt  = plot_dt/2

## MANUALLY CHANGE CASE
t_end, plot_dt = 60*60+1, 30*60
f = 1.0
n_0 = (f**3)*2**27 + 2**26 + 2**25
R_0 = 10.e-6/f
X_0 = (4.*np.pi/3.)*(R_0**3)
M_0 = X_0*RHO_WATER
m_tot_ana = 1.0
delta_V /= 1.

print """
CASE - {casename:s}
      t_end = {t_end:d} seconds
    plot_dt = {plot_dt:d} seconds

INITIAL DISTRIBUTION
    n_0 = {n_0:3.1e} m^-3
    R_0 = {R_0:3.2e} m
    X_0 = {X_0:2.1e} m^3
    M_0 = {M_0:2.1e} kg
      m = {m_tot_ana:2.2f} g m^-3
""".format(casename=casename, **settings.__dict__)

def mom_dist(x, l=0., x0=X_0):
    return (x**l)*(n_0/x0)*np.exp(-x/x0)

dist = sse(scale=X_0)
size_dist = lambda x: dist.pdf(x) # m^-3
number_dens = lambda x: (n_0/RHO_WATER) * size_dist(x) # m^-6
mass_dens = lambda x: (x*RHO_WATER)*(n_0/RHO_WATER)*size_dist(x) # is volume

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

Rs = np.logspace(0., np.log10(5e3), 250)

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
        if isinstance(xi, (int, long)):
            print "converting xi -> array"
            xi = np.ones_like(r_grid, dtype=float)*xi
        else:
            print "xi_i type issue", type(xi), xi

    sigma_0 = 0.62
    sigma = sigma_0*(len(r_grid)**(-1./5.))

    xx = Rs
    print len(Rs), len(r_grid), len(xi)
    yy = gtilde_jit(Rs*1e-6, r_grid, xi, sigma) * 1e3 # g m^-3

    return xx, yy

def bin_mass_density(sd_list, div_lnr=False,
                     log_r_min=-7, log_r_max=-2, n_bins=20):
    """ Given a list of superdroplets, compute a binned estimate of 
    the mass density function. 
    
    Parameters
    ----------
    sd_list : list of Superdroplet objects
    div_lnr : bool, optionl
        Divide through by lnr to make comparison with Shima et al, 2009
    log_r_min, log_r_max : floats, optional
    n_bins : int, optional
    
    Returns
    -------
    rs, gyt : array of length `n_bins`
        arrays containing the bin centers in radius space and
        the estimate of the mass density function at those points
    
    """
    
    # Read the superdroplet atributes
    rs = np.array([sd.rcubed**(1./3.) for sd in sd_list])
    ns = np.array([sd.multi for sd in sd_list])
    xs = np.array([sd.mass for sd in sd_list])*1e-3

    df = pd.DataFrame({'r': rs, 'n': ns, 'x': xs})
    
    # Create the bin space
    r_bins = np.logspace(log_r_min, log_r_max, n_bins)
    x_bins = (4.*np.pi/3.)*(r_bins**3)*RHO_WATER
    df['x_bins'] = pd.cut(df['x'], x_bins)

    n_binned = df.groupby('x_bins').sum()['n']
    
    # Compute the number density function
    x_widths = x_bins[1:] - x_bins[:-1]
    x_centers = np.sqrt(x_bins[1:] * x_bins[:-1])

    r_widths = r_bins[1:] - r_bins[:-1]
    r_centers = np.sqrt(r_bins[1:] * r_bins[:-1])

    nxt = n_binned / delta_V / x_widths
    
    # Compute the mass density funtion.
    # Note we convert `x` from kg/m3 -> g/m3
    gyt = 3.*((x_centers*1e3)**2) * nxt
    
    return r_centers, gyt

if DIAG_PLOTS:
    print
    print "SMOOTHING PARTICLES"
    xx, yy = kde_plot(r_grid=r_grid, xi=xi_i)

    fig = plt.figure(1)
    plt.clf()
    plt.plot(xx, yy, '-k', lw=2, label='original')
    plt.semilogx()
    plt.xlim(1, 5000)
    # plt.ylim(0, 2.7)

    plt.xlabel("r ($\mu$m)")
    plt.ylabel("g(ln r) (g/m$^3$/unit ln r)")


def to_sd_array(sds):
    return sds
    # return np.asarray(sds, dtype=Superdroplet)

def sort_sds(sds):
    cmpfun = attrgetter('multi')
    sds = sorted(sds, key=cmpfun)
    return to_sd_array(sds)

def main(profile=False):

    if not PROFILE:
        end = raw_input("Begin simulation? ('n' to break)") == 'n'
        if end: sys.exit() 


    print
    print "BEGINNING MAIN ROUTINE"

    # Generate list of Superdroplets
    sds = [Superdroplet(xi_i, r**3., 0.) for r in r_grid]
    sds = to_sd_array(sds)
    wm0 = np.sum([s.multi*s.mass for s in sds])/1e3
    sdss = [sort_sds(sds), ]
    print "Initial water mass = ", wm0

    results = []

    c_step = partial(cython_step, t_c=t_c, delta_V=delta_V)

    t, ti = 0., 0
    n_drops = len(sds)
    n_init = n_drops*1
    wms = [np.sum([s.multi*s.mass for s in sds])/1e3, ]
    xi_s = [np.sum([s.multi for s in sds]), ]

    while t < t_end:
        t += t_c
        ti += 1

        if profile and ti > 10:
            break

        print "STEP %d (%5.1f s)" % (ti, t)

        sds = c_step(sds)
        # sds = recycle(sds)
        # sds = to_sd_array(sds)

        print len(sds)
        if len(sds) < n_drops:
            print "   Lost %d superdrops" % (n_drops - len(sds), )
            n_drops = len(sds)
        else:
            print "   No change in superdrop number"

        print

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

        if t % out_dt == 0:
            wms.append(np.sum([s.multi*s.mass for s in sds])/1e3)
            xi_s.append(np.sum([s.multi for s in sds]))
            sdss.append(sort_sds(sds))
        print "--"*40

    if DIAG_PLOTS:
        xx, yy = kde_plot(sds)
        plt.figure(1)
        plt.plot(xx, yy, '-b', label='final')
        plt.legend()
        plt.draw()

        plt.figure(2)
        plt.plot(wms, color='k')
        plt.twinx()
        plt.plot(xi_s, color='r')
        plt.legend()

        plt.show()

    return wms, xi_s, sdss

if __name__ == "__main__":

    if DEBUG:
        from math import floor
        ifloor = lambda x: int(floor(x))
        n_part = len(r_grid)
        scaling = (n_part*(n_part - 1)/2.)/ifloor(n_part/2)

        # Generate list of Superdroplets
        sds = [Superdroplet(xi_i, r**3., 0.) for r in r_grid]
        wm0 = np.sum([s.mass for s in sds])/1e3
        print "Initial water mass = ", wm0

        results = []
        n_steps = int(t_end / t_c)
        c_step = partial(cython_step, t_c=t_c, delta_V=delta_V)

    elif PROFILE:
        cProfile.runctx("main(profile=True)", globals(), locals(),
                         "sd_profile.prof")
        s = pstats.Stats("sd_profile.prof")
        s.strip_dirs().sort_stats("time").print_stats()

    else:
        out = main(profile=PROFILE)

        if DIAG_PLOTS:
            plt.figure(3)
            plt.clf()

            sdss = out[-1]
            for sds in sdss[:-1:2]:
                r_centers, gyt = \
                    bin_mass_density(sds, True, n_bins=51, 
                                     log_r_min=-7, log_r_max=np.log(5e-2))
                r_centers *= 1e6 * 10
                plt.plot(r_centers, gyt)

            plt.semilogx()
            plt.xlabel('droplet radius (micron)')
            plt.ylabel('mass density function, g(g m$^{-3}$)')