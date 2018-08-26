#!/usr/bin/env python
import sys

import pandas as pd
import numpy as np
from numpy import random as npr
from scipy.integrate import quad
from scipy.stats import expon as sse

import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="ticks", rc={ 'axes.labelsize': 14,
                            'ytick.labelsize': 12,
                            'xtick.labelsize': 12,
                            'legend.fontsize': 13, } )

PYXIMPORT = False
if PYXIMPORT:
    import pyximport
    pyximport.install(
        setup_args = { 'include_dirs': np.get_include() },
    )

# from sd_cy import *
from sd_cy import step as cython_step, recycle
from sd_cy import GOLOVIN, LONG, HYDRO, HALL
from superdroplet import Superdroplet

from cases import get_case
from utils import ( exp_dist_moments, ifloor, estimate_n0,
                    s_to_min_s, fmt_time, colortext_legend )
from plot import *

from numba import jit
from collections import OrderedDict
from functools import partial
from operator import attrgetter

import pstats, cProfile

RHO_WATER = 1e3 # kg/m^3

WRITE_OUT  = True
DIAG_PLOTS = False
DEBUG      = False
PROFILE    = False

# Don't do diag plots and profile simultaneously
if PROFILE and DIAG_PLOTS: DIAG_PLOTS = False

## Cell/experiment setup
delta_V = 1e6  # Cell volume, m^3
t_c     = 1.0  # timestep, seconds
n_part  = 2**17 # number of superdroplets to use in simulation
casename = "shima_golo"
kernel = HALL

settings = get_case(casename)
t_end, plot_dt, n_0, R_0, X_0, M_0, m_tot_ana = settings

## MANUALLY CHANGE CASE
# t_end, plot_dt = 60*60+1, 30*60
# m_tot_ana = 1.0
# f = 1.0
# R_0 = 10.e-6/f
# X_0 = (4.*np.pi/3.)*(R_0**3)
# M_0 = X_0*RHO_WATER
# # n_0 = (f**3)*2**27 + 2**26 + 2**25
# n_0 = estimate_n0(R_0, m_tot_ana)

print("""
CASE - {casename:s}
      t_end = {t_end:d} seconds
    plot_dt = {plot_dt:d} seconds

INITIAL DISTRIBUTION
    n_0 = {n_0:3.1e} m^-3
    R_0 = {R_0:3.2e} m
    X_0 = {X_0:2.1e} m^3
    M_0 = {M_0:2.1e} kg
      m = {m_tot_ana:2.2f} g m^-3
""".format(casename=casename, **settings._asdict()))


# Helper function for checking if we're at an output time
out_dt  = 1200. # plot_dt/1
is_output_step = lambda ti: ti % int(out_dt / t_c) == 0
is_plot_step = lambda ti: ti % int(plot_dt / t_c) == 0

mom_dist = exp_dist_moments
dist = sse(scale=X_0)
size_dist = lambda x: dist.pdf(x) # m^-3
number_dens = lambda x: (n_0/RHO_WATER) * size_dist(x) # m^-6
mass_dens = lambda x: (x*RHO_WATER)*(n_0/RHO_WATER)*size_dist(x) # is volume

x_grid = np.sort( dist.rvs(size=n_part) ) # m^3
r_grid = np.power( x_grid*3./np.pi/4., 1./3. ) # m
m_grid = x_grid*RHO_WATER # kg

print("GRID SETUP")
print("   radii: %1.3e - %1.3e m" % (r_grid[0], r_grid[-1]))
print("  volume: %1.3e - %1.3e m^3" % (x_grid[0], x_grid[-1]))
print("    mass: %1.3e - %1.3e kg" % (m_grid[0], m_grid[-1]))
print()
m_tot_cal = quad(mom_dist, m_grid[0], m_grid[-1],
                 args=(n_0, M_0, 1.))[0]*1e3

print("TOTAL MASS")
print("   estimate:", m_tot_ana)
print("    numeric:", m_tot_cal)
print()
total_droplets = delta_V*n_0 # unitless

## MULTIPLICITY BASED ON TOTAL NUMBER OF PARTICLES
xi_i = ifloor(total_droplets / float(n_part))

## MULTIPLICITY BASED ON SIZE DISTRIBUTION

## MULTIPLICITY BASED ON MASS
# xi_i = ifloor(m_tot_cal*delta_V/np.sum((4.*np.pi/3.)*(RHO_WATER*1e3)*x_grid))
# sd_mass = np.sum(xi_i*(4.*np.pi/3.)*(RHO_WATER*1e3)*x_grid)/delta_V
# print "sd_mass:", sd_mass, m_tot_cal

print("SD SETUP")
print("   N_s:", n_part)
print("  xi_i:", xi_i)
total_xi = xi_i*n_part # unitless
N_per_SD = total_droplets / total_xi
print(" N per SD_xi: ", N_per_SD)

if DIAG_PLOTS:
    print()
    print("SMOOTHING PARTICLES")
    xx, yy = kde_plot(r_grid=r_grid, xi=xi_i, delta_V=delta_V)

    fig = plt.figure(1)
    plt.clf()
    plt.plot(xx, yy, '-k', lw=2, label='initial')
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

    print()
    print("BEGINNING MAIN ROUTINE")

    # Generate list of Superdroplets
    sds = [Superdroplet(xi_i, r**3., 0.) for r in r_grid]
    sds = to_sd_array(sds)
    wm0 = np.sum([s.multi*s.mass for s in sds])/1e3
    sdss = [sort_sds(sds), ]
    print("Initial water mass = ", wm0)

    if not PROFILE:
        end = input("Begin simulation? ('n' to break) ") == 'n'
        if end: sys.exit()

    results = []

    c_step = partial(cython_step, t_c=t_c, delta_V=delta_V,
                     kern=kernel)

    t, ti = 0., 0
    n_drops = len(sds)
    n_init = n_drops*1
    wms = [np.sum([s.multi*s.mass for s in sds])/1e3, ]
    xi_s = [np.sum([s.multi for s in sds]), ]

    color_map = OrderedDict(initial='k')

    while t < t_end:

        minutes, seconds = s_to_min_s(t)

        if profile and ti > 10:
            break

        print("STEP %d (%s)" % (ti, fmt_time(minutes, seconds)))
        print(ti, out_dt, ti % out_dt, is_output_step(ti))
        print(minutes, seconds)

        sds = c_step(sds)
        # sds = recycle(sds)
        sds = to_sd_array(sds)

        print(len(sds))
        if len(sds) < n_drops:
            print("   Lost %d superdrops" % (n_drops - len(sds), ))
            n_drops = len(sds)
        else:
            print("   No change in superdrop number")

        print()

        if (is_output_step(ti) or is_plot_step(ti)) and \
           (DIAG_PLOTS or WRITE_OUT):

            # pre-compute KDE data for plotting/saving
            xx, yy = kde_plot(sds, delta_V=delta_V)

        if is_output_step(ti) and WRITE_OUT and (ti > 0):
            # save KDE-estimate of droplet size distribution
            output_df = pd.DataFrame({'size': xx, 'wm': yy})
            output_df.to_csv("out_{:04d}.csv".format(int(t)),
                             index=False)

            wms.append(np.sum([s.multi*s.mass for s in sds])/1e3)
            xi_s.append(np.sum([s.multi for s in sds]))
            sdss.append(sort_sds(sds))

        if DIAG_PLOTS and is_plot_step(ti) and (ti > 0):

            minutes_str = fmt_time(minutes, seconds, True)

            plt.figure(1)
            lines = plt.plot(xx, yy, lw=2.5, label=minutes_str)
            l = lines[0]
            results.append(yy)

            r_centers, gyt = \
                bin_mass_density(sds, True, n_bins=41,
                                 log_r_min=-7, log_r_max=np.log(5e-2),
                                 delta_V=delta_V)

            r_centers *= 1e6 * 10
            plt.plot(r_centers, gyt, ls='--', color=l.get_color(),
                    # label='%s [binned]' % minutes_str
            )

            color_map[minutes_str] = l.get_color()

            plt.draw()

        if n_drops < n_init/(2**3):
            if DIAG_PLOTS:
                xx, yy = kde_plot(sds, delta_V=delta_V)
                plt.figure(1)
                plt.plot(xx, yy, label='crash')
                plt.draw()

            print("Superdroplet number dropped too low")
            break

        print("--"*40)

        ti += 1
        t = ti * t_c # diagnose time from step to avoid accumulating
                     # floating point errors which screw things up

    if DIAG_PLOTS:
        xx, yy = kde_plot(sds, delta_V=delta_V)
        plt.figure(1)
        # plt.plot(xx, yy, '-b', label='final')
        plt.legend(loc='best')

        ax = plt.gca()
        colortext_legend(color_map, ax, fontsize=14,
                         loc='upper left',
                         bbox_to_anchor=(-0.1, 1.0))

        plt.draw()

        ## Additional diagnostics -
        ## 1) total water mass, black (should be flat)
        ## 2) total SD number, red (should decrease over time)
        # plt.figure(2)
        # plt.clf()
        # plt.plot(wms, color='k')
        # plt.twinx()
        # plt.plot(xi_s, color='r')
        # plt.legend()

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
        print("Initial water mass = ", wm0)

        results = []
        n_steps = int(t_end / t_c)
        c_step = partial(cython_step, t_c=t_c, delta_V=delta_V,
                         kern=kernel)

        print("DEBUG STEP")

        sds = c_step(sds)
        # sds = recycle(sds)
        # sds = to_sd_array(sds)

    elif PROFILE:
        cProfile.runctx("main(profile=True)", globals(), locals(),
                         "sd_profile.prof")
        s = pstats.Stats("sd_profile.prof")
        s.strip_dirs().sort_stats("time").print_stats()

    else:
        out = main(profile=PROFILE)

        # if DIAG_PLOTS:
        #     plt.figure(3)
        #     plt.clf()

        #     sdss = out[-1]
        #     for sds in sdss:
        #         r_centers, gyt = \
        #             bin_mass_density(sds, True, n_bins=61,
        #                              log_r_min=-7, log_r_max=np.log(5e-2)
        #                              delta_V=delta_V)
        #         r_centers *= 1e6 * 10
        #         plt.plot(r_centers, gyt)

        #     plt.semilogx()
        #     plt.xlabel('droplet radius (micron)')
        #     plt.ylabel('mass density function, g(g m$^{-3}$)')
