
import pandas as pd
from numba import jit
import numpy as np

Rs = np.logspace(0., np.log10(5e3), 250)

RHO_WATER = 1e3

@jit
def gtilde_jit(R, r_grid, xi, sigma, delta_V=1e6):
    nr = len(R)
    nrg = len(r_grid)
    assert len(r_grid) == len(xi)

    ln_R = np.log(R)
    ln_r_grid = np.log(r_grid)

    A = (1./np.sqrt(2.*np.pi)/sigma)

    result = np.zeros_like(Rs)
    for i in range(nr):
        s = 0.
        for j in range(nrg):
            Y = ln_R[i] - ln_r_grid[j]
            expon = (Y**2.)/2./(sigma**2.)
            W = A*np.exp(-expon)

            s += xi[j]*(r_grid[j]**3.)*W
        s *= RHO_WATER*4.*np.pi/3./delta_V
        result[i] = s

    return result


def kde_plot(sds_list=None, r_grid=None, xi=None, delta_V=1e6):

    if r_grid is None:
        r_grid = np.array([sd.rcubed for sd in sds_list])
        r_grid = np.power( r_grid, 1./3. )

        xi = np.array([sd.multi for sd in sds_list])

    else:
        assert xi is not None
        if isinstance(xi, (int, )):
            print("converting xi -> array")
            xi = np.ones_like(r_grid, dtype=float)*xi
        else:
            print("xi_i type issue", type(xi), xi)

    # Use the larger value because it generally smooths things more nicely
    sigma_0 = 1.5 # 0.62
    sigma = sigma_0*(len(r_grid)**(-1./5.))

    xx = Rs
    yy = gtilde_jit(Rs*1e-6, r_grid, xi, sigma, delta_V) * 1e3 # g m^-3

    return xx, yy


def bin_mass_density(sd_list, div_lnr=False,
                     log_r_min=-7, log_r_max=-2, n_bins=20,
                     delta_V=1e6):
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
