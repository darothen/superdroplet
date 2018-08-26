import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

RHO_WATER = 1e3 # kg / m3

def ifloor(x):
    return int(np.floor(x))

def exp_dist_moments(x, n0, x0, l=0.):
    """ Exponential distribution PDF multiplied by one
    of its ordinate moments. """
    return (x**l)*(n0/x0)*np.exp(-x/x0)

def s_to_min_s(seconds):
    """ Convert a quantity in seconds to minutes and seconds. """
    minutes, seconds = divmod(seconds, 60)
    return minutes, seconds

def fmt_time(minutes, seconds, minutes_only=False):
    """ Yield a print-suitable format for the time in minutes and
    seconds. """

    minutes_str = "{minutes:>4d} min".format(minutes=int(minutes))
    seconds_str = "{seconds:>2.1f} sec".format(seconds=seconds)

    if minutes_only:
        return minutes_str
    else:
        return " ".join([minutes_str, seconds_str])

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

def colortext_legend(text_color_map, ax, text_labels=None, **kwargs):
    """ Add a custom-built legend to a plot, where all the items in
    the legend have colored text corresponding to colored elements
    in the figure.

    Parameters:
    -----------
    text_color_map : dict
        A mapping from labels -> colors which will be used to construct
        the patch elements to include in the legend.
    ax : Axes
        The axes instance on which to add the legend.
    text_labels : dict
        A mapping from labels -> longer descriptions to be used in their
        place in the legend.
    **kwargs : dict-like
        Additional arguments to pass to the legend function.

    Returns:
    --------
    leg : legend object
        The legend, for further customization

    """

    legend_elems = []

    for text, color in text_color_map.items():
        legend_elems.append(
            ( mpatches.Rectangle((0, 0), 1, 1,
                                 facecolor='none',
                                 edgecolor='none'),
              text )
        )

    # Set up a legend with colored-text elements
    elems, texts = zip(*legend_elems)
    ax.legend(elems, texts, **kwargs)
    leg = ax.get_legend()

    # Change the label color
    for label in leg.get_texts():
        old_label = label.get_text()
        if text_labels is None:
            new_label = old_label
        else:
            try:
                new_label = text_labels[old_label]
            except KeyError:
                new_label = old_label

        plt.setp(label, color=text_color_map[old_label])
        label.set_text(new_label)

    return leg