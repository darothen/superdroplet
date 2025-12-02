#!/usr/bin/env python
"""
Analyze a set of SCE model outputs to make a movie of the autoconversion
process
"""

from glob import glob
import os
import pandas as pd
import re

import matplotlib.pyplot as plt
plt.ioff()
import seaborn as sns
sns.set(style='ticks', context='talk')

PATH_TO_DATA = os.getcwd()
PATH_TO_SAVE = os.path.join(os.getcwd(), 'anim')

def plot_dist(radii, wms, time):

    aspect = 16./9.
    size = 4.
    width = size * aspect
    height = size

    fig = plt.figure(figsize=(width, height))
    ax = fig.add_subplot(111)

    ax.plot(radii, wms, lw=2, color='k')
    plt.semilogx()

    ax.set_xlim(5., 5000.)
    ax.set_ylim(0., 1.8)

    ax.set_xlabel("r, $\mu$m")
    ax.set_ylabel("Mass Density Distribution"
                  "\n"
                  "g(ln r), g/m$^3$/unit ln( r)" )
    ax.set_title("%d seconds" % time, loc='right')
    sns.despine()

    return fig, ax

# Read in the output timesteps
output_files = glob(os.path.join(PATH_TO_DATA, "out_*.csv"))

for i, output_fn in enumerate(output_files):

    fn = os.path.basename(output_fn)
    print fn

    # Extract timestep from filename
    m = re.match("out_(\d{4}).csv", fn)
    timestep = int(m.group(1))

    data = pd.read_csv(output_fn)

    fig, ax = plot_dist(data['size'], data['wm'], timestep)
    plt.savefig(os.path.join(PATH_TO_SAVE, "out_{:04d}.png".format(timestep)),
                dpi=150, transparent=False, bbox_inches='tight')
    plt.close(fig)

    # if i > 10: break
