#!/usr/bin/env python
from __future__ import print_function

import glob
import os

import numpy as np
import pandas as pd

from plot import *
from superdroplet import Superdroplet

import matplotlib.pyplot as plt

OUT_DIR = ""

if __name__ == "__main__":

    out_files = glob.glob(os.path.join(OUT_DIR, "*output.txt"))
    print(out_files)

    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    plt.semilogx()
    plt.xlim(1, 5000)

    output = []
    for i, fn in enumerate(out_files):
        print("Reading %s" % fn)
        data = pd.read_csv(fn, names=['rcubed', 'multi'])
        sds = []
        for _, row in data.iterrows():
            droplet = Superdroplet(row.multi, row.rcubed, 0.)
            sds.append(droplet)
        output.append((i, sds))

        xx, yy = kde_plot(sds)
        ax.plot(xx, yy, lw=2.5)
        plt.draw()

    plt.show()
