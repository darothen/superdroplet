A Cython implementation of the [Shima et al, 2009](http://dx.doi.org/10.1002/qj.441) Superdroplet algorithm for simulating stochastic collision/coalescence by a nascent droplet population, with modifications based on [Arabas et al, 2015](http://10.5194/gmd-8-1677-2015).

## Installation

None required; you can run **sce.py** directly (preferablly in an IPython terminal).

## Notes

- The slowest part of the algorithm is consistently in Python realm, when we iterate through a list of Superdroplets and compute properties. Very clearly, this should be deferred to a Cython method - potentially even a Cython datastructure which contains a list of Superdroplets ("Superpopulation") and provides easy access to ensemble-level properties.

- Is there a faster way to do np.random.shuffle?

- For some reason, collisions aren't simulated when there are 2**17 superdroplets. The most I can have and still get it to work is 2**15 + (2**15 / 4)

- Increasing delta_V to 1e6 m^3 (like in the Shima et al paper) works fine for Golovin's kernel. Works for the hydrodynamic kernel, too, but the distribution doesn't seem to narrow and grow as much as I'd expect it to.
