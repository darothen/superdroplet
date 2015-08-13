A Cython implementation of the [Shima et al, 2009](http://dx.doi.org/10.1002/qj.441) Superdroplet algorithm for simulating stochastic collision/coalescence by a nascent droplet population, with modifications based on [Arabas et al, 2015](http://10.5194/gmd-8-1677-2015).

## Installation

None required; you can run **sce.py** directly with the variable `PYXIMPORT` set to `True`. Else, you can build an extension in place:

```bash
$ python setup.py build_ext --inplace
```

## Line profiling

First, build the module using **setup.py** to be sure that the line tracing macros were injected into the Cython module. From an IPython terminal,

```python
>>> run sce.py
>>> %load_ext line_profiler
>>> %lprun -f main -f step main(profile=True)

```


## Notes

- The slowest part of the algorithm is consistently in Python realm, when we iterate through a list of Superdroplets and compute properties. Very clearly, this should be deferred to a Cython method - potentially even a Cython datastructure which contains a list of Superdroplets ("Superpopulation") and provides easy access to ensemble-level properties.

- Is there a faster way to do np.random.shuffle?

- For some reason, collisions aren't simulated when there are 2**17 superdroplets. The most I can have and still get it to work is 2**15 + (2**15 / 4)

- Increasing delta_V to 1e6 m^3 (like in the Shima et al paper) works fine for Golovin's kernel. Works for the hydrodynamic kernel, too, but the distribution doesn't seem to narrow and grow as much as I'd expect it to.

- **hydro** and **golovin** are all elemental operations, but are the most expensive part of **step** because they get called so many times - I think the overhead for function calling is just totalling eating them up. The shuffle operation is a distant second - nearly ten times faster, even when `n_part` is very large. However, inlining the **golovin** operations didn't reduce total run time.