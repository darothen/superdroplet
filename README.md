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
    + probably not necessary, since that's not the leading-order bottleneck each step.
     
- For some reason, collisions aren't simulated when there are 2**17 superdroplets. The most I can have and still get it to work is 2**15 + (2**15 / 4)
    + **SOLVED**: issue was an integer buffer overflow with the variable `n_part` in the Cython file. Computing `n_part*(n_part - 1)` overflowed and reset it to 0 when `n_part` was sufficiently large. To fix this, I cast the length-inspection operation to a *long and forced `scaling` to be a *double*.

- **hydro** and **golovin** are all elemental operations, but are the most expensive part of **step** because they get called so many times - I think the overhead for function calling is just totalling eating them up. The shuffle operation is a distant second - nearly ten times faster, even when `n_part` is very large. However, inlining the **golovin** operations didn't reduce total run time.
    + These have all been encapsulated in one **kernel** function

- Implemented additional cases in `cases.py`. There are two major issues:
    + Using the cases where `n_0 = 3e8` yield integer-overflows when passing `xi` to the `Superdroplet` constructor. It should be made a long.
    + ~~Reducing `delta_V` tends to result in floating point exceptions within ~600 timesteps.~~ **FIXED** following move to `long`-type multiplicities

- Implemented the Hall collection efficiencies as done by Bott (1998). It's not *really* slow, but there is a good amount of overhead associated with all of its calls. It might be a good idea to obfuscate the computation and compute them for the entire list of Superdroplet pairs in one fell swoop, kind of as implemented by Bott.

- re-added some manual case setting configuration. With a Bott-like initial population and a re-scaling of number and box volume, can reproduce bi-modal droplet growth!!!!

- Biggest sticking point seems to be that the larger droplets grow too fast and congeal into a larger mode too quickly. This seems related to the terminal velocity; using the parameterization from Simmel et al (2002), decrease the Tv by about a factor of 10, which slows down droplet growth, and fixes the problem. It also matches the figure from Lamb and Verlinde (9.5; p 392) really well.

- I disagree with the `simmel_long` cases, because I think they use an incorrect version of the Long kernel (I use Bott's, theirs is very weird)
