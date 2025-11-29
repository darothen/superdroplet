# Modern Fortran Superdroplet Model

An implementation of the superdroplet collision-coalescence model in modern Fortran.

By "modern" Fortran we refer to both the language standard and coding style
(specifically, Fortran 2023 - the latest standard) as well as the build system
(the [Fortran Package Manager](https://fpm.fortran-lang.org/index.html)). We
have only tested this using `gfortran v15.2.0`, which only offers partial
support for the complete Fortran 2023 standard (most notably - lack of the
new enumeration type).

## Quick Start

```bash
# Build with all dependencies and optimizations
fpm build --flag "-std=f2023" --profile release

# Run
fpm run --profile release

# Visualie output from a run
../plot_output.py collision_output.csv

# Build documentation
uvx ford ford.md
```