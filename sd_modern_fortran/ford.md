---
project: Superdroplet Model (Modern Fortran)
summary: Stochastic collision-coalescence simulation using the superdroplet method
author: Daniel Rothenberg
email: daniel@danielrothenberg.com
github: https://github.com/darothen
license: by-nc
src_dir: ./src
         ./app
output_dir: ./doc
project_github: https://github.com/darothen/superdroplet
project_download: https://github.com/darothen/superdroplet/releases
docmark: <
display: public
         protected
         private
source: true
graph: true
search: true
macro: TEST
       LOGIC=.true.
extra_mods: json_module: http://jacobwilliams.github.io/json-fortran/
preprocess: true
predocmark: >
predocmark_alt: #
md_extensions: markdown.extensions.toc
               markdown.extensions.smarty

---

# Superdroplet Model - Modern Fortran Implementation

## Overview

This is a modern Fortran implementation of the stochastic collision-coalescence superdroplet method
originally described by Shima et al. (2009). The model simulates cloud droplet collisions and
coalescence processes using a computationally efficient superdroplet representation.

## Key Features

- **Modern Fortran Design**: Utilizes Fortran 2008/2018 features including:
  - Type-bound procedures for object-oriented programming
  - Abstract interfaces
  - Procedure pointers
  - Modern module structure

- **Performance Optimizations**:
  - Cached terminal velocity calculations
  - Index-based shuffling to avoid copying large structures
  - Precomputed collision kernel constants
  - Optimized for compiler inlining and vectorization

- **Multiple Collision Kernels**:
  - Golovin kernel (analytical test case)
  - Hydrodynamic kernel
  - Long kernel with collection efficiency

## Algorithm

The superdroplet method represents the cloud droplet population using a smaller number of
"superdroplets", each representing many identical physical droplets. The algorithm:

1. Initialize superdroplets with random masses from an exponential distribution
2. For each timestep:
   - Randomly shuffle droplets to form collision pairs
   - Compute collision probabilities using the selected kernel
   - Perform stochastic collision/coalescence based on probabilities
   - Update droplet properties

## References

Shima, S., Kusano, K., Kawano, A., Sugiyama, T., & Kawahara, S. (2009).
The super-droplet method for the numerical simulation of clouds and precipitation:
A particle-based and probabilistic microphysics model coupled with a non-hydrostatic model.
*Quarterly Journal of the Royal Meteorological Society*, 135(642), 1307-1320.

## Building

This project uses the Fortran Package Manager (fpm):

```bash
# Build the project
fpm build

# Run the simulation
fpm run

# Generate documentation
ford ford.md

# Run tests (if available)
fpm test
```

## Project Structure

- `src/` - Core library modules
  - `constants.f90` - Physical and mathematical constants
  - `droplet.f90` - Droplet type and operations
  - `kernels.f90` - Collision kernel implementations
  - `collisions.f90` - Main collision algorithm
  - `time.f90` - Stopwatch utility for timing
  - `util.f90` - Utility functions

- `app/` - Application executables
  - `main.f90` - Main simulation driver

- `test/` - Unit tests

## License

Copyright 2025, Daniel Rothenberg. All rights reserved.