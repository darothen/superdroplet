# SuperdropletModel.jl

A modern Julia implementation of the Shima et al. (2009) superdroplet method for simulating stochastic collision-coalescence in clouds.

This implementation follows the refined architecture of the reference implementations in Rust, Python, and Fortran, with proper module structure and modern Julia best practices.

## Installation

This package uses standard Julia package management with `Project.toml`. To set up:

```bash
cd sd_julia-v2
```

Then in Julia:

```julia
# Activate the project environment
using Pkg
Pkg.activate(".")

# Install dependencies
Pkg.instantiate()
```

## Testing

Run the test suite to verify the installation:

```julia
using Pkg
Pkg.activate(".")
Pkg.test()
```

Or from the command line:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Usage

### Basic Usage

```julia
using SuperdropletModel

# Run simulation with default parameters
run_simulation()
```

This will:
- Simulate 131,072 superdroplets (2^17)
- Run for 3600 seconds (1 hour) with 1-second timesteps
- Use the Long collision kernel
- Output results to `collision_output.csv` every 600 seconds

### Custom Parameters

To enable debug output:

```julia
run_simulation(debug=true)
```

To disable CSV output:

```julia
run_simulation(plot=false)
```

### Output

The simulation produces a CSV file `collision_output.csv` with the following format:
- First row: `-9999` followed by bin midpoints (radius in microns)
- Subsequent rows: Time (seconds) followed by total water mass in each bin (kg)

The output is compatible with visualization tools used for the other reference implementations.

## Architecture

The implementation follows a modular structure:

```
src/
├── SuperdropletModel.jl     # Main module
├── core/
│   ├── constants.jl         # Physical and mathematical constants
│   ├── config.jl            # Model configuration
│   └── droplet.jl           # Droplet data structure
├── physics/
│   ├── kernels.jl           # Collision kernels
│   └── collision.jl         # Collision-coalescence algorithm
├── utils/
│   ├── math.jl              # Mathematical utilities
│   ├── time.jl              # Time tracking
│   └── binning.jl           # Droplet size binning
└── sce.jl                   # Main simulation driver
```

## Reference

Shima, S., Kusano, K., Kawano, A., Sugiyama, T., & Kawahara, S. (2009). The super-droplet method for the numerical simulation of clouds and precipitation: A particle-based and probabilistic microphysics model coupled with a non-hydrostatic model. *Quarterly Journal of the Royal Meteorological Society*, 135(642), 1307-1320. [doi:10.1002/qj.441](http://dx.doi.org/10.1002/qj.441)
