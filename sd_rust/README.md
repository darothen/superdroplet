# Rust Superdroplet Model

An implementation of the superdroplet collision-coalescence model in Rust.

> [!NOTE]
> I experimented with using `rayon` to implement some parallelization
> but for the majority of use cases, it doesn't help much at all (you
> need very large simulations to clear the overhead). By default, 
> paralleliation is disabled but you can re-enable it by setting the
> value of `PARALLEL_THRESHOLD` in [****src/lib/collisions.rs****]()
> to something lower than the current value.


## Quick Start

```bash
# Build and run with optimizations
cargo run --release

# OPTIONAL - Control thread count
RAYON_NUM_THREADS=8 cargo run --release

# Benchmark
time cargo run --release
```

## Project Structure

```
sd_rust/
├── src/
│   ├── main.rs          # Main simulation driver
│   └── lib.rs           # Core modules:
│       ├── droplets     # Droplet data structure
│       ├── constants    # Physical constants
│       ├── math         # Mathematical utilities
│       ├── model        # Collision-coalescence (parallelized)
│       ├── time         # Time tracking
│       └── io           # Output and binning
├── Cargo.toml           # Dependencies and build config
└── plot.png             # Output visualization
```

## Key Implementation Details

### Parallelization Strategy

The collision step can optionally process droplet pairs in parallel, thanks to the way we split the droplet vector for pairing collision candidates:

1. **Shuffle** droplets randomly
2. **Split** into two non-overlapping halves
3. **Pre-generate** random numbers for all pairs
4. **Process pairs in parallel** using Rayon
5. **Aggregate** statistics after completion

### Adaptive Execution

(see note above)

```rust
const PARALLEL_THRESHOLD: usize = 1000;

if n_pairs > PARALLEL_THRESHOLD {
    // Use parallel processing
} else {
    // Use sequential processing (lower overhead)
}
```

### Cached Properties

```rust
pub struct Droplet {
    pub multi: usize,
    pub rcubed: f64,
    pub solute: f64,
    // Cached values (computed once, reused many times)
    volume: f64,
    mass: f64,
}
```

## Configuration

Edit `src/main.rs` to configure:

```rust
const DEBUG: bool = false;           // Enable/disable debug output
const PLOT: bool = true;             // Enable/disable plotting

let n_part = math::pow2(17);         // Number of superdroplets
let t_end: u32 = 3600;               // Simulation time (seconds)
let plot_dt: u32 = 600;              // Plotting interval
let
```

## Dependencies

- **rand** (0.9.2): Random number generation
- **rand_distr** (0.5.1): Exponential distribution
- **rayon** (1.11.0): Data parallelism
- **plotters** (0.3.7): Visualization

## Building

```bash
# Debug build (slow, for development)
cargo build

# Release build (fast, for production)
cargo build --release

# Run tests
cargo test

# Clean build artifacts
cargo clean
```

## Benchmarking

```bash
# Simple timing
time cargo run --release

# Detailed profiling (requires cargo-flamegraph)
cargo install flamegraph
cargo flamegraph --release
open flamegraph.svg
```
