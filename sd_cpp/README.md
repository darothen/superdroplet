# Superdroplet Collision-Coalescence Model - C++ Implementation

A high-performance, modern C++ implementation of the super-droplet method for simulating collision-coalescence in clouds, following Shima et al. (2009).

## Requirements

- C++20 compatible compiler (GCC 10+, Clang 10+, MSVC 2019+)
- CMake 3.20 or later
- GoogleTest (automatically fetched during build)

## Performance

The C++ implementation includes several performance optimizations:
- ✅ Compile-time evaluation (constexpr pow2)
- ✅ FMA (fused multiply-add) instructions
- ✅ Link-time optimization (LTO)
- ✅ Aggressive inlining of hot functions
- ✅ Pre-computed scaling factors
- ✅ Optimized compiler flags (`-O3 -march=native -ffast-math -flto`)

**Phase 1 optimizations implemented** - Estimated 15-30% speedup over baseline.

See [PERFORMANCE_OPTIMIZATION.md](PERFORMANCE_OPTIMIZATION.md) for detailed analysis and [OPTIMIZATIONS_IMPLEMENTED.md](OPTIMIZATIONS_IMPLEMENTED.md) for what's been done.

## Testing

### Running Tests

Using Justfile (recommended):
```bash
# Run all tests
just test

# Run with verbose output
just test-verbose

Or using CMake/CTest directly:
```bash
cd build
ctest --output-on-failure
```

See [TEST_SUITE_SUMMARY.md](TEST_SUITE_SUMMARY.md) for detailed test documentation.

## Building

### Quick Build

```bash
# Create build directory
mkdir build && cd build

# Configure and build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build .

# Run the simulation
./sd_cpp
```

### Using the provided script

```bash
./build.sh
```

## Project Structure

```
sd_cpp/
├── CMakeLists.txt          # Build configuration
├── README.md               # This file
├── build.sh                # Quick build script
├── include/                # Public headers
│   ├── core/
│   │   ├── constants.hpp   # Physical constants
│   │   └── droplet.hpp     # Droplet data structure
│   ├── physics/
│   │   ├── collision.hpp   # Collision algorithm
│   │   └── kernels.hpp     # Collision kernels
│   ├── utils/
│   │   ├── math.hpp        # Math utilities
│   │   └── time.hpp        # Time tracking
│   └── io/
│       └── binning.hpp     # Binning for output
├── src/                    # Implementation files
│   ├── main.cpp            # Main executable
│   ├── core/
│   │   └── droplet.cpp
│   ├── physics/
│   │   ├── collision.cpp
│   │   └── kernels.cpp
│   ├── utils/
│   │   ├── math.cpp
│   │   └── time.cpp
│   └── io/
│       └── binning.cpp
└── tests/                  # Unit tests (135 tests)
    ├── CMakeLists.txt
    ├── test_constants.cpp
    ├── test_droplet.cpp
    ├── test_math.cpp
    ├── test_time.cpp
    ├── test_kernels.cpp
    ├── test_collision.cpp
    └── test_binning.cpp
```
