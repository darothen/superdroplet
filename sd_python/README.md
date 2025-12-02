# Python Superdroplet Model 

A particle-based cloud microphysics simulation library.

## Overview

## Project Structure

```
sd_python/
├── src/sd_python/
│   ├── __init__.py          # Package initialization
│   ├── __main__.py          # Module entry point (python -m sd_python)
│   ├── cli.py               # Command-line interface
│   ├── core/
│   │   ├── constants.py     # Physical constants
│   │   └── droplet.py       # Droplet data structure and physics
│   ├── physics/
│   │   ├── collision.py     # Collision/coalescence logic
│   │   └── kernels.py       # Collision kernels
│   └── utils/
│       ├── math.py          # Mathematical utilities
│       └── time.py          # Time tracking utilities
├── tests/                   # Comprehensive test suite
├── Justfile                 # Common tasks
├── pyproject.toml           # Project configuration
└── README.md                # This file
```

## Development

### Running Tests

```bash
# Run all tests
uv run pytest

# Run with coverage
uv run pytest --cov=sd_python --cov-report=term-missing

# Run specific test file
uv run pytest tests/test_core_droplet.py -v
```

### Performance Profiling

In general, the shuffling is the slowest part of this implementation, although kernel computation isn't far behind it.

```bash
# Record performance cProfile
uv run python profile_cprofile.py

# Visualize with flameprof
uvx flameprof profile_output.prof > profile_output.svg

# Visualize with snakeviz
uvx snakeviz profile_output.prof
```