#!/usr/bin/env python3
"""
Profile the superdroplet simulation using cProfile.

This is a fallback profiler that doesn't require special permissions
like py-spy does. It's less detailed but always works.

Usage:
    uv run python profile_cprofile.py
"""

import cProfile
import pstats
import sys
from pathlib import Path

# Add src to path so we can import sd_python
sys.path.insert(0, str(Path(__file__).parent / "src"))

from sd_python.sce import sce


def main():
    profiler = cProfile.Profile()

    print("Starting profiling with cProfile...")
    print("This will be slower than py-spy but doesn't require special permissions.")
    print()

    profiler.enable()

    # Run the simulation
    sce()

    profiler.disable()

    print("\n" + "=" * 80)
    print("PROFILING RESULTS")
    print("=" * 80 + "\n")

    # Save detailed stats to file
    stats = pstats.Stats(profiler)
    stats.dump_stats("profile_output.prof")
    print("✓ Saved detailed profile to: profile_output.prof")
    print("  View with: python -m pstats profile_output.prof")
    print()

    # Print top functions by cumulative time
    print("Top 30 functions by cumulative time:")
    print("-" * 80)
    stats.sort_stats("cumulative")
    stats.print_stats(30)

    print("\n" + "=" * 80)
    print("Top 20 functions by internal time (time excluding calls):")
    print("-" * 80)
    stats.sort_stats("time")
    stats.print_stats(20)

    print("\n" + "=" * 80)
    print("Callers of collision_step (the hot loop):")
    print("-" * 80)
    stats.print_callers("collision_step")

    print("\n" + "=" * 80)
    print("What collision_step calls:")
    print("-" * 80)
    stats.print_callees("collision_step")

    print("\n✓ Profiling complete!")
    print("\nTo visualize with SnakeViz (if installed):")
    print("  pip install snakeviz")
    print("  snakeviz profile_output.prof")


if __name__ == "__main__":
    main()
