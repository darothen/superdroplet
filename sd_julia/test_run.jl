#!/usr/bin/env julia

"""
Simple test script to run the SuperdropletModel simulation.

Usage:
    julia test_run.jl
"""

# Add the src directory to the load path
push!(LOAD_PATH, joinpath(@__DIR__, "src"))

# Load the module
using SuperdropletModel

# Run the simulation
println("Starting superdroplet simulation...")
println("=" ^ 70)

run_simulation(debug=false, plot=true)

println("\n" * "=" ^ 70)
println("Simulation complete. Output written to collision_output.csv")
