"""
Stochastic collision/coalescence program using the superdroplet model.

Main simulation driver for the superdroplet collision-coalescence model.
"""

using CSV
using Printf
using ProgressMeter
using Random
using Statistics

"""
    run_simulation(; debug::Bool=false, plot::Bool=true)

Run the stochastic collision/coalescence simulation.

# Arguments
- `debug`: Enable debug output (default: false)
- `plot`: Enable CSV output for plotting (default: true)

# Returns
- Nothing (writes output to CSV file)
"""
function run_simulation(; debug::Bool=false, plot::Bool=true)
    
    # Main simulation parameters
    n_0 = 27.0 * 2.0^23  # Total number concentration
    r_0 = 30.531e-6 / 3.0  # Total droplet radius
    x_0 = FOUR_THIRD * PI * r_0^3  # Total droplet volume
    m_0 = x_0 * RHO_WATER  # Total droplet water mass
    delta_v = 1.0e6  # Total parcel volume
    t_c = 1  # Model timestep (seconds)
    kernel = Long  # Collision kernel
    
    n_part = 2^17  # Total number of superdroplets
    t_end = 3600  # Total simulation time (seconds)
    plot_dt = 600  # Output interval time
    smooth_window = 9  # Smoothing window for the median
    
    # CSV output configuration
    output_fn = "collision_output.csv"
    nr = 250  # Number of bins for the output binning grid
    r_max_bin = 5.0e3  # Maximum radius for the output binning grid (microns)
    
    # Simulation setup
    model_config = ModelConfig(
        t_c,
        delta_v,
        n_part,
        kernel,
        debug
    )
    println(model_config)
    
    # Initialize RNG with seed for reproducibility
    Random.seed!(0)
    
    # Initialize the droplet array
    total_droplets = delta_v * n_0
    xi_i = floor(Int, total_droplets / n_part)
    
    # Generate exponentially distributed masses
    # Sample from Exp(λ) where λ = 1/x_0 using inverse transform: X = -ln(U)/λ = -x_0*ln(U)
    x_grid = [-x_0 * log(rand()) for _ in 1:n_part]
    sort!(x_grid)
    
    # Compute the radius grid
    r_grid = [cbrt(x * THREE_FOURTH / PI) for x in x_grid]
    
    # Create the initial droplet array
    droplets = [Droplet(xi_i, radius) for radius in r_grid]
    
    println("\nGRID SETUP")
    @printf("   radii: %.3e - %.3e m\n", droplets[1].radius, droplets[end].radius)
    @printf("  volume: %.3e - %.3e m^3\n", droplets[1].volume, droplets[end].volume)
    @printf("    mass: %.3e - %.3e kg\n", droplets[1].mass, droplets[end].mass)
    
    println("\nSD SETUP")
    @printf("        N_s: %d\n", n_part)
    @printf("       xi_i: %d\n", xi_i)
    @printf("N per SD_xi: %.1f\n", total_droplets / xi_i / n_part)
    
    wm_0 = total_water(droplets)
    @printf("\nInitial total water mass = %.3e kg\n\n", wm_0)
    
    # Set up the output binning grid
    bin_edge_exponents = range(0.0, log10(r_max_bin), length=nr+1)
    bin_edges = [10.0^exp for exp in bin_edge_exponents]
    bin_mids = [0.5 * (left + right) for (left, right) in zip(bin_edges[1:end-1], bin_edges[2:end])]
    
    # Initialize CSV output
    csv_io = nothing
    if plot
        csv_io = open(output_fn, "w")
        # Write header row with bin midpoints
        write(csv_io, "-9999")
        for mid in bin_mids
            @printf(csv_io, ",%.5e", mid)
        end
        write(csv_io, "\n")
    end
    
    # Main loop
    stopwatch = Stopwatch(0)
    step = 0
    
    println("BEGINNING MAIN SIMULATION LOOP\n")
    
    # Create progress bar
    prog = Progress(t_end, dt=1.0, desc="Running simulation: ", 
                    barglyphs=BarGlyphs("[=> ]"), barlen=50)
    
    # Main simulation loop
    while total_seconds(stopwatch) <= t_end
        # Output binned data at specified intervals
        if plot && (step % plot_dt == 0)
            bin_values = bin_droplets(droplets, bin_edges)
            smoothed_bin_values = rolling_median(bin_values, smooth_window)
            
            # Write to CSV
            @printf(csv_io, "%d", total_seconds(stopwatch))
            for value in smoothed_bin_values
                @printf(csv_io, ",%.5e", value)
            end
            write(csv_io, "\n")
            flush(csv_io)
        end
        
        # Perform collision step
        collision_result = collision_step!(droplets, model_config)
        
        # Update progress bar with collision statistics
        prog_msg = @sprintf(
            "Collisions: %d | Probs: %.2f - %.2f [%d] | Water: %.2e kg",
            collision_result.counter,
            collision_result.min_prob,
            collision_result.max_prob,
            collision_result.big_probs,
            collision_result.total_water
        )
        next!(prog, showvalues=[(:Status, prog_msg)])
        
        # Increment time
        step += 1
        increment!(stopwatch, t_c)
    end
    
    # Close CSV file
    if plot
        close(csv_io)
    end
    
    # Final statistics
    println("\n\nSimulation completed successfully.")
    wm_final = total_water(droplets)
    @printf("Remaining water mass: %.3e kg (%.1f%% of initial)\n",
            wm_final, wm_final / wm_0 * 100.0)
    
    return nothing
end
