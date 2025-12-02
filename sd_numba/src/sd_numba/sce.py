"""Stochastic collision/coalescence program using the superdroplet model."""

import csv

import numpy as np
from numba.typed import List as TypedList

from sd_numba.core.config import ModelConfig
from sd_numba.core.constants import FOUR_THIRD, PI, RHO_WATER, THREE_FOURTH
from sd_numba.core.droplet import Droplet, compute_total_water
from sd_numba.physics.collision import collision_step
from sd_numba.physics.kernels import Kernel
from sd_numba.utils.binning import bin_droplets
from sd_numba.utils.math import rolling_median
from sd_numba.utils.time import Stopwatch

#: Enable debug mode - print detailed information
DEBUG = False
#: Enable outputting plottable data to a file
PLOT = True


def sce():
    """Run the stochastic collision/coalescence program."""

    # Main simulation parameters
    n_0 = 27.0 * 2.0**23  # Total number concentration
    r_0 = 30.531e-6 / 3.0  # Total droplet radius
    x_0 = FOUR_THIRD * PI * r_0**3  # Total droplet volume
    m_0 = x_0 * RHO_WATER  # Total droplet water mass
    delta_v = 1e6  # Total parcel volume
    t_c = 1.0  # Model timestep (seconds)
    kernel = Kernel.GOLOVIN  # Collision kernel

    n_part = 2**17  # Total number of superdroplets
    # t_end = 3600  # Total simulation time (seconds)
    t_end = 601
    plot_dt = 600  # Output interval time
    smooth_window = 9  # Smoothing window for the median

    # CSV output configuration
    output_fn = "collision_output.csv"
    nr = 250  # Number of bins for the output binning grid
    r_max_bin = 5e3  # Maximum radius for the output binning grid (microns)
    csv_val_fmt = "{val:1.5e}"

    # Simulation setup
    model_config = ModelConfig(
        step_seconds=int(t_c),
        delta_v=delta_v,
        num_droplets=n_part,
        kernel=kernel,
        debug=DEBUG,
    )
    print(model_config)

    # Initialize the droplet array
    total_droplets = delta_v * n_0
    xi_i = int(total_droplets / n_part)
    x_grid = np.random.exponential(scale=x_0, size=n_part)
    x_grid.sort()
    r_grid = np.cbrt(x_grid * THREE_FOURTH / PI)

    # Create the initial droplet array
    droplets = TypedList([Droplet.new(xi_i, radius) for radius in r_grid])
    m_grid = np.asarray([d.mass for d in droplets])

    print("\nGRID SETUP")
    print(f"   radii: {r_grid[0]:.3e} - {r_grid[-1]:.3e} m")
    print(f"  volume: {x_grid[0]:.3e} - {x_grid[-1]:.3e} m^3")
    print(f"    mass: {m_grid[0]:.3e} - {m_grid[-1]:.3e} kg")

    print("\nSD SETUP")
    print(f"        N_s: {n_part}")
    print(f"       xi_i: {xi_i}")
    print(f"N per SD_xi: {total_droplets / xi_i / n_part}")

    wm_0 = compute_total_water(droplets)
    print(f"\nInitial total water mass = {wm_0:.3e} kg\n")

    # Set up the output binning grid
    bin_edge_exponents = np.linspace(0.0, np.log10(r_max_bin), nr + 1)
    bin_edges = np.power(10.0, bin_edge_exponents)
    bin_mids = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    if PLOT:
        bin_values = bin_droplets(droplets, bin_edges)

        with open(output_fn, "w") as csv_file:
            csv_writer = csv.writer(csv_file)
            row_vals = [-9999] + [csv_val_fmt.format(val=mid) for mid in bin_mids]
            csv_writer.writerow(row_vals)

    # Main loop
    masses = [
        wm_0,
    ]
    stopwatch = Stopwatch()
    step = 0
    print("\nBEGINNING MAIN SIMULATION LOOP\n")

    while stopwatch.total_seconds() <= t_end:
        if PLOT and step % plot_dt == 0:
            # print(f"   PLOTTING ({stopwatch}) ... ")
            bin_values = bin_droplets(droplets, bin_edges)
            smoothed_bin_values = rolling_median(bin_values, window_size=smooth_window)
            with open(output_fn, "a") as csv_file:
                csv_writer = csv.writer(csv_file)
                row_vals = [stopwatch.total_seconds()] + [
                    csv_val_fmt.format(val=value) for value in smoothed_bin_values
                ]
                csv_writer.writerow(row_vals)

        collision_step_result = collision_step(droplets, model_config)

        counter = collision_step_result.counter
        min_prob = collision_step_result.min_prob
        max_prob = collision_step_result.max_prob
        big_probs = collision_step_result.big_probs
        total_water = compute_total_water(droplets)
        collision_message = (
            f"Collisions: {collision_step_result.counter} | "
            f"Probabilities: {collision_step_result.min_prob:.2f} - {collision_step_result.max_prob:.2f} [{collision_step_result.big_probs}] | "
            f"Total water: {total_water:.2e} kg"
        )
        print(f"Step {step:5d} | {stopwatch} | {collision_message}")
        # main_pbar.update(t_c)
        # desc_pbar.set_postfix_str(collision_message)
        step += 1
        stopwatch.increment(int(t_c))

    ## Clean up
    # main_pbar.close()


if __name__ == "__main__":
    sce()
