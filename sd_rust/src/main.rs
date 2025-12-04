use indicatif::{ProgressBar, ProgressStyle};
use rand::SeedableRng;
use rand::rngs::SmallRng;
use rand_distr::{Distribution, Exp};
use sd_rust::{Droplet, Kernel, constants, io, physics, utils};
use std::fs::File;
use std::io::Write;
use std::time::Instant;

const DEBUG: bool = false;
const PLOT: bool = true;

fn main() {
    let n0: f64 = 27.0 * utils::math::pow2(23); // Total number concentration
    let r0: f64 = 30.531e-6 / 3.0; // Total droplet radius
    let x0: f64 = (4.0 / 3.0) * constants::PI * r0.powi(3); // Total droplet volume
    let _m0: f64 = x0 * constants::RHO_WATER; // Total droplet water mass
    let delta_v: f64 = 1e6; // Total parcel volume
    let t_c: u32 = 1; // Model timestep (seconds)
    let kernel: Kernel = Kernel::Golovin; // Collision kernel

    let n_part: usize = utils::math::pow2(19 as usize) as usize; // Total number of superdroplets
    let t_end: u32 = 3600; // Total simulation time (seconds)
    let plot_dt: u32 = 600; // Output interval time
    let smooth_window: usize = 9;

    // CSV output configuration
    let output_fn = String::from("collision_output.csv");
    let mut csv_writer: Option<csv::Writer<std::fs::File>> = None;

    // Use the Xoshiro256++ RNG. The default RNG uses the ChaCha12 algorithm and is
    // cryptographically secure -- way overkill for our Monte Carlo simulations.
    let mut rng = SmallRng::seed_from_u64(0);

    let pb = ProgressBar::new(t_end as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{elapsed} -> {eta} [{bar:60.cyan/blue}] {pos:>7}/{len:7}\n{wide_msg}")
            .unwrap()
            .progress_chars("=> "),
    );

    // SIMULATION SETUP
    let model_config = physics::ModelConfig {
        step_seconds: t_c,
        delta_v: delta_v,
        num_droplets: n_part as u32,
        kernel: kernel,
        debug: DEBUG,
    };
    println!("\nMODEL SETUP\n{}", model_config);

    // Initialize the droplet array;
    // First, pre-compute the superdroplet multiplicity
    let total_droplets = delta_v * n0;
    let xi_i = (total_droplets / (n_part as f64)).floor() as usize;

    // Compute droplet masses according to an exponential distribution
    let dist = Exp::new(1. / x0).unwrap();
    // TODO: re-write this initialization to populate a fixed-size vector
    let mut x_grid: Vec<f64> = dist.sample_iter(&mut rng).take(n_part).collect();
    // NOTE: We have to manually sort with our own comparator since one isn't implemented
    // for float types. This should be safe because we know by construction we will
    // not have any NaNs.
    x_grid.sort_by(|a, b| a.partial_cmp(b).unwrap());
    // let x_min = &x_grid[0];
    // let x_max = &x_grid[n_part - 1];
    // println!("`x_grid`: [{:.2e}, {:.2e}]", x_min, x_max);

    // Compute the radius grid
    let r_grid: Vec<f64> = x_grid
        .into_iter()
        .map(|x| (x * constants::THREE_FOURTH / constants::PI).cbrt())
        .collect();
    // let r_min = &r_grid[0];
    // let r_max = &r_grid[n_part - 1];
    // println!("`rgrid`: [{:.2e}, {:.2e}]", r_min, r_max);

    // Create the initial droplet vector
    let mut droplets: Vec<Droplet> = r_grid
        .into_iter()
        .map(|radius| Droplet::new(xi_i, radius))
        .collect();

    println!("\nGRID SETUP");
    println!(
        "   radii: {:5.3e} - {:5.3e} m",
        &droplets[0].radius(),
        &droplets[n_part - 1].radius()
    );
    println!(
        "  volume: {:5.3e} - {:5.3e} m^3",
        &droplets[0].volume(),
        &droplets[n_part - 1].volume()
    );
    println!(
        "    mass: {:5.3e} - {:5.3e} kg",
        &droplets[0].mass(),
        &droplets[n_part - 1].mass()
    );

    println!("\nSD SETUP");
    println!("        N_s: {:}", n_part);
    println!("       xi_i: {:}", xi_i);
    println!("N per SD_xi: {:}", total_droplets as usize / xi_i / n_part);

    let wm0 = sd_rust::core::total_water(&droplets);
    println!("\nInitial total water mass = {:5.3e} kg\n", wm0);

    // Initialize CSV output
    if PLOT {
        let bin_grid = io::bin_droplets(&droplets);
        let mut writer = csv::Writer::from_path(output_fn).unwrap();
        writer
            .write_record(
                std::iter::once("-9999".to_string())
                    .chain(bin_grid.mids.iter().map(|x| format!("{:e}", x)))
                    .collect::<Vec<String>>(),
            )
            .unwrap();
        csv_writer = Some(writer);
    }

    // Main loop
    let mut masses: Vec<f64> = vec![wm0];
    let mut stopwatch = utils::time::Stopwatch::new(0);
    let mut step: u32 = 0;

    println!("\nBEGINNING MAIN SIMULATION LOOP\n");
    
    // Start timing the main simulation loop
    let start_time = Instant::now();
    
    while stopwatch.total_seconds() <= t_end {
        if PLOT && step % plot_dt == 0 {
            let bin_grid = io::bin_droplets(&droplets);

            // Compute rolling median of the values of the bin grid
            // Pad the edges to avoid truncation
            let half_window = smooth_window / 2;
            let smoothed_values: Vec<f64> = (0..bin_grid.values.len())
                .map(|i| {
                    let start = i.saturating_sub(half_window);
                    let end = (i + half_window + 1).min(bin_grid.values.len());
                    utils::math::median(&bin_grid.values[start..end])
                })
                .collect();
            if let Some(ref mut writer) = csv_writer {
                writer
                    .write_record(
                        std::iter::once(format!("{:}", stopwatch.total_seconds()))
                            .chain(smoothed_values.iter().map(|x| format!("{:e}", x)))
                            .collect::<Vec<String>>(),
                    )
                    .unwrap();
            }
        }

        let collision_step_result = physics::collision_step(&mut droplets, &model_config, &mut rng)
            .expect("Something went wrong in collision step.");
        masses.push(sd_rust::core::total_water(&droplets));

        utils::time::increment(&mut stopwatch, t_c);
        pb.inc(t_c as u64);
        pb.set_message(format!(
            "Collisions: {} | Probabilities: {:.2} - {:.2} [{:}] | Total water: {:.2e} kg",
            collision_step_result.counter,
            collision_step_result.min_prob,
            collision_step_result.max_prob,
            collision_step_result.big_probs,
            collision_step_result.total_water
        ));
        step += 1;
    }

    // End timing
    let elapsed_time = start_time.elapsed();
    let elapsed_seconds = elapsed_time.as_secs_f64();

    if PLOT {
        if let Some(ref mut writer) = csv_writer {
            writer.flush().unwrap();
        }
    }

    println!("\n\nSimulation completed successfully.");
    println!("Runtime: {:.6} seconds", elapsed_seconds);
    let final_mass = masses.last().unwrap();
    println!(
        "Remaining water mass: {:5.3e} kg ({:.1}% of initial)",
        final_mass,
        final_mass / wm0 * 100.0
    );

    // Write timing to file
    let mut time_file = File::create("time.out").expect("Failed to create time.out");
    writeln!(time_file, "{:.6}", elapsed_seconds).expect("Failed to write to time.out");
}
