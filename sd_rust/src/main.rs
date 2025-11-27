use indicatif::{ProgressBar, ProgressStyle};
use plotters::prelude::*;
use rand::SeedableRng;
use rand::rngs::SmallRng;
use rand_distr::{Distribution, Exp};
use sd_rust::{Droplet, Kernel, constants, io, physics, utils};

const DEBUG: bool = false;
const PLOT: bool = true;

fn main() {
    let n0: f64 = 27.0 * utils::math::pow2(23); // Total number concentration
    let r0: f64 = 30.531e-6 / 3.0; // Total droplet radius
    let x0: f64 = (4.0 / 3.0) * constants::PI * r0.powi(3); // Total droplet volume
    let _m0: f64 = x0 * constants::RHO_WATER; // Total droplet water mass
    let delta_v: f64 = 1e6; // Total parcel volume
    let t_c: u32 = 1; // Model timestep (seconds)
    let kernel: Kernel = Kernel::Long; // Collision kernel

    let n_part: usize = utils::math::pow2(17 as usize) as usize; // Total number of superdroplets
    let t_end: u32 = 3600; // Total simulation time (seconds)
    let plot_dt: u32 = 600; // Output interval time
    let smooth_window: usize = 9;

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

    // Main loop
    let mut masses: Vec<f64> = vec![wm0];
    let mut stopwatch = utils::time::Stopwatch::new(0);
    let mut step: u32 = 0;

    println!("\nBEGINNING MAIN SIMULATION LOOP\n");
    let mut chart = None;
    if PLOT {
        print!("   INITIALIZING PLOT ... ");
        let root_drawing_area = BitMapBackend::new("plot.png", (1024, 768)).into_drawing_area();
        root_drawing_area.fill(&WHITE).unwrap();
        let bin_grid = io::bin_droplets(&droplets);
        let smoothed_values: Vec<f64> = bin_grid
            .values
            .windows(smooth_window)
            .map(|window| utils::math::median(&window))
            .collect();
        let (x_min, x_max) = (bin_grid.mids[0], bin_grid.mids[bin_grid.mids.len() - 1]);
        let (y_min, y_max) = (
            *smoothed_values
                .iter()
                .min_by(|a, b| a.total_cmp(b))
                .unwrap(),
            *smoothed_values
                .iter()
                .max_by(|a, b| a.total_cmp(b))
                .unwrap(),
        );

        // println!(
        //     "   PLOT RANGE: [{:.2e}, {:.2e}] x [{:.2e}, {:.2e}]",
        //     x_min, x_max, y_min, y_max
        // );

        let mut chart_obj = ChartBuilder::on(&root_drawing_area)
            // .set_label_area_size(LabelAreaPosition::Left, 40)
            // .set_label_area_size(LabelAreaPosition::Bottom, 40)
            .x_label_area_size(40)
            .y_label_area_size(40)
            .build_cartesian_2d((x_min..x_max).log_scale(), y_min..y_max)
            .unwrap();
        chart_obj.configure_mesh().draw().unwrap();
        chart_obj
            .draw_series(LineSeries::new(
                bin_grid
                    .mids
                    .iter()
                    .zip(smoothed_values.iter())
                    .map(|(x, y)| (*x, *y)),
                &RED,
            ))
            .unwrap()
            .label("t=0 (initial)")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));
        chart = Some(chart_obj);
        println!(" done.\n");
    }

    while stopwatch.total_seconds() <= t_end {
        if PLOT && step % plot_dt == 0 && step > 0 {
            let bin_grid = io::bin_droplets(&droplets);

            // Compute rolling median of the values of the bin grid
            let smoothed_values: Vec<f64> = bin_grid
                .values
                .windows(smooth_window)
                .map(|window| utils::math::median(&window))
                .collect();

            chart
                .as_mut()
                .unwrap()
                .draw_series(LineSeries::new(
                    bin_grid
                        .mids
                        .iter()
                        .zip(smoothed_values.iter())
                        .map(|(x, y)| (*x, *y)),
                    &RED,
                ))
                .unwrap()
                .label(&format!("t={} s", stopwatch.total_seconds()))
                .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));
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

    if PLOT {
        chart
            .as_mut()
            .unwrap()
            .configure_series_labels()
            .background_style(&WHITE.mix(0.8))
            .border_style(&BLACK)
            .draw()
            .unwrap();
    }

    println!("\n\nSimulation completed successfully.");
    let final_mass = masses.last().unwrap();
    println!(
        "Remaining water mass: {:5.3e} kg ({:.1}% of initial)",
        final_mass,
        final_mass / wm0 * 100.0
    );
}
