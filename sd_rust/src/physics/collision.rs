//! Collision-coalescence algorithm for superdroplets.
//!
//! Implements the Shima et al. (2009) super-droplet method for simulating
//! collision-coalescence in clouds. The algorithm can run sequentially or
//! in parallel using Rayon.

use rand::Rng;
use rand::distr::StandardUniform;
use rayon::prelude::*;
use std::fmt::{self, Display};

use crate::core::{Droplet, total_water};
use crate::physics::Kernel;
use crate::utils::math::knuth_shuffle;

/// Parallelization threshold - only use parallel processing when
/// the number of pairs exceeds this value.
///
/// Set to `usize::MAX` to disable parallelization entirely.
/// Set to `0` to always use parallelization.
// const PARALLEL_THRESHOLD: usize = 0;
const PARALLEL_THRESHOLD: usize = usize::MAX;

/// Configuration for the collision-coalescence model.
pub struct ModelConfig {
    /// Timestep duration in seconds
    pub step_seconds: u32,

    /// Total parcel volume (m³)
    pub delta_v: f64,

    /// Number of superdroplets
    pub num_droplets: u32,

    /// Collision kernel
    pub kernel: Kernel,

    /// Enable debug output
    pub debug: bool,
}

impl Display for ModelConfig {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "ModelConfig {{
            step_seconds: {},
            delta_v: {},
            num_droplets: {},
            kernel: {:?},
            debug: {}
        }}",
            self.step_seconds, self.delta_v, self.num_droplets, self.kernel, self.debug
        )
    }
}

pub struct CollisionStepResult {
    pub counter: u32,
    pub big_probs: u32,
    pub max_prob: f64,
    pub min_prob: f64,
    pub total_xi: f64,
    pub total_water: f64,
}

/// Performs one collision-coalescence timestep.
///
/// This function automatically chooses between parallel and sequential
/// execution based on the problem size and `PARALLEL_THRESHOLD`.
///
/// # Arguments
/// * `droplets` - Mutable slice of droplets to process
/// * `model_config` - Model configuration
/// * `rng` - Random number generator
///
/// # Returns
/// `Ok(())` on success, `Err` if something goes wrong
///
/// # Example
/// ```no_run
/// use sd_rust::{Droplet, Kernel, physics::{collision_step, ModelConfig}};
///
/// let mut droplets = vec![Droplet::new(1000, 1e-6); 2048];
/// let config = ModelConfig {
///     step_seconds: 1,
///     delta_v: 1e6,
///     num_droplets: 2048,
///     debug: false,
///     kernel: Kernel::Golovin,
/// };
/// let mut rng = rand::rng();
///
/// collision_step(&mut droplets, &config, &mut rng).unwrap();
/// ```
pub fn collision_step(
    droplets: &mut Vec<Droplet>,
    model_config: &ModelConfig,
    rng: &mut impl Rng,
) -> Result<CollisionStepResult, std::io::Error> {
    let half_n_part = (model_config.num_droplets as f64 / 2.0).floor() as usize;

    let collision_step_result: Result<CollisionStepResult, std::io::Error>;

    // Use adaptive parallelization based on problem size
    if half_n_part > PARALLEL_THRESHOLD {
        collision_step_result = collision_step_parallel(droplets, model_config, rng)
    } else {
        collision_step_result = collision_step_sequential(droplets, model_config, rng)
    }

    collision_step_result
}

/// Parallelized collision step using Rayon.
///
/// Processes droplet pairs in parallel, using a map-reduce pattern
/// to avoid intermediate allocations.
fn collision_step_parallel(
    droplets: &mut Vec<Droplet>,
    model_config: &ModelConfig,
    rng: &mut impl Rng,
) -> Result<CollisionStepResult, std::io::Error> {
    let t_c = model_config.step_seconds;
    let delta_v = model_config.delta_v;
    let kernel_fn = model_config.kernel.to_function();

    // Permute the droplet list
    knuth_shuffle(droplets, rng);

    // Generate candidate pairs
    let n_part = model_config.num_droplets as f64;
    let half_n_part = (n_part / 2.0).floor() as usize;
    let scaling = (n_part * (n_part - 1.0) / 2.0) / (half_n_part as f64);

    // Pre-generate all random numbers for parallel iteration
    let random_numbers: Vec<f64> = (0..half_n_part)
        .map(|_| rng.sample(StandardUniform))
        .collect();

    // Split droplets into two halves
    let (left_half, right_half) = droplets.split_at_mut(half_n_part);

    // Process pairs in parallel using map-reduce pattern
    let (counter, big_probs, max_prob, min_prob) = left_half
        .par_iter_mut()
        .zip(right_half.par_iter_mut())
        .zip(random_numbers.par_iter())
        .map(|((sd_j, sd_k), &phi)| {
            // Compute collision kernel
            let k_ij = kernel_fn(sd_j, sd_k);

            let max_xi = std::cmp::max(sd_j.multi, sd_k.multi);
            let min_xi = std::cmp::min(sd_j.multi, sd_k.multi);

            if min_xi == 0 {
                return (0u32, 0u32, 0.0f64, 1.0f64);
            }

            let prob = scaling * max_xi as f64 * (t_c as f64 / delta_v) * k_ij;
            let big_prob = if prob > 1.0 { 1u32 } else { 0u32 };

            // Check for collision and coalesce if necessary
            let collision_occurred = if (prob - prob.floor()) > phi {
                let gamma = prob.floor() + 1.0;
                if sd_j.multi < sd_k.multi {
                    multi_coalesce_inline(sd_k, sd_j, gamma);
                } else {
                    multi_coalesce_inline(sd_j, sd_k, gamma);
                }
                1u32
            } else {
                0u32
            };

            (collision_occurred, big_prob, prob, prob)
        })
        .reduce(
            || (0u32, 0u32, 0.0f64, 1.0f64),
            |(c1, bp1, max1, min1), (c2, bp2, max2, min2)| {
                (
                    c1 + c2,
                    bp1 + bp2,
                    max1.max(max2),
                    if min2 > 0.0 { min1.min(min2) } else { min1 },
                )
            },
        );

    let total_xi = droplets.iter().map(|d| d.multi as f64).sum::<f64>();
    let total_water = total_water(droplets);
    if model_config.debug {
        println!("{} collisions simulated (parallel)", counter);
        println!(
            "Max/min probabilities (count): {:.2} {:.2} {}",
            min_prob, max_prob, big_probs
        );
        println!("Total ξ: {}", total_xi);
    }

    Ok(CollisionStepResult {
        counter,
        big_probs,
        max_prob,
        min_prob,
        total_xi,
        total_water,
    })
}

/// Sequential collision step for small problems.
///
/// Uses a simple loop to process pairs sequentially. This avoids
/// the overhead of parallelization for small problem sizes.
fn collision_step_sequential(
    droplets: &mut Vec<Droplet>,
    model_config: &ModelConfig,
    rng: &mut impl Rng,
) -> Result<CollisionStepResult, std::io::Error> {
    let t_c = model_config.step_seconds;
    let delta_v = model_config.delta_v;
    let kernel_fn = model_config.kernel.to_function();

    // Permute the droplet list
    knuth_shuffle(droplets, rng);

    // Generate candidate pairs
    let n_part = model_config.num_droplets as f64;
    let half_n_part = (n_part / 2.0).floor() as usize;
    let scaling = (n_part * (n_part - 1.0) / 2.0) / (half_n_part as f64);

    let mut counter: u32 = 0;
    let mut big_probs: u32 = 0;
    let mut max_prob: f64 = 0.0;
    let mut min_prob: f64 = 1.0;

    for i in 0..half_n_part {
        // Split the vector to get mutable references to two different elements
        let (left, right) = droplets.split_at_mut(i + half_n_part);
        let sd_j = &mut left[i];
        let sd_k = &mut right[0];

        let phi: f64 = rng.sample(StandardUniform);
        let k_ij = kernel_fn(sd_j, sd_k);

        let max_xi = std::cmp::max(sd_j.multi, sd_k.multi);
        let min_xi = std::cmp::min(sd_j.multi, sd_k.multi);
        if min_xi == 0 {
            continue;
        }

        let prob = scaling * max_xi as f64 * (t_c as f64 / delta_v) * k_ij;
        max_prob = if prob > max_prob { prob } else { max_prob };
        min_prob = if prob < min_prob { prob } else { min_prob };
        if prob > 1.0 {
            big_probs += 1;
        }

        // Check for collision and coalesce if necessary
        if (prob - prob.floor()) > phi {
            let gamma = prob.floor() + 1.0;
            if sd_j.multi < sd_k.multi {
                multi_coalesce_inline(sd_k, sd_j, gamma);
            } else {
                multi_coalesce_inline(sd_j, sd_k, gamma);
            }
            counter += 1;
        }
    }

    let total_xi = droplets.iter().map(|d| d.multi as f64).sum::<f64>();
    let total_water = total_water(droplets);
    if model_config.debug {
        println!("{} collisions simulated (sequential)", counter);
        println!(
            "Max/min probabilities (count): {:.2} {:.2} {}",
            min_prob, max_prob, big_probs
        );
        println!("Total ξ: {}", total_xi);
    }

    Ok(CollisionStepResult {
        counter,
        big_probs,
        max_prob,
        min_prob,
        total_xi,
        total_water,
    })
}

/// Inlined coalescence function for maximum performance.
///
/// Merges two droplets according to the multi-coalescence algorithm.
/// This function is marked `#[inline]` to allow the compiler to optimize
/// it into the calling loop.
///
/// # Arguments
/// * `sd_j` - First droplet (must have larger or equal multiplicity)
/// * `sd_k` - Second droplet
/// * `gamma` - Number of coalescence events
#[inline(always)]
fn multi_coalesce_inline(sd_j: &mut Droplet, sd_k: &mut Droplet, gamma: f64) {
    let ratio = sd_j.multi as f64 / sd_k.multi as f64;
    let gamma_t = if gamma < ratio { gamma } else { ratio };
    let excess = sd_j.multi - ((gamma_t * sd_k.multi as f64).floor() as usize);

    if excess > 0 {
        // Case 1: Some droplets from sd_j remain unpaired
        sd_j.multi = excess;
        // sd_k.multi stays the same

        // sd_j.rcubed stays the same
        let rcubed_k_p = gamma_t.mul_add(sd_j.rcubed, sd_k.rcubed);
        sd_k.update_rcubed(rcubed_k_p);

        // sd_j.solute stays the same
        sd_k.solute = gamma_t.mul_add(sd_j.solute, sd_k.solute);
    } else {
        // Case 2: All droplets from sd_j are paired, split sd_k
        let multi_j_p = (sd_k.multi as f64 / 2.0).floor() as usize;
        let multi_k_p = sd_k.multi - multi_j_p;

        let rcubed_new = gamma_t.mul_add(sd_j.rcubed, sd_k.rcubed);
        let solute_new = gamma_t.mul_add(sd_j.solute, sd_k.solute);

        sd_j.multi = multi_j_p;
        sd_k.multi = multi_k_p;

        sd_j.update_rcubed(rcubed_new);
        sd_k.update_rcubed(rcubed_new);

        sd_j.solute = solute_new;
        sd_k.solute = solute_new;
    }
}

/// Public API for multi-coalescence (kept for compatibility).
///
/// This wraps the inlined version for use outside the collision step.
pub fn multi_coalesce(
    sd_j: &mut Droplet,
    sd_k: &mut Droplet,
    gamma: f64,
) -> Result<(), std::io::Error> {
    multi_coalesce_inline(sd_j, sd_k, gamma);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_multi_coalesce() {
        let mut d1 = Droplet::new(1000, 10e-6);
        let mut d2 = Droplet::new(500, 20e-6);

        let old_mass_1 = d1.mass();
        let old_mass_2 = d2.mass();

        multi_coalesce(&mut d1, &mut d2, 1.0).unwrap();

        // At least one droplet should have changed
        assert!(d1.mass() != old_mass_1 || d2.mass() != old_mass_2);

        // Droplets should still have positive multiplicity
        assert!(d1.multi > 0 || d2.multi > 0);
    }
}
