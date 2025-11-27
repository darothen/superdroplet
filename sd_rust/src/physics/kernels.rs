//! Collision kernel implementations.
//!
//! Collision kernels determine the rate at which droplets collide
//! based on their physical properties.

use crate::core::Droplet;
use crate::core::constants::{FOUR_THIRD, PI};
use crate::utils::math::{max_f64, min_f64};

/// Function pointer type for collision kernels.
pub type KernelFn = fn(&Droplet, &Droplet) -> f64;

#[derive(Debug)]
pub enum Kernel {
    Golovin,
    Hydro,
    Long,
}

impl Kernel {
    /// Computes the collision kernel for two droplets.
    ///
    /// Note: this is just a reference for completion, we use a function
    /// pointer `to_function` instead in the actual implementation.
    ///
    /// # Arguments
    /// * `sd_j` - First droplet
    /// * `sd_k` - Second droplet
    ///
    /// # Returns
    /// Collision kernel value (m³/s)
    #[inline(always)]
    pub fn compute(&self, sd_j: &Droplet, sd_k: &Droplet) -> f64 {
        match self {
            Kernel::Golovin => golovin_kernel(sd_j, sd_k),
            Kernel::Hydro => hydro_kernel(sd_j, sd_k),
            Kernel::Long => long_kernel(sd_j, sd_k),
        }
    }

    /// Converts the kernel enum to a function pointer.
    ///
    /// # Returns
    /// Function pointer to the kernel function.
    #[inline(always)]
    pub fn to_function(&self) -> KernelFn {
        match self {
            Kernel::Golovin => golovin_kernel,
            Kernel::Hydro => hydro_kernel,
            Kernel::Long => long_kernel,
        }
    }
}

/// Golovin collision kernel constant (m³/s)
const GOLOVIN_B: f64 = 1500.0;

/// Pre-computed constant for Golovin kernel
pub const GOLOVIN_CONSTANT: f64 = GOLOVIN_B * FOUR_THIRD * PI;

/// Computes the Golovin collision kernel for two droplets.
///
/// The Golovin kernel is a simple analytical kernel often used for testing:
/// K(i,j) = B * (V_i + V_j)
///
/// # Arguments
/// * `sd_j` - First droplet
/// * `sd_k` - Second droplet
///
/// # Returns
/// Collision kernel value (m³/s)
#[inline(always)]
pub fn golovin_kernel(sd_j: &Droplet, sd_k: &Droplet) -> f64 {
    GOLOVIN_CONSTANT * (sd_j.rcubed + sd_k.rcubed)
}

#[inline(always)]
fn calc_hydro_kernel(e_coal: f64, e_coll: f64, r_sum: f64, tv_diff: f64) -> f64 {
    e_coal * e_coll * PI * r_sum * r_sum * tv_diff.abs()
}

/// Computes the hydrodynamic collision kernel for two droplets.
///
/// The hydrodynamic collision kernel is a simple analytical kernel often used for testing:
/// K(i,j) = E_coal * E_coll * PI * r_sum * r_sum * tv_diff
///
/// # Arguments
/// * `e_coal` - Coalescence efficiency
/// * `e_coll` - Collection efficiency
/// * `r_sum` - Sum of the radii of the two droplets
/// * `tv_diff` - Difference in the terminal velocities of the two droplets
///
/// # Returns
/// Collision kernel value (m³/s)
#[inline(always)]
pub fn hydro_kernel(sd_j: &Droplet, sd_k: &Droplet) -> f64 {
    let r_j = sd_j.radius();
    let r_k = sd_k.radius();
    let tv_j = sd_j.terminal_velocity();
    let tv_k = sd_k.terminal_velocity();

    let tv_diff = tv_j - tv_k;
    let r_sum = r_j + r_k;

    calc_hydro_kernel(1.0, 1.0, r_sum, tv_diff)
}

/// Computes the Long collision kernel for two droplets.
///
/// The Long collision kernel is a simple analytical kernel often used for testing:
/// K(i,j) = E_coal * E_coll * PI * r_sum * r_sum * tv_diff
///
/// # Arguments
/// * `sd_j` - First droplet
/// * `sd_k` - Second droplet
///
/// # Returns
/// Collision kernel value (m³/s)
#[inline(always)]
pub fn long_kernel(sd_j: &Droplet, sd_k: &Droplet) -> f64 {
    let r_j = sd_j.radius();
    let r_k = sd_k.radius();
    let tv_j = sd_j.terminal_velocity();
    let tv_k = sd_k.terminal_velocity();

    let tv_diff = tv_j - tv_k;
    let r_sum = r_j + r_k;

    // Convert radii to microns for collection efficiency calculation
    let (r_small, r_large) = if r_j < r_k {
        (r_j * 1e6, r_k * 1e6)
    } else {
        (r_k * 1e6, r_j * 1e6)
    };

    // Collection efficiency cut-off in limit of very large drops
    let mut e_coll = if r_large >= 50.0 {
        1.0
    } else {
        4.5e-4 * r_large * r_large * (1.0 - 3.0 / (max_f64(3.0, r_small) + 1e-2))
    };

    // Limit collection efficiency to 0 <= E_coll <= 1.0
    e_coll = min_f64(e_coll, 1.0);
    e_coll = max_f64(e_coll, 0.0);

    calc_hydro_kernel(1.0, e_coll, r_sum, tv_diff)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_golovin_kernel() {
        let d1 = Droplet::new(1000, 10e-6);
        let d2 = Droplet::new(2000, 20e-6);
        let k = golovin_kernel(&d1, &d2);
        assert!(k > 0.0);
    }

    #[test]
    fn test_golovin_symmetric() {
        let d1 = Droplet::new(1000, 10e-6);
        let d2 = Droplet::new(2000, 20e-6);
        let k1 = golovin_kernel(&d1, &d2);
        let k2 = golovin_kernel(&d2, &d1);
        assert_eq!(k1, k2);
    }
}
