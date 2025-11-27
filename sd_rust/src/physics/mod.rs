//! Physics algorithms for the superdroplet model.
//!
//! This module contains the scientific algorithms that drive the simulation,
//! including collision-coalescence and collision kernels.

pub mod collision;
pub mod kernels;

// Re-export commonly used items
pub use collision::{ModelConfig, collision_step, multi_coalesce};
pub use kernels::{Kernel, golovin_kernel, hydro_kernel, long_kernel};
