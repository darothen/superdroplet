//! Superdroplet Collision-Coalescence Model
//!
//! A high-performance implementation of the super-droplet method for simulating
//! collision-coalescence in clouds, following Shima et al. (2009).
//!
//! # Quick Start
//!
//! ```no_run
//! use sd_rust::{Droplet, Kernel, physics, utils};
//!
//! // Create droplets
//! let mut droplets = vec![Droplet::new(1000, 1e-6); 2048];
//!
//! // Configure model
//! let config = physics::ModelConfig {
//!     step_seconds: 1,
//!     delta_v: 1e6,
//!     num_droplets: 2048,
//!     debug: false,
//!     kernel: Kernel::Long,
//! };
//!
//! // Run collision step
//! let mut rng = rand::rng();
//! physics::collision_step(&mut droplets, &config, &mut rng).unwrap();
//! ```
//!
//! # Module Organization
//!
//! - [`core`] - Core data structures (Droplet, constants)
//! - [`physics`] - Physics algorithms (collision, kernels)
//! - [`utils`] - Utility functions (math, time)
//! - [`io`] - Input/output (binning)

pub mod core;
pub mod io;
pub mod physics;
pub mod utils;

// Re-export commonly used items for convenience
pub use core::{Droplet, constants};
pub use physics::{Kernel, ModelConfig, collision_step};
