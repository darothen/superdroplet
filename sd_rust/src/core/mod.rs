//! Core data structures and domain concepts for the superdroplet model.
//!
//! This module contains the fundamental types and constants used throughout
//! the simulation, including the `Droplet` struct and physical constants.

pub mod constants;
pub mod droplet;

// Re-export commonly used items
pub use droplet::{total_water, Droplet};

