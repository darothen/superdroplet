//! Input/Output functionality for the superdroplet model.
//!
//! This module handles data transformation for external consumption,
//! including binning droplets for analysis and visualization.

pub mod binning;

// Re-export commonly used items
pub use binning::{bin_droplets, BinGrid, NR};

