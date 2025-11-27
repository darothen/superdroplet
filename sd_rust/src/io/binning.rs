//! Droplet binning for analysis and output.
//!
//! Bins droplets by size for visualization and analysis purposes.

use std::sync::OnceLock;

use crate::core::Droplet;
use crate::utils::math::generate_linear_grid;

/// Number of bins for droplet size distribution
pub const NR: usize = 250;

/// Maximum radius for binning (microns)
const R_MAX: f64 = 5e3;

/// Grid structure for binned droplet data.
///
/// Contains bin edges, midpoints, and aggregated values (e.g., water mass).
pub struct BinGrid<const N: usize> {
    /// Left edges of bins (microns)
    pub lefts: [f64; N],

    /// Right edges of bins (microns)
    pub rights: [f64; N],

    /// Midpoints of bins (microns)
    pub mids: [f64; N],

    /// Aggregated values in each bin (e.g., total water mass)
    pub values: [f64; N],
}

static RADIUS_GRID: OnceLock<Vec<f64>> = OnceLock::new();

/// Returns the radius grid for binning output.
#[inline]
pub fn get_radius_grid() -> &'static Vec<f64> {
    RADIUS_GRID.get_or_init(|| {
        generate_linear_grid(0.0, R_MAX.log10(), NR + 1)
            .iter()
            // Opting to not optimie this exponentiation for clarity.
            .map(|x| 10.0_f64.powf(*x))
            .collect()
    })
}

/// Bins droplets by size and computes total water mass in each bin.
///
/// Uses logarithmically-spaced bins and sorts droplets by size for
/// efficient binning. Avoids cloning the droplet vector by sorting
/// indices instead.
///
/// # Arguments
/// * `droplets` - Slice of droplets to bin
///
/// # Returns
/// A `BinGrid` containing bin edges and water mass in each bin
///
/// # Example
/// ```no_run
/// use sd_rust::io::bin_droplets;
/// use sd_rust::core::Droplet;
///
/// let droplets = vec![Droplet::new(1000, 10e-6); 100];
/// let grid = bin_droplets(&droplets);
/// println!("First bin: {} to {} microns", grid.lefts[0], grid.rights[0]);
/// ```
pub fn bin_droplets(droplets: &[Droplet]) -> BinGrid<NR> {
    let mut bin_grid = BinGrid::<NR> {
        lefts: [0.0; NR],
        rights: [0.0; NR],
        mids: [0.0; NR],
        values: [0.0; NR],
    };

    // Sort indices instead of cloning droplets (much faster!)
    let mut sorted_indices: Vec<usize> = (0..droplets.len()).collect();
    sorted_indices.sort_by(|&a, &b| droplets[a].rcubed.partial_cmp(&droplets[b].rcubed).unwrap());

    let mut bin_idx = 0;
    let mut droplet_idx = 0;

    // Generate logarithmically-spaced radius grid
    let r_grid: &Vec<f64> = get_radius_grid();

    // Iterate over bin edges
    for (left, right) in r_grid.iter().zip(r_grid.iter().skip(1)) {
        let mut value = 0.0; // Total water mass in the bin
        let rcubed_max = (right * 1e-6).powi(3); // Convert microns to meters, then cube

        // Add all droplets that fall in this bin
        for &idx in sorted_indices.iter().skip(droplet_idx) {
            let droplet = &droplets[idx];
            if droplet.rcubed < rcubed_max {
                value += droplet.mass() * (droplet.multi as f64);
                droplet_idx += 1;
            } else {
                break;
            }
        }

        bin_grid.lefts[bin_idx] = *left;
        bin_grid.rights[bin_idx] = *right;
        bin_grid.mids[bin_idx] = (left + right) / 2.0;
        bin_grid.values[bin_idx] = value;
        bin_idx += 1;
    }

    bin_grid
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bin_droplets() {
        let droplets = vec![
            Droplet::new(1000, 10e-6),
            Droplet::new(2000, 20e-6),
            Droplet::new(3000, 30e-6),
        ];

        let grid = bin_droplets(&droplets);

        // Check that bins are created
        assert_eq!(grid.lefts.len(), NR);
        assert_eq!(grid.rights.len(), NR);
        assert_eq!(grid.mids.len(), NR);
        assert_eq!(grid.values.len(), NR);

        // Check that some bins have non-zero values
        let total_mass: f64 = grid.values.iter().sum();
        assert!(total_mass > 0.0);
    }

    #[test]
    fn test_bin_edges_increasing() {
        let droplets = vec![Droplet::new(1000, 10e-6)];
        let grid = bin_droplets(&droplets);

        // Check that bin edges are monotonically increasing
        for i in 0..NR - 1 {
            assert!(grid.lefts[i] < grid.lefts[i + 1]);
            assert!(grid.rights[i] < grid.rights[i + 1]);
        }
    }
}
