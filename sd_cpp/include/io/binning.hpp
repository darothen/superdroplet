#pragma once

#include "core/droplet.hpp"
#include <vector>
#include <array>

namespace sd_cpp {
namespace io {

/// Number of bins for droplet size distribution
constexpr size_t NR = 250;

/// Grid structure for binned droplet data.
///
/// Contains bin edges, midpoints, and aggregated values (e.g., water mass).
struct BinGrid {
    /// Left edges of bins (microns)
    std::array<double, NR> lefts;

    /// Right edges of bins (microns)
    std::array<double, NR> rights;

    /// Midpoints of bins (microns)
    std::array<double, NR> mids;

    /// Aggregated values in each bin (e.g., total water mass)
    std::array<double, NR> values;
};

/// Returns the radius grid for binning output.
///
/// Uses a lazy-initialized static variable for efficiency.
const std::vector<double>& get_radius_grid();

/// Bins droplets by size and computes total water mass in each bin.
///
/// Uses logarithmically-spaced bins and sorts droplets by size for
/// efficient binning.
///
/// @param droplets Vector of droplets to bin
/// @return A BinGrid containing bin edges and water mass in each bin
BinGrid bin_droplets(const std::vector<Droplet>& droplets);

} // namespace io
} // namespace sd_cpp
