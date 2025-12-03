#include "io/binning.hpp"
#include "utils/math.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>

namespace sd_cpp {
namespace io {

/// Maximum radius for binning (microns)
constexpr double R_MAX = 5e3;

const std::vector<double>& get_radius_grid() {
    static std::vector<double> radius_grid = []() {
        auto linear_grid = utils::generate_linear_grid(0.0, std::log10(R_MAX), NR + 1);
        std::vector<double> result;
        result.reserve(linear_grid.size());
        
        for (double x : linear_grid) {
            result.push_back(std::pow(10.0, x));
        }
        
        return result;
    }();
    
    return radius_grid;
}

BinGrid bin_droplets(const std::vector<Droplet>& droplets) {
    BinGrid bin_grid;
    
    // Initialize all values to zero
    bin_grid.lefts.fill(0.0);
    bin_grid.rights.fill(0.0);
    bin_grid.mids.fill(0.0);
    bin_grid.values.fill(0.0);
    
    // Sort indices instead of cloning droplets (much faster!)
    std::vector<size_t> sorted_indices(droplets.size());
    std::iota(sorted_indices.begin(), sorted_indices.end(), 0);
    
    std::sort(sorted_indices.begin(), sorted_indices.end(),
        [&droplets](size_t a, size_t b) {
            return droplets[a].rcubed < droplets[b].rcubed;
        });
    
    size_t bin_idx = 0;
    size_t droplet_idx = 0;
    
    // Get the logarithmically-spaced radius grid
    const auto& r_grid = get_radius_grid();
    
    // Iterate over bin edges
    for (size_t i = 0; i < r_grid.size() - 1; ++i) {
        const double left = r_grid[i];
        const double right = r_grid[i + 1];
        
        double value = 0.0; // Total water mass in the bin
        const double rcubed_max = std::pow(right * 1e-6, 3.0); // Convert microns to meters, then cube
        
        // Add all droplets that fall in this bin
        while (droplet_idx < sorted_indices.size()) {
            const size_t idx = sorted_indices[droplet_idx];
            const Droplet& droplet = droplets[idx];
            
            if (droplet.rcubed < rcubed_max) {
                value += droplet.mass() * static_cast<double>(droplet.multi);
                droplet_idx++;
            } else {
                break;
            }
        }
        
        bin_grid.lefts[bin_idx] = left;
        bin_grid.rights[bin_idx] = right;
        bin_grid.mids[bin_idx] = (left + right) / 2.0;
        bin_grid.values[bin_idx] = value;
        bin_idx++;
    }
    
    return bin_grid;
}

} // namespace io
} // namespace sd_cpp
