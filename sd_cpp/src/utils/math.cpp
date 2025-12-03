#include "utils/math.hpp"
#include <algorithm>
#include <cmath>

namespace sd_cpp {
namespace utils {

double median(const std::vector<double>& values) {
    const size_t n = values.size();
    
    // Short circuit some base cases
    if (n == 0) return 0.0;
    if (n == 1) return values[0];
    if (n == 2) return (values[0] + values[1]) / 2.0;
    
    // Create a copy and sort it
    std::vector<double> sorted = values;
    std::sort(sorted.begin(), sorted.end());
    
    // Return median based on even/odd length
    if (n & 1) {
        return sorted[n / 2];
    } else {
        return (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0;
    }
}

std::vector<double> generate_linear_grid(double start, double stop, size_t n) {
    if (n == 0) {
        return std::vector<double>();
    }
    if (n == 1) {
        return std::vector<double>{start};
    }
    
    std::vector<double> result;
    result.reserve(n);
    const double step = (stop - start) / static_cast<double>(n - 1);
    
    for (size_t i = 0; i < n; ++i) {
        result.push_back(start + static_cast<double>(i) * step);
    }
    
    return result;
}

} // namespace utils
} // namespace sd_cpp
