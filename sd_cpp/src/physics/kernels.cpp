#include "physics/kernels.hpp"
#include "core/constants.hpp"
#include "utils/math.hpp"
#include <algorithm>
#include <cmath>

namespace sd_cpp {
namespace physics {

using namespace constants;

// Golovin collision kernel constant (m³/s)
constexpr double GOLOVIN_B = 1500.0;
constexpr double GOLOVIN_CONSTANT = GOLOVIN_B * FOUR_THIRD * PI;

double golovin_kernel(const Droplet& sd_j, const Droplet& sd_k) {
    return GOLOVIN_CONSTANT * (sd_j.rcubed + sd_k.rcubed);
}

static inline double calc_hydro_kernel(double e_coal, double e_coll, 
                                       double r_sum, double tv_diff) {
    return e_coal * e_coll * PI * r_sum * r_sum * std::abs(tv_diff);
}

double hydro_kernel(const Droplet& sd_j, const Droplet& sd_k) {
    const double r_j = sd_j.radius();
    const double r_k = sd_k.radius();
    const double tv_j = sd_j.terminal_velocity();
    const double tv_k = sd_k.terminal_velocity();

    const double tv_diff = tv_j - tv_k;
    const double r_sum = r_j + r_k;

    return calc_hydro_kernel(1.0, 1.0, r_sum, tv_diff);
}

double long_kernel(const Droplet& sd_j, const Droplet& sd_k) {
    const double r_j = sd_j.radius();
    const double r_k = sd_k.radius();
    const double tv_j = sd_j.terminal_velocity();
    const double tv_k = sd_k.terminal_velocity();

    const double tv_diff = tv_j - tv_k;
    const double r_sum = r_j + r_k;

    // Convert radii to microns for collection efficiency calculation
    const double r_small = std::min(r_j, r_k) * 1e6;
    const double r_large = std::max(r_j, r_k) * 1e6;

    // Collection efficiency cut-off in limit of very large drops
    double e_coll;
    if (r_large >= 50.0) {
        e_coll = 1.0;
    } else {
        e_coll = 4.5e-4 * r_large * r_large * 
                 (1.0 - 3.0 / (utils::max_f64(3.0, r_small) + 1e-2));
    }

    // Limit collection efficiency to 0 <= E_coll <= 1.0
    e_coll = std::clamp(e_coll, 0.0, 1.0);

    return calc_hydro_kernel(1.0, e_coll, r_sum, tv_diff);
}

KernelFn get_kernel_function(Kernel kernel) {
    switch (kernel) {
        case Kernel::Golovin:
            return golovin_kernel;
        case Kernel::Hydro:
            return hydro_kernel;
        case Kernel::Long:
            return long_kernel;
        default:
            return golovin_kernel;
    }
}

const char* kernel_name(Kernel kernel) {
    switch (kernel) {
        case Kernel::Golovin:
            return "Golovin";
        case Kernel::Hydro:
            return "Hydro";
        case Kernel::Long:
            return "Long";
        default:
            return "Unknown";
    }
}

} // namespace physics
} // namespace sd_cpp

