#pragma once

#include "core/droplet.hpp"

namespace sd_cpp {
namespace physics {

/// Function pointer type for collision kernels.
using KernelFn = double (*)(const Droplet&, const Droplet&);

/// Collision kernel types
enum class Kernel {
    Golovin,
    Hydro,
    Long
};

/// Computes the Golovin collision kernel for two droplets.
///
/// The Golovin kernel is a simple analytical kernel often used for testing:
/// K(i,j) = B * (V_i + V_j)
///
/// @param sd_j First droplet
/// @param sd_k Second droplet
/// @return Collision kernel value (m³/s)
double golovin_kernel(const Droplet& sd_j, const Droplet& sd_k);

/// Computes the hydrodynamic collision kernel for two droplets.
///
/// @param sd_j First droplet
/// @param sd_k Second droplet
/// @return Collision kernel value (m³/s)
double hydro_kernel(const Droplet& sd_j, const Droplet& sd_k);

/// Computes the Long collision kernel for two droplets.
///
/// @param sd_j First droplet
/// @param sd_k Second droplet
/// @return Collision kernel value (m³/s)
double long_kernel(const Droplet& sd_j, const Droplet& sd_k);

/// Converts a kernel enum to a function pointer.
///
/// @param kernel The kernel type
/// @return Function pointer to the kernel function
KernelFn get_kernel_function(Kernel kernel);

/// Returns the name of a kernel as a string.
///
/// @param kernel The kernel type
/// @return Kernel name
const char* kernel_name(Kernel kernel);

} // namespace physics
} // namespace sd_cpp

