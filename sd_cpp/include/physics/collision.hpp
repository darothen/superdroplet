#pragma once

#include "core/droplet.hpp"
#include "physics/kernels.hpp"
#include <random>
#include <vector>

namespace sd_cpp {
namespace physics {

/// Configuration for the collision-coalescence model.
struct ModelConfig {
  /// Timestep duration in seconds
  unsigned int step_seconds;

  /// Total parcel volume (m³)
  double delta_v;

  /// Number of superdroplets
  unsigned int num_droplets;

  /// Collision kernel
  Kernel kernel;

  /// Enable debug output
  bool debug;
};

/// Result of a collision step.
struct CollisionStepResult {
  /// Number of collisions that occurred
  unsigned int counter;

  /// Number of probabilities > 1
  unsigned int big_probs;

  /// Maximum probability encountered
  double max_prob;

  /// Minimum probability encountered
  double min_prob;

  /// Total multiplicity
  double total_xi;

  /// Total water mass
  double total_water;
};

/// Performs one collision-coalescence timestep.
///
/// @param droplets Mutable vector of droplets to process
/// @param config Model configuration
/// @param rng Random number generator
/// @return Collision step result
CollisionStepResult collision_step(std::vector<Droplet> &droplets,
                                   const ModelConfig &config,
                                   std::mt19937_64 &rng);

/// Inlined coalescence function for maximum performance.
///
/// Merges two droplets according to the multi-coalescence algorithm.
///
/// @param sd_j First droplet (must have larger or equal multiplicity)
/// @param sd_k Second droplet
/// @param gamma Number of coalescence events
__attribute__((always_inline)) inline void
multi_coalesce(Droplet &sd_j, Droplet &sd_k, double gamma) {
  const double ratio =
      static_cast<double>(sd_j.multi) / static_cast<double>(sd_k.multi);
  const double gamma_t = std::min(gamma, ratio);
  const size_t excess =
      sd_j.multi - static_cast<size_t>(std::floor(gamma_t * sd_k.multi));

  if (excess > 0) {
    // Case 1: Some droplets from sd_j remain unpaired
    sd_j.multi = excess;
    // sd_k.multi stays the same

    // sd_j.rcubed stays the same
    const double rcubed_k_p = std::fma(gamma_t, sd_j.rcubed, sd_k.rcubed);
    sd_k.update_rcubed(rcubed_k_p);

    // sd_j.solute stays the same
    sd_k.solute = std::fma(gamma_t, sd_j.solute, sd_k.solute);
  } else {
    // Case 2: All droplets from sd_j are paired, split sd_k
    const size_t multi_j_p = static_cast<size_t>(std::floor(sd_k.multi / 2.0));
    const size_t multi_k_p = sd_k.multi - multi_j_p;

    const double rcubed_new = std::fma(gamma_t, sd_j.rcubed, sd_k.rcubed);
    const double solute_new = std::fma(gamma_t, sd_j.solute, sd_k.solute);

    sd_j.multi = multi_j_p;
    sd_k.multi = multi_k_p;

    sd_j.update_rcubed(rcubed_new);
    sd_k.update_rcubed(rcubed_new);

    sd_j.solute = solute_new;
    sd_k.solute = solute_new;
  }
}

} // namespace physics
} // namespace sd_cpp
