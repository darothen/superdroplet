#include "physics/collision.hpp"
#include "core/constants.hpp"
#include "utils/math.hpp"
#include <algorithm>
#include <random>

namespace sd_cpp {
namespace physics {

CollisionStepResult collision_step(std::vector<Droplet> &droplets,
                                   const ModelConfig &config,
                                   std::mt19937_64 &rng) {
  const unsigned int t_c = config.step_seconds;
  const double delta_v = config.delta_v;
  const KernelFn kernel_fn = get_kernel_function(config.kernel);

  // Permute the droplet list
  utils::knuth_shuffle(droplets, rng);

  // Generate candidate pairs
  const double n_part = static_cast<double>(config.num_droplets);
  const size_t half_n_part = static_cast<size_t>(std::floor(n_part / 2.0));
  const double scaling =
      (n_part * (n_part - 1.0) / 2.0) / static_cast<double>(half_n_part);

  // Pre-compute scaling factor to reduce multiplications in loop
  const double scaling_factor = scaling * (static_cast<double>(t_c) / delta_v);

  unsigned int counter = 0;
  unsigned int big_probs = 0;
  double max_prob = 0.0;
  double min_prob = 1.0;

  std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);

  for (size_t i = 0; i < half_n_part; ++i) {
    Droplet &sd_j = droplets[i];
    Droplet &sd_k = droplets[i + half_n_part];

    const double phi = uniform_dist(rng);
    const double k_ij = kernel_fn(sd_j, sd_k);

    // Cache droplet multiplicities
    const size_t multi_j = sd_j.multi;
    const size_t multi_k = sd_k.multi;
    const size_t max_xi = std::max(multi_j, multi_k);
    const size_t min_xi = std::min(multi_j, multi_k);

    if (min_xi == 0) {
      continue;
    }

    const double prob = scaling_factor * static_cast<double>(max_xi) * k_ij;

    max_prob = std::max(prob, max_prob);
    min_prob = std::min(prob, min_prob);

    if (prob > 1.0) {
      big_probs++;
    }

    // Check for collision and coalesce if necessary
    if ((prob - std::floor(prob)) > phi) {
      const double gamma = std::floor(prob) + 1.0;
      if (sd_j.multi < sd_k.multi) {
        multi_coalesce(sd_k, sd_j, gamma);
      } else {
        multi_coalesce(sd_j, sd_k, gamma);
      }
      counter++;
    }
  }

  // Calculate total multiplicity and water mass
  double total_xi = 0.0;
  for (const auto &d : droplets) {
    total_xi += static_cast<double>(d.multi);
  }

  const double total_water_mass = total_water(droplets);

  if (config.debug) {
    // Debug output could be added here if needed
  }

  return CollisionStepResult{counter,  big_probs, max_prob,
                             min_prob, total_xi,  total_water_mass};
}

} // namespace physics
} // namespace sd_cpp
