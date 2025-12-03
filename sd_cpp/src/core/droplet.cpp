#include "core/droplet.hpp"
#include "core/constants.hpp"
#include <cmath>
#include <numeric>

namespace sd_cpp {

using namespace constants;

Droplet::Droplet(size_t multi, double radius)
    : multi(multi), rcubed(radius * radius * radius), solute(0.0),
      radius_(radius) {
  volume_ = rcubed * PI * FOUR_THIRD;
  mass_ = volume_ * RHO_WATER;
  terminal_velocity_ = compute_terminal_velocity(radius_, mass_);
}

void Droplet::update_rcubed(double new_rcubed) {
  rcubed = new_rcubed;
  radius_ = std::cbrt(new_rcubed);
  volume_ = new_rcubed * PI * FOUR_THIRD;
  mass_ = volume_ * RHO_WATER;
  terminal_velocity_ = compute_terminal_velocity(radius_, mass_);
}

double Droplet::compute_terminal_velocity(double radius, double mass) {
  const double d = 2.0 * radius * 1e6; // diameter, m -> μm
  const double x = mass * 1e3;         // mass, kg -> g

  double alpha, x_to_beta;

  if (d <= 134.43) {
    // x^(2/3) = cbrt(x)^2
    const double cbrt_x = std::cbrt(x);
    alpha = 4.5795e5;
    x_to_beta = cbrt_x * cbrt_x;
  } else if (d <= 1511.64) {
    // x^(1/3) = cbrt(x)
    alpha = 4962.0;
    x_to_beta = std::cbrt(x);
  } else if (d <= 3477.84) {
    // x^(1/6) = sqrt(cbrt(x))
    alpha = 1732.0;
    x_to_beta = std::sqrt(std::cbrt(x));
  } else {
    alpha = 917.0;
    x_to_beta = 1.0;
  }

  return 1e-2 * alpha * x_to_beta; // cm/s -> m/s
}

double total_water(const std::vector<Droplet> &droplets) {
  return std::accumulate(droplets.begin(), droplets.end(), 0.0,
                         [](double sum, const Droplet &d) {
                           return sum + d.mass() * static_cast<double>(d.multi);
                         });
}

} // namespace sd_cpp
