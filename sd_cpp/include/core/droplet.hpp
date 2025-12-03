#pragma once

#include <cstddef>
#include <vector>

namespace sd_cpp {

/// Represents a superdroplet in the simulation.
///
/// A superdroplet is a computational particle that represents many
/// identical physical droplets. The `multi` field indicates how many
/// real droplets this superdroplet represents.
class Droplet {
public:
  /// Multiplicity - number of real droplets represented
  size_t multi;

  /// Radius cubed (m³) - stored for efficiency
  double rcubed;

  /// Solute mass (kg)
  double solute;

  /// Creates a new droplet with the given multiplicity and radius.
  ///
  /// @param multi Number of real droplets this superdroplet represents
  /// @param radius Radius of the droplet in meters
  Droplet(size_t multi, double radius);

  /// Updates the radius cubed and recalculates cached properties.
  ///
  /// This method should be called whenever the droplet size changes
  /// (e.g., after coalescence).
  ///
  /// @param new_rcubed New value for radius cubed (m³)
  void update_rcubed(double new_rcubed);

  /// Returns the volume of the droplet (m³).
  [[nodiscard]] inline double volume() const { return volume_; }

  /// Returns the radius of the droplet (m).
  [[nodiscard]] inline double radius() const { return radius_; }

  /// Returns the mass of the droplet (kg).
  [[nodiscard]] inline double mass() const { return mass_; }

  /// Returns the terminal velocity of the droplet (m/s).
  [[nodiscard]] inline double terminal_velocity() const {
    return terminal_velocity_;
  }

private:
  /// Cached radius (m)
  double radius_;

  /// Cached volume (m³)
  double volume_;

  /// Cached mass (kg)
  double mass_;

  /// Cached terminal velocity (m/s)
  double terminal_velocity_;

  /// Computes terminal velocity using Beard (1976) parameterization.
  ///
  /// @param radius Droplet radius in meters
  /// @param mass Droplet mass in kilograms
  /// @return Terminal velocity in m/s
  static double compute_terminal_velocity(double radius, double mass);
};

/// Calculates the total water mass in a collection of droplets.
///
/// @param droplets Vector of droplets to sum
/// @return Total water mass in kilograms
double total_water(const std::vector<Droplet> &droplets);

} // namespace sd_cpp
