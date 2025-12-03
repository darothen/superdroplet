#pragma once

#include <cmath>

namespace sd_cpp {
namespace constants {

/// Mathematical constant π
constexpr double PI = 3.1415926535897932384626433832;

/// Density of water (kg/m³)
constexpr double RHO_WATER = 1e3;

/// Density of air (kg/m³)
constexpr double RHO_AIR = 1.0;

/// One third (1/3)
constexpr double THIRD = 1.0 / 3.0;

/// Three fourths (3/4)
constexpr double THREE_FOURTH = 3.0 / 4.0;

/// Four thirds (4/3)
constexpr double FOUR_THIRD = 4.0 / 3.0;

/// Multiplicity threshold for special handling
constexpr double MULTI_THRESH = 1e4;

} // namespace constants
} // namespace sd_cpp

