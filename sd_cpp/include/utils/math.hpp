#pragma once

#include <random>
#include <vector>

namespace sd_cpp {
namespace utils {

/// Computes 2^n for integer exponents.
///
/// @param n The exponent
/// @return The result of 2^n as double
inline constexpr double pow2(size_t n) {
  return static_cast<double>(1ULL << n);
}

/// Computes the median of a vector of double values.
///
/// @param values The vector of values to compute the median of
/// @return The median of the values
double median(const std::vector<double> &values);

/// Performs an in-place Knuth shuffle (Fisher-Yates shuffle) on a vector.
///
/// @param v The vector to shuffle
/// @param rng Random number generator
template <typename T>
void knuth_shuffle(std::vector<T> &v, std::mt19937_64 &rng) {
  const size_t l = v.size();
  std::uniform_int_distribution<size_t> dist;

  for (size_t n = 0; n < l; ++n) {
    const size_t max_idx = l - n - 1;
    const size_t i = dist(
        rng, std::uniform_int_distribution<size_t>::param_type{0, max_idx});
    std::swap(v[i], v[max_idx]);
  }
}

/// Generates a linear grid of n points between start and stop.
///
/// @param start The start of the grid
/// @param stop The end of the grid
/// @param n The number of points in the grid
/// @return A vector of n evenly-spaced points between start and stop
std::vector<double> generate_linear_grid(double start, double stop, size_t n);

/// Returns the minimum of two double values.
///
/// @param a The first value
/// @param b The second value
/// @return The minimum of the two values
inline double min_f64(double a, double b) { return (a < b) ? a : b; }

/// Returns the maximum of two double values.
///
/// @param a The first value
/// @param b The second value
/// @return The maximum of the two values
inline double max_f64(double a, double b) { return (a > b) ? a : b; }

} // namespace utils
} // namespace sd_cpp
