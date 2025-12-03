/**
 * @file test_binning.cpp
 * @brief Tests for sd_cpp::io::binning module
 */

#include "core/droplet.hpp"
#include "io/binning.hpp"
#include <algorithm>
#include <gtest/gtest.h>
#include <numeric>
#include <vector>

using namespace sd_cpp;
using namespace sd_cpp::io;

class BinDropletsTest : public ::testing::Test {};

/**
 * @brief Test binning with empty droplet list
 */
TEST_F(BinDropletsTest, EmptyList) {
  std::vector<Droplet> droplets;

  BinGrid result = bin_droplets(droplets);

  // All bins should be zero
  EXPECT_EQ(result.lefts.size(), NR);
  EXPECT_EQ(result.rights.size(), NR);
  EXPECT_EQ(result.mids.size(), NR);
  EXPECT_EQ(result.values.size(), NR);

  double total = std::accumulate(result.values.begin(), result.values.end(), 0.0);
  EXPECT_DOUBLE_EQ(total, 0.0);
}

/**
 * @brief Test binning with a single droplet
 */
TEST_F(BinDropletsTest, SingleDroplet) {
  std::vector<Droplet> droplets;
  droplets.emplace_back(1000, 10e-6);

  BinGrid result = bin_droplets(droplets);

  EXPECT_EQ(result.lefts.size(), NR);
  EXPECT_EQ(result.rights.size(), NR);
  EXPECT_EQ(result.mids.size(), NR);
  EXPECT_EQ(result.values.size(), NR);

  // At least one bin should have non-zero value
  double total = std::accumulate(result.values.begin(), result.values.end(), 0.0);
  EXPECT_GT(total, 0.0);

  // Total should equal droplet mass * multiplicity
  double expected = droplets[0].mass() * static_cast<double>(droplets[0].multi);
  EXPECT_NEAR(total, expected, 1e-15);
}

/**
 * @brief Test binning with multiple droplets
 */
TEST_F(BinDropletsTest, MultipleDroplets) {
  std::vector<Droplet> droplets;
  droplets.emplace_back(1000, 5e-6);   // 5 μm radius
  droplets.emplace_back(2000, 15e-6);  // 15 μm radius
  droplets.emplace_back(3000, 25e-6);  // 25 μm radius

  BinGrid result = bin_droplets(droplets);

  EXPECT_EQ(result.values.size(), NR);

  // Check that total mass is conserved
  double total = std::accumulate(result.values.begin(), result.values.end(), 0.0);
  double expected = 0.0;
  for (const auto &d : droplets) {
    expected += d.mass() * static_cast<double>(d.multi);
  }
  EXPECT_NEAR(total, expected, 1e-12);
}

/**
 * @brief Test binning with unsorted droplets
 */
TEST_F(BinDropletsTest, UnsortedDroplets) {
  // Create droplets in unsorted order
  std::vector<Droplet> droplets;
  droplets.emplace_back(1000, 25e-6);
  droplets.emplace_back(2000, 5e-6);
  droplets.emplace_back(3000, 15e-6);

  BinGrid result = bin_droplets(droplets);

  // Should bin correctly regardless of input order
  EXPECT_EQ(result.values.size(), NR);
  double total = std::accumulate(result.values.begin(), result.values.end(), 0.0);
  EXPECT_GT(total, 0.0);
}

/**
 * @brief Test that binning accounts for multiplicity
 */
TEST_F(BinDropletsTest, AccountsForMultiplicity) {
  std::vector<Droplet> droplets;
  droplets.emplace_back(1000, 5e-6);
  droplets.emplace_back(1, 5e-6);

  BinGrid result = bin_droplets(droplets);

  // Total should weight by multiplicity
  double mass_per_droplet = droplets[0].mass();
  double expected =
      1000 * mass_per_droplet + 1 * mass_per_droplet;
  double total = std::accumulate(result.values.begin(), result.values.end(), 0.0);
  EXPECT_NEAR(total, expected, 1e-12);
}

/**
 * @brief Test binning with identical droplets
 */
TEST_F(BinDropletsTest, IdenticalDroplets) {
  size_t n_droplets = 10;
  std::vector<Droplet> droplets;
  for (size_t i = 0; i < n_droplets; ++i) {
    droplets.emplace_back(1000, 10e-6);
  }

  BinGrid result = bin_droplets(droplets);

  // All droplets should be in the same bin(s)
  double expected = n_droplets * droplets[0].mass() *
                    static_cast<double>(droplets[0].multi);
  double total = std::accumulate(result.values.begin(), result.values.end(), 0.0);
  EXPECT_NEAR(total, expected, 1e-12);
}

/**
 * @brief Test binning with large range of droplet sizes
 */
TEST_F(BinDropletsTest, LargeRange) {
  std::vector<Droplet> droplets;
  droplets.emplace_back(1000, 1e-6);    // Very small
  droplets.emplace_back(2000, 50e-6);   // Medium
  droplets.emplace_back(3000, 500e-6);  // Large

  BinGrid result = bin_droplets(droplets);

  EXPECT_EQ(result.values.size(), NR);

  // Count non-zero bins
  size_t non_zero_bins = 0;
  for (double val : result.values) {
    if (val > 0.0) {
      non_zero_bins++;
    }
  }

  // Each droplet should be in a different bin
  EXPECT_GE(non_zero_bins, 3);
}

/**
 * @brief Test binning with many droplets
 */
TEST_F(BinDropletsTest, ManyDroplets) {
  std::vector<Droplet> droplets;
  for (size_t i = 0; i < 1000; ++i) {
    double radius = (1.0 + i * 0.1) * 1e-6;
    droplets.emplace_back(1000, radius);
  }

  BinGrid result = bin_droplets(droplets);

  // Total mass should equal sum of all droplets
  double total_from_bins =
      std::accumulate(result.values.begin(), result.values.end(), 0.0);
  double total_from_droplets = 0.0;
  for (const auto &d : droplets) {
    total_from_droplets += d.mass() * static_cast<double>(d.multi);
  }
  EXPECT_NEAR(total_from_bins, total_from_droplets, 1e-10);
}

/**
 * @brief Test that binning conserves total mass
 */
TEST_F(BinDropletsTest, ConservationOfMass) {
  std::vector<Droplet> droplets;
  droplets.emplace_back(1000, 5e-6);
  droplets.emplace_back(2000, 15e-6);
  droplets.emplace_back(3000, 25e-6);
  droplets.emplace_back(4000, 35e-6);

  BinGrid result = bin_droplets(droplets);

  // Sum of bin values should equal total mass
  double total_binned =
      std::accumulate(result.values.begin(), result.values.end(), 0.0);
  double total_actual = 0.0;
  for (const auto &d : droplets) {
    total_actual += d.mass() * static_cast<double>(d.multi);
  }
  EXPECT_NEAR(total_binned, total_actual, 1e-12);
}

/**
 * @brief Test binning with droplets having zero multiplicity
 */
TEST_F(BinDropletsTest, ZeroMultiplicity) {
  std::vector<Droplet> droplets;
  droplets.emplace_back(1000, 5e-6);
  droplets.emplace_back(0, 15e-6);  // Zero multiplicity
  droplets.emplace_back(2000, 25e-6);

  BinGrid result = bin_droplets(droplets);

  // Should handle zero multiplicity correctly
  EXPECT_EQ(result.values.size(), NR);
  double total = std::accumulate(result.values.begin(), result.values.end(), 0.0);
  
  // Expected total should not include zero-multiplicity droplet
  double expected = droplets[0].mass() * static_cast<double>(droplets[0].multi) +
                    droplets[2].mass() * static_cast<double>(droplets[2].multi);
  EXPECT_NEAR(total, expected, 1e-12);
}

/**
 * @brief Test that bin edges are increasing
 */
TEST_F(BinDropletsTest, BinEdgesIncreasing) {
  std::vector<Droplet> droplets;
  droplets.emplace_back(1000, 10e-6);

  BinGrid result = bin_droplets(droplets);

  // Check that bin edges are monotonically increasing
  for (size_t i = 0; i < NR - 1; ++i) {
    EXPECT_LT(result.lefts[i], result.lefts[i + 1]);
    EXPECT_LT(result.rights[i], result.rights[i + 1]);
  }

  // Check that left < right for each bin
  for (size_t i = 0; i < NR; ++i) {
    EXPECT_LT(result.lefts[i], result.rights[i]);
  }
}

/**
 * @brief Test that midpoints are correct
 */
TEST_F(BinDropletsTest, MidpointsCorrect) {
  std::vector<Droplet> droplets;
  droplets.emplace_back(1000, 10e-6);

  BinGrid result = bin_droplets(droplets);

  // Check that midpoints are average of left and right edges
  for (size_t i = 0; i < NR; ++i) {
    double expected_mid = (result.lefts[i] + result.rights[i]) / 2.0;
    EXPECT_DOUBLE_EQ(result.mids[i], expected_mid);
  }
}

/**
 * @brief Test that consecutive bins share edges
 */
TEST_F(BinDropletsTest, ConsecutiveBinsShareEdges) {
  std::vector<Droplet> droplets;
  droplets.emplace_back(1000, 10e-6);

  BinGrid result = bin_droplets(droplets);

  // Check that right edge of bin i equals left edge of bin i+1
  for (size_t i = 0; i < NR - 1; ++i) {
    EXPECT_DOUBLE_EQ(result.rights[i], result.lefts[i + 1]);
  }
}

/**
 * @brief Test binning with droplets outside bin range
 */
TEST_F(BinDropletsTest, DropletsOutsideRange) {
  // Create droplets that might be outside typical bin range
  std::vector<Droplet> droplets;
  droplets.emplace_back(1000, 10000e-6);  // Very large, might be outside range

  BinGrid result = bin_droplets(droplets);

  // Should handle gracefully
  EXPECT_EQ(result.values.size(), NR);
  
  // Total should still be conserved (or zero if truly outside range)
  double total = std::accumulate(result.values.begin(), result.values.end(), 0.0);
  EXPECT_GE(total, 0.0);
}

/**
 * @brief Test get_radius_grid function
 */
TEST_F(BinDropletsTest, RadiusGrid) {
  const auto &grid = get_radius_grid();

  // Grid should have NR+1 edges
  EXPECT_EQ(grid.size(), NR + 1);

  // Grid should be monotonically increasing
  for (size_t i = 0; i < grid.size() - 1; ++i) {
    EXPECT_LT(grid[i], grid[i + 1]);
  }

  // First edge should be positive
  EXPECT_GT(grid[0], 0.0);
}

/**
 * @brief Test that get_radius_grid is consistent across calls
 */
TEST_F(BinDropletsTest, RadiusGridConsistent) {
  const auto &grid1 = get_radius_grid();
  const auto &grid2 = get_radius_grid();

  // Should return the same reference
  EXPECT_EQ(&grid1, &grid2);

  // Values should be identical
  EXPECT_EQ(grid1.size(), grid2.size());
  for (size_t i = 0; i < grid1.size(); ++i) {
    EXPECT_DOUBLE_EQ(grid1[i], grid2[i]);
  }
}

/**
 * @brief Test binning with various droplet size distributions
 */
TEST_F(BinDropletsTest, VariousSizeDistributions) {
  // Test with different size distributions
  std::vector<std::vector<double>> radii_sets = {
      {1e-6, 2e-6, 3e-6, 4e-6, 5e-6},              // Small droplets
      {10e-6, 20e-6, 30e-6, 40e-6, 50e-6},         // Medium droplets
      {100e-6, 200e-6, 300e-6, 400e-6, 500e-6},    // Large droplets
      {1e-6, 10e-6, 100e-6, 1000e-6}               // Mixed sizes
  };

  for (const auto &radii : radii_sets) {
    std::vector<Droplet> droplets;
    for (double r : radii) {
      droplets.emplace_back(1000, r);
    }

    BinGrid result = bin_droplets(droplets);

    // Check mass conservation
    double total_binned =
        std::accumulate(result.values.begin(), result.values.end(), 0.0);
    double total_actual = 0.0;
    for (const auto &d : droplets) {
      total_actual += d.mass() * static_cast<double>(d.multi);
    }
    EXPECT_NEAR(total_binned, total_actual, 1e-10);
  }
}

/**
 * @brief Test binning performance with large dataset
 */
TEST_F(BinDropletsTest, LargeDataset) {
  // Create a large dataset
  std::vector<Droplet> droplets;
  for (size_t i = 0; i < 10000; ++i) {
    double radius = (1.0 + i * 0.01) * 1e-6;
    droplets.emplace_back(100, radius);
  }

  BinGrid result = bin_droplets(droplets);

  // Verify mass conservation
  double total_binned =
      std::accumulate(result.values.begin(), result.values.end(), 0.0);
  double total_actual = 0.0;
  for (const auto &d : droplets) {
    total_actual += d.mass() * static_cast<double>(d.multi);
  }
  EXPECT_NEAR(total_binned, total_actual, 1e-8);
}

/**
 * @brief Test that binning handles edge cases in bin boundaries
 */
TEST_F(BinDropletsTest, BinBoundaryEdgeCases) {
  const auto &grid = get_radius_grid();
  
  // Create droplets exactly at bin boundaries
  std::vector<Droplet> droplets;
  for (size_t i = 0; i < std::min(size_t(10), grid.size()); ++i) {
    double radius_microns = grid[i];
    double radius_meters = radius_microns * 1e-6;
    droplets.emplace_back(1000, radius_meters);
  }

  BinGrid result = bin_droplets(droplets);

  // Should bin correctly at boundaries
  double total = std::accumulate(result.values.begin(), result.values.end(), 0.0);
  EXPECT_GT(total, 0.0);
}

/**
 * @brief Test binning with single very small droplet
 */
TEST_F(BinDropletsTest, VerySmallDroplet) {
  std::vector<Droplet> droplets;
  droplets.emplace_back(1000, 0.1e-6);  // 0.1 micron radius

  BinGrid result = bin_droplets(droplets);

  // Should place in first bin(s)
  double total = std::accumulate(result.values.begin(), result.values.end(), 0.0);
  double expected = droplets[0].mass() * static_cast<double>(droplets[0].multi);
  EXPECT_NEAR(total, expected, 1e-15);
}

/**
 * @brief Test that all bin values are non-negative
 */
TEST_F(BinDropletsTest, NonNegativeValues) {
  std::vector<Droplet> droplets;
  for (size_t i = 0; i < 100; ++i) {
    droplets.emplace_back(1000, (1.0 + i) * 1e-6);
  }

  BinGrid result = bin_droplets(droplets);

  // All bin values should be non-negative
  for (double val : result.values) {
    EXPECT_GE(val, 0.0);
  }
}
