/**
 * @file test_droplet.cpp
 * @brief Tests for sd_cpp::Droplet class
 */

#include "core/constants.hpp"
#include "core/droplet.hpp"
#include <gtest/gtest.h>
#include <vector>

using namespace sd_cpp;
using namespace sd_cpp::constants;

class DropletCreationTest : public ::testing::Test {};

/**
 * @brief Test creating a droplet with basic properties
 */
TEST_F(DropletCreationTest, NewDropletBasicProperties) {
  size_t multi = 1000;
  double radius = 10e-6; // 10 microns

  Droplet droplet(multi, radius);

  EXPECT_EQ(droplet.multi, multi);
  EXPECT_DOUBLE_EQ(droplet.radius(), radius);
  EXPECT_DOUBLE_EQ(droplet.solute, 0.0);
}

/**
 * @brief Test that computed properties are calculated correctly
 */
TEST_F(DropletCreationTest, NewDropletComputedProperties) {
  double radius = 10e-6; // 10 microns
  Droplet droplet(1, radius);

  double expected_rcubed = radius * radius * radius;
  double expected_volume = FOUR_THIRD * PI * expected_rcubed;
  double expected_mass = expected_volume * RHO_WATER;

  EXPECT_DOUBLE_EQ(droplet.rcubed, expected_rcubed);
  EXPECT_NEAR(droplet.volume(), expected_volume, 1e-25);
  EXPECT_NEAR(droplet.mass(), expected_mass, 1e-15);
}

/**
 * @brief Test that terminal velocity is computed
 */
TEST_F(DropletCreationTest, NewDropletTerminalVelocity) {
  Droplet droplet(1, 10e-6);
  EXPECT_GT(droplet.terminal_velocity(), 0.0);
}

class TerminalVelocityTest : public ::testing::Test {};

/**
 * @brief Test terminal velocity for small droplets (d <= 134.43 μm)
 */
TEST(TerminalVelocityTest, SmallDroplet) {
  double radius = 50e-6; // 50 microns, diameter = 100 μm
  Droplet droplet(1, radius);

  EXPECT_GT(droplet.terminal_velocity(), 0.0);
  EXPECT_LT(droplet.terminal_velocity(), 10.0); // Should be reasonable
}

/**
 * @brief Test terminal velocity for medium droplets (134.43 < d <= 1511.64 μm)
 */
TEST(TerminalVelocityTest, MediumDroplet) {
  double radius = 500e-6; // 500 microns, diameter = 1000 μm
  Droplet droplet(1, radius);

  EXPECT_GT(droplet.terminal_velocity(), 0.0);
  EXPECT_LT(droplet.terminal_velocity(), 10.0);
}

/**
 * @brief Test terminal velocity for large droplets (1511.64 < d <= 3477.84 μm)
 */
TEST(TerminalVelocityTest, LargeDroplet) {
  double radius = 1500e-6; // 1500 microns, diameter = 3000 μm
  Droplet droplet(1, radius);

  EXPECT_GT(droplet.terminal_velocity(), 0.0);
  EXPECT_LT(droplet.terminal_velocity(), 10.0);
}

/**
 * @brief Test terminal velocity for very large droplets (d > 3477.84 μm)
 */
TEST(TerminalVelocityTest, VeryLargeDroplet) {
  double radius = 2000e-6; // 2000 microns, diameter = 4000 μm
  Droplet droplet(1, radius);

  EXPECT_GT(droplet.terminal_velocity(), 0.0);
  EXPECT_LT(droplet.terminal_velocity(), 10.0);
}

/**
 * @brief Test that terminal velocity increases with size
 */
TEST(TerminalVelocityTest, IncreasesWithSize) {
  std::vector<double> radii = {10e-6, 50e-6, 100e-6, 500e-6};
  std::vector<double> terminal_velocities;

  for (double radius : radii) {
    Droplet droplet(1, radius);
    terminal_velocities.push_back(droplet.terminal_velocity());
  }

  // Check that terminal velocities generally increase
  for (size_t i = 0; i < terminal_velocities.size() - 1; ++i) {
    EXPECT_GE(terminal_velocities[i + 1], terminal_velocities[i]);
  }
}

class TotalWaterTest : public ::testing::Test {};

/**
 * @brief Test total water for empty vector
 */
TEST_F(TotalWaterTest, EmptyVector) {
  std::vector<Droplet> droplets;
  double total = total_water(droplets);
  EXPECT_DOUBLE_EQ(total, 0.0);
}

/**
 * @brief Test total water for a single droplet
 */
TEST_F(TotalWaterTest, SingleDroplet) {
  std::vector<Droplet> droplets;
  droplets.emplace_back(1000, 10e-6);

  double total = total_water(droplets);
  double expected = droplets[0].mass() * droplets[0].multi;
  EXPECT_DOUBLE_EQ(total, expected);
}

/**
 * @brief Test total water for multiple droplets
 */
TEST_F(TotalWaterTest, MultipleDroplets) {
  std::vector<Droplet> droplets;
  droplets.emplace_back(1000, 10e-6);
  droplets.emplace_back(2000, 20e-6);
  droplets.emplace_back(3000, 30e-6);

  double total = total_water(droplets);

  double expected = 0.0;
  for (const auto &d : droplets) {
    expected += d.mass() * d.multi;
  }
  EXPECT_NEAR(total, expected, 1e-15);
}

/**
 * @brief Test total water with high multiplicity droplets
 */
TEST_F(TotalWaterTest, HighMultiplicity) {
  std::vector<Droplet> droplets;
  droplets.emplace_back(10000, 5e-6);
  droplets.emplace_back(5000, 20e-6);
  droplets.emplace_back(1000, 100e-6);

  double total = total_water(droplets);

  double expected = 0.0;
  for (const auto &d : droplets) {
    expected += d.mass() * d.multi;
  }
  EXPECT_NEAR(total, expected, 1e-12);
}

/**
 * @brief Test total water with a large list of droplets
 */
TEST_F(TotalWaterTest, LargeList) {
  std::vector<Droplet> droplets;
  for (size_t i = 0; i < 1000; ++i) {
    droplets.emplace_back(100, (1.0 + i * 0.01) * 1e-6);
  }

  double total = total_water(droplets);

  double expected = 0.0;
  for (const auto &d : droplets) {
    expected += d.mass() * d.multi;
  }
  EXPECT_NEAR(total, expected, 1e-12);
}

class DropletUpdateRcubedTest : public ::testing::Test {};

/**
 * @brief Test that update_rcubed increases derived properties
 */
TEST_F(DropletUpdateRcubedTest, IncreasesSize) {
  Droplet droplet(1000, 10e-6);
  double original_radius = droplet.radius();
  double original_volume = droplet.volume();
  double original_mass = droplet.mass();

  // Double the rcubed
  double new_rcubed = droplet.rcubed * 2.0;
  droplet.update_rcubed(new_rcubed);

  EXPECT_DOUBLE_EQ(droplet.rcubed, new_rcubed);
  EXPECT_GT(droplet.radius(), original_radius);
  EXPECT_GT(droplet.volume(), original_volume);
  EXPECT_GT(droplet.mass(), original_mass);
}

/**
 * @brief Test that update_rcubed updates terminal velocity
 */
TEST_F(DropletUpdateRcubedTest, UpdatesTerminalVelocity) {
  Droplet droplet(1000, 10e-6);
  double original_tv = droplet.terminal_velocity();

  // Increase the size significantly
  double new_rcubed = droplet.rcubed * 8.0; // Double the radius cubed
  droplet.update_rcubed(new_rcubed);

  // Terminal velocity should change
  EXPECT_NE(droplet.terminal_velocity(), original_tv);
}

/**
 * @brief Test that properties remain consistent after update
 */
TEST_F(DropletUpdateRcubedTest, ConsistentProperties) {
  Droplet droplet(1000, 10e-6);
  double new_radius = 20e-6; // New radius = 20 microns
  double new_rcubed = new_radius * new_radius * new_radius;

  droplet.update_rcubed(new_rcubed);

  // Check consistency
  EXPECT_DOUBLE_EQ(droplet.rcubed, new_rcubed);
  EXPECT_NEAR(droplet.volume(), FOUR_THIRD * PI * new_rcubed, 1e-20);
  EXPECT_NEAR(droplet.mass(), droplet.volume() * RHO_WATER, 1e-15);
}
