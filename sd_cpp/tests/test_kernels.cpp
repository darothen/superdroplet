/**
 * @file test_kernels.cpp
 * @brief Tests for sd_cpp::physics::kernels module
 */

#include "core/constants.hpp"
#include "core/droplet.hpp"
#include "physics/kernels.hpp"
#include <cmath>
#include <gtest/gtest.h>

using namespace sd_cpp;
using namespace sd_cpp::physics;
using namespace sd_cpp::constants;

class GolovinKernelTest : public ::testing::Test {
protected:
  Droplet small_droplet{1, 10e-6};
  Droplet medium_droplet{1, 50e-6};
  Droplet large_droplet{1, 100e-6};
};

/**
 * @brief Test Golovin kernel basic computation
 */
TEST_F(GolovinKernelTest, BasicComputation) {
  double result = golovin_kernel(small_droplet, medium_droplet);

  // Expected: GOLOVIN_B * (4/3) * PI * (rcubed_1 + rcubed_2)
  constexpr double GOLOVIN_B = 1500.0;
  double expected = GOLOVIN_B * FOUR_THIRD * PI *
                    (small_droplet.rcubed + medium_droplet.rcubed);

  EXPECT_NEAR(result, expected, 1e-15);
}

/**
 * @brief Test that Golovin kernel is symmetric
 */
TEST_F(GolovinKernelTest, Symmetric) {
  double result1 = golovin_kernel(small_droplet, large_droplet);
  double result2 = golovin_kernel(large_droplet, small_droplet);
  EXPECT_DOUBLE_EQ(result1, result2);
}

/**
 * @brief Test that Golovin kernel returns positive values
 */
TEST_F(GolovinKernelTest, PositiveValues) {
  double result = golovin_kernel(small_droplet, medium_droplet);
  EXPECT_GT(result, 0.0);
}

/**
 * @brief Test that Golovin kernel increases with droplet size
 */
TEST_F(GolovinKernelTest, IncreasesWithSize) {
  Droplet d1(1, 10e-6);
  Droplet d2(1, 50e-6);
  Droplet d3(1, 100e-6);

  double k12 = golovin_kernel(d1, d2);
  double k23 = golovin_kernel(d2, d3);

  EXPECT_GT(k23, k12);
}

class HydroKernelTest : public ::testing::Test {
protected:
  Droplet small_droplet{1, 10e-6};
  Droplet medium_droplet{1, 50e-6};
  Droplet large_droplet{1, 200e-6};
};

/**
 * @brief Test basic hydrodynamic kernel computation
 */
TEST_F(HydroKernelTest, BasicComputation) {
  double result = hydro_kernel(small_droplet, large_droplet);
  EXPECT_GT(result, 0.0);
}

/**
 * @brief Test that hydro kernel is symmetric in droplet order
 */
TEST_F(HydroKernelTest, Symmetric) {
  double result1 = hydro_kernel(small_droplet, large_droplet);
  double result2 = hydro_kernel(large_droplet, small_droplet);
  // Due to abs(tv_diff), results should be the same
  EXPECT_DOUBLE_EQ(result1, result2);
}

/**
 * @brief Test hydro kernel for droplets of the same size
 */
TEST_F(HydroKernelTest, SameSizeDroplets) {
  Droplet d1(1, 50e-6);
  Droplet d2(1, 50e-6);

  double result = hydro_kernel(d1, d2);
  // Same size droplets should have same terminal velocity, so kernel ~0
  EXPECT_NEAR(result, 0.0, 1e-10);
}

/**
 * @brief Test that hydro kernel increases with size difference
 */
TEST_F(HydroKernelTest, IncreasesWithSizeDifference) {
  Droplet d_small(1, 10e-6);
  Droplet d_medium(1, 50e-6);
  Droplet d_large(1, 200e-6);

  double k_small_medium = hydro_kernel(d_small, d_medium);
  double k_small_large = hydro_kernel(d_small, d_large);

  EXPECT_GT(k_small_large, k_small_medium);
}

class LongKernelTest : public ::testing::Test {
protected:
  Droplet small_droplet{1, 10e-6};
  Droplet medium_droplet{1, 30e-6};
  Droplet large_droplet{1, 200e-6};
};

/**
 * @brief Test basic Long kernel computation
 */
TEST_F(LongKernelTest, BasicComputation) {
  double result = long_kernel(small_droplet, large_droplet);
  EXPECT_GT(result, 0.0);
}

/**
 * @brief Test that Long kernel is symmetric in droplet order
 */
TEST_F(LongKernelTest, Symmetric) {
  double result1 = long_kernel(small_droplet, large_droplet);
  double result2 = long_kernel(large_droplet, small_droplet);
  EXPECT_DOUBLE_EQ(result1, result2);
}

/**
 * @brief Test Long kernel for droplets of the same size
 */
TEST_F(LongKernelTest, SameSizeDroplets) {
  Droplet d1(1, 50e-6);
  Droplet d2(1, 50e-6);

  double result = long_kernel(d1, d2);
  // Same size droplets should have same terminal velocity, so kernel ~0
  EXPECT_NEAR(result, 0.0, 1e-10);
}

/**
 * @brief Test collection efficiency for large droplets (>= 50 μm)
 */
TEST_F(LongKernelTest, CollectionEfficiencyLargeDroplet) {
  Droplet d_small(1, 10e-6);
  Droplet d_large(1, 60e-6); // > 50 μm

  double result = long_kernel(d_small, d_large);
  EXPECT_GT(result, 0.0);

  // For large droplets, e_coll should be 1.0
  // So Long kernel should equal hydro kernel
  double hydro_result = hydro_kernel(d_small, d_large);
  EXPECT_NEAR(result, hydro_result, 1e-15);
}

/**
 * @brief Test collection efficiency for small droplets (< 50 μm)
 */
TEST_F(LongKernelTest, CollectionEfficiencySmallDroplets) {
  Droplet d1(1, 10e-6);
  Droplet d2(1, 30e-6);

  double result = long_kernel(d1, d2);

  // For small droplets, e_coll is calculated
  // Result should be less than or equal to hydro kernel
  double hydro_result = hydro_kernel(d1, d2);
  EXPECT_GE(result, 0.0);
  EXPECT_LE(result, hydro_result + 1e-15);
}

/**
 * @brief Test that Long kernel is generally <= hydro kernel
 */
TEST_F(LongKernelTest, LessThanOrEqualHydroKernel) {
  std::vector<double> radii = {10e-6, 20e-6, 30e-6, 40e-6, 60e-6, 100e-6};

  for (double r1 : radii) {
    for (double r2 : radii) {
      if (r1 == r2)
        continue;

      Droplet d1(1, r1);
      Droplet d2(1, r2);

      double long_result = long_kernel(d1, d2);
      double hydro_result = hydro_kernel(d1, d2);

      // Long kernel should be <= hydro kernel due to collection efficiency
      EXPECT_LE(long_result, hydro_result + 1e-15);
    }
  }
}

/**
 * @brief Test r_small and r_large assignment (bug fix verification)
 */
TEST_F(LongKernelTest, RSmallRLargeAssignment) {
  // This test verifies the bug fix where r_small and r_large
  // were correctly assigned using min/max

  Droplet d_small(1, 10e-6);  // 10 μm
  Droplet d_large(1, 100e-6); // 100 μm

  // Both orderings should give the same result
  double result1 = long_kernel(d_small, d_large);
  double result2 = long_kernel(d_large, d_small);

  // Results should be identical (symmetry)
  EXPECT_DOUBLE_EQ(result1, result2);

  // For large droplets (>= 50 μm), e_coll should be 1.0
  // So Long kernel should equal hydro kernel
  double hydro_result = hydro_kernel(d_small, d_large);
  EXPECT_NEAR(result1, hydro_result, 1e-15);
}

/**
 * @brief Test Long kernel with small collector, large collected
 */
TEST_F(LongKernelTest, SmallCollectorLargeCollected) {
  // r_small < 50 μm, r_large >= 50 μm case
  Droplet d1(1, 20e-6); // 20 μm (small)
  Droplet d2(1, 60e-6); // 60 μm (large)

  double result = long_kernel(d1, d2);

  // Since r_large >= 50, e_coll = 1.0
  double hydro_result = hydro_kernel(d1, d2);
  EXPECT_NEAR(result, hydro_result, 1e-15);
}

class KernelFunctionTest : public ::testing::Test {};

/**
 * @brief Test get_kernel_function returns correct function for Golovin
 */
TEST_F(KernelFunctionTest, GetGolovinKernel) {
  KernelFn fn = get_kernel_function(Kernel::Golovin);

  Droplet d1(1, 10e-6);
  Droplet d2(1, 20e-6);

  double result = fn(d1, d2);
  double expected = golovin_kernel(d1, d2);

  EXPECT_DOUBLE_EQ(result, expected);
}

/**
 * @brief Test get_kernel_function returns correct function for Hydro
 */
TEST_F(KernelFunctionTest, GetHydroKernel) {
  KernelFn fn = get_kernel_function(Kernel::Hydro);

  Droplet d1(1, 10e-6);
  Droplet d2(1, 20e-6);

  double result = fn(d1, d2);
  double expected = hydro_kernel(d1, d2);

  EXPECT_DOUBLE_EQ(result, expected);
}

/**
 * @brief Test get_kernel_function returns correct function for Long
 */
TEST_F(KernelFunctionTest, GetLongKernel) {
  KernelFn fn = get_kernel_function(Kernel::Long);

  Droplet d1(1, 10e-6);
  Droplet d2(1, 20e-6);

  double result = fn(d1, d2);
  double expected = long_kernel(d1, d2);

  EXPECT_DOUBLE_EQ(result, expected);
}

/**
 * @brief Test kernel_name returns correct names
 */
TEST_F(KernelFunctionTest, KernelNames) {
  EXPECT_STREQ(kernel_name(Kernel::Golovin), "Golovin");
  EXPECT_STREQ(kernel_name(Kernel::Hydro), "Hydro");
  EXPECT_STREQ(kernel_name(Kernel::Long), "Long");
}
