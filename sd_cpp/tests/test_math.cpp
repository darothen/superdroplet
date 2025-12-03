/**
 * @file test_math.cpp
 * @brief Tests for sd_cpp::utils::math module
 */

#include "utils/math.hpp"
#include <algorithm>
#include <gtest/gtest.h>
#include <random>
#include <vector>

using namespace sd_cpp::utils;

class Pow2Test : public ::testing::Test {};

/**
 * @brief Test pow2 function with various inputs
 */
TEST_F(Pow2Test, BasicPowers) {
  EXPECT_DOUBLE_EQ(pow2(0), 1.0);
  EXPECT_DOUBLE_EQ(pow2(1), 2.0);
  EXPECT_DOUBLE_EQ(pow2(2), 4.0);
  EXPECT_DOUBLE_EQ(pow2(3), 8.0);
  EXPECT_DOUBLE_EQ(pow2(10), 1024.0);
}

class MedianTest : public ::testing::Test {};

/**
 * @brief Test median with odd number of elements
 */
TEST_F(MedianTest, OddElements) {
  std::vector<double> values = {5.0, 1.0, 3.0, 2.0, 4.0};
  double result = median(values);
  EXPECT_DOUBLE_EQ(result, 3.0);
}

/**
 * @brief Test median with even number of elements
 */
TEST_F(MedianTest, EvenElements) {
  std::vector<double> values = {5.0, 1.0, 3.0, 2.0, 4.0, 6.0};
  double result = median(values);
  EXPECT_DOUBLE_EQ(result, 3.5); // Average of 3 and 4
}

/**
 * @brief Test median with single element
 */
TEST_F(MedianTest, SingleElement) {
  std::vector<double> values = {5.0};
  double result = median(values);
  EXPECT_DOUBLE_EQ(result, 5.0);
}

/**
 * @brief Test median with two elements
 */
TEST_F(MedianTest, TwoElements) {
  std::vector<double> values = {3.0, 7.0};
  double result = median(values);
  EXPECT_DOUBLE_EQ(result, 5.0); // Average of 3 and 7
}

/**
 * @brief Test median with identical elements
 */
TEST_F(MedianTest, IdenticalElements) {
  std::vector<double> values = {5.0, 5.0, 5.0, 5.0, 5.0};
  double result = median(values);
  EXPECT_DOUBLE_EQ(result, 5.0);
}

/**
 * @brief Test median with already sorted elements
 */
TEST_F(MedianTest, SortedElements) {
  std::vector<double> values = {1.0, 2.0, 3.0, 4.0, 5.0};
  double result = median(values);
  EXPECT_DOUBLE_EQ(result, 3.0);
}

class KnuthShuffleTest : public ::testing::Test {};

/**
 * @brief Test that shuffle changes order
 */
TEST_F(KnuthShuffleTest, ChangesOrder) {
  std::vector<int> values(100);
  for (size_t i = 0; i < values.size(); ++i) {
    values[i] = static_cast<int>(i);
  }

  std::vector<int> original = values;
  std::mt19937_64 rng(42);
  knuth_shuffle(values, rng);

  // With high probability, shuffled should differ from original
  bool different = false;
  for (size_t i = 0; i < values.size(); ++i) {
    if (values[i] != original[i]) {
      different = true;
      break;
    }
  }
  EXPECT_TRUE(different);
}

/**
 * @brief Test that shuffle preserves all elements
 */
TEST_F(KnuthShuffleTest, PreservesElements) {
  std::vector<int> values = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  std::vector<int> original = values;

  std::mt19937_64 rng(42);
  knuth_shuffle(values, rng);

  std::sort(values.begin(), values.end());
  std::sort(original.begin(), original.end());

  EXPECT_EQ(values, original);
}

/**
 * @brief Test shuffle with single element
 */
TEST_F(KnuthShuffleTest, SingleElement) {
  std::vector<int> values = {42};
  std::mt19937_64 rng(42);
  knuth_shuffle(values, rng);

  EXPECT_EQ(values.size(), 1);
  EXPECT_EQ(values[0], 42);
}

/**
 * @brief Test shuffle is deterministic with same seed
 */
TEST_F(KnuthShuffleTest, Deterministic) {
  std::vector<int> values1(20);
  std::vector<int> values2(20);
  for (size_t i = 0; i < 20; ++i) {
    values1[i] = static_cast<int>(i);
    values2[i] = static_cast<int>(i);
  }

  std::mt19937_64 rng1(12345);
  std::mt19937_64 rng2(12345);

  knuth_shuffle(values1, rng1);
  knuth_shuffle(values2, rng2);

  EXPECT_EQ(values1, values2);
}

class GenerateLinearGridTest : public ::testing::Test {};

/**
 * @brief Test generating a linear grid
 */
TEST_F(GenerateLinearGridTest, BasicGrid) {
  auto grid = generate_linear_grid(0.0, 10.0, 11);

  EXPECT_EQ(grid.size(), 11);
  EXPECT_DOUBLE_EQ(grid[0], 0.0);
  EXPECT_DOUBLE_EQ(grid[10], 10.0);

  // Check spacing
  for (size_t i = 0; i < grid.size() - 1; ++i) {
    EXPECT_NEAR(grid[i + 1] - grid[i], 1.0, 1e-10);
  }
}

/**
 * @brief Test grid with two points
 */
TEST_F(GenerateLinearGridTest, TwoPoints) {
  auto grid = generate_linear_grid(0.0, 10.0, 2);

  EXPECT_EQ(grid.size(), 2);
  EXPECT_DOUBLE_EQ(grid[0], 0.0);
  EXPECT_DOUBLE_EQ(grid[1], 10.0);
}

/**
 * @brief Test grid with single point
 */
TEST_F(GenerateLinearGridTest, SinglePoint) {
  auto grid = generate_linear_grid(5.0, 5.0, 1);

  EXPECT_EQ(grid.size(), 1);
  EXPECT_DOUBLE_EQ(grid[0], 5.0);
}

/**
 * @brief Test grid with negative range
 */
TEST_F(GenerateLinearGridTest, NegativeRange) {
  auto grid = generate_linear_grid(-5.0, 5.0, 11);

  EXPECT_EQ(grid.size(), 11);
  EXPECT_DOUBLE_EQ(grid[0], -5.0);
  EXPECT_DOUBLE_EQ(grid[10], 5.0);
}

class MinMaxTest : public ::testing::Test {};

/**
 * @brief Test min_f64 function
 */
TEST_F(MinMaxTest, MinF64) {
  EXPECT_DOUBLE_EQ(min_f64(3.0, 5.0), 3.0);
  EXPECT_DOUBLE_EQ(min_f64(5.0, 3.0), 3.0);
  EXPECT_DOUBLE_EQ(min_f64(-1.0, 1.0), -1.0);
  EXPECT_DOUBLE_EQ(min_f64(0.0, 0.0), 0.0);
}

/**
 * @brief Test max_f64 function
 */
TEST_F(MinMaxTest, MaxF64) {
  EXPECT_DOUBLE_EQ(max_f64(3.0, 5.0), 5.0);
  EXPECT_DOUBLE_EQ(max_f64(5.0, 3.0), 5.0);
  EXPECT_DOUBLE_EQ(max_f64(-1.0, 1.0), 1.0);
  EXPECT_DOUBLE_EQ(max_f64(0.0, 0.0), 0.0);
}
