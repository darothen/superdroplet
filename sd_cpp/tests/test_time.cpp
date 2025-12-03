/**
 * @file test_time.cpp
 * @brief Tests for sd_cpp::utils::Stopwatch class
 */

#include "utils/time.hpp"
#include <gtest/gtest.h>

using namespace sd_cpp::utils;

class StopwatchCreationTest : public ::testing::Test {};

/**
 * @brief Test creating stopwatch with 0 seconds
 */
TEST_F(StopwatchCreationTest, ZeroSeconds) {
  Stopwatch sw(0);
  EXPECT_EQ(sw.minutes, 0);
  EXPECT_EQ(sw.seconds, 0);
}

/**
 * @brief Test creating stopwatch with less than a minute
 */
TEST_F(StopwatchCreationTest, LessThanMinute) {
  Stopwatch sw(45);
  EXPECT_EQ(sw.minutes, 0);
  EXPECT_EQ(sw.seconds, 45);
}

/**
 * @brief Test creating stopwatch with exactly 60 seconds
 */
TEST_F(StopwatchCreationTest, ExactlyOneMinute) {
  Stopwatch sw(60);
  EXPECT_EQ(sw.minutes, 1);
  EXPECT_EQ(sw.seconds, 0);
}

/**
 * @brief Test creating stopwatch with minutes and seconds
 */
TEST_F(StopwatchCreationTest, MinutesAndSeconds) {
  Stopwatch sw(90);
  EXPECT_EQ(sw.minutes, 1);
  EXPECT_EQ(sw.seconds, 30);
}

/**
 * @brief Test creating stopwatch with multiple minutes
 */
TEST_F(StopwatchCreationTest, MultipleMinutes) {
  Stopwatch sw(185); // 3 minutes and 5 seconds
  EXPECT_EQ(sw.minutes, 3);
  EXPECT_EQ(sw.seconds, 5);
}

/**
 * @brief Test creating stopwatch with large value
 */
TEST_F(StopwatchCreationTest, LargeValue) {
  Stopwatch sw(3661); // 61 minutes and 1 second
  EXPECT_EQ(sw.minutes, 61);
  EXPECT_EQ(sw.seconds, 1);
}

class StopwatchIncrementTest : public ::testing::Test {};

/**
 * @brief Test basic increment operation
 */
TEST_F(StopwatchIncrementTest, BasicIncrement) {
  Stopwatch sw(0);
  sw.increment(10);
  EXPECT_EQ(sw.seconds, 10);
  EXPECT_EQ(sw.minutes, 0);
}

/**
 * @brief Test increment without overflow to minutes
 */
TEST_F(StopwatchIncrementTest, NoOverflow) {
  Stopwatch sw(30);
  sw.increment(15);
  EXPECT_EQ(sw.seconds, 45);
  EXPECT_EQ(sw.minutes, 0);
}

/**
 * @brief Test increment with overflow to minutes
 */
TEST_F(StopwatchIncrementTest, WithOverflow) {
  Stopwatch sw(50);
  sw.increment(20);
  EXPECT_EQ(sw.seconds, 10);
  EXPECT_EQ(sw.minutes, 1);
}

/**
 * @brief Test increment that results in exactly a minute boundary
 */
TEST_F(StopwatchIncrementTest, ExactlyToMinute) {
  Stopwatch sw(40);
  sw.increment(20);
  EXPECT_EQ(sw.seconds, 0);
  EXPECT_EQ(sw.minutes, 1);
}

/**
 * @brief Test increment that adds multiple minutes
 */
TEST_F(StopwatchIncrementTest, MultipleMinutes) {
  Stopwatch sw(30);
  sw.increment(150); // 2 minutes and 30 seconds
  EXPECT_EQ(sw.seconds, 0);
  EXPECT_EQ(sw.minutes, 3);
}

/**
 * @brief Test increment from stopwatch already having minutes
 */
TEST_F(StopwatchIncrementTest, FromExistingMinutes) {
  Stopwatch sw(120); // 2 minutes
  sw.increment(45);
  EXPECT_EQ(sw.seconds, 45);
  EXPECT_EQ(sw.minutes, 2);
}

/**
 * @brief Test multiple increments
 */
TEST_F(StopwatchIncrementTest, MultipleIncrements) {
  Stopwatch sw(0);
  sw.increment(30);
  sw.increment(30);
  sw.increment(30);
  EXPECT_EQ(sw.seconds, 30);
  EXPECT_EQ(sw.minutes, 1);
}

/**
 * @brief Test increment by zero
 */
TEST_F(StopwatchIncrementTest, IncrementZero) {
  Stopwatch sw(45);
  sw.increment(0);
  EXPECT_EQ(sw.seconds, 45);
  EXPECT_EQ(sw.minutes, 0);
}

class StopwatchTotalSecondsTest : public ::testing::Test {};

/**
 * @brief Test total_seconds for zero time
 */
TEST_F(StopwatchTotalSecondsTest, ZeroTime) {
  Stopwatch sw(0);
  EXPECT_EQ(sw.total_seconds(), 0);
}

/**
 * @brief Test total_seconds with only seconds
 */
TEST_F(StopwatchTotalSecondsTest, OnlySeconds) {
  Stopwatch sw(45);
  EXPECT_EQ(sw.total_seconds(), 45);
}

/**
 * @brief Test total_seconds with only minutes
 */
TEST_F(StopwatchTotalSecondsTest, OnlyMinutes) {
  Stopwatch sw(120);
  EXPECT_EQ(sw.total_seconds(), 120);
}

/**
 * @brief Test total_seconds with minutes and seconds
 */
TEST_F(StopwatchTotalSecondsTest, MixedTime) {
  Stopwatch sw(185); // 3 minutes and 5 seconds
  EXPECT_EQ(sw.total_seconds(), 185);
}

/**
 * @brief Test total_seconds after increment
 */
TEST_F(StopwatchTotalSecondsTest, AfterIncrement) {
  Stopwatch sw(60);
  sw.increment(45);
  EXPECT_EQ(sw.total_seconds(), 105);
}

class StopwatchStringTest : public ::testing::Test {};

/**
 * @brief Test string representation of zero time
 */
TEST_F(StopwatchStringTest, ZeroTime) {
  Stopwatch sw(0);
  std::string result = sw.to_string();
  EXPECT_FALSE(result.empty());
  // Should contain "0" for both minutes and seconds
}

/**
 * @brief Test string representation with only seconds
 */
TEST_F(StopwatchStringTest, OnlySeconds) {
  Stopwatch sw(45);
  std::string result = sw.to_string();
  EXPECT_FALSE(result.empty());
}

/**
 * @brief Test string representation with only minutes
 */
TEST_F(StopwatchStringTest, OnlyMinutes) {
  Stopwatch sw(120);
  std::string result = sw.to_string();
  EXPECT_FALSE(result.empty());
}

/**
 * @brief Test string representation with minutes and seconds
 */
TEST_F(StopwatchStringTest, MixedTime) {
  Stopwatch sw(185);
  std::string result = sw.to_string();
  EXPECT_FALSE(result.empty());
}

class StopwatchComparisonTest : public ::testing::Test {};

/**
 * @brief Test equality for same time
 */
TEST_F(StopwatchComparisonTest, EqualitySameTime) {
  Stopwatch sw1(90);
  Stopwatch sw2(90);
  EXPECT_EQ(sw1, sw2);
}

/**
 * @brief Test inequality for different times
 */
TEST_F(StopwatchComparisonTest, InequalityDifferentTime) {
  Stopwatch sw1(90);
  Stopwatch sw2(120);
  EXPECT_NE(sw1, sw2);
}

/**
 * @brief Test less than comparison
 */
TEST_F(StopwatchComparisonTest, LessThan) {
  Stopwatch sw1(60);
  Stopwatch sw2(120);
  EXPECT_LT(sw1, sw2);
  EXPECT_FALSE(sw2 < sw1);
}

/**
 * @brief Test less than or equal comparison
 */
TEST_F(StopwatchComparisonTest, LessThanOrEqual) {
  Stopwatch sw1(60);
  Stopwatch sw2(120);
  Stopwatch sw3(120);
  EXPECT_LE(sw1, sw2);
  EXPECT_LE(sw2, sw3);
  EXPECT_FALSE(sw2 <= sw1);
}

/**
 * @brief Test greater than comparison
 */
TEST_F(StopwatchComparisonTest, GreaterThan) {
  Stopwatch sw1(120);
  Stopwatch sw2(60);
  EXPECT_GT(sw1, sw2);
  EXPECT_FALSE(sw2 > sw1);
}

/**
 * @brief Test greater than or equal comparison
 */
TEST_F(StopwatchComparisonTest, GreaterThanOrEqual) {
  Stopwatch sw1(120);
  Stopwatch sw2(60);
  Stopwatch sw3(120);
  EXPECT_GE(sw1, sw2);
  EXPECT_GE(sw1, sw3);
  EXPECT_FALSE(sw2 >= sw1);
}

/**
 * @brief Test comparison with zero time
 */
TEST_F(StopwatchComparisonTest, ComparisonWithZero) {
  Stopwatch sw_zero(0);
  Stopwatch sw_pos(10);
  EXPECT_LT(sw_zero, sw_pos);
  EXPECT_LE(sw_zero, sw_pos);
  EXPECT_GT(sw_pos, sw_zero);
  EXPECT_GE(sw_pos, sw_zero);
}

/**
 * @brief Test comparison after increment
 */
TEST_F(StopwatchComparisonTest, AfterIncrement) {
  Stopwatch sw1(60);
  Stopwatch sw2(30);
  sw2.increment(30);
  EXPECT_EQ(sw1, sw2);
}

class StopwatchEdgeCasesTest : public ::testing::Test {};

/**
 * @brief Test increment with very large value
 */
TEST_F(StopwatchEdgeCasesTest, LargeIncrement) {
  Stopwatch sw(0);
  sw.increment(7200); // 2 hours
  EXPECT_EQ(sw.minutes, 120);
  EXPECT_EQ(sw.seconds, 0);
}

/**
 * @brief Test multiple overflows
 */
TEST_F(StopwatchEdgeCasesTest, MultipleOverflows) {
  Stopwatch sw(59);
  sw.increment(61); // Should overflow twice
  EXPECT_EQ(sw.minutes, 2);
  EXPECT_EQ(sw.seconds, 0);
}

/**
 * @brief Test that increment maintains consistency with total_seconds
 */
TEST_F(StopwatchEdgeCasesTest, IncrementConsistency) {
  Stopwatch sw(0);
  std::vector<unsigned int> increments = {15, 30, 45, 60, 90, 120};
  unsigned int expected_total = 0;

  for (unsigned int inc : increments) {
    sw.increment(inc);
    expected_total += inc;
    EXPECT_EQ(sw.total_seconds(), expected_total);
  }
}
