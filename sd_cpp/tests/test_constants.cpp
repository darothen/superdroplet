/**
 * @file test_constants.cpp
 * @brief Tests for sd_cpp::constants module
 */

#include "core/constants.hpp"
#include <cmath>
#include <gtest/gtest.h>

using namespace sd_cpp::constants;

/**
 * @brief Test that PI constant is correct
 */
TEST(ConstantsTest, PiConstant) { EXPECT_DOUBLE_EQ(PI, M_PI); }

/**
 * @brief Test that water density is correct
 */
TEST(ConstantsTest, RhoWaterConstant) { EXPECT_DOUBLE_EQ(RHO_WATER, 1000.0); }

/**
 * @brief Test that air density is correct
 */
TEST(ConstantsTest, RhoAirConstant) { EXPECT_DOUBLE_EQ(RHO_AIR, 1.0); }

/**
 * @brief Test that THIRD constant is 1/3
 */
TEST(ConstantsTest, ThirdConstant) { EXPECT_DOUBLE_EQ(THIRD, 1.0 / 3.0); }

/**
 * @brief Test that THREE_FOURTH constant is 3/4
 */
TEST(ConstantsTest, ThreeFourthConstant) {
  EXPECT_DOUBLE_EQ(THREE_FOURTH, 3.0 / 4.0);
}

/**
 * @brief Test that FOUR_THIRD constant is 4/3
 */
TEST(ConstantsTest, FourThirdConstant) {
  EXPECT_DOUBLE_EQ(FOUR_THIRD, 4.0 / 3.0);
}

/**
 * @brief Test that MULTI_THRESH constant exists and is positive
 */
TEST(ConstantsTest, MultiThreshConstant) {
  EXPECT_GT(MULTI_THRESH, 0.0);
  EXPECT_DOUBLE_EQ(MULTI_THRESH, 1e4);
}
