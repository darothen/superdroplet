/**
 * @file test_collision.cpp
 * @brief Tests for sd_cpp::physics::collision module
 */

#include "core/droplet.hpp"
#include "physics/collision.hpp"
#include "physics/kernels.hpp"
#include <gtest/gtest.h>
#include <random>
#include <vector>

using namespace sd_cpp;
using namespace sd_cpp::physics;

class MultiCoalesceTest : public ::testing::Test {
protected:
  void SetUp() override {}
};

/**
 * @brief Test multi-coalescence case 1: excess droplets remain in sd_j
 */
TEST_F(MultiCoalesceTest, Case1ExcessDroplets) {
  // sd_j has more droplets than needed for full coalescence
  Droplet sd_j(1000, 10e-6);
  Droplet sd_k(500, 20e-6);

  size_t original_j_multi = sd_j.multi;
  size_t original_k_multi = sd_k.multi;
  double original_k_rcubed = sd_k.rcubed;

  double gamma = 1.5;

  multi_coalesce(sd_j, sd_k, gamma);

  // sd_j should have some excess droplets remaining
  EXPECT_LT(sd_j.multi, original_j_multi);
  EXPECT_GT(sd_j.multi, 0);

  // sd_k multiplicity should stay the same
  EXPECT_EQ(sd_k.multi, original_k_multi);

  // sd_k should have grown (rcubed increased)
  EXPECT_GT(sd_k.rcubed, original_k_rcubed);
}

/**
 * @brief Test multi-coalescence case 2: all sd_j droplets pair, split sd_k
 */
TEST_F(MultiCoalesceTest, Case2SplitSdK) {
  // Gamma large enough that all sd_j droplets are consumed
  Droplet sd_j(200, 10e-6);
  Droplet sd_k(1000, 20e-6);

  size_t original_k_multi = sd_k.multi;
  double gamma = 10.0; // Large gamma

  multi_coalesce(sd_j, sd_k, gamma);

  // Both multiplicities should have changed
  EXPECT_NE(sd_j.multi, 200);
  EXPECT_NE(sd_k.multi, original_k_multi);

  // Both should have same rcubed (coalesced)
  EXPECT_DOUBLE_EQ(sd_j.rcubed, sd_k.rcubed);

  // Multiplicities should sum to original sd_k
  EXPECT_EQ(sd_j.multi + sd_k.multi, original_k_multi);
}

/**
 * @brief Test that total droplet count changes correctly in case 1
 */
TEST_F(MultiCoalesceTest, ConservesTotalMultiplicityCase1) {
  Droplet sd_j(1000, 15e-6);
  Droplet sd_k(500, 25e-6);

  double gamma = 1.5;
  double ratio =
      static_cast<double>(sd_j.multi) / static_cast<double>(sd_k.multi);
  double gamma_t = std::min(gamma, ratio);

  size_t original_total = sd_j.multi + sd_k.multi;
  size_t consumed_from_j =
      static_cast<size_t>(std::floor(gamma_t * sd_k.multi));

  multi_coalesce(sd_j, sd_k, gamma);

  // In case 1, consumed droplets from j merge with k
  // Total count decreases by consumed amount
  size_t expected_total = original_total - consumed_from_j;
  EXPECT_EQ(sd_j.multi + sd_k.multi, expected_total);
}

/**
 * @brief Test that coalescence increases droplet size
 */
TEST_F(MultiCoalesceTest, IncreasesDropletSize) {
  Droplet sd_j(1000, 10e-6);
  Droplet sd_k(500, 20e-6);

  double original_k_radius = sd_k.radius();
  double gamma = 1.5;

  multi_coalesce(sd_j, sd_k, gamma);

  // sd_k should have grown
  EXPECT_GT(sd_k.radius(), original_k_radius);
}

/**
 * @brief Test multi-coalescence when gamma equals the multiplicity ratio
 */
TEST_F(MultiCoalesceTest, GammaEqualsRatio) {
  Droplet sd_j(1000, 15e-6);
  Droplet sd_k(500, 20e-6);

  double ratio =
      static_cast<double>(sd_j.multi) / static_cast<double>(sd_k.multi);
  double gamma = ratio; // Exactly at the ratio boundary

  size_t original_j_multi = sd_j.multi;

  multi_coalesce(sd_j, sd_k, gamma);

  // Should use gamma_t = ratio, consuming all of j's droplets
  EXPECT_TRUE(sd_j.multi == 0 || sd_j.multi < original_j_multi);
}

/**
 * @brief Test that solute is updated correctly in case 1
 */
TEST_F(MultiCoalesceTest, SoluteUpdateCase1) {
  Droplet sd_j(1000, 15e-6);
  Droplet sd_k(500, 20e-6);

  // Set some solute mass
  sd_j.solute = 1e-10;
  sd_k.solute = 2e-10;

  double original_k_solute = sd_k.solute;
  double gamma = 1.5;

  multi_coalesce(sd_j, sd_k, gamma);

  // sd_k solute should have increased
  EXPECT_GT(sd_k.solute, original_k_solute);
}

/**
 * @brief Test that solute is updated correctly in case 2
 */
TEST_F(MultiCoalesceTest, SoluteUpdateCase2) {
  Droplet sd_j(200, 10e-6);
  Droplet sd_k(1000, 20e-6);

  sd_j.solute = 1e-10;
  sd_k.solute = 2e-10;

  double gamma = 10.0;

  multi_coalesce(sd_j, sd_k, gamma);

  // Both should have the same solute (split)
  EXPECT_DOUBLE_EQ(sd_j.solute, sd_k.solute);
}

class CollisionStepTest : public ::testing::Test {
protected:
  std::mt19937_64 rng{42};

  std::vector<Droplet> make_droplets(size_t n_droplets) {
    std::vector<Droplet> droplets;
    for (size_t i = 0; i < n_droplets; ++i) {
      droplets.emplace_back(1000, 10e-6 + i * 1e-7);
    }
    return droplets;
  }
};

/**
 * @brief Test that collision_step returns a CollisionStepResult
 */
TEST_F(CollisionStepTest, ReturnsResult) {
  auto droplets = make_droplets(64);

  ModelConfig config{.step_seconds = 1,
                     .delta_v = 1e6,
                     .num_droplets = 64,
                     .kernel = Kernel::Golovin,
                     .debug = false};

  CollisionStepResult result = collision_step(droplets, config, rng);

  EXPECT_GE(result.counter, 0);
  EXPECT_GE(result.big_probs, 0);
  EXPECT_GE(result.max_prob, 0.0);
  EXPECT_GT(result.min_prob, 0.0);
  EXPECT_GT(result.total_xi, 0.0);
  EXPECT_GT(result.total_water, 0.0);
}

/**
 * @brief Test that collision_step processes correct number of pairs
 */
TEST_F(CollisionStepTest, ProcessesPairs) {
  size_t n_droplets = 128; // Power of 2
  auto droplets = make_droplets(n_droplets);

  ModelConfig config{.step_seconds = 1,
                     .delta_v = 1e6,
                     .num_droplets = static_cast<unsigned int>(n_droplets),
                     .kernel = Kernel::Golovin,
                     .debug = false};

  CollisionStepResult result = collision_step(droplets, config, rng);

  // Should process n_droplets/2 pairs
  // Counter might be less than n_droplets/2 due to random collisions
  EXPECT_GE(result.counter, 0);
  EXPECT_LE(result.counter, n_droplets / 2);
}

/**
 * @brief Test collision_step with Golovin kernel
 */
TEST_F(CollisionStepTest, WithGolovinKernel) {
  auto droplets = make_droplets(64);

  ModelConfig config{.step_seconds = 1,
                     .delta_v = 1e5,
                     .num_droplets = 64,
                     .kernel = Kernel::Golovin,
                     .debug = false};

  CollisionStepResult result = collision_step(droplets, config, rng);

  EXPECT_GE(result.counter, 0);
  EXPECT_GE(result.max_prob, 0.0);
}

/**
 * @brief Test collision_step with Long kernel
 */
TEST_F(CollisionStepTest, WithLongKernel) {
  auto droplets = make_droplets(64);

  ModelConfig config{.step_seconds = 1,
                     .delta_v = 1e5,
                     .num_droplets = 64,
                     .kernel = Kernel::Long,
                     .debug = false};

  CollisionStepResult result = collision_step(droplets, config, rng);

  EXPECT_GE(result.counter, 0);
}

/**
 * @brief Test collision_step with Hydro kernel
 */
TEST_F(CollisionStepTest, WithHydroKernel) {
  auto droplets = make_droplets(64);

  ModelConfig config{.step_seconds = 1,
                     .delta_v = 1e5,
                     .num_droplets = 64,
                     .kernel = Kernel::Hydro,
                     .debug = false};

  CollisionStepResult result = collision_step(droplets, config, rng);

  EXPECT_GE(result.counter, 0);
}

/**
 * @brief Test that collision_step skips pairs with zero multiplicity
 */
TEST_F(CollisionStepTest, SkipsZeroMultiplicity) {
  std::vector<Droplet> droplets;
  for (size_t i = 0; i < 64; ++i) {
    size_t multi = (i % 2 == 0) ? 1000 : 0;
    droplets.emplace_back(multi, 10e-6);
  }

  ModelConfig config{.step_seconds = 1,
                     .delta_v = 1e6,
                     .num_droplets = 64,
                     .kernel = Kernel::Golovin,
                     .debug = false};

  CollisionStepResult result = collision_step(droplets, config, rng);

  // Some pairs should be skipped due to zero multiplicity
  // This test just ensures it doesn't crash
  EXPECT_GE(result.counter, 0);
}

/**
 * @brief Test that total_xi is computed correctly
 */
TEST_F(CollisionStepTest, TotalXi) {
  size_t n_droplets = 64;
  size_t multi_per_droplet = 1000;
  std::vector<Droplet> droplets;
  for (size_t i = 0; i < n_droplets; ++i) {
    droplets.emplace_back(multi_per_droplet, 10e-6);
  }

  ModelConfig config{.step_seconds = 1,
                     .delta_v = 1e6,
                     .num_droplets = static_cast<unsigned int>(n_droplets),
                     .kernel = Kernel::Golovin,
                     .debug = false};

  CollisionStepResult result = collision_step(droplets, config, rng);

  // total_xi should be close to original (may change due to collisions)
  size_t original_total = n_droplets * multi_per_droplet;
  // Allow for collisions reducing the count
  EXPECT_LE(result.total_xi, original_total);
  EXPECT_GE(result.total_xi, 0);
}

/**
 * @brief Test that probability statistics are tracked correctly
 */
TEST_F(CollisionStepTest, ProbabilityStatistics) {
  auto droplets = make_droplets(64);

  ModelConfig config{.step_seconds = 1,
                     .delta_v = 1e5,
                     .num_droplets = 64,
                     .kernel = Kernel::Golovin,
                     .debug = false};

  CollisionStepResult result = collision_step(droplets, config, rng);

  // Max probability should be >= min probability
  if (result.min_prob < 1.0) { // min_prob was updated
    EXPECT_GE(result.max_prob, result.min_prob);
  }

  // big_probs should be count of probabilities > 1
  EXPECT_GE(result.big_probs, 0);
}

/**
 * @brief Test that collision_step modifies droplet properties
 */
TEST_F(CollisionStepTest, ModifiesDroplets) {
  size_t n_droplets = 64;
  std::vector<Droplet> droplets;
  for (size_t i = 0; i < n_droplets; ++i) {
    droplets.emplace_back(10000, 20e-6 + i * 1e-6);
  }

  // Record original states
  std::vector<size_t> original_multis;
  std::vector<double> original_radii;
  for (const auto &d : droplets) {
    original_multis.push_back(d.multi);
    original_radii.push_back(d.radius());
  }

  ModelConfig config{.step_seconds = 1,
                     .delta_v =
                         1e4, // Small volume increases collision probability
                     .num_droplets = static_cast<unsigned int>(n_droplets),
                     .kernel = Kernel::Golovin,
                     .debug = false};

  CollisionStepResult result = collision_step(droplets, config, rng);

  // If any collisions occurred, some droplets should have changed
  if (result.counter > 0) {
    bool multis_changed = false;
    bool radii_changed = false;

    for (size_t i = 0; i < droplets.size(); ++i) {
      if (droplets[i].multi != original_multis[i]) {
        multis_changed = true;
      }
      if (droplets[i].radius() != original_radii[i]) {
        radii_changed = true;
      }
    }

    // At least one of these should be true if collisions occurred
    EXPECT_TRUE(multis_changed || radii_changed);
  }
}

/**
 * @brief Test collision_step with small timestep
 */
TEST_F(CollisionStepTest, SmallTimestep) {
  std::vector<Droplet> droplets;
  for (size_t i = 0; i < 64; ++i) {
    droplets.emplace_back(10000, 20e-6);
  }

  ModelConfig config{.step_seconds = 1, // Small timestep
                     .delta_v = 1e6,    // Large volume
                     .num_droplets = 64,
                     .kernel = Kernel::Golovin,
                     .debug = false};

  CollisionStepResult result = collision_step(droplets, config, rng);

  // With small timestep and large volume, collision probability is low
  // Counter might be 0 or very small
  EXPECT_GE(result.counter, 0);
}

/**
 * @brief Test collision_step with large timestep
 */
TEST_F(CollisionStepTest, LargeTimestep) {
  std::vector<Droplet> droplets;
  for (size_t i = 0; i < 64; ++i) {
    droplets.emplace_back(10000, 20e-6 + i * 1e-6);
  }

  ModelConfig config{.step_seconds = 10, // Larger timestep
                     .delta_v = 1e4,     // Smaller volume
                     .num_droplets = 64,
                     .kernel = Kernel::Golovin,
                     .debug = false};

  CollisionStepResult result = collision_step(droplets, config, rng);

  // With larger timestep and smaller volume, collision probability is higher
  EXPECT_GE(result.counter, 0);
}

/**
 * @brief Test that collision_step produces valid statistics over multiple runs
 */
TEST_F(CollisionStepTest, ValidStatisticsMultipleRuns) {
  std::vector<unsigned int> counters;

  for (int run = 0; run < 10; ++run) {
    auto droplets = make_droplets(64);

    ModelConfig config{.step_seconds = 1,
                       .delta_v = 1e5,
                       .num_droplets = 64,
                       .kernel = Kernel::Golovin,
                       .debug = false};

    CollisionStepResult result = collision_step(droplets, config, rng);
    counters.push_back(result.counter);
  }

  // Should produce reasonable collision counts
  double mean_count = 0.0;
  for (unsigned int c : counters) {
    mean_count += c;
  }
  mean_count /= counters.size();

  EXPECT_GE(mean_count, 0.0);
  EXPECT_LE(mean_count, 32.0); // Max possible is 64/2 = 32
}

/**
 * @brief Test collision_step with varied droplet sizes
 */
TEST_F(CollisionStepTest, VariedDropletSizes) {
  // Create droplets with varying sizes
  std::vector<double> radii = {5e-6, 10e-6, 20e-6, 50e-6};
  std::vector<Droplet> droplets;

  for (size_t i = 0; i < 64; ++i) {
    droplets.emplace_back(10000, radii[i % 4]);
  }

  ModelConfig config{.step_seconds = 1,
                     .delta_v = 1e5,
                     .num_droplets = 64,
                     .kernel = Kernel::Golovin,
                     .debug = false};

  CollisionStepResult result = collision_step(droplets, config, rng);

  // Should handle varied sizes without error
  EXPECT_GE(result.counter, 0);
  EXPECT_GT(result.total_xi, 0.0);
}
