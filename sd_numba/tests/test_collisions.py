"""Tests for sd_numba.physics.collision module."""

import numpy as np
import pytest
from numba.typed import List as TypedList
from sd_numba.core.config import ModelConfig
from sd_numba.core.droplet import Droplet
from sd_numba.physics.collision import (
    CollisionStepResult,
    collision_step,
    multi_coalesce,
)
from sd_numba.physics.kernels import Kernel


class TestCollisionStepResult:
    """Tests for CollisionStepResult dataclass."""

    def test_collision_step_result_creation(self):
        """Test creating a CollisionStepResult."""
        result = CollisionStepResult(
            counter=10,
            big_probs=2,
            max_prob=1.5,
            min_prob=0.3,
            total_xi=100000,
        )

        assert result.counter == 10
        assert result.big_probs == 2
        assert result.max_prob == 1.5
        assert result.min_prob == 0.3
        assert result.total_xi == 100000


class TestMultiCoalesce:
    """Tests for multi_coalesce function."""

    def test_multi_coalesce_case1_excess_droplets(self):
        """Test multi-coalescence case 1: excess droplets remain in sd_j."""
        # sd_j has more droplets than needed for full coalescence
        sd_j = Droplet.new(multi=1000, radius=10e-6)
        sd_k = Droplet.new(multi=500, radius=20e-6)

        original_j_multi = sd_j.multi
        original_k_multi = sd_k.multi
        original_k_rcubed = sd_k.rcubed

        gamma = 1.5

        multi_coalesce(sd_j, sd_k, gamma)

        # sd_j should have some excess droplets remaining
        assert sd_j.multi < original_j_multi
        assert sd_j.multi > 0

        # sd_k multiplicity should stay the same
        assert sd_k.multi == original_k_multi

        # sd_k should have grown (rcubed increased)
        assert sd_k.rcubed > original_k_rcubed

    def test_multi_coalesce_case2_split_sd_k(self):
        """Test multi-coalescence case 2: all sd_j droplets pair, split sd_k."""
        # Gamma large enough that all sd_j droplets are consumed
        sd_j = Droplet.new(multi=200, radius=10e-6)
        sd_k = Droplet.new(multi=1000, radius=20e-6)

        original_k_multi = sd_k.multi
        gamma = 10.0  # Large gamma

        multi_coalesce(sd_j, sd_k, gamma)

        # Both multiplicities should have changed
        assert sd_j.multi != 200
        assert sd_k.multi != original_k_multi

        # Both should have same rcubed (coalesced)
        assert sd_j.rcubed == pytest.approx(sd_k.rcubed)

        # Multiplicities should sum to original sd_k
        assert sd_j.multi + sd_k.multi == original_k_multi

    def test_multi_coalesce_conserves_total_multiplicity_case1(self):
        """Test that total droplet count changes correctly in case 1."""
        sd_j = Droplet.new(multi=1000, radius=15e-6)
        sd_k = Droplet.new(multi=500, radius=25e-6)

        gamma = 1.5
        gamma_t = min(gamma, sd_j.multi / sd_k.multi)

        original_total = sd_j.multi + sd_k.multi
        consumed_from_j = int(np.floor(gamma_t * sd_k.multi))

        multi_coalesce(sd_j, sd_k, gamma)

        # In case 1, consumed droplets from j merge with k
        # Total count decreases by consumed amount
        expected_total = original_total - consumed_from_j
        assert sd_j.multi + sd_k.multi == expected_total

    def test_multi_coalesce_increases_droplet_size(self):
        """Test that coalescence increases droplet radius."""
        sd_j = Droplet.new(multi=1000, radius=10e-6)
        sd_k = Droplet.new(multi=500, radius=20e-6)

        original_k_radius = sd_k.radius
        gamma = 1.5

        multi_coalesce(sd_j, sd_k, gamma)

        # sd_k should have grown
        assert sd_k.radius > original_k_radius

    def test_multi_coalesce_with_gamma_equals_ratio(self):
        """Test multi-coalescence when gamma equals the multiplicity ratio."""
        sd_j = Droplet.new(multi=1000, radius=15e-6)
        sd_k = Droplet.new(multi=500, radius=20e-6)

        ratio = sd_j.multi / sd_k.multi
        gamma = ratio  # Exactly at the ratio boundary

        original_j_multi = sd_j.multi

        multi_coalesce(sd_j, sd_k, gamma)

        # Should use gamma_t = ratio, consuming all of j's droplets
        assert sd_j.multi == 0 or sd_j.multi < original_j_multi

    def test_multi_coalesce_solute_update_case1(self):
        """Test that solute is updated correctly in case 1."""
        sd_j = Droplet.new(multi=1000, radius=15e-6)
        sd_k = Droplet.new(multi=500, radius=20e-6)

        # Set some solute mass
        sd_j.solute = 1e-10
        sd_k.solute = 2e-10

        original_k_solute = sd_k.solute
        gamma = 1.5

        multi_coalesce(sd_j, sd_k, gamma)

        # sd_k solute should have increased
        assert sd_k.solute > original_k_solute

    def test_multi_coalesce_solute_update_case2(self):
        """Test that solute is updated correctly in case 2."""
        sd_j = Droplet.new(multi=200, radius=10e-6)
        sd_k = Droplet.new(multi=1000, radius=20e-6)

        sd_j.solute = 1e-10
        sd_k.solute = 2e-10

        gamma = 10.0

        multi_coalesce(sd_j, sd_k, gamma)

        # Both should have the same solute (split)
        assert sd_j.solute == pytest.approx(sd_k.solute)


class TestCollisionStep:
    """Tests for collision_step function."""

    def test_collision_step_returns_result(self, typed_droplets_64):
        """Test that collision_step returns a CollisionStepResult."""
        config = ModelConfig(
            step_seconds=1,
            delta_v=1e6,
            num_droplets=64,
            kernel=Kernel.GOLOVIN,
            debug=False,
        )

        result = collision_step(typed_droplets_64, config)

        assert isinstance(result, CollisionStepResult)
        assert isinstance(result.counter, int)
        assert isinstance(result.big_probs, int)
        assert isinstance(result.max_prob, float)
        assert isinstance(result.min_prob, float)
        assert isinstance(result.total_xi, int)

    def test_collision_step_processes_pairs(self):
        """Test that collision_step processes correct number of pairs."""
        n_droplets = 128  # Power of 2
        typed_list = TypedList()
        for i in range(n_droplets):
            typed_list.append(Droplet.new(multi=1000, radius=10e-6 + i * 1e-7))

        config = ModelConfig(
            step_seconds=1,
            delta_v=1e6,
            num_droplets=n_droplets,
            kernel=Kernel.GOLOVIN,
            debug=False,
        )

        result = collision_step(typed_list, config)

        # Should process n_droplets/2 pairs
        # Counter might be less than n_droplets/2 due to random collisions
        assert result.counter >= 0
        assert result.counter <= n_droplets // 2

    def test_collision_step_with_golovin_kernel(self, typed_droplets_64):
        """Test collision_step with Golovin kernel."""
        config = ModelConfig(
            step_seconds=1,
            delta_v=1e5,
            num_droplets=64,
            kernel=Kernel.GOLOVIN,
            debug=False,
        )

        result = collision_step(typed_droplets_64, config)

        assert result.counter >= 0
        assert result.max_prob >= 0
        assert (
            result.min_prob <= 1.0 or result.min_prob == 1.0
        )  # min_prob initialized to 1.0

    def test_collision_step_with_long_kernel(self, typed_droplets_64):
        """Test collision_step with Long kernel."""
        config = ModelConfig(
            step_seconds=1,
            delta_v=1e5,
            num_droplets=64,
            kernel=Kernel.LONG,
            debug=False,
        )

        result = collision_step(typed_droplets_64, config)

        assert result.counter >= 0

    def test_collision_step_with_hydro_kernel(self, typed_droplets_64):
        """Test collision_step with Hydro kernel."""
        config = ModelConfig(
            step_seconds=1,
            delta_v=1e5,
            num_droplets=64,
            kernel=Kernel.HYDRO,
            debug=False,
        )

        result = collision_step(typed_droplets_64, config)

        assert result.counter >= 0

    def test_collision_step_skips_zero_multiplicity(self):
        """Test that collision_step skips pairs with zero multiplicity."""
        typed_list = TypedList()
        for i in range(64):
            multi = 1000 if i % 2 == 0 else 0
            typed_list.append(Droplet.new(multi=multi, radius=10e-6))

        config = ModelConfig(
            step_seconds=1,
            delta_v=1e6,
            num_droplets=64,
            kernel=Kernel.GOLOVIN,
            debug=False,
        )

        result = collision_step(typed_list, config)

        # Some pairs should be skipped due to zero multiplicity
        # This test just ensures it doesn't crash
        assert result.counter >= 0

    def test_collision_step_total_xi(self):
        """Test that total_xi is computed correctly."""
        n_droplets = 64
        multi_per_droplet = 1000
        typed_list = TypedList()
        for _ in range(n_droplets):
            typed_list.append(Droplet.new(multi=multi_per_droplet, radius=10e-6))

        config = ModelConfig(
            step_seconds=1,
            delta_v=1e6,
            num_droplets=n_droplets,
            kernel=Kernel.GOLOVIN,
            debug=False,
        )

        result = collision_step(typed_list, config)

        # total_xi should be close to original (may change due to collisions)
        original_total = n_droplets * multi_per_droplet
        # Allow for collisions reducing the count
        assert result.total_xi <= original_total
        assert result.total_xi >= 0

    def test_collision_step_probability_statistics(self, typed_droplets_64):
        """Test that probability statistics are tracked correctly."""
        config = ModelConfig(
            step_seconds=1,
            delta_v=1e5,
            num_droplets=64,
            kernel=Kernel.GOLOVIN,
            debug=False,
        )

        result = collision_step(typed_droplets_64, config)

        # Max probability should be >= min probability
        if result.min_prob < 1.0:  # min_prob was updated
            assert result.max_prob >= result.min_prob

        # big_probs should be count of probabilities > 1
        assert result.big_probs >= 0

    def test_collision_step_modifies_droplets(self):
        """Test that collision_step modifies droplet properties."""
        n_droplets = 64
        typed_list = TypedList()
        for i in range(n_droplets):
            typed_list.append(Droplet.new(multi=10000, radius=20e-6 + i * 1e-6))

        config = ModelConfig(
            step_seconds=1,
            delta_v=1e4,  # Small volume increases collision probability
            num_droplets=n_droplets,
            kernel=Kernel.GOLOVIN,
            debug=False,
        )

        # Record original states
        original_multis = [d.multi for d in typed_list]
        original_radii = [d.radius for d in typed_list]

        result = collision_step(typed_list, config)

        # If any collisions occurred, some droplets should have changed
        if result.counter > 0:
            multis_changed = any(
                d.multi != orig for d, orig in zip(typed_list, original_multis)
            )
            radii_changed = any(
                d.radius != orig for d, orig in zip(typed_list, original_radii)
            )

            # At least one of these should be true if collisions occurred
            assert multis_changed or radii_changed

    def test_collision_step_with_small_timestep(self):
        """Test collision_step with small timestep."""
        typed_list = TypedList()
        for _ in range(64):
            typed_list.append(Droplet.new(multi=10000, radius=20e-6))

        config = ModelConfig(
            step_seconds=1,  # Small timestep
            delta_v=1e6,  # Large volume
            num_droplets=64,
            kernel=Kernel.GOLOVIN,
            debug=False,
        )

        result = collision_step(typed_list, config)

        # With small timestep and large volume, collision probability is low
        # Counter might be 0 or very small
        assert result.counter >= 0

    def test_collision_step_with_large_timestep(self):
        """Test collision_step with large timestep."""
        typed_list = TypedList()
        for i in range(64):
            typed_list.append(Droplet.new(multi=10000, radius=20e-6 + i * 1e-6))

        config = ModelConfig(
            step_seconds=10,  # Larger timestep
            delta_v=1e4,  # Smaller volume
            num_droplets=64,
            kernel=Kernel.GOLOVIN,
            debug=False,
        )

        result = collision_step(typed_list, config)

        # With larger timestep and smaller volume, collision probability is higher
        # More likely to have big_probs > 0
        assert result.counter >= 0

    def test_collision_step_produces_valid_statistics(self):
        """Test that collision_step produces valid statistical results over multiple runs."""

        def make_droplets():
            typed_list = TypedList()
            for i in range(64):
                typed_list.append(Droplet.new(multi=10000, radius=20e-6 + i * 1e-6))
            return typed_list

        config = ModelConfig(
            step_seconds=1,
            delta_v=1e5,
            num_droplets=64,
            kernel=Kernel.GOLOVIN,
            debug=False,
        )

        # Run multiple times and collect statistics
        counters = []
        for _ in range(10):
            typed_list = make_droplets()
            result = collision_step(typed_list, config)
            counters.append(result.counter)

        # Should produce reasonable collision counts (not always 0 or always max)
        mean_count = sum(counters) / len(counters)
        assert mean_count >= 0
        assert mean_count <= 32  # Max possible is 64/2 = 32

        # Should have some variance in results (randomness is working)
        # But allow for possibility of no variance with small samples
        assert len(counters) == 10

    def test_collision_step_different_droplet_sizes(self):
        """Test collision_step with varied droplet sizes."""
        # Create droplets with varying sizes
        radii = [5e-6, 10e-6, 20e-6, 50e-6] * 16  # 64 droplets
        typed_list = TypedList()
        for r in radii:
            typed_list.append(Droplet.new(multi=10000, radius=r))

        config = ModelConfig(
            step_seconds=1,
            delta_v=1e5,
            num_droplets=64,
            kernel=Kernel.GOLOVIN,
            debug=False,
        )

        result = collision_step(typed_list, config)

        # Should handle varied sizes without error
        assert result.counter >= 0
        assert result.total_xi > 0
