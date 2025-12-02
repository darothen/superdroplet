"""Tests for sd_python.physics.collision module."""

import pytest
from sd_python.core.config import ModelConfig
from sd_python.core.droplet import Droplet
from sd_python.physics.collision import (
    CollisionStepResult,
    collision_step,
    multi_coalesce,
)
from sd_python.physics.kernels import Kernel


class TestCollisionStepResult:
    """Tests for CollisionStepResult dataclass."""

    def test_collision_step_result_creation(self):
        """Test creating a CollisionStepResult."""
        result = CollisionStepResult(
            counter=10,
            big_probs=2,
            max_prob=1.5,
            min_prob=0.3,
            total_xi=100000.0,
        )

        assert result.counter == 10
        assert result.big_probs == 2
        assert result.max_prob == 1.5
        assert result.min_prob == 0.3
        assert result.total_xi == 100000.0


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
        consumed_from_j = int(gamma_t * sd_k.multi)

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
        original_k_multi = sd_k.multi

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

    def test_collision_step_returns_result(self):
        """Test that collision_step returns a CollisionStepResult."""
        droplets = [
            Droplet.new(multi=1000, radius=10e-6 + i * 1e-6) for i in range(100)
        ]
        config = ModelConfig(
            step_seconds=1,
            delta_v=1e6,
            num_droplets=100,
            kernel=Kernel.GOLOVIN,
            debug=False,
        )

        result = collision_step(droplets, config)

        assert isinstance(result, CollisionStepResult)
        assert isinstance(result.counter, int)
        assert isinstance(result.big_probs, int)
        assert isinstance(result.max_prob, float)
        assert isinstance(result.min_prob, float)
        assert isinstance(result.total_xi, int)

    def test_collision_step_processes_pairs(self):
        """Test that collision_step processes correct number of pairs."""
        n_droplets = 128  # Power of 2
        droplets = [
            Droplet.new(multi=1000, radius=10e-6 + i * 1e-7) for i in range(n_droplets)
        ]
        config = ModelConfig(
            step_seconds=1,
            delta_v=1e6,
            num_droplets=n_droplets,
            kernel=Kernel.GOLOVIN,
            debug=False,
        )

        result = collision_step(droplets, config)

        # Should process n_droplets/2 pairs
        # Counter might be less than n_droplets/2 due to random collisions
        assert result.counter >= 0
        assert result.counter <= n_droplets // 2

    def test_collision_step_with_golovin_kernel(self):
        """Test collision_step with Golovin kernel."""
        droplets = [
            Droplet.new(multi=10000, radius=20e-6 + i * 1e-6) for i in range(64)
        ]
        config = ModelConfig(
            step_seconds=1,
            delta_v=1e5,
            num_droplets=64,
            kernel=Kernel.GOLOVIN,
            debug=False,
        )

        result = collision_step(droplets, config)

        assert result.counter >= 0
        assert result.max_prob >= 0
        assert (
            result.min_prob <= 1.0 or result.min_prob == 1.0
        )  # min_prob initialized to 1.0

    def test_collision_step_with_long_kernel(self):
        """Test collision_step with Long kernel."""
        droplets = [
            Droplet.new(multi=10000, radius=20e-6 + i * 1e-6) for i in range(64)
        ]
        config = ModelConfig(
            step_seconds=1,
            delta_v=1e5,
            num_droplets=64,
            kernel=Kernel.LONG,
            debug=False,
        )

        result = collision_step(droplets, config)

        assert result.counter >= 0

    def test_collision_step_with_hydro_kernel(self):
        """Test collision_step with Hydro kernel."""
        droplets = [
            Droplet.new(multi=10000, radius=20e-6 + i * 1e-6) for i in range(64)
        ]
        config = ModelConfig(
            step_seconds=1,
            delta_v=1e5,
            num_droplets=64,
            kernel=Kernel.HYDRO,
            debug=False,
        )

        result = collision_step(droplets, config)

        assert result.counter >= 0

    def test_collision_step_skips_zero_multiplicity(self):
        """Test that collision_step skips pairs with zero multiplicity."""
        droplets = [
            Droplet.new(multi=1000 if i % 2 == 0 else 0, radius=10e-6)
            for i in range(64)
        ]
        config = ModelConfig(
            step_seconds=1,
            delta_v=1e6,
            num_droplets=64,
            kernel=Kernel.GOLOVIN,
            debug=False,
        )

        result = collision_step(droplets, config)

        # Some pairs should be skipped due to zero multiplicity
        # This test just ensures it doesn't crash
        assert result.counter >= 0

    def test_collision_step_total_xi(self):
        """Test that total_xi is computed correctly."""
        n_droplets = 64
        multi_per_droplet = 1000
        droplets = [
            Droplet.new(multi=multi_per_droplet, radius=10e-6)
            for _ in range(n_droplets)
        ]
        config = ModelConfig(
            step_seconds=1,
            delta_v=1e6,
            num_droplets=n_droplets,
            kernel=Kernel.GOLOVIN,
            debug=False,
        )

        result = collision_step(droplets, config)

        # total_xi should be close to original (may change due to collisions)
        original_total = n_droplets * multi_per_droplet
        # Allow for collisions reducing the count
        assert result.total_xi <= original_total
        assert result.total_xi >= 0

    def test_collision_step_probability_statistics(self):
        """Test that probability statistics are tracked correctly."""
        droplets = [
            Droplet.new(multi=10000, radius=20e-6 + i * 1e-6) for i in range(64)
        ]
        config = ModelConfig(
            step_seconds=1,
            delta_v=1e5,
            num_droplets=64,
            kernel=Kernel.GOLOVIN,
            debug=False,
        )

        result = collision_step(droplets, config)

        # Max probability should be >= min probability
        if result.min_prob < 1.0:  # min_prob was updated
            assert result.max_prob >= result.min_prob

        # big_probs should be count of probabilities > 1
        assert result.big_probs >= 0

    def test_collision_step_modifies_droplets(self):
        """Test that collision_step modifies droplet properties."""
        n_droplets = 64
        droplets = [
            Droplet.new(multi=10000, radius=20e-6 + i * 1e-6) for i in range(n_droplets)
        ]
        config = ModelConfig(
            step_seconds=1,
            delta_v=1e4,  # Small volume increases collision probability
            num_droplets=n_droplets,
            kernel=Kernel.GOLOVIN,
            debug=False,
        )

        # Record original states
        original_multis = [d.multi for d in droplets]
        original_radii = [d.radius for d in droplets]

        result = collision_step(droplets, config)

        # If any collisions occurred, some droplets should have changed
        if result.counter > 0:
            multis_changed = any(
                d.multi != orig for d, orig in zip(droplets, original_multis)
            )
            radii_changed = any(
                d.radius != orig for d, orig in zip(droplets, original_radii)
            )

            # At least one of these should be true if collisions occurred
            assert multis_changed or radii_changed

    def test_collision_step_with_small_timestep(self):
        """Test collision_step with small timestep."""
        droplets = [Droplet.new(multi=10000, radius=20e-6) for _ in range(64)]
        config = ModelConfig(
            step_seconds=1,  # Small timestep
            delta_v=1e6,  # Large volume
            num_droplets=64,
            kernel=Kernel.GOLOVIN,
            debug=False,
        )

        result = collision_step(droplets, config)

        # With small timestep and large volume, collision probability is low
        # Counter might be 0 or very small
        assert result.counter >= 0

    def test_collision_step_with_large_timestep(self):
        """Test collision_step with large timestep."""
        droplets = [
            Droplet.new(multi=10000, radius=20e-6 + i * 1e-6) for i in range(64)
        ]
        config = ModelConfig(
            step_seconds=10,  # Larger timestep
            delta_v=1e4,  # Smaller volume
            num_droplets=64,
            kernel=Kernel.GOLOVIN,
            debug=False,
        )

        result = collision_step(droplets, config)

        # With larger timestep and smaller volume, collision probability is higher
        # More likely to have big_probs > 0
        assert result.counter >= 0

    def test_collision_step_reproducibility_with_seed(self):
        """Test that collision_step gives consistent results with same random seed."""
        import random

        droplets1 = [
            Droplet.new(multi=10000, radius=20e-6 + i * 1e-6) for i in range(64)
        ]
        droplets2 = [
            Droplet.new(multi=10000, radius=20e-6 + i * 1e-6) for i in range(64)
        ]

        config = ModelConfig(
            step_seconds=1,
            delta_v=1e5,
            num_droplets=64,
            kernel=Kernel.GOLOVIN,
            debug=False,
        )

        random.seed(42)
        result1 = collision_step(droplets1, config)

        random.seed(42)
        result2 = collision_step(droplets2, config)

        # Results should be identical with same seed
        assert result1.counter == result2.counter
        assert result1.big_probs == result2.big_probs
        assert result1.max_prob == pytest.approx(result2.max_prob)
        assert result1.min_prob == pytest.approx(result2.min_prob)

    def test_collision_step_different_droplet_sizes(self):
        """Test collision_step with varied droplet sizes."""
        # Create droplets with varying sizes
        radii = [5e-6, 10e-6, 20e-6, 50e-6] * 16  # 64 droplets
        droplets = [Droplet.new(multi=10000, radius=r) for r in radii]

        config = ModelConfig(
            step_seconds=1,
            delta_v=1e5,
            num_droplets=64,
            kernel=Kernel.GOLOVIN,
            debug=False,
        )

        result = collision_step(droplets, config)

        # Should handle varied sizes without error
        assert result.counter >= 0
        assert result.total_xi > 0
