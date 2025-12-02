"""Tests for sd_numba.physics.kernels module."""

import pytest
from sd_numba.core.constants import FOUR_THIRD, PI
from sd_numba.core.droplet import Droplet
from sd_numba.physics.kernels import (
    GOLOVIN_CONSTANT,
    Kernel,
    calc_hydro_kernel,
    compute_kernel,
    golovin_kernel,
    hydro_kernel,
    long_kernel,
)


class TestKernelEnum:
    """Tests for Kernel enumeration."""

    def test_kernel_enum_values(self):
        """Test that kernel enum has correct values."""
        assert Kernel.GOLOVIN == 1
        assert Kernel.HYDRO == 2
        assert Kernel.LONG == 3

    def test_compute_kernel_golovin(self, small_droplet, medium_droplet):
        """Test compute_kernel with GOLOVIN kernel."""
        result = compute_kernel(int(Kernel.GOLOVIN), small_droplet, medium_droplet)
        expected = golovin_kernel(small_droplet, medium_droplet)
        assert result == pytest.approx(expected)

    def test_compute_kernel_hydro(self, small_droplet, medium_droplet):
        """Test compute_kernel with HYDRO kernel."""
        result = compute_kernel(int(Kernel.HYDRO), small_droplet, medium_droplet)
        expected = hydro_kernel(small_droplet, medium_droplet)
        assert result == pytest.approx(expected)

    def test_compute_kernel_long(self, small_droplet, medium_droplet):
        """Test compute_kernel with LONG kernel."""
        result = compute_kernel(int(Kernel.LONG), small_droplet, medium_droplet)
        expected = long_kernel(small_droplet, medium_droplet)
        assert result == pytest.approx(expected)

    def test_compute_kernel_invalid(self, small_droplet, medium_droplet):
        """Test that invalid kernel returns 0.0."""
        # In nopython mode, we return 0.0 for invalid kernels (can't raise)
        result = compute_kernel(999, small_droplet, medium_droplet)
        assert result == 0.0


class TestGolovinKernel:
    """Tests for Golovin collision kernel."""

    def test_golovin_constant(self):
        """Test that Golovin constant is correct."""
        expected = 1500.0 * FOUR_THIRD * PI
        assert GOLOVIN_CONSTANT == pytest.approx(expected)

    def test_golovin_kernel_basic(self, small_droplet, medium_droplet):
        """Test Golovin kernel basic computation."""
        result = golovin_kernel(small_droplet, medium_droplet)
        expected = GOLOVIN_CONSTANT * (small_droplet.rcubed + medium_droplet.rcubed)
        assert result == pytest.approx(expected)

    def test_golovin_kernel_symmetric(self, small_droplet, large_droplet):
        """Test that Golovin kernel is symmetric."""
        result1 = golovin_kernel(small_droplet, large_droplet)
        result2 = golovin_kernel(large_droplet, small_droplet)
        assert result1 == pytest.approx(result2)

    def test_golovin_kernel_positive(self, small_droplet, medium_droplet):
        """Test that Golovin kernel returns positive values."""
        result = golovin_kernel(small_droplet, medium_droplet)
        assert result > 0

    def test_golovin_kernel_increases_with_size(self):
        """Test that Golovin kernel increases with droplet size."""
        d1 = Droplet.new(multi=1, radius=10e-6)
        d2 = Droplet.new(multi=1, radius=50e-6)
        d3 = Droplet.new(multi=1, radius=100e-6)

        k12 = golovin_kernel(d1, d2)
        k23 = golovin_kernel(d2, d3)

        assert k23 > k12


class TestCalcHydroKernel:
    """Tests for hydrodynamic kernel calculation helper."""

    def test_calc_hydro_kernel_basic(self):
        """Test basic hydrodynamic kernel calculation."""
        e_coal = 1.0
        e_coll = 1.0
        r_sum = 100e-6  # 100 microns
        tv_diff = 1.0  # 1 m/s

        result = calc_hydro_kernel(e_coal, e_coll, r_sum, tv_diff)
        expected = e_coal * e_coll * PI * r_sum * r_sum * abs(tv_diff)
        assert result == pytest.approx(expected)

    def test_calc_hydro_kernel_negative_velocity_diff(self):
        """Test that negative velocity differences are handled correctly."""
        e_coal = 1.0
        e_coll = 1.0
        r_sum = 100e-6
        tv_diff = -1.0

        result = calc_hydro_kernel(e_coal, e_coll, r_sum, tv_diff)
        assert result > 0  # Should be positive due to abs()

    def test_calc_hydro_kernel_zero_velocity_diff(self):
        """Test hydrodynamic kernel with zero velocity difference."""
        e_coal = 1.0
        e_coll = 1.0
        r_sum = 100e-6
        tv_diff = 0.0

        result = calc_hydro_kernel(e_coal, e_coll, r_sum, tv_diff)
        assert result == pytest.approx(0.0)

    def test_calc_hydro_kernel_efficiency_factors(self):
        """Test hydrodynamic kernel with different efficiency factors."""
        r_sum = 100e-6
        tv_diff = 1.0

        result_full = calc_hydro_kernel(1.0, 1.0, r_sum, tv_diff)
        result_half = calc_hydro_kernel(0.5, 1.0, r_sum, tv_diff)

        assert result_half == pytest.approx(result_full * 0.5)


class TestHydroKernel:
    """Tests for hydrodynamic collision kernel."""

    def test_hydro_kernel_basic(self, small_droplet, large_droplet):
        """Test basic hydrodynamic kernel computation."""
        result = hydro_kernel(small_droplet, large_droplet)
        assert result > 0

    def test_hydro_kernel_symmetric(self, small_droplet, large_droplet):
        """Test that hydro kernel is symmetric in droplet order."""
        result1 = hydro_kernel(small_droplet, large_droplet)
        result2 = hydro_kernel(large_droplet, small_droplet)
        # Note: Due to terminal velocity differences, this may not be exactly symmetric
        # but the absolute value should be the same due to abs(tv_diff)
        assert result1 == pytest.approx(result2)

    def test_hydro_kernel_same_size_droplets(self):
        """Test hydro kernel for droplets of the same size."""
        d1 = Droplet.new(multi=1, radius=50e-6)
        d2 = Droplet.new(multi=1, radius=50e-6)

        result = hydro_kernel(d1, d2)
        # Same size droplets should have same terminal velocity, so kernel ~0
        assert result == pytest.approx(0.0, abs=1e-10)

    def test_hydro_kernel_increases_with_size_difference(self):
        """Test that hydro kernel increases with size difference."""
        d_small = Droplet.new(multi=1, radius=10e-6)
        d_medium = Droplet.new(multi=1, radius=50e-6)
        d_large = Droplet.new(multi=1, radius=200e-6)

        k_small_medium = hydro_kernel(d_small, d_medium)
        k_small_large = hydro_kernel(d_small, d_large)

        assert k_small_large > k_small_medium


class TestLongKernel:
    """Tests for Long collision kernel."""

    def test_long_kernel_basic(self, small_droplet, large_droplet):
        """Test basic Long kernel computation."""
        result = long_kernel(small_droplet, large_droplet)
        assert result > 0

    def test_long_kernel_symmetric(self, small_droplet, large_droplet):
        """Test that Long kernel is symmetric in droplet order."""
        result1 = long_kernel(small_droplet, large_droplet)
        result2 = long_kernel(large_droplet, small_droplet)
        assert result1 == pytest.approx(result2)

    def test_long_kernel_same_size_droplets(self):
        """Test Long kernel for droplets of the same size."""
        d1 = Droplet.new(multi=1, radius=50e-6)
        d2 = Droplet.new(multi=1, radius=50e-6)

        result = long_kernel(d1, d2)
        # Same size droplets should have same terminal velocity, so kernel ~0
        assert result == pytest.approx(0.0, abs=1e-10)

    def test_long_kernel_collection_efficiency_large_droplet(self):
        """Test that collection efficiency is 1.0 for large droplets (>= 50 μm)."""
        d_small = Droplet.new(multi=1, radius=10e-6)
        d_large = Droplet.new(multi=1, radius=60e-6)  # > 50 μm

        result = long_kernel(d_small, d_large)

        # For large droplets, e_coll should be 1.0
        assert result > 0

        # Verify it's less than or equal to hydro kernel (as expected)
        hydro_result = hydro_kernel(d_small, d_large)
        assert result <= hydro_result

    def test_long_kernel_collection_efficiency_small_droplets(self):
        """Test collection efficiency calculation for small droplets (< 50 μm)."""
        d1 = Droplet.new(multi=1, radius=10e-6)
        d2 = Droplet.new(multi=1, radius=30e-6)

        result = long_kernel(d1, d2)

        # For small droplets, e_coll is calculated
        # The result should be positive but less than hydro kernel
        hydro_result = hydro_kernel(d1, d2)
        assert 0 <= result <= hydro_result

    def test_long_kernel_collection_efficiency_bounds(self):
        """Test that collection efficiency is bounded between 0 and 1."""
        # Test with various droplet sizes
        radii = [5e-6, 10e-6, 20e-6, 30e-6, 40e-6]

        for r1 in radii:
            for r2 in radii:
                if r1 == r2:
                    continue
                d1 = Droplet.new(multi=1, radius=r1)
                d2 = Droplet.new(multi=1, radius=r2)

                result = long_kernel(d1, d2)
                # Result should be non-negative
                assert result >= 0

    def test_long_kernel_vs_hydro_kernel(self):
        """Test that Long kernel is generally less than or equal to hydro kernel."""
        radii = [10e-6, 20e-6, 30e-6, 40e-6, 60e-6, 100e-6]

        for r1 in radii:
            for r2 in radii:
                if r1 == r2:
                    continue
                d1 = Droplet.new(multi=1, radius=r1)
                d2 = Droplet.new(multi=1, radius=r2)

                long_result = long_kernel(d1, d2)
                hydro_result = hydro_kernel(d1, d2)

                # Long kernel should be <= hydro kernel due to collection efficiency
                assert (
                    long_result <= hydro_result + 1e-15
                )  # Small tolerance for numerical errors

    def test_long_kernel_r_small_r_large_assignment(self):
        """Test that r_small and r_large are correctly assigned (bug fix verification)."""
        # This test specifically verifies the bug fix where r_small and r_large
        # were incorrectly swapped. The fix ensures:
        # - r_small = min(r_j, r_k)
        # - r_large = max(r_j, r_k)

        d_small = Droplet.new(multi=1, radius=10e-6)  # 10 μm
        d_large = Droplet.new(multi=1, radius=100e-6)  # 100 μm

        # Both orderings should give the same result
        result1 = long_kernel(d_small, d_large)
        result2 = long_kernel(d_large, d_small)

        # Results should be identical (symmetry)
        assert result1 == pytest.approx(result2)

        # For large droplets (>= 50 μm), e_coll should be 1.0
        # So Long kernel should equal hydro kernel
        hydro_result = hydro_kernel(d_small, d_large)
        assert result1 == pytest.approx(hydro_result)

    def test_long_kernel_small_collector_large_collected(self):
        """Test Long kernel with specific small and large droplets."""
        # r_small < 50 μm, r_large >= 50 μm case
        d1 = Droplet.new(multi=1, radius=20e-6)  # 20 μm (small)
        d2 = Droplet.new(multi=1, radius=60e-6)  # 60 μm (large)

        result = long_kernel(d1, d2)

        # Since r_large >= 50, e_coll = 1.0
        hydro_result = hydro_kernel(d1, d2)
        assert result == pytest.approx(hydro_result)
