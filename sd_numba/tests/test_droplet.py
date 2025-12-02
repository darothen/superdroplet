"""Tests for sd_python.core.droplet module."""

import math

import pytest
from sd_python.core.constants import FOUR_THIRD, PI, RHO_WATER
from sd_python.core.droplet import (
    Droplet,
    compute_terminal_velocity,
    compute_total_water,
)


class TestDropletCreation:
    """Tests for Droplet creation."""

    def test_new_droplet_basic_properties(self):
        """Test creating a droplet with basic properties."""
        multi = 1000
        radius = 10e-6  # 10 microns

        droplet = Droplet.new(multi=multi, radius=radius)

        assert droplet.multi == multi
        assert droplet.radius == radius
        assert droplet.solute == 0.0
        assert droplet.density == RHO_WATER

    def test_new_droplet_computed_properties(self):
        """Test that computed properties are calculated correctly."""
        radius = 10e-6  # 10 microns
        droplet = Droplet.new(multi=1, radius=radius)

        expected_rcubed = radius**3
        expected_volume = FOUR_THIRD * PI * expected_rcubed
        expected_mass = expected_volume * RHO_WATER

        assert droplet.rcubed == pytest.approx(expected_rcubed)
        assert droplet.volume == pytest.approx(expected_volume)
        assert droplet.mass == pytest.approx(expected_mass)

    def test_new_droplet_terminal_velocity(self):
        """Test that terminal velocity is computed."""
        droplet = Droplet.new(multi=1, radius=10e-6)
        assert droplet.terminal_velocity > 0


class TestTerminalVelocity:
    """Tests for terminal velocity computation."""

    def test_compute_terminal_velocity_small_droplet(self):
        """Test terminal velocity for small droplets (d <= 134.43 μm)."""
        radius = 50e-6  # 50 microns, diameter = 100 μm
        volume = FOUR_THIRD * PI * radius**3
        mass = volume * RHO_WATER

        tv = compute_terminal_velocity(radius, mass)

        assert tv > 0
        assert tv < 10.0  # Should be reasonable for small droplets

    def test_compute_terminal_velocity_medium_droplet(self):
        """Test terminal velocity for medium droplets (134.43 < d <= 1511.64 μm)."""
        radius = 500e-6  # 500 microns, diameter = 1000 μm
        volume = FOUR_THIRD * PI * radius**3
        mass = volume * RHO_WATER

        tv = compute_terminal_velocity(radius, mass)

        assert tv > 0
        assert tv < 10.0  # Should be reasonable for medium droplets

    def test_compute_terminal_velocity_large_droplet(self):
        """Test terminal velocity for large droplets (1511.64 < d <= 3477.84 μm)."""
        radius = 1500e-6  # 1500 microns, diameter = 3000 μm
        volume = FOUR_THIRD * PI * radius**3
        mass = volume * RHO_WATER

        tv = compute_terminal_velocity(radius, mass)

        assert tv > 0
        assert tv < 10.0  # Should be reasonable for large droplets

    def test_compute_terminal_velocity_very_large_droplet(self):
        """Test terminal velocity for very large droplets (d > 3477.84 μm)."""
        radius = 2000e-6  # 2000 microns, diameter = 4000 μm
        volume = FOUR_THIRD * PI * radius**3
        mass = volume * RHO_WATER

        tv = compute_terminal_velocity(radius, mass)

        assert tv > 0
        assert tv < 10.0  # Should be reasonable for very large droplets

    def test_terminal_velocity_increases_with_size(self):
        """Test that terminal velocity increases with droplet size."""
        radii = [10e-6, 50e-6, 100e-6, 500e-6]
        terminal_velocities = []

        for radius in radii:
            volume = FOUR_THIRD * PI * radius**3
            mass = volume * RHO_WATER
            tv = compute_terminal_velocity(radius, mass)
            terminal_velocities.append(tv)

        # Check that terminal velocities generally increase
        for i in range(len(terminal_velocities) - 1):
            assert terminal_velocities[i + 1] >= terminal_velocities[i]


class TestComputeTotalWater:
    """Tests for total water computation."""

    def test_compute_total_water_empty_list(self):
        """Test total water for empty list of droplets."""
        total = compute_total_water([])
        assert total == 0.0

    def test_compute_total_water_single_droplet(self):
        """Test total water for a single droplet."""
        droplet = Droplet.new(multi=1000, radius=10e-6)
        total = compute_total_water([droplet])
        expected = droplet.mass * droplet.multi
        assert total == pytest.approx(expected)

    def test_compute_total_water_multiple_droplets(self, droplet_collection):
        """Test total water for multiple droplets."""
        total = compute_total_water(droplet_collection)

        expected = sum(d.mass * d.multi for d in droplet_collection)
        assert total == pytest.approx(expected)

    def test_compute_total_water_high_multiplicity(self):
        """Test total water with high multiplicity droplets."""
        droplets = [
            Droplet.new(multi=10000, radius=5e-6),
            Droplet.new(multi=5000, radius=20e-6),
            Droplet.new(multi=1000, radius=100e-6),
        ]
        total = compute_total_water(droplets)
        expected = sum(d.mass * d.multi for d in droplets)
        assert total == pytest.approx(expected)
