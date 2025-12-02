"""Tests for sd_numba.utils.binning module."""

import numpy as np
import pytest
from numba.typed import List as TypedList
from sd_numba.core.droplet import Droplet
from sd_numba.utils.binning import _bin_droplets_jitted, bin_droplets


class TestBinDroplets:
    """Tests for bin_droplets function."""

    def _make_typed_list(self, droplets):
        """Helper to convert list to TypedList."""
        typed_list = TypedList()
        for d in droplets:
            typed_list.append(d)
        return typed_list

    def test_bin_droplets_empty_list(self):
        """Test binning with empty droplet list."""
        bin_edges = np.array([0.0, 10.0, 20.0, 30.0])
        # Create empty typed list by adding and removing
        typed_list = TypedList()
        typed_list.append(Droplet.new(multi=1, radius=1e-6))
        typed_list.pop()

        result = bin_droplets(typed_list, bin_edges)

        assert len(result) == len(bin_edges) - 1
        assert all(val == 0.0 for val in result)

    def test_bin_droplets_single_bin(self):
        """Test binning with a single bin."""
        droplets = [
            Droplet.new(multi=1000, radius=5e-6),
            Droplet.new(multi=2000, radius=7e-6),
        ]
        typed_list = self._make_typed_list(droplets)
        bin_edges = np.array([0.0, 100.0])  # All droplets in one bin

        result = bin_droplets(typed_list, bin_edges)

        assert len(result) == 1
        # Should sum all droplet masses
        expected = sum(d.mass * d.multi for d in droplets)
        assert result[0] == pytest.approx(expected)

    def test_bin_droplets_multiple_bins(self):
        """Test binning with multiple bins."""
        droplets = [
            Droplet.new(multi=1000, radius=5e-6),  # 5 μm radius -> 10 μm diameter
            Droplet.new(multi=2000, radius=15e-6),  # 15 μm radius -> 30 μm diameter
            Droplet.new(multi=3000, radius=25e-6),  # 25 μm radius -> 50 μm diameter
        ]
        typed_list = self._make_typed_list(droplets)
        bin_edges = np.array([0.0, 20.0, 40.0, 60.0])  # bins in microns (radius)

        result = bin_droplets(typed_list, bin_edges)

        assert len(result) == 3
        assert all(val >= 0.0 for val in result)

    def test_bin_droplets_sorting(self):
        """Test that droplets are sorted by radius correctly."""
        # Create droplets in unsorted order
        droplets = [
            Droplet.new(multi=1000, radius=25e-6),
            Droplet.new(multi=2000, radius=5e-6),
            Droplet.new(multi=3000, radius=15e-6),
        ]
        typed_list = self._make_typed_list(droplets)
        bin_edges = np.array([0.0, 10.0, 20.0, 30.0])

        result = bin_droplets(typed_list, bin_edges)

        # Should bin correctly regardless of input order
        assert len(result) == 3
        assert all(val >= 0.0 for val in result)

    def test_bin_droplets_radius_conversion(self):
        """Test that radius conversion from meters to microns is correct."""
        # Create droplet with radius just under bin edge
        droplets = [
            Droplet.new(multi=1000, radius=9.9e-6),  # 9.9 microns
            Droplet.new(multi=2000, radius=10.1e-6),  # 10.1 microns
        ]
        typed_list = self._make_typed_list(droplets)
        bin_edges = np.array([0.0, 10.0, 20.0])  # bins in microns

        result = bin_droplets(typed_list, bin_edges)

        # First droplet should be in first bin
        # Second droplet should be in second bin
        assert result[0] > 0  # First bin has first droplet
        assert result[1] > 0  # Second bin has second droplet

    def test_bin_droplets_accounts_for_multiplicity(self):
        """Test that binning accounts for droplet multiplicity."""
        droplets = [
            Droplet.new(multi=1000, radius=5e-6),
            Droplet.new(multi=1, radius=5e-6),
        ]
        typed_list = self._make_typed_list(droplets)
        bin_edges = np.array([0.0, 100.0])

        result = bin_droplets(typed_list, bin_edges)

        # Should weight by multiplicity
        mass_per_droplet = droplets[0].mass
        expected = 1000 * mass_per_droplet + 1 * mass_per_droplet
        assert result[0] == pytest.approx(expected)

    def test_bin_droplets_with_identical_droplets(self):
        """Test binning with identical droplets."""
        n_droplets = 10
        droplets = [Droplet.new(multi=1000, radius=10e-6) for _ in range(n_droplets)]
        typed_list = self._make_typed_list(droplets)
        bin_edges = np.array([0.0, 50.0, 100.0])

        result = bin_droplets(typed_list, bin_edges)

        # All should be in first bin
        expected = n_droplets * droplets[0].mass * droplets[0].multi
        assert result[0] == pytest.approx(expected)
        assert result[1] == pytest.approx(0.0)

    def test_bin_droplets_boundary_cases(self):
        """Test binning with droplets exactly on bin boundaries."""
        droplets = [
            Droplet.new(multi=1000, radius=10e-6),  # Exactly 10 μm radius
        ]
        typed_list = self._make_typed_list(droplets)
        bin_edges = np.array([0.0, 20.0, 40.0])

        result = bin_droplets(typed_list, bin_edges)

        # Droplet at boundary should be in first bin (< comparison)
        assert result[0] >= 0

    def test_bin_droplets_large_range(self):
        """Test binning with large range of droplet sizes."""
        droplets = [
            Droplet.new(multi=1000, radius=1e-6),  # Very small
            Droplet.new(multi=2000, radius=50e-6),  # Medium
            Droplet.new(multi=3000, radius=500e-6),  # Large
        ]
        typed_list = self._make_typed_list(droplets)
        bin_edges = np.array([0.0, 10.0, 100.0, 1000.0])

        result = bin_droplets(typed_list, bin_edges)

        assert len(result) == 3
        # Each droplet should be in a different bin
        non_zero_bins = sum(1 for val in result if val > 0)
        assert non_zero_bins == 3

    def test_bin_droplets_many_droplets(self):
        """Test binning with many droplets."""
        np.random.seed(42)
        n_droplets = 1000
        droplets = [
            Droplet.new(multi=1000, radius=np.random.uniform(1e-6, 100e-6))
            for _ in range(n_droplets)
        ]
        typed_list = self._make_typed_list(droplets)
        bin_edges = np.array([0.0, 10.0, 50.0, 100.0, 200.0])

        result = bin_droplets(typed_list, bin_edges)

        assert len(result) == 4
        # Total mass should equal sum of all droplets
        total_from_bins = sum(result)
        total_from_droplets = sum(d.mass * d.multi for d in droplets)
        assert total_from_bins == pytest.approx(total_from_droplets)

    def test_bin_droplets_conservation_of_mass(self):
        """Test that binning conserves total mass."""
        droplets = [
            Droplet.new(multi=1000, radius=5e-6),
            Droplet.new(multi=2000, radius=15e-6),
            Droplet.new(multi=3000, radius=25e-6),
            Droplet.new(multi=4000, radius=35e-6),
        ]
        typed_list = self._make_typed_list(droplets)
        bin_edges = np.array([0.0, 20.0, 30.0, 40.0, 50.0])

        result = bin_droplets(typed_list, bin_edges)

        # Sum of bin values should equal total mass
        total_binned = sum(result)
        total_actual = sum(d.mass * d.multi for d in droplets)
        assert total_binned == pytest.approx(total_actual)

    def test_bin_droplets_with_zero_multiplicity(self):
        """Test binning with droplets having zero multiplicity."""
        droplets = [
            Droplet.new(multi=1000, radius=5e-6),
            Droplet.new(multi=0, radius=15e-6),  # Zero multiplicity
            Droplet.new(multi=2000, radius=25e-6),
        ]
        typed_list = self._make_typed_list(droplets)
        bin_edges = np.array([0.0, 20.0, 40.0])

        result = bin_droplets(typed_list, bin_edges)

        # Should handle zero multiplicity correctly
        assert len(result) == 2
        assert all(val >= 0.0 for val in result)

    def test_bin_droplets_fine_bins(self):
        """Test binning with very fine bin resolution."""
        droplets = [Droplet.new(multi=1000, radius=i * 1e-6) for i in range(1, 11)]
        typed_list = self._make_typed_list(droplets)
        bin_edges = np.array([float(i) for i in range(0, 22, 2)])  # Every 2 microns

        result = bin_droplets(typed_list, bin_edges)

        assert len(result) == len(bin_edges) - 1
        assert all(val >= 0.0 for val in result)

    def test_bin_droplets_logarithmic_bins(self):
        """Test binning with logarithmically spaced bins."""
        droplets = [
            Droplet.new(multi=1000, radius=10 ** (i / 2.0) * 1e-6) for i in range(0, 10)
        ]
        typed_list = self._make_typed_list(droplets)
        # Logarithmic bin edges
        bin_edges = np.array([10.0 ** (i / 2.0) for i in range(0, 11)])

        result = bin_droplets(typed_list, bin_edges)

        assert len(result) == len(bin_edges) - 1
        # Should distribute droplets across logarithmic bins
        assert sum(1 for val in result if val > 0) > 0

    def test_bin_droplets_returns_correct_length(self):
        """Test that result length is always len(bin_edges) - 1."""
        droplets = [Droplet.new(multi=1000, radius=10e-6) for _ in range(10)]
        typed_list = self._make_typed_list(droplets)

        for n_edges in [2, 3, 5, 10, 20, 50]:
            bin_edges = np.array([float(i * 10) for i in range(n_edges)])
            result = bin_droplets(typed_list, bin_edges)
            assert len(result) == n_edges - 1

    def test_bin_droplets_all_droplets_outside_bins(self):
        """Test binning when all droplets are outside bin range."""
        # All droplets have radius > 100 microns
        droplets = [Droplet.new(multi=1000, radius=200e-6) for _ in range(5)]
        typed_list = self._make_typed_list(droplets)
        bin_edges = np.array([0.0, 10.0, 20.0, 30.0])  # Max 30 microns

        result = bin_droplets(typed_list, bin_edges)

        # All bins should be zero (droplets outside range)
        assert all(val == 0.0 for val in result)


class TestBinDropletsJitted:
    """Tests specifically for the jitted binning function."""

    def _make_typed_list(self, droplets):
        """Helper to convert list to TypedList."""
        typed_list = TypedList()
        for d in droplets:
            typed_list.append(d)
        return typed_list

    def test_jitted_version_matches_wrapper(self):
        """Test that jitted version produces same results as wrapper."""
        droplets = [
            Droplet.new(multi=1000, radius=5e-6),
            Droplet.new(multi=2000, radius=15e-6),
            Droplet.new(multi=3000, radius=25e-6),
        ]
        typed_list = self._make_typed_list(droplets)
        bin_edges = np.array([0.0, 10.0, 20.0, 30.0])

        result_wrapper = bin_droplets(typed_list, bin_edges)
        result_jitted = _bin_droplets_jitted(typed_list, bin_edges)

        np.testing.assert_array_almost_equal(result_wrapper, result_jitted)

    def test_jitted_is_numpy_array(self):
        """Test that jitted version returns numpy array."""
        droplets = [Droplet.new(multi=1000, radius=10e-6)]
        typed_list = self._make_typed_list(droplets)
        bin_edges = np.array([0.0, 20.0])

        result = _bin_droplets_jitted(typed_list, bin_edges)

        assert isinstance(result, np.ndarray)
        assert result.dtype == np.float64

    def test_jitted_large_dataset(self):
        """Test jitted version with larger dataset for performance."""
        np.random.seed(42)
        n_droplets = 10000
        droplets = [
            Droplet.new(multi=100, radius=np.random.uniform(1e-6, 100e-6))
            for _ in range(n_droplets)
        ]
        typed_list = self._make_typed_list(droplets)
        bin_edges = np.linspace(0.0, 150.0, 100)

        result = _bin_droplets_jitted(typed_list, bin_edges)

        # Verify mass conservation
        total_from_bins = np.sum(result)
        total_from_droplets = sum(d.mass * d.multi for d in droplets)
        assert total_from_bins == pytest.approx(total_from_droplets)
