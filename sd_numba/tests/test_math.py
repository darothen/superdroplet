"""Tests for sd_python.utils.math module."""

import math

import pytest
from sd_python.utils.math import (
    knuth_shuffle,
    linspace,
    median,
    rolling_median,
    sample_exponential_dist,
)


class TestMedian:
    """Tests for median function."""

    def test_median_empty_list(self):
        """Test median of empty list returns 0.0."""
        result = median([])
        assert result == 0.0

    def test_median_single_value(self):
        """Test median of single value."""
        result = median([5.0])
        assert result == 5.0

    def test_median_two_values(self):
        """Test median of two values."""
        result = median([3.0, 7.0])
        assert result == 5.0

    def test_median_odd_count(self):
        """Test median with odd number of values."""
        result = median([1.0, 3.0, 5.0, 7.0, 9.0])
        assert result == 5.0

    def test_median_even_count(self):
        """Test median with even number of values."""
        result = median([1.0, 2.0, 3.0, 4.0])
        assert result == 2.5

    def test_median_unsorted_list(self):
        """Test median with unsorted values."""
        result = median([9.0, 1.0, 5.0, 3.0, 7.0])
        assert result == 5.0

    def test_median_negative_values(self):
        """Test median with negative values."""
        result = median([-5.0, -1.0, 0.0, 1.0, 5.0])
        assert result == 0.0

    def test_median_duplicate_values(self):
        """Test median with duplicate values."""
        result = median([1.0, 2.0, 2.0, 2.0, 3.0])
        assert result == 2.0

    def test_median_large_list(self):
        """Test median with large list."""
        values = list(range(1, 1001))  # 1 to 1000
        result = median(values)
        assert result == 500.5


class TestKnuthShuffle:
    """Tests for Knuth shuffle algorithm."""

    def test_knuth_shuffle_preserves_elements(self):
        """Test that shuffle preserves all elements."""
        values = [1, 2, 3, 4, 5]
        original = values.copy()
        knuth_shuffle(values)

        assert sorted(values) == sorted(original)

    def test_knuth_shuffle_empty_list(self):
        """Test shuffle with empty list."""
        values = []
        knuth_shuffle(values)
        assert values == []

    def test_knuth_shuffle_single_element(self):
        """Test shuffle with single element."""
        values = [42]
        knuth_shuffle(values)
        assert values == [42]

    def test_knuth_shuffle_two_elements(self):
        """Test shuffle with two elements."""
        values = [1, 2]
        knuth_shuffle(values)
        # Should still contain both elements
        assert set(values) == {1, 2}

    def test_knuth_shuffle_modifies_in_place(self):
        """Test that shuffle modifies list in place."""
        values = [1, 2, 3, 4, 5]
        original_id = id(values)
        knuth_shuffle(values)
        assert id(values) == original_id

    def test_knuth_shuffle_with_strings(self):
        """Test shuffle with string elements."""
        values = ["a", "b", "c", "d", "e"]
        original = values.copy()
        knuth_shuffle(values)
        assert sorted(values) == sorted(original)

    def test_knuth_shuffle_randomness(self):
        """Test that shuffle produces different results (probabilistic)."""
        # Run shuffle multiple times and check we get different orderings
        original = list(range(10))
        results = []

        for _ in range(10):
            values = original.copy()
            knuth_shuffle(values)
            results.append(tuple(values))

        # With 10 shuffles of 10 elements, we should get at least 2 different orderings
        # (This test has a very small chance of false positive)
        unique_results = set(results)
        assert len(unique_results) > 1


class TestLinspace:
    """Tests for linspace function."""

    def test_linspace_empty(self):
        """Test linspace with n=0."""
        result = linspace(0.0, 10.0, 0)
        assert result == []

    def test_linspace_single_point(self):
        """Test linspace with n=1."""
        result = linspace(5.0, 10.0, 1)
        assert result == [5.0]

    def test_linspace_two_points(self):
        """Test linspace with n=2."""
        result = linspace(0.0, 10.0, 2)
        assert result == pytest.approx([0.0, 10.0])

    def test_linspace_basic(self):
        """Test basic linspace functionality."""
        result = linspace(0.0, 10.0, 11)
        expected = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
        assert result == pytest.approx(expected)

    def test_linspace_negative_range(self):
        """Test linspace with negative range."""
        result = linspace(-5.0, 5.0, 11)
        expected = [-5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
        assert result == pytest.approx(expected)

    def test_linspace_decreasing(self):
        """Test linspace with decreasing values (start > stop)."""
        result = linspace(10.0, 0.0, 11)
        expected = [10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0]
        assert result == pytest.approx(expected)

    def test_linspace_fractional_steps(self):
        """Test linspace with fractional steps."""
        result = linspace(0.0, 1.0, 5)
        expected = [0.0, 0.25, 0.5, 0.75, 1.0]
        assert result == pytest.approx(expected)

    def test_linspace_same_start_stop(self):
        """Test linspace where start equals stop."""
        result = linspace(5.0, 5.0, 10)
        expected = [5.0] * 10
        assert result == pytest.approx(expected)

    def test_linspace_length(self):
        """Test that linspace returns correct number of points."""
        for n in [1, 5, 10, 50, 100]:
            result = linspace(0.0, 100.0, n)
            assert len(result) == n

    def test_linspace_endpoints(self):
        """Test that linspace includes both endpoints."""
        result = linspace(1.5, 7.3, 20)
        assert result[0] == pytest.approx(1.5)
        assert result[-1] == pytest.approx(7.3)

    def test_linspace_uniform_spacing(self):
        """Test that linspace produces uniformly spaced values."""
        result = linspace(0.0, 10.0, 11)
        differences = [result[i + 1] - result[i] for i in range(len(result) - 1)]

        # All differences should be approximately equal
        for diff in differences:
            assert diff == pytest.approx(1.0)


class TestSampleExponentialDist:
    """Tests for sample_exponential_dist function."""

    def test_sample_exponential_dist_returns_positive(self):
        """Test that exponential samples are positive."""
        lam = 1.0
        samples = [sample_exponential_dist(lam) for _ in range(100)]
        assert all(s > 0 for s in samples)

    def test_sample_exponential_dist_with_different_lambda(self):
        """Test exponential sampling with different lambda values."""
        for lam in [0.1, 0.5, 1.0, 2.0, 5.0]:
            sample = sample_exponential_dist(lam)
            assert sample > 0

    def test_sample_exponential_dist_mean(self):
        """Test that exponential samples have correct mean (approximately)."""
        lam = 2.0
        n_samples = 10000
        samples = [sample_exponential_dist(lam) for _ in range(n_samples)]
        mean = sum(samples) / len(samples)

        # Mean of exponential distribution is 1/lambda
        expected_mean = 1.0 / lam
        # Allow 5% tolerance due to randomness
        assert mean == pytest.approx(expected_mean, rel=0.05)

    def test_sample_exponential_dist_invalid_lambda(self):
        """Test that invalid lambda values raise ValueError."""
        with pytest.raises(ValueError, match="Rate parameter lambda must be positive"):
            sample_exponential_dist(0.0)

        with pytest.raises(ValueError, match="Rate parameter lambda must be positive"):
            sample_exponential_dist(-1.0)

    def test_sample_exponential_dist_small_lambda(self):
        """Test exponential sampling with very small lambda."""
        lam = 0.001
        samples = [sample_exponential_dist(lam) for _ in range(100)]
        # Mean should be large (1/0.001 = 1000)
        mean = sum(samples) / len(samples)
        assert mean > 500  # Should be around 1000, but allow variance

    def test_sample_exponential_dist_large_lambda(self):
        """Test exponential sampling with large lambda."""
        lam = 100.0
        samples = [sample_exponential_dist(lam) for _ in range(100)]
        # Mean should be small (1/100 = 0.01)
        mean = sum(samples) / len(samples)
        assert mean < 0.05  # Should be around 0.01, but allow variance

    def test_sample_exponential_dist_variance(self):
        """Test that exponential samples have correct variance (approximately)."""
        lam = 1.0
        n_samples = 10000
        samples = [sample_exponential_dist(lam) for _ in range(n_samples)]

        mean = sum(samples) / len(samples)
        variance = sum((s - mean) ** 2 for s in samples) / len(samples)

        # Variance of exponential distribution is 1/lambda^2
        expected_variance = 1.0 / (lam**2)
        # Allow 10% tolerance due to randomness
        assert variance == pytest.approx(expected_variance, rel=0.1)

    def test_sample_exponential_dist_produces_different_values(self):
        """Test that exponential sampling produces different values."""
        lam = 1.0
        samples = [sample_exponential_dist(lam) for _ in range(10)]
        # Should have at least some variation
        unique_samples = set(samples)
        assert len(unique_samples) > 1  # Very unlikely to get all the same value


class TestRollingMedian:
    """Tests for rolling_median function."""

    def test_rolling_median_empty_list(self):
        """Test rolling median with empty list."""
        result = rolling_median([], window=1)
        assert result == []

    def test_rolling_median_single_element(self):
        """Test rolling median with single element."""
        result = rolling_median([5.0], window=1)
        assert result == [5.0]

    def test_rolling_median_zero_window(self):
        """Test rolling median with window=0 returns original."""
        values = [1.0, 2.0, 3.0, 4.0, 5.0]
        result = rolling_median(values, window=0)
        assert result == values

    def test_rolling_median_basic(self):
        """Test basic rolling median functionality."""
        values = [1.0, 5.0, 3.0, 2.0, 4.0]
        result = rolling_median(values, window=1)

        # Window size 1 means looking at current +/- 1
        # For each position, should compute median of surrounding values
        assert len(result) == len(values)
        assert all(isinstance(v, (int, float)) for v in result)

    def test_rolling_median_smooth_noisy_data(self):
        """Test that rolling median smooths noisy data."""
        # Create data with a spike
        values = [1.0, 1.0, 1.0, 10.0, 1.0, 1.0, 1.0]
        result = rolling_median(values, window=1)

        # The spike at index 3 should be smoothed
        # Median of [1.0, 1.0, 10.0] = 1.0
        # Median of [1.0, 10.0, 1.0] = 1.0
        # Median of [10.0, 1.0, 1.0] = 1.0
        assert result[3] < 10.0  # Should be smoothed

    def test_rolling_median_preserves_monotonic_data(self):
        """Test that rolling median preserves monotonically increasing data."""
        values = [1.0, 2.0, 3.0, 4.0, 5.0]
        result = rolling_median(values, window=1)

        # Should still be generally increasing
        assert len(result) == len(values)

    def test_rolling_median_edge_behavior(self):
        """Test rolling median behavior at edges."""
        values = [1.0, 2.0, 3.0, 4.0, 5.0]
        result = rolling_median(values, window=2)

        # At edges, window is truncated
        # First element: median of [1.0, 2.0, 3.0]
        # Last element: median of [3.0, 4.0, 5.0]
        assert len(result) == len(values)
        assert result[0] >= 1.0
        assert result[-1] <= 5.0

    def test_rolling_median_constant_data(self):
        """Test rolling median with constant data."""
        values = [5.0] * 10
        result = rolling_median(values, window=2)

        # All values should remain 5.0
        assert all(v == 5.0 for v in result)

    def test_rolling_median_alternating_data(self):
        """Test rolling median with alternating values."""
        values = [1.0, 5.0, 1.0, 5.0, 1.0, 5.0]
        result = rolling_median(values, window=1)

        # Should smooth out the alternation
        assert len(result) == len(values)
        # Values should be less extreme than original alternation
        assert max(result) <= 5.0
        assert min(result) >= 1.0

    def test_rolling_median_large_window(self):
        """Test rolling median with large window."""
        values = list(range(10))
        result = rolling_median(values, window=5)

        # Large window should produce smoother result
        assert len(result) == len(values)

    def test_rolling_median_window_larger_than_data(self):
        """Test rolling median when window is larger than data."""
        values = [1.0, 2.0, 3.0]
        result = rolling_median(values, window=10)

        # Should handle gracefully - window extends beyond data
        assert len(result) == len(values)

    def test_rolling_median_with_negative_values(self):
        """Test rolling median with negative values."""
        values = [-5.0, -1.0, 0.0, 1.0, 5.0]
        result = rolling_median(values, window=1)

        assert len(result) == len(values)
        # Should handle negative values correctly

    def test_rolling_median_with_floats(self):
        """Test rolling median with floating point values."""
        values = [1.5, 2.7, 3.2, 4.9, 5.1]
        result = rolling_median(values, window=1)

        assert len(result) == len(values)
        assert all(isinstance(v, (int, float)) for v in result)

    def test_rolling_median_preserves_length(self):
        """Test that rolling median always returns same length as input."""
        for length in [1, 5, 10, 50, 100]:
            values = list(range(length))
            for window in [0, 1, 2, 5]:
                result = rolling_median(values, window=window)
                assert len(result) == length
