"""Tests for sd_numba.utils.math module."""

import numpy as np
import pytest
from sd_numba.utils.math import rolling_median


class TestRollingMedian:
    """Tests for rolling_median function (bottleneck-based)."""

    def test_rolling_median_basic(self):
        """Test basic rolling median functionality."""
        values = np.array([1.0, 5.0, 3.0, 2.0, 4.0])
        result = rolling_median(values, window_size=3)

        assert len(result) == len(values)
        assert all(isinstance(v, (int, float, np.floating)) for v in result)

    def test_rolling_median_single_element(self):
        """Test rolling median with single element."""
        values = np.array([5.0])
        result = rolling_median(values, window_size=1)
        assert len(result) == 1
        assert result[0] == pytest.approx(5.0)

    def test_rolling_median_smooth_noisy_data(self):
        """Test that rolling median smooths noisy data."""
        # Create data with a spike
        values = np.array([1.0, 1.0, 1.0, 10.0, 1.0, 1.0, 1.0])
        result = rolling_median(values, window_size=3)

        # The spike at index 3 should be smoothed
        assert result[3] < 10.0  # Should be smoothed

    def test_rolling_median_constant_data(self):
        """Test rolling median with constant data."""
        values = np.array([5.0] * 10)
        result = rolling_median(values, window_size=3)

        # All values should remain approximately 5.0
        assert all(v == pytest.approx(5.0) for v in result)

    def test_rolling_median_preserves_length(self):
        """Test that rolling median always returns same length as input."""
        for length in [1, 5, 10, 50, 100]:
            values = np.arange(length, dtype=np.float64)
            for window_size in [1, 3, 5]:
                result = rolling_median(values, window_size=window_size)
                assert len(result) == length

    def test_rolling_median_odd_window_required(self):
        """Test that even window size raises error."""
        values = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        with pytest.raises(ValueError, match="window_size must be an odd number"):
            rolling_median(values, window_size=4)

    def test_rolling_median_with_numpy_array(self):
        """Test rolling median returns numpy array."""
        values = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        result = rolling_median(values, window_size=3)

        assert isinstance(result, np.ndarray)
        assert result.dtype == np.float64

    def test_rolling_median_larger_window(self):
        """Test rolling median with larger window."""
        values = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0])
        result = rolling_median(values, window_size=5)

        assert len(result) == len(values)
        # Center values should be the median of 5 values
        assert result[4] == pytest.approx(5.0)  # median of [3,4,5,6,7]

    def test_rolling_median_edge_behavior(self):
        """Test rolling median behavior at edges (edge padding)."""
        values = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        result = rolling_median(values, window_size=3)

        # At edges, should use edge padding
        assert len(result) == len(values)
        # First element padded with edge value, median of [1,1,2] = 1
        # Last element padded with edge value, median of [4,5,5] = 5
