"""Tests for sd_python.utils.time module."""

import pytest
from sd_python.utils.time import Stopwatch


class TestStopwatchCreation:
    """Tests for Stopwatch creation."""

    def test_stopwatch_zero_seconds(self):
        """Test creating stopwatch with 0 seconds."""
        sw = Stopwatch(0)
        assert sw.minutes == 0
        assert sw.seconds == 0

    def test_stopwatch_less_than_minute(self):
        """Test creating stopwatch with less than a minute."""
        sw = Stopwatch(45)
        assert sw.minutes == 0
        assert sw.seconds == 45

    def test_stopwatch_exactly_one_minute(self):
        """Test creating stopwatch with exactly 60 seconds."""
        sw = Stopwatch(60)
        assert sw.minutes == 1
        assert sw.seconds == 0

    def test_stopwatch_minutes_and_seconds(self):
        """Test creating stopwatch with minutes and seconds."""
        sw = Stopwatch(90)
        assert sw.minutes == 1
        assert sw.seconds == 30

    def test_stopwatch_multiple_minutes(self):
        """Test creating stopwatch with multiple minutes."""
        sw = Stopwatch(185)  # 3 minutes and 5 seconds
        assert sw.minutes == 3
        assert sw.seconds == 5

    def test_stopwatch_large_value(self):
        """Test creating stopwatch with large number of seconds."""
        sw = Stopwatch(3661)  # 61 minutes and 1 second
        assert sw.minutes == 61
        assert sw.seconds == 1


class TestStopwatchIncrement:
    """Tests for Stopwatch increment method."""

    def test_increment_basic(self):
        """Test basic increment operation."""
        sw = Stopwatch(0)
        sw.increment(10)
        assert sw.seconds == 10
        assert sw.minutes == 0

    def test_increment_no_overflow(self):
        """Test increment without overflow to minutes."""
        sw = Stopwatch(30)
        sw.increment(15)
        assert sw.seconds == 45
        assert sw.minutes == 0

    def test_increment_with_overflow(self):
        """Test increment with overflow to minutes."""
        sw = Stopwatch(50)
        sw.increment(20)
        assert sw.seconds == 10
        assert sw.minutes == 1

    def test_increment_exactly_to_minute(self):
        """Test increment that results in exactly a minute boundary."""
        sw = Stopwatch(40)
        sw.increment(20)
        assert sw.seconds == 0
        assert sw.minutes == 1

    def test_increment_multiple_minutes(self):
        """Test increment that adds multiple minutes."""
        sw = Stopwatch(30)
        sw.increment(150)  # 2 minutes and 30 seconds
        assert sw.seconds == 0
        assert sw.minutes == 3

    def test_increment_from_existing_minutes(self):
        """Test increment from stopwatch already having minutes."""
        sw = Stopwatch(120)  # 2 minutes
        sw.increment(45)
        assert sw.seconds == 45
        assert sw.minutes == 2

    def test_increment_multiple_times(self):
        """Test multiple increments."""
        sw = Stopwatch(0)
        sw.increment(30)
        sw.increment(30)
        sw.increment(30)
        assert sw.seconds == 30
        assert sw.minutes == 1

    def test_increment_zero(self):
        """Test increment by zero."""
        sw = Stopwatch(45)
        sw.increment(0)
        assert sw.seconds == 45
        assert sw.minutes == 0


class TestStopwatchTotalSeconds:
    """Tests for Stopwatch total_seconds method."""

    def test_total_seconds_zero(self):
        """Test total_seconds for zero time."""
        sw = Stopwatch(0)
        assert sw.total_seconds() == 0

    def test_total_seconds_only_seconds(self):
        """Test total_seconds with only seconds."""
        sw = Stopwatch(45)
        assert sw.total_seconds() == 45

    def test_total_seconds_only_minutes(self):
        """Test total_seconds with only minutes."""
        sw = Stopwatch(120)
        assert sw.total_seconds() == 120

    def test_total_seconds_mixed(self):
        """Test total_seconds with minutes and seconds."""
        sw = Stopwatch(185)  # 3 minutes and 5 seconds
        assert sw.total_seconds() == 185

    def test_total_seconds_after_increment(self):
        """Test total_seconds after increment."""
        sw = Stopwatch(60)
        sw.increment(45)
        assert sw.total_seconds() == 105


class TestStopwatchStringRepresentation:
    """Tests for Stopwatch string representations."""

    def test_str_zero(self):
        """Test string representation of zero time."""
        sw = Stopwatch(0)
        assert str(sw) == "   0 min  0 sec"

    def test_str_only_seconds(self):
        """Test string representation with only seconds."""
        sw = Stopwatch(45)
        assert str(sw) == "   0 min 45 sec"

    def test_str_only_minutes(self):
        """Test string representation with only minutes."""
        sw = Stopwatch(120)
        assert str(sw) == "   2 min  0 sec"

    def test_str_mixed(self):
        """Test string representation with minutes and seconds."""
        sw = Stopwatch(185)
        assert str(sw) == "   3 min  5 sec"

    def test_repr_basic(self):
        """Test repr representation."""
        sw = Stopwatch(90)
        assert repr(sw) == "Stopwatch(minutes=1, seconds=30)"


class TestStopwatchComparison:
    """Tests for Stopwatch comparison operations."""

    def test_equality_same_time(self):
        """Test equality for same time."""
        sw1 = Stopwatch(90)
        sw2 = Stopwatch(90)
        assert sw1 == sw2

    def test_equality_different_time(self):
        """Test inequality for different times."""
        sw1 = Stopwatch(90)
        sw2 = Stopwatch(120)
        assert sw1 != sw2

    def test_less_than(self):
        """Test less than comparison."""
        sw1 = Stopwatch(60)
        sw2 = Stopwatch(120)
        assert sw1 < sw2
        assert not sw2 < sw1

    def test_less_than_or_equal(self):
        """Test less than or equal comparison."""
        sw1 = Stopwatch(60)
        sw2 = Stopwatch(120)
        sw3 = Stopwatch(120)
        assert sw1 <= sw2
        assert sw2 <= sw3
        assert not sw2 <= sw1

    def test_greater_than(self):
        """Test greater than comparison."""
        sw1 = Stopwatch(120)
        sw2 = Stopwatch(60)
        assert sw1 > sw2
        assert not sw2 > sw1

    def test_greater_than_or_equal(self):
        """Test greater than or equal comparison."""
        sw1 = Stopwatch(120)
        sw2 = Stopwatch(60)
        sw3 = Stopwatch(120)
        assert sw1 >= sw2
        assert sw1 >= sw3
        assert not sw2 >= sw1

    def test_comparison_zero(self):
        """Test comparison with zero time."""
        sw_zero = Stopwatch(0)
        sw_pos = Stopwatch(10)
        assert sw_zero < sw_pos
        assert sw_zero <= sw_pos
        assert sw_pos > sw_zero
        assert sw_pos >= sw_zero

    def test_comparison_after_increment(self):
        """Test comparison after increment operations."""
        sw1 = Stopwatch(60)
        sw2 = Stopwatch(30)
        sw2.increment(30)
        assert sw1 == sw2

    def test_sorting(self):
        """Test that stopwatches can be sorted."""
        stopwatches = [
            Stopwatch(120),
            Stopwatch(30),
            Stopwatch(90),
            Stopwatch(0),
            Stopwatch(180),
        ]
        sorted_stopwatches = sorted(stopwatches)

        assert sorted_stopwatches[0].total_seconds() == 0
        assert sorted_stopwatches[1].total_seconds() == 30
        assert sorted_stopwatches[2].total_seconds() == 90
        assert sorted_stopwatches[3].total_seconds() == 120
        assert sorted_stopwatches[4].total_seconds() == 180


class TestStopwatchEdgeCases:
    """Tests for Stopwatch edge cases."""

    def test_large_increment(self):
        """Test increment with very large value."""
        sw = Stopwatch(0)
        sw.increment(7200)  # 2 hours
        assert sw.minutes == 120
        assert sw.seconds == 0

    def test_multiple_overflows(self):
        """Test multiple overflow scenarios."""
        sw = Stopwatch(59)
        sw.increment(61)  # Should overflow twice
        assert sw.minutes == 2
        assert sw.seconds == 0

    def test_increment_preserves_consistency(self):
        """Test that increment maintains consistency with total_seconds."""
        sw = Stopwatch(0)
        increments = [15, 30, 45, 60, 90, 120]
        expected_total = 0

        for inc in increments:
            sw.increment(inc)
            expected_total += inc
            assert sw.total_seconds() == expected_total
