"""Tests for sd_numba.core.constants module."""

import math

import pytest
from sd_numba.core import constants


def test_pi_constant():
    """Test that PI constant is correct."""
    assert constants.PI == math.pi


def test_rho_water_constant():
    """Test that water density is correct."""
    assert constants.RHO_WATER == 1000.0


def test_rho_air_constant():
    """Test that air density is correct."""
    assert constants.RHO_AIR == 1.0


def test_third_constant():
    """Test that THIRD constant is 1/3."""
    assert constants.THIRD == pytest.approx(1.0 / 3.0)


def test_three_fourth_constant():
    """Test that THREE_FOURTH constant is 3/4."""
    assert constants.THREE_FOURTH == pytest.approx(3.0 / 4.0)


def test_four_third_constant():
    """Test that FOUR_THIRD constant is 4/3."""
    assert constants.FOUR_THIRD == pytest.approx(4.0 / 3.0)
