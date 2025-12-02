"""Shared test fixtures for the sd_numba test suite."""

import numpy as np
import pytest
from numba.typed import List as TypedList
from sd_numba.core.droplet import Droplet


@pytest.fixture
def small_droplet():
    """Create a small droplet (10 micron radius)."""
    return Droplet.new(multi=1000, radius=10e-6)


@pytest.fixture
def medium_droplet():
    """Create a medium droplet (50 micron radius)."""
    return Droplet.new(multi=500, radius=50e-6)


@pytest.fixture
def large_droplet():
    """Create a large droplet (100 micron radius)."""
    return Droplet.new(multi=100, radius=100e-6)


@pytest.fixture
def droplet_collection(small_droplet, medium_droplet, large_droplet):
    """Create a collection of droplets for testing (as Python list)."""
    return [small_droplet, medium_droplet, large_droplet]


@pytest.fixture
def typed_droplet_collection(small_droplet, medium_droplet, large_droplet):
    """Create a typed list collection of droplets for testing."""
    typed_list = TypedList()
    typed_list.append(small_droplet)
    typed_list.append(medium_droplet)
    typed_list.append(large_droplet)
    return typed_list


@pytest.fixture
def typed_droplets_64():
    """Create a typed list of 64 droplets for collision tests."""
    typed_list = TypedList()
    for i in range(64):
        typed_list.append(Droplet.new(multi=10000, radius=20e-6 + i * 1e-6))
    return typed_list


@pytest.fixture
def bin_edges():
    """Create bin edges for binning tests."""
    return np.array([0.0, 10.0, 20.0, 30.0, 100.0])
