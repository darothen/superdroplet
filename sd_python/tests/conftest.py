"""Shared test fixtures for the sd_python test suite."""

import pytest
from sd_python.core.droplet import Droplet


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
    """Create a collection of droplets for testing."""
    return [small_droplet, medium_droplet, large_droplet]
