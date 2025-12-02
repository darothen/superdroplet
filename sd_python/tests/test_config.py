"""Tests for sd_python.core.config module."""

import pytest
from sd_python.core.config import ModelConfig
from sd_python.physics.kernels import Kernel


class TestModelConfig:
    """Tests for ModelConfig dataclass."""

    def test_model_config_creation(self):
        """Test creating a ModelConfig with all required fields."""
        config = ModelConfig(
            step_seconds=1,
            delta_v=1e6,
            num_droplets=1024,
            kernel=Kernel.GOLOVIN,
            debug=False,
        )

        assert config.step_seconds == 1
        assert config.delta_v == 1e6
        assert config.num_droplets == 1024
        assert config.kernel == Kernel.GOLOVIN
        assert config.debug is False

    def test_model_config_with_different_kernels(self):
        """Test ModelConfig with different kernel types."""
        for kernel in [Kernel.GOLOVIN, Kernel.HYDRO, Kernel.LONG]:
            config = ModelConfig(
                step_seconds=1,
                delta_v=1e6,
                num_droplets=512,
                kernel=kernel,
                debug=False,
            )
            assert config.kernel == kernel

    def test_model_config_with_debug_enabled(self):
        """Test ModelConfig with debug flag enabled."""
        config = ModelConfig(
            step_seconds=1,
            delta_v=1e6,
            num_droplets=256,
            kernel=Kernel.GOLOVIN,
            debug=True,
        )
        assert config.debug is True

    def test_model_config_with_large_timestep(self):
        """Test ModelConfig with large timestep value."""
        config = ModelConfig(
            step_seconds=10,
            delta_v=1e6,
            num_droplets=1024,
            kernel=Kernel.GOLOVIN,
            debug=False,
        )
        assert config.step_seconds == 10

    def test_model_config_with_small_volume(self):
        """Test ModelConfig with small parcel volume."""
        config = ModelConfig(
            step_seconds=1,
            delta_v=1e3,
            num_droplets=128,
            kernel=Kernel.LONG,
            debug=False,
        )
        assert config.delta_v == 1e3

    def test_model_config_with_large_volume(self):
        """Test ModelConfig with large parcel volume."""
        config = ModelConfig(
            step_seconds=1,
            delta_v=1e9,
            num_droplets=2048,
            kernel=Kernel.HYDRO,
            debug=False,
        )
        assert config.delta_v == 1e9

    def test_model_config_with_power_of_two_droplets(self):
        """Test ModelConfig with power-of-2 number of droplets."""
        for power in [6, 8, 10, 12, 14, 16]:
            num_droplets = 2**power
            config = ModelConfig(
                step_seconds=1,
                delta_v=1e6,
                num_droplets=num_droplets,
                kernel=Kernel.GOLOVIN,
                debug=False,
            )
            assert config.num_droplets == num_droplets

    def test_model_config_string_representation(self):
        """Test ModelConfig string representation."""
        config = ModelConfig(
            step_seconds=1,
            delta_v=1e6,
            num_droplets=1024,
            kernel=Kernel.GOLOVIN,
            debug=False,
        )

        config_str = str(config)

        # Check that string contains key information
        assert "ModelConfig" in config_str
        assert "step_seconds=1" in config_str
        assert "delta_v=1000000.0" in config_str
        assert "num_droplets=1024" in config_str
        assert "kernel" in config_str
        assert "debug=False" in config_str

    def test_model_config_immutability(self):
        """Test that ModelConfig fields can be modified (dataclass is mutable by default)."""
        config = ModelConfig(
            step_seconds=1,
            delta_v=1e6,
            num_droplets=1024,
            kernel=Kernel.GOLOVIN,
            debug=False,
        )

        # Dataclass is mutable by default, so this should work
        config.step_seconds = 2
        assert config.step_seconds == 2

    def test_model_config_equality(self):
        """Test ModelConfig equality comparison."""
        config1 = ModelConfig(
            step_seconds=1,
            delta_v=1e6,
            num_droplets=1024,
            kernel=Kernel.GOLOVIN,
            debug=False,
        )
        config2 = ModelConfig(
            step_seconds=1,
            delta_v=1e6,
            num_droplets=1024,
            kernel=Kernel.GOLOVIN,
            debug=False,
        )
        config3 = ModelConfig(
            step_seconds=2,  # Different
            delta_v=1e6,
            num_droplets=1024,
            kernel=Kernel.GOLOVIN,
            debug=False,
        )

        assert config1 == config2
        assert config1 != config3

    def test_model_config_with_zero_timestep(self):
        """Test ModelConfig with zero timestep (edge case)."""
        config = ModelConfig(
            step_seconds=0,
            delta_v=1e6,
            num_droplets=1024,
            kernel=Kernel.GOLOVIN,
            debug=False,
        )
        assert config.step_seconds == 0

    def test_model_config_realistic_simulation_parameters(self):
        """Test ModelConfig with realistic simulation parameters."""
        config = ModelConfig(
            step_seconds=1,
            delta_v=1e6,  # 1 cubic meter
            num_droplets=2**17,  # 131,072 droplets
            kernel=Kernel.GOLOVIN,
            debug=False,
        )

        assert config.step_seconds == 1
        assert config.delta_v == 1e6
        assert config.num_droplets == 131072
        assert config.kernel == Kernel.GOLOVIN
        assert config.debug is False

    def test_model_config_fields_are_correct_types(self):
        """Test that ModelConfig fields have correct types."""
        config = ModelConfig(
            step_seconds=1,
            delta_v=1e6,
            num_droplets=1024,
            kernel=Kernel.GOLOVIN,
            debug=False,
        )

        assert isinstance(config.step_seconds, int)
        assert isinstance(config.delta_v, float)
        assert isinstance(config.num_droplets, int)
        assert isinstance(config.kernel, Kernel)
        assert isinstance(config.debug, bool)
