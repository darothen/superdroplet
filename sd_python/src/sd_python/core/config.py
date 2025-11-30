"""Configuration for the stochastic collision/coalescence program."""

import dataclasses
import textwrap

from sd_python.physics.kernels import Kernel


@dataclasses.dataclass
class ModelConfig:
    """Configuration for the stochastic collision/coalescence program.

    Args:
        step_seconds: Timestep duration in seconds
        delta_v: Total parcel volume in m^3
        num_droplets: Total number of superdroplets
        kernel: Collision kernel
        debug: Enable debug output
    """

    step_seconds: int
    delta_v: float
    num_droplets: int
    kernel: Kernel
    debug: bool

    def __str__(self):
        return textwrap.dedent(f"""
            ModelConfig(
                step_seconds={self.step_seconds},
                delta_v={self.delta_v},
                num_droplets={self.num_droplets},
                kernel={self.kernel},
                debug={self.debug}
            )
        """)
