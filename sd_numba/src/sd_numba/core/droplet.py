"""
Superdroplet data structure and related methods.
"""

import numpy as np
from numba import jit, types
from numba.experimental import jitclass

from sd_numba.core.constants import FOUR_THIRD, PI, RHO_WATER


@jitclass(
    spec=[
        ("multi", types.int64),
        ("rcubed", types.float64),
        ("solute", types.float64),
        ("density", types.float64),
        ("radius", types.float64),
        ("volume", types.float64),
        ("mass", types.float64),
        ("terminal_velocity", types.float64),
    ]
)
class Droplet:
    """Represents a superdroplet in the simulation.

    A superdroplet is a computational particle that represents many
    identical physical droplets. The `multi` field indicates how many
    real droplets this superdroplet represents.

    In this implementation, we eagerly compute and cache important derived properties
    to avoid expensive recalculations.

    Attributes:
        multi (int): The multiplicity of the superdroplet.
        rcubed (float): The cube of the droplet radius (m^3).
        solute (float): The solute mass of the droplet (kg).
        density (float): The density of the droplet (kg/m^3).
        radius (float): The radius of the droplet (m).
        volume (float): The volume of the droplet (m^3).
        mass (float): The mass of the droplet (kg).
        terminal_velocity (float): The terminal velocity of the droplet (m/s).
    """

    multi: int
    rcubed: float
    solute: float
    density: float
    radius: float
    volume: float
    mass: float
    terminal_velocity: float

    def __init__(
        self,
        multi: int,
        rcubed: float,
        solute: float,
        density: float,
        radius: float,
        volume: float,
        mass: float,
        terminal_velocity: float,
    ):
        self.multi = multi
        self.rcubed = rcubed
        self.solute = solute
        self.density = density
        self.radius = radius
        self.volume = volume
        self.mass = mass
        self.terminal_velocity = terminal_velocity

    @staticmethod
    def new(multi: int, radius: float) -> "Droplet":
        """Create a new Droplet instance.

        Args:
            multi (int): The multiplicity of the superdroplet.
            radius (float): The radius of the droplet in meters.

        Returns:
            Droplet: A new Droplet instance.
        """
        rcubed = radius * radius * radius
        volume = rcubed * FOUR_THIRD * PI
        mass = volume * RHO_WATER
        terminal_velocity = compute_terminal_velocity(radius, mass)

        return Droplet(
            multi,
            rcubed,
            0.0,
            RHO_WATER,
            radius,
            volume,
            mass,
            terminal_velocity,
        )

    def update_rcubed(self, new_rcubed: float):
        """Update the rcubed value and recalculate the derived properties."""
        self.rcubed = new_rcubed
        self.radius = np.cbrt(new_rcubed)
        self.volume = new_rcubed * FOUR_THIRD * PI
        self.mass = self.volume * RHO_WATER
        self.terminal_velocity = compute_terminal_velocity(self.radius, self.mass)


@jit(nopython=True, fastmath=True)
def compute_terminal_velocity(radius: float, mass: float) -> float:
    """Compute the terminal velocity of a droplet using the Beard (1976) parameterization.

    Args:
        radius: The radius of the droplet in meters.
        mass: The mass of the droplet in kilograms.

    Returns:
        The terminal velocity of the droplet in meters per second.
    """
    d = 2.0 * radius * 1e6  # diameter, m -> μm
    x = mass * 1e3  # mass, kg -> g

    if d <= 134.43:
        alpha = 4.5795e5
        cbrt_x = np.cbrt(x)
        # x**(2/3) = cbrt(x)^2 = cbrt(x) * cbrt(x)
        x_to_beta = cbrt_x * cbrt_x
    elif d <= 1511.64:
        alpha = 4962.0
        x_to_beta = np.cbrt(x)
    elif d <= 3477.84:
        alpha = 1732.0
        # x**(1/6) = sqrt(cbrt(x))
        x_to_beta = np.sqrt(np.cbrt(x))
    else:
        alpha = 917.0
        x_to_beta = 1.0

    return 1e-2 * alpha * x_to_beta  # cm/s -> m/s


@jit(nopython=True, fastmath=True)
def _compute_total_water_jitted(droplets) -> float:
    """Compute the total water mass using jitted loop.

    Args:
        droplets: TypedList of Droplet instances

    Returns:
        The total water mass in kilograms.
    """
    total = 0.0
    for droplet in droplets:
        total += droplet.mass * droplet.multi
    return total


def compute_total_water(droplets) -> float:
    """Compute the total water mass of a list of droplets.

    Args:
        droplets: The list of Droplet instances (TypedList or list)

    Returns:
        The total water mass in kilograms.
    """
    return _compute_total_water_jitted(droplets)
