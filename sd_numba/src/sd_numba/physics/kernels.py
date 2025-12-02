import enum

from numba import jit

from sd_numba.core.constants import FOUR_THIRD, PI
from sd_numba.core.droplet import Droplet


class Kernel(enum.IntEnum):
    """Collision kernel enumeration."""

    GOLOVIN = 1
    HYDRO = 2
    LONG = 3


@jit(nopython=True, fastmath=True)
def compute_kernel(kernel_type: int, sd_j: Droplet, sd_k: Droplet) -> float:
    """Calculate the collision kernel for two droplets."""
    if kernel_type == Kernel.GOLOVIN:
        return golovin_kernel(sd_j, sd_k)
    elif kernel_type == Kernel.HYDRO:
        return hydro_kernel(sd_j, sd_k)
    elif kernel_type == Kernel.LONG:
        return long_kernel(sd_j, sd_k)
    else:
        return 0.0


#: Golovin collision kernel constant (m³/s)
GOLOVIN_CONSTANT = 1500.0 * FOUR_THIRD * PI


@jit(nopython=True, fastmath=True)
def golovin_kernel(sd_j: Droplet, sd_k: Droplet) -> float:
    """Calculate the Golovin collision kernel for two droplets.

    The Golovin collision kernel is a simple analytical kernel often used for testing:
    K(i,j) = B * (V_i + V_j)

    Args:
        sd_j (Droplet): First droplet
        sd_k (Droplet): Second droplet

    Returns:
        Collision kernel value (m³/s)
    """
    return GOLOVIN_CONSTANT * (sd_j.rcubed + sd_k.rcubed)


@jit(nopython=True, fastmath=True)
def calc_hydro_kernel(
    e_coal: float, e_coll: float, r_sum: float, tv_diff: float
) -> float:
    """Calculate the hydrodynamic collision kernel."""
    return e_coal * e_coll * PI * r_sum * r_sum * abs(tv_diff)


@jit(nopython=True, fastmath=True)
def hydro_kernel(sd_j: Droplet, sd_k: Droplet) -> float:
    """Calculate the hydrodynamic collision kernel for two droplets.

    The hydrodynamic collision kernel is a simple analytical kernel often used for testing:
    K(i,j) = E_coal * E_coll * PI * r_sum * r_sum * abs(tv_diff)

    Args:
        sd_j (Droplet): First droplet
        sd_k (Droplet): Second droplet

    Returns:
        Collision kernel value (m³/s)
    """
    r_j = sd_j.radius
    r_k = sd_k.radius
    tv_j = sd_j.terminal_velocity
    tv_k = sd_k.terminal_velocity
    tv_diff = tv_j - tv_k
    r_sum = r_j + r_k
    return calc_hydro_kernel(1.0, 1.0, r_sum, tv_diff)


@jit(nopython=True, fastmath=True)
def long_kernel(sd_j: Droplet, sd_k: Droplet) -> float:
    """Calculate the long collision kernel for two droplets.

    The long collision kernel is a simple analytical kernel often used for testing:
    K(i,j) = PI * r_sum * r_sum * tv_diff

    Args:
        sd_j (Droplet): First droplet
        sd_k (Droplet): Second droplet

    Returns:
        Collision kernel value (m³/s)
    """
    r_j = sd_j.radius
    r_k = sd_k.radius
    tv_j = sd_j.terminal_velocity
    tv_k = sd_k.terminal_velocity
    tv_diff = tv_j - tv_k
    r_sum = r_j + r_k

    # Convert radii to microns for collection efficiency calculation
    (r_small, r_large) = (r_j, r_k) if r_j < r_k else (r_k, r_j)
    r_small = r_small * 1e6  # convert to microns
    r_large = r_large * 1e6  # convert to microns

    if r_large >= 50.0:
        e_coll = 1.0
    else:
        e_coll = 4.5e-4 * r_large * r_large * (1.0 - 3.0 / (max(3.0, r_small) + 1e-2))

    e_coll = max(0.0, min(e_coll, 1.0))
    return calc_hydro_kernel(1.0, e_coll, r_sum, tv_diff)
