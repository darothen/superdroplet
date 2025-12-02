from typing import Sequence

import numpy as np
import numpy.typing as npt
from numba import jit

from sd_numba.core.droplet import Droplet


@jit(nopython=True, fastmath=True)
def _bin_droplets_jitted(
    droplets, bin_edges: npt.NDArray[np.float64]
) -> npt.NDArray[np.float64]:
    """Jitted version of bin_droplets for performance.

    Args:
        droplets: TypedList of Droplet objects to bin.
        bin_edges: NumPy array of bin edges.

    Returns:
        NumPy array of bin values.
    """
    n_droplets = len(droplets)

    # Extract radii, masses, and multis into arrays for efficient sorting
    radii = np.empty(n_droplets, dtype=np.float64)
    masses = np.empty(n_droplets, dtype=np.float64)
    multis = np.empty(n_droplets, dtype=np.int32)

    for i in range(n_droplets):
        radii[i] = droplets[i].radius
        masses[i] = droplets[i].mass
        multis[i] = droplets[i].multi

    sorted_indices = np.argsort(radii)

    bin_values = np.zeros(len(bin_edges) - 1, dtype=np.float64)
    droplet_idx = 0

    for bin_idx in range(len(bin_edges) - 1):
        bin_r_max = bin_edges[bin_idx + 1] * 1e-6
        bin_value = 0.0

        while droplet_idx < n_droplets:
            idx = sorted_indices[droplet_idx]
            if radii[idx] < bin_r_max:
                bin_value += masses[idx] * multis[idx]
                droplet_idx += 1
            else:
                break

        bin_values[bin_idx] = bin_value

    return bin_values


def bin_droplets(
    droplets: Sequence[Droplet], bin_edges: npt.NDArray[np.float64]
) -> npt.NDArray[np.float64]:
    """Bin droplets into bins defined by `bin_edges`.

    Hard-coded to sort `droplets` by the `radius` attribute in each Droplet, and to
    tally up the total mass of droplets in each bin.

    Args:
        droplets: Sequence of Droplet objects to bin.
        bin_edges: NumPy array of bin edges.

    Returns:
        NumPy array of bin values.
    """
    return _bin_droplets_jitted(droplets, bin_edges)
