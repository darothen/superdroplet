from typing import Sequence

import numpy as np
import numpy.typing as npt

from sd_numba.core.droplet import Droplet


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
        List of bin values.
    """

    # Compute a list of the indices corresponding to the sorted droplets
    sorted_indices = np.argsort([d.radius for d in droplets])

    bin_values = np.zeros(len(bin_edges) - 1)
    bin_idx = 0
    droplet_idx = 0
    for right in bin_edges[1:]:
        bin_value = 0.0
        bin_r_max = right * 1e-6

        # We won't vectorize this loop because it already uses an optimization to
        # minimize the number of comparisons we perform.
        for idx in sorted_indices[droplet_idx:]:
            droplet = droplets[idx]
            if droplet.radius < bin_r_max:
                bin_value += droplet.mass * droplet.multi
                droplet_idx += 1
            else:
                break

        bin_values[bin_idx] = bin_value
        bin_idx += 1

    return bin_values
