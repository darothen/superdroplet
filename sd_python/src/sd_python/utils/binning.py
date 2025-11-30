from typing import Sequence

from sd_python.core.droplet import Droplet


def bin_droplets(
    droplets: Sequence[Droplet], bin_edges: Sequence[float]
) -> list[float]:
    """Bin droplets into bins defined by `bin_edges`.

    Hard-coded to sort `droplets` by the `radius` attribute in each Droplet, and to
    tally up the total mass of droplets in each bin.

    Args:
        droplets: Sequence of droplets to bin.
        bin_edges: Sequence of bin edges.

    Returns:
        List of bin values.
    """

    # Compute a list of the indices corresponding to the sorted droplets
    sorted_indices = sorted(range(len(droplets)), key=lambda i: droplets[i].radius)

    bin_values = [0.0] * (len(bin_edges) - 1)
    bin_idx = 0
    droplet_idx = 0
    for left, right in zip(bin_edges[:-1], bin_edges[1:]):
        bin_value = 0.0
        bin_r_max = right * 1e-6

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
