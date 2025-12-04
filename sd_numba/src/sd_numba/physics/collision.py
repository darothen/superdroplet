import dataclasses

import numpy as np
from numba import jit
from numba.typed import List as TypedList

from sd_numba.core.config import ModelConfig
from sd_numba.core.droplet import Droplet
from sd_numba.physics.kernels import Kernel, compute_kernel


@dataclasses.dataclass
class CollisionStepResult:
    counter: int
    big_probs: int
    max_prob: float
    min_prob: float
    total_xi: int


@jit(nopython=True, fastmath=True)
def multi_coalesce(sd_j: Droplet, sd_k: Droplet, gamma: float):
    """Performs a multi-coalescence event.

    Args:
        sd_j: Droplet to coalesce.
        sd_k: Droplet to coalesce.
        gamma: Multi-coalescence factor.
    """
    # Cache attribute lookups
    multi_j = sd_j.multi
    multi_k = sd_k.multi

    ratio = multi_j / multi_k
    gamma_t = gamma if gamma < ratio else ratio

    # Pre-compute floor once
    gamma_t_multi_k = gamma_t * multi_k
    excess = multi_j - int(np.floor(gamma_t_multi_k))

    if excess > 0:
        # Case 1: Some droplets from sd_j remain unpaired
        sd_j.multi = excess
        # sd_k.multi stays the same

        # sd_j.rcubed stays the same
        rcubed_k_p = gamma_t * sd_j.rcubed + sd_k.rcubed
        sd_k.update_rcubed(rcubed_k_p)

        # sd_j.solute stays the same
        sd_k.solute = gamma_t * sd_j.solute + sd_k.solute
    else:
        # Case 2: All droplets from sd_j are paired, split sd_k
        multi_j_p = int(np.floor(multi_k / 2.0))
        multi_k_p = multi_k - multi_j_p

        rcubed_new = gamma_t * sd_j.rcubed + sd_k.rcubed
        solute_new = gamma_t * sd_j.solute + sd_k.solute

        sd_j.multi = multi_j_p
        sd_k.multi = multi_k_p

        sd_j.update_rcubed(rcubed_new)
        sd_k.update_rcubed(rcubed_new)
        sd_j.solute = solute_new
        sd_k.solute = solute_new


def collision_step(
    typed_droplets: TypedList[Droplet], model_config: ModelConfig
) -> CollisionStepResult:
    """Performs one collision-coalescence timestep."""
    result_tup = _collision_step(
        typed_droplets,
        model_config.step_seconds,
        model_config.delta_v,
        int(model_config.kernel),
    )
    return CollisionStepResult(*result_tup)


@jit(nopython=True, fastmath=True)
def _collision_step(
    droplets: TypedList[Droplet], t_c: float, delta_v: float, kernel: Kernel
) -> tuple[int, int, float, float, int]:
    """Performs one collision-coalescence timestep.

    Args:
        droplets: List of droplets to process.
        t_c: Timestep duration in seconds.
        delta_v: Total parcel volume in m^3.
        kernel: Collision kernel.

    Returns:
        Tuple containing:
            counter: Number of collisions.
            big_probs: Number of probabilities greater than 1.
            max_prob: Maximum probability.
            min_prob: Minimum probability.
            total_xi: Total multiplicity.
    """

    # Generate candidate pairs
    n_part = len(droplets)
    half_n_part = n_part // 2  # n_part should always be a power of 2, so this is safe

    # Pre-compute combined scaling factor to reduce multiplications in loop
    n_part_float = float(n_part)
    scaling = (n_part_float * (n_part_float - 1.0) / 2.0) / float(half_n_part)
    scaling_factor = scaling * float(t_c) / delta_v

    counter = 0
    big_probs = 0
    max_prob = 0.0
    min_prob = 1.0

    # Shuffle the droplet list by permuting the indices
    sampled_indices = np.random.permutation(n_part)

    for pair_idx in range(half_n_part):
        # Pair up: split into two halves and pair corresponding elements
        # (0, half_n_part), (1, half_n_part+1), (2, half_n_part+2), etc.
        idx_j = sampled_indices[pair_idx]
        idx_k = sampled_indices[pair_idx + half_n_part]

        sd_j = droplets[idx_j]
        sd_k = droplets[idx_k]

        # Hoist attribute lookups to local variables
        multi_j = sd_j.multi
        multi_k = sd_k.multi

        # Skip if either droplet has zero multiplicity
        if multi_j == 0 or multi_k == 0:
            continue

        phi = np.random.rand()
        k_ij = compute_kernel(kernel, sd_j, sd_k)

        # Determine max/min multiplicity
        if multi_j < multi_k:
            max_xi = multi_k
            min_xi = multi_j
        else:
            max_xi = multi_j
            min_xi = multi_k

        prob = scaling_factor * float(max_xi) * k_ij

        # Update statistics
        if prob > max_prob:
            max_prob = prob
        if prob < min_prob:
            min_prob = prob
        if prob > 1:
            big_probs += 1

        # Check for collision
        prob_floor = np.floor(prob)
        if (prob - prob_floor) > phi:
            gamma = prob_floor + 1.0
            if multi_j < multi_k:
                multi_coalesce(sd_k, sd_j, gamma)
            else:
                multi_coalesce(sd_j, sd_k, gamma)
            counter += 1

    total_xi = 0
    for droplet in droplets:
        total_xi += droplet.multi

    return (counter, big_probs, max_prob, min_prob, total_xi)
