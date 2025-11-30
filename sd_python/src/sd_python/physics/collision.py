import dataclasses
import math
import random

from sd_python.core.config import ModelConfig
from sd_python.core.constants import FOUR_THIRD, PI, RHO_WATER
from sd_python.core.droplet import Droplet, compute_total_water
from sd_python.utils.math import knuth_shuffle


@dataclasses.dataclass
class CollisionStepResult:
    counter: int
    big_probs: int
    max_prob: float
    min_prob: float
    total_xi: float


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
    excess = multi_j - int(math.floor(gamma_t_multi_k))

    if excess > 0:
        # Case 1: Some droplets from sd_j remain unpaired
        sd_j.multi = excess
        # sd_k.multi stays the same

        # sd_j.rcubed stays the same
        # rcubed_k_p = gamma_t * sd_j.rcubed + sd_k.rcubed
        # Apparently there's a fused-multiply-add op in the Python stdlib
        rcubed_k_p = math.fma(gamma_t, sd_j.rcubed, sd_k.rcubed)
        sd_k.update_rcubed(rcubed_k_p)

        # sd_j.solute stays the same
        sd_k.solute = math.fma(gamma_t, sd_j.solute, sd_k.solute)
    else:
        # Case 2: All droplets from sd_j are paired, split sd_k
        multi_j_p = int(math.floor(multi_k / 2.0))
        multi_k_p = multi_k - multi_j_p

        rcubed_new = math.fma(gamma_t, sd_j.rcubed, sd_k.rcubed)
        solute_new = math.fma(gamma_t, sd_j.solute, sd_k.solute)

        sd_j.multi = multi_j_p
        sd_k.multi = multi_k_p

        sd_j.update_rcubed(rcubed_new)
        sd_k.update_rcubed(rcubed_new)
        sd_j.solute = solute_new
        sd_k.solute = solute_new


def collision_step(
    droplets: list[Droplet], model_config: ModelConfig
) -> CollisionStepResult:
    """Performs one collision-coalescence timestep.

    Args:
        droplets: List of droplets to process.
        model_config: Model configuration.

    Returns:
        CollisionStepResult: Result of the collision-coalescence timestep.
    """

    t_c = model_config.step_seconds
    delta_v = model_config.delta_v
    kernel = model_config.kernel

    # Permute the droplet list
    knuth_shuffle(droplets)

    # Generate candidate pairs
    n_part = float(len(droplets))
    half_n_part = int(
        n_part / 2.0
    )  # n_part should always be a power of 2, so this is safe

    # Pre-compute combined scaling factor to reduce multiplications in loop
    scaling = (n_part * (n_part - 1.0) / 2.0) / float(half_n_part)
    # Pre-generate all random numbers for parallel iteration
    # TODO: Do we need this optimization? I doubt it does anything here.
    # random_numbers = [random.random() for _ in range(half_n_part)]
    scaling_factor = scaling * float(t_c) / delta_v

    counter = 0
    big_probs = 0
    max_prob = 0.0
    min_prob = 1.0

    # Cache kernel method reference to avoid repeated attribute lookups
    kernel_compute = kernel.compute

    for i in range(half_n_part):
        sd_j = droplets[i]
        sd_k = droplets[i + half_n_part]

        # Hoist attribute lookups to local variables
        multi_j = sd_j.multi
        multi_k = sd_k.multi

        # Skip if either droplet has zero multiplicity
        if multi_j == 0 or multi_k == 0:
            continue

        phi = random.random()
        k_ij = kernel_compute(sd_j, sd_k)

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
        prob_floor = math.floor(prob)
        if (prob - prob_floor) > phi:
            gamma = prob_floor + 1.0
            if multi_j < multi_k:
                multi_coalesce(sd_k, sd_j, gamma)
            else:
                multi_coalesce(sd_j, sd_k, gamma)
            counter += 1

    total_xi = sum(d.multi for d in droplets)

    return CollisionStepResult(
        counter=counter,
        big_probs=big_probs,
        max_prob=max_prob,
        min_prob=min_prob,
        total_xi=total_xi,
    )
