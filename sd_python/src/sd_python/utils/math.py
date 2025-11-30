"""Mathematical utilities."""

import math
import random
from typing import Any, Sequence


def median(values: Sequence[float]) -> float:
    """Compute the median of a sequence of numbers."""

    n = len(values)
    match n:
        case 0:
            return 0.0
        case 1:
            return values[0]
        case 2:
            return (values[0] + values[1]) / 2.0
        case _:
            # In all other cases, we actuall need to sort the values.
            sorted_values = sorted(values)
            if n & 1:  # Odd number of values, using bitwise AND to check
                return sorted_values[n // 2]
            else:
                return (sorted_values[n // 2 - 1] + sorted_values[n // 2]) / 2.0


def knuth_shuffle(values: Sequence[Any]):
    """Shuffle a sequence of values in-place using the Knuth shuffle algorithm."""
    n = len(values)
    for i in range(n - 1, 0, -1):
        j = random.randint(0, i)
        values[i], values[j] = values[j], values[i]


def linspace(start: float, stop: float, n: int) -> list[float]:
    """Generate a linear space of `n` points between `start` and `stop`."""
    if n == 0:
        return []
    if n == 1:
        return [start]
    step = (stop - start) / (n - 1)
    return [start + i * step for i in range(n)]


def sample_exponential_dist(lam: float) -> float:
    """Generate a sample from an exponential distribution.

    Uses the inverse transform sampling method: if U ~ Uniform(0,1),
    then X = -ln(U)/λ follows an exponential distribution with rate λ.

    Args:
        lam: The rate parameter (lambda) of the exponential distribution.
             Must be positive.

    Returns:
        A sample from the exponential distribution.

    Raises:
        ValueError: If lam <= 0 or n < 0.

    Example:
        >>> samples = [sample_exponential_dist(lam=0.5) for _ in range(1000)]
        >>> # Mean should be approximately 1/lam = 2.0
        >>> mean = sum(samples) / len(samples)
        >>> assert mean == pytest.approx(1/lam, rel=1e-2)
    """
    if lam <= 0:
        raise ValueError(f"Rate parameter lambda must be positive, got {lam}")

    # Use inverse transform sampling: X = -ln(U)/λ where U ~ Uniform(0,1)
    return -math.log(random.random()) / lam


def rolling_median(values: Sequence[float], window: int) -> list[float]:
    """Smooth a sequence of values using a windowed median filter.

    The `window` parameter is the sie of the window looking in one direction -
    so a window size of 3 will look at the current value and the two values before and after.
    """
    n = len(values)
    if n == 0:
        return []
    if n == 1:
        return values

    if window == 0:
        return values

    smoothed_values = []
    for i in range(n):
        start = max(0, i - window)
        end = min(n, i + window + 1)
        smoothed_values.append(median(values[start:end]))
    return smoothed_values
