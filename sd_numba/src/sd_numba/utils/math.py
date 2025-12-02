"""Mathematical utilities."""

import bottleneck as bn
import numpy as np
import numpy.typing as npt


def rolling_median(
    values: npt.NDArray[np.float64],
    window_size: int,
) -> npt.NDArray[np.float64]:
    """Compute the centered rolling median of a 1D NumPy array.

    Args:
        values: The input 1D array.
        window_size: The size of the rolling window.

    Returns:
        The array containing the rolling medians.
    """
    if window_size % 2 == 0:
        raise ValueError("window_size must be an odd number for symmetrical padding.")

    # Pad the array to center the window
    # For a centered window, we need to pad by window_size // 2 on each side
    half_window = window_size // 2
    padded_values = np.pad(values, pad_width=half_window, mode="edge")

    # Apply the trailing median to the padded array
    result = bn.move_median(padded_values, window_size, min_count=1, axis=-1)

    # Extract the centered portion (remove the padding effect)
    # The result has extra elements at the start due to padding, so we slice appropriately
    return result[window_size - 1 :]
