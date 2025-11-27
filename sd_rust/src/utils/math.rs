//! Mathematical utility functions.

/// Computes 2^n for integer exponents.
///
/// # Arguments
/// * `n` - The exponent
///
/// # Returns
/// The result of 2^n as f64
///
/// # Example
/// ```
/// use sd_rust::utils::math::pow2;
/// assert_eq!(pow2(10), 1024.0);
/// ```
pub fn pow2(n: usize) -> f64 {
    let base: f64 = 2.0;
    base.powi(n as i32)
}

/// Computes the median of a slice of f64 values.
///
/// # Arguments
/// * `values` - The slice of values to compute the median of
///
/// # Returns
/// The median of the values
///
/// # Example
/// ```
/// use sd_rust::utils::math::median;
/// let values = vec![1.0, 2.0, 3.0, 4.0, 5.0];
/// assert_eq!(median(&values), 3.0);
/// ```
pub fn median(values: &[f64]) -> f64 {
    let n = values.len();
    // Short circuit some base cases
    match n {
        0 => return 0.0,
        1 => return values[0],
        2 => return (values[0] + values[1]) / 2.0,
        _ => (),
    }
    // This isn't going to be the fastest, but since our window size is very small,
    // we can easily sort the values in the window and then return the median
    // through inspection of the even / odd length cases.
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.total_cmp(b));
    if (n & 1) == 1 {
        sorted[n / 2]
    } else {
        (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
    }
}

/// Performs an in-place Knuth shuffle (Fisher-Yates shuffle) on a slice.
///
/// # Arguments
/// * `v` - The slice to shuffle
/// * `rng` - Random number generator
///
/// # Example
/// ```
/// use sd_rust::utils::math::knuth_shuffle;
/// let mut v = vec![1, 2, 3, 4, 5];
/// let mut rng = rand::rng();
/// knuth_shuffle(&mut v, &mut rng);
/// ```
pub fn knuth_shuffle<T>(v: &mut [T], rng: &mut impl rand::Rng) {
    let l = v.len();
    for n in 0..l {
        let i = rng.random_range(0..=(l - n - 1));
        v.swap(i, l - n - 1);
    }
}

/// Generates a linear grid of n points between start and stop.
///
/// # Arguments
/// * `start` - The start of the grid
/// * `stop` - The end of the grid
/// * `n` - The number of points in the grid
///
/// # Returns
/// A vector of n evenly-spaced points between start and stop
///
/// # Example
/// ```
/// use sd_rust::utils::math::generate_linear_grid;
/// let grid = generate_linear_grid(0.0, 1.0, 11);
/// assert_eq!(grid.len(), 11);
/// assert!((grid[0] - 0.0).abs() < 1e-10);
/// assert!((grid[10] - 1.0).abs() < 1e-10);
/// ```
pub fn generate_linear_grid(start: f64, stop: f64, n: usize) -> Vec<f64> {
    if n == 0 {
        return Vec::new();
    }
    if n == 1 {
        return vec![start];
    }
    let mut result = Vec::with_capacity(n);
    let step = (stop - start) / (n - 1) as f64;

    for i in 0..n {
        result.push(start + (i as f64) * step);
    }

    result
}

/// Returns the minimum of two f64 values.
///
/// # Arguments
/// * `a` - The first value
/// * `b` - The second value
///
/// # Returns
/// The minimum of the two values
///
/// # Example
/// ```
/// use sd_rust::utils::math::min_f64;
/// assert_eq!(min_f64(1.0, 2.0), 1.0);
/// assert_eq!(min_f64(2.0, 1.0), 1.0);
/// ```
pub fn min_f64(a: f64, b: f64) -> f64 {
    if a < b { a } else { b }
}

/// Returns the maximum of two f64 values.
///
/// # Arguments
/// * `a` - The first value
/// * `b` - The second value
///
/// # Returns
/// The maximum of the two values
///
/// # Example
/// ```
/// use sd_rust::utils::math::max_f64;
/// assert_eq!(max_f64(1.0, 2.0), 2.0);
/// assert_eq!(max_f64(2.0, 1.0), 2.0);
/// ```
pub fn max_f64(a: f64, b: f64) -> f64 {
    if a > b { a } else { b }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pow2() {
        assert_eq!(pow2(0), 1.0);
        assert_eq!(pow2(1), 2.0);
        assert_eq!(pow2(10), 1024.0);
    }

    #[test]
    fn test_median_odd() {
        let values = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        assert_eq!(median(&values), 3.0);
    }

    #[test]
    fn test_median_even() {
        let values = vec![1.0, 2.0, 3.0, 4.0];
        assert_eq!(median(&values), 2.5);
    }

    #[test]
    fn test_generate_linear_grid() {
        let grid = generate_linear_grid(0.0, 1.0, 11);
        assert_eq!(grid.len(), 11);
        assert!((grid[0] - 0.0).abs() < 1e-10);
        assert!((grid[5] - 0.5).abs() < 1e-10);
        assert!((grid[10] - 1.0).abs() < 1e-10);
    }
}
