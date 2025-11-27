//! Time tracking utilities for the simulation.

/// A stopwatch to measure elapsed time in minutes and seconds.
///
/// # Example
/// ```
/// use sd_rust::utils::time::{Stopwatch, increment};
/// let mut stopwatch = Stopwatch::new(0);
/// increment(&mut stopwatch, 90);
/// assert_eq!(stopwatch.minutes, 1);
/// assert_eq!(stopwatch.seconds, 30);
/// ```
pub struct Stopwatch {
    /// Minutes component
    pub minutes: u32,
    /// Seconds component (0-59)
    pub seconds: u32,
}

impl Stopwatch {
    /// Creates a new stopwatch initialized to the given number of seconds.
    ///
    /// # Arguments
    /// * `seconds` - Initial time in seconds
    ///
    /// # Example
    /// ```
    /// use sd_rust::utils::time::Stopwatch;
    /// let stopwatch = Stopwatch::new(120);
    /// assert_eq!(stopwatch.minutes, 2);
    /// assert_eq!(stopwatch.seconds, 0);
    /// ```
    pub fn new(seconds: u32) -> Stopwatch {
        let seconds_over = seconds % 60;
        let minutes = (seconds - seconds_over) / 60;
        Stopwatch {
            minutes,
            seconds: seconds_over,
        }
    }

    /// Returns the total time in seconds.
    pub fn total_seconds(&self) -> u32 {
        self.minutes * 60 + self.seconds
    }
}

/// Increments a stopwatch by the given number of seconds.
///
/// # Arguments
/// * `stopwatch` - The stopwatch to increment
/// * `seconds` - Number of seconds to add
pub fn increment(stopwatch: &mut Stopwatch, seconds: u32) {
    stopwatch.seconds += seconds;
    if stopwatch.seconds >= 60 {
        let seconds_overflow = stopwatch.seconds % 60;
        let minutes_overflow = (stopwatch.seconds - seconds_overflow) / 60;
        stopwatch.minutes += minutes_overflow;
        stopwatch.seconds = seconds_overflow;
    }
}

impl std::fmt::Display for Stopwatch {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:4} minutes, {:2} seconds", self.minutes, self.seconds)
    }
}

impl PartialEq for Stopwatch {
    fn eq(&self, other: &Self) -> bool {
        self.total_seconds() == other.total_seconds()
    }
}

impl Eq for Stopwatch {}

impl PartialOrd for Stopwatch {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Stopwatch {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.total_seconds().cmp(&other.total_seconds())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_stopwatch_new() {
        let sw = Stopwatch::new(125);
        assert_eq!(sw.minutes, 2);
        assert_eq!(sw.seconds, 5);
    }

    #[test]
    fn test_increment() {
        let mut sw = Stopwatch::new(50);
        increment(&mut sw, 15);
        assert_eq!(sw.total_seconds(), 65);
        assert_eq!(sw.minutes, 1);
        assert_eq!(sw.seconds, 5);
    }

    #[test]
    fn test_display() {
        let sw = Stopwatch::new(125);
        let s = format!("{}", sw);
        assert!(s.contains("2"));
        assert!(s.contains("5"));
    }
}
