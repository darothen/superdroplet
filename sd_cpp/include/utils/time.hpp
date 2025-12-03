#pragma once

#include <string>
#include <compare>

namespace sd_cpp {
namespace utils {

/// A stopwatch to measure elapsed time in minutes and seconds.
class Stopwatch {
public:
    /// Minutes component
    unsigned int minutes;
    
    /// Seconds component (0-59)
    unsigned int seconds;

    /// Creates a new stopwatch initialized to the given number of seconds.
    ///
    /// @param seconds Initial time in seconds
    explicit Stopwatch(unsigned int seconds = 0);

    /// Returns the total time in seconds.
    [[nodiscard]] unsigned int total_seconds() const;

    /// Increments the stopwatch by the given number of seconds.
    ///
    /// @param seconds_to_add Number of seconds to add
    void increment(unsigned int seconds_to_add);

    /// Formats the stopwatch as a string.
    [[nodiscard]] std::string to_string() const;

    /// Comparison operators
    auto operator<=>(const Stopwatch& other) const = default;
};

} // namespace utils
} // namespace sd_cpp
