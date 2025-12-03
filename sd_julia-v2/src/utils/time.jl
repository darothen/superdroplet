"""
Time tracking utilities for the simulation.
"""

using Printf

"""
    Stopwatch

A stopwatch to measure elapsed time in minutes and seconds.

# Fields
- `minutes::Int`: Minutes component
- `seconds::Int`: Seconds component (0-59)
"""
mutable struct Stopwatch
    minutes::Int
    seconds::Int
end

"""
    Stopwatch(seconds::Int=0)

Create a new stopwatch initialized to the given number of seconds.

# Arguments
- `seconds`: Initial time in seconds (default: 0)

# Example
```julia
stopwatch = Stopwatch(120)  # 2 minutes, 0 seconds
```
"""
function Stopwatch(seconds::Int=0)
    seconds_over = mod(seconds, 60)
    minutes = div(seconds - seconds_over, 60)
    return Stopwatch(minutes, seconds_over)
end

"""
    total_seconds(stopwatch::Stopwatch) -> Int

Return the total time in seconds.

# Arguments
- `stopwatch`: The stopwatch to query

# Returns
- Total time in seconds
"""
function total_seconds(stopwatch::Stopwatch)::Int
    return stopwatch.minutes * 60 + stopwatch.seconds
end

"""
    increment!(stopwatch::Stopwatch, seconds::Int)

Increment a stopwatch by the given number of seconds.

# Arguments
- `stopwatch`: The stopwatch to increment (modified in place)
- `seconds`: Number of seconds to add
"""
function increment!(stopwatch::Stopwatch, seconds::Int)
    stopwatch.seconds += seconds
    if stopwatch.seconds >= 60
        seconds_overflow = mod(stopwatch.seconds, 60)
        minutes_overflow = div(stopwatch.seconds - seconds_overflow, 60)
        stopwatch.minutes += minutes_overflow
        stopwatch.seconds = seconds_overflow
    end
end

"""
Custom display for Stopwatch.
"""
function Base.show(io::IO, sw::Stopwatch)
    @printf(io, "%4d minutes, %2d seconds", sw.minutes, sw.seconds)
end
