"""
Mathematical utility functions.
"""

using Random
using Statistics

"""
    knuth_shuffle!(v::AbstractVector)

Perform an in-place Knuth shuffle (Fisher-Yates shuffle) on a vector.

# Arguments
- `v`: The vector to shuffle (modified in place)

# Example
```julia
v = [1, 2, 3, 4, 5]
knuth_shuffle!(v)
```
"""
function knuth_shuffle!(v::AbstractVector)
    l = length(v)
    for n in 0:(l-2)
        i = rand(1:(l-n))
        v[i], v[l-n] = v[l-n], v[i]
    end
    return v
end

"""
    median(values::AbstractVector{T}) -> T where T<:Real

Compute the median of a vector of values.

# Arguments
- `values`: Vector of values

# Returns
- The median value
"""
function median(values::AbstractVector{T})::T where T<:Real
    n = length(values)
    
    # Handle edge cases
    if n == 0
        return zero(T)
    elseif n == 1
        return values[1]
    elseif n == 2
        return (values[1] + values[2]) / 2
    end
    
    # Sort and return median
    sorted = sort(values)
    if isodd(n)
        return sorted[div(n, 2) + 1]
    else
        return (sorted[div(n, 2)] + sorted[div(n, 2) + 1]) / 2
    end
end

"""
    rolling_median(values::Vector{Float64}, window_size::Int) -> Vector{Float64}

Apply rolling median smoothing to a vector of values.

# Arguments
- `values`: Vector of values to smooth
- `window_size`: Size of the rolling window

# Returns
- Smoothed vector of the same length as input
"""
function rolling_median(values::Vector{Float64}, window_size::Int)::Vector{Float64}
    n = length(values)
    smoothed = zeros(Float64, n)
    half_window = div(window_size, 2)
    for i in 1:n
        start_idx = max(1, i - half_window)
        end_idx = min(n, i + half_window)
        window_vals = values[start_idx:end_idx]
        smoothed[i] = median(window_vals)
    end
    
    return smoothed
end
