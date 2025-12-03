"""
Droplet binning for analysis and output.

Bins droplets by size for visualization and analysis purposes.
"""

"""
    bin_droplets(droplets::Vector{Droplet}, bin_edges::Vector{Float64}) -> Vector{Float64}

Bin droplets by size and compute total water mass in each bin.

Uses logarithmically-spaced bins and sorts droplets by size for
efficient binning.

# Arguments
- `droplets`: Vector of droplets to bin
- `bin_edges`: Bin edges in microns (μm)

# Returns
- Vector of total water mass (kg) in each bin

# Example
```julia
# Create logarithmically-spaced bins from 1 to 1000 microns
bin_edges = [10.0^exp for exp in range(0.0, 3.0, length=101)]
bin_values = bin_droplets(droplets, bin_edges)
```
"""
function bin_droplets(droplets::Vector{Droplet}, bin_edges::Vector{Float64})::Vector{Float64}
    n_bins = length(bin_edges) - 1
    bin_values = zeros(Float64, n_bins)
    
    # Sort indices by radius (not rcubed) for correct binning
    sorted_indices = sortperm([d.radius for d in droplets])
    
    droplet_idx = 1
    n_droplets = length(droplets)
    
    # Iterate over bin edges
    for bin_idx in 1:n_bins
        bin_value = 0.0
        r_max = bin_edges[bin_idx + 1] * 1.0e-6  # Convert microns to meters
        
        # Add all droplets that fall in this bin
        while droplet_idx <= n_droplets
            idx = sorted_indices[droplet_idx]
            droplet = droplets[idx]
            
            if droplet.radius < r_max
                bin_value += droplet.mass * Float64(droplet.multi)
                droplet_idx += 1
            else
                break
            end
        end
        
        bin_values[bin_idx] = bin_value
    end
    
    return bin_values
end
