"""
Collision kernel implementations.

Collision kernels determine the rate at which droplets collide
based on their physical properties.
"""

"""
    Kernel

Enumeration of available collision kernels.

# Values
- `Golovin`: Simple analytical kernel for testing
- `Hydro`: Hydrodynamic kernel with unit efficiencies
- `Long`: Long kernel with collection efficiency parameterization
"""
@enum Kernel begin
    Golovin
    Hydro
    Long
end

# Golovin kernel constant
const GOLOVIN_B = 1500.0
const GOLOVIN_CONSTANT = GOLOVIN_B * FOUR_THIRD * PI

"""
    golovin_kernel(sd_j::Droplet, sd_k::Droplet) -> Float64

Compute the Golovin collision kernel for two droplets.

The Golovin kernel is a simple analytical kernel often used for testing:
K(i,j) = B * (V_i + V_j)

# Arguments
- `sd_j`: First droplet
- `sd_k`: Second droplet

# Returns
- Collision kernel value (m³/s)
"""
function golovin_kernel(sd_j::Droplet, sd_k::Droplet)::Float64
    return GOLOVIN_CONSTANT * (sd_j.rcubed + sd_k.rcubed)
end

"""
    calc_hydro_kernel(e_coal::Float64, e_coll::Float64, r_sum::Float64, tv_diff::Float64) -> Float64

Calculate hydrodynamic kernel value.

# Arguments
- `e_coal`: Coalescence efficiency
- `e_coll`: Collection efficiency
- `r_sum`: Sum of radii
- `tv_diff`: Difference in terminal velocities

# Returns
- Collision kernel value (m³/s)
"""
function calc_hydro_kernel(e_coal::Float64, e_coll::Float64, r_sum::Float64, tv_diff::Float64)::Float64
    return e_coal * e_coll * PI * r_sum * r_sum * abs(tv_diff)
end

"""
    hydro_kernel(sd_j::Droplet, sd_k::Droplet) -> Float64

Compute the hydrodynamic collision kernel for two droplets.

Assumes unit coalescence and collection efficiencies.

# Arguments
- `sd_j`: First droplet
- `sd_k`: Second droplet

# Returns
- Collision kernel value (m³/s)
"""
function hydro_kernel(sd_j::Droplet, sd_k::Droplet)::Float64
    r_j = sd_j.radius
    r_k = sd_k.radius
    tv_j = sd_j.terminal_velocity
    tv_k = sd_k.terminal_velocity
    
    tv_diff = tv_j - tv_k
    r_sum = r_j + r_k
    
    return calc_hydro_kernel(1.0, 1.0, r_sum, tv_diff)
end

"""
    long_kernel(sd_j::Droplet, sd_k::Droplet) -> Float64

Compute the Long collision kernel for two droplets.

Uses a parameterization of collection efficiency based on droplet sizes.

# Arguments
- `sd_j`: First droplet
- `sd_k`: Second droplet

# Returns
- Collision kernel value (m³/s)
"""
function long_kernel(sd_j::Droplet, sd_k::Droplet)::Float64
    r_j = sd_j.radius
    r_k = sd_k.radius
    tv_j = sd_j.terminal_velocity
    tv_k = sd_k.terminal_velocity
    
    tv_diff = tv_j - tv_k
    r_sum = r_j + r_k
    
    # Convert radii to microns for collection efficiency calculation
    r_small, r_large = if r_j < r_k
        (r_j * 1.0e6, r_k * 1.0e6)
    else
        (r_k * 1.0e6, r_j * 1.0e6)
    end
    
    # Collection efficiency cut-off in limit of very large drops
    e_coll = if r_large >= 50.0
        1.0
    else
        4.5e-4 * r_large * r_large * (1.0 - 3.0 / (max(3.0, r_small) + 1.0e-2))
    end
    
    # Limit collection efficiency to 0 <= E_coll <= 1.0
    e_coll = min(e_coll, 1.0)
    e_coll = max(e_coll, 0.0)
    
    return calc_hydro_kernel(1.0, e_coll, r_sum, tv_diff)
end

"""
    compute_kernel(kernel::Kernel, sd_j::Droplet, sd_k::Droplet) -> Float64

Dispatch function to compute collision kernel based on kernel type.

# Arguments
- `kernel`: Kernel type to use
- `sd_j`: First droplet
- `sd_k`: Second droplet

# Returns
- Collision kernel value (m³/s)
"""
function compute_kernel(kernel::Kernel, sd_j::Droplet, sd_k::Droplet)::Float64
    if kernel == Golovin
        return golovin_kernel(sd_j, sd_k)
    elseif kernel == Hydro
        return hydro_kernel(sd_j, sd_k)
    elseif kernel == Long
        return long_kernel(sd_j, sd_k)
    else
        error("Unknown kernel type: $kernel")
    end
end
