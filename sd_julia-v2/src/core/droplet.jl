"""
Superdroplet data structure and methods.

A superdroplet represents a collection of identical cloud droplets,
characterized by their multiplicity and physical properties.
"""

"""
    Droplet

Represents a superdroplet in the simulation.

A superdroplet is a computational particle that represents many
identical physical droplets. The `multi` field indicates how many
real droplets this superdroplet represents.

# Fields
- `multi::Int`: Multiplicity - number of real droplets represented
- `rcubed::Float64`: Radius cubed (m³) - stored for efficiency
- `solute::Float64`: Solute mass (kg)
- `density::Float64`: Density (kg/m³)
- `radius::Float64`: Cached radius (m)
- `volume::Float64`: Cached volume (m³)
- `mass::Float64`: Cached mass (kg)
- `terminal_velocity::Float64`: Cached terminal velocity (m/s)
"""
mutable struct Droplet
    multi::Int
    rcubed::Float64
    solute::Float64
    density::Float64
    radius::Float64
    volume::Float64
    mass::Float64
    terminal_velocity::Float64
end

"""
    Droplet(multi::Int, radius::Float64)

Create a new droplet with the given multiplicity and radius.

# Arguments
- `multi`: Number of real droplets this superdroplet represents
- `radius`: Radius of the droplet in meters

# Returns
- A new `Droplet` instance with pre-computed cached properties

# Example
```julia
droplet = Droplet(1000, 1.0e-6)  # 1000 droplets of 1 micron radius
```
"""
function Droplet(multi::Int, radius::Float64)
    rcubed = radius^3
    volume = rcubed * FOUR_THIRD * PI
    mass = volume * RHO_WATER
    terminal_velocity = compute_terminal_velocity(radius, mass)
    
    return Droplet(
        multi,
        rcubed,
        0.0,  # solute
        RHO_WATER,
        radius,
        volume,
        mass,
        terminal_velocity
    )
end

"""
    update_rcubed!(droplet::Droplet, new_rcubed::Float64)

Update the radius cubed and recalculate cached properties.

This method should be called whenever the droplet size changes
(e.g., after coalescence).

# Arguments
- `droplet`: The droplet to update
- `new_rcubed`: New value for radius cubed (m³)
"""
function update_rcubed!(droplet::Droplet, new_rcubed::Float64)
    droplet.rcubed = new_rcubed
    droplet.radius = cbrt(new_rcubed)
    droplet.volume = new_rcubed * FOUR_THIRD * PI
    droplet.mass = droplet.volume * RHO_WATER
    droplet.terminal_velocity = compute_terminal_velocity(droplet.radius, droplet.mass)
end

"""
    compute_terminal_velocity(radius::Float64, mass::Float64) -> Float64

Compute terminal velocity using Beard (1976) parameterization.

# Arguments
- `radius`: Droplet radius in meters
- `mass`: Droplet mass in kilograms

# Returns
- Terminal velocity in m/s
"""
function compute_terminal_velocity(radius::Float64, mass::Float64)::Float64
    d = 2.0 * radius * 1.0e6  # diameter, m -> μm
    x = mass * 1.0e3  # mass, kg -> g
    
    if d <= 134.43
        alpha = 4.5795e5
        cbrt_x = cbrt(x)
        # x^(2/3) = cbrt(x)^2
        x_to_beta = cbrt_x * cbrt_x
    elseif d <= 1511.64
        alpha = 4962.0
        x_to_beta = cbrt(x)
    elseif d <= 3477.84
        alpha = 1732.0
        # x^(1/6) = sqrt(cbrt(x))
        x_to_beta = sqrt(cbrt(x))
    else
        alpha = 917.0
        x_to_beta = 1.0
    end
    
    return 1.0e-2 * alpha * x_to_beta  # cm/s -> m/s
end

"""
    total_water(droplets::Vector{Droplet}) -> Float64

Calculate the total water mass in a collection of droplets.

# Arguments
- `droplets`: Vector of droplets to sum

# Returns
- Total water mass in kilograms
"""
function total_water(droplets::Vector{Droplet})::Float64
    return sum(d.mass * Float64(d.multi) for d in droplets)
end
