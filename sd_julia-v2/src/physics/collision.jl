"""
Collision-coalescence algorithm for superdroplets.

Implements the Shima et al. (2009) super-droplet method for simulating
collision-coalescence in clouds.
"""

using Random

"""
    CollisionStepResult

Result of a collision-coalescence timestep.

# Fields
- `counter::Int`: Number of collisions that occurred
- `big_probs::Int`: Number of probabilities > 1.0
- `max_prob::Float64`: Maximum collision probability
- `min_prob::Float64`: Minimum collision probability
- `total_xi::Float64`: Total multiplicity across all droplets
- `total_water::Float64`: Total water mass (kg)
"""
struct CollisionStepResult
    counter::Int
    big_probs::Int
    max_prob::Float64
    min_prob::Float64
    total_xi::Float64
    total_water::Float64
end

"""
    multi_coalesce!(sd_j::Droplet, sd_k::Droplet, gamma::Float64)

Perform a multi-coalescence event between two droplets.

This function implements the coalescence algorithm, updating both droplets
in place according to the multi-coalescence rules. The first droplet (sd_j)
must have multiplicity >= the second droplet (sd_k).

# Arguments
- `sd_j`: First droplet (must have larger or equal multiplicity)
- `sd_k`: Second droplet
- `gamma`: Number of coalescence events
"""
function multi_coalesce!(sd_j::Droplet, sd_k::Droplet, gamma::Float64)
    ratio = Float64(sd_j.multi) / Float64(sd_k.multi)
    gamma_t = gamma < ratio ? gamma : ratio
    excess = sd_j.multi - floor(Int, gamma_t * Float64(sd_k.multi))
    
    if excess > 0
        # Case 1: Some droplets from sd_j remain unpaired
        sd_j.multi = excess
        # sd_k.multi stays the same
        
        # sd_j.rcubed stays the same
        rcubed_k_p = gamma_t * sd_j.rcubed + sd_k.rcubed
        update_rcubed!(sd_k, rcubed_k_p)
        
        # sd_j.solute stays the same
        sd_k.solute = gamma_t * sd_j.solute + sd_k.solute
    else
        # Case 2: All droplets from sd_j are paired, split sd_k
        multi_j_p = floor(Int, Float64(sd_k.multi) / 2.0)
        multi_k_p = sd_k.multi - multi_j_p
        
        rcubed_new = gamma_t * sd_j.rcubed + sd_k.rcubed
        solute_new = gamma_t * sd_j.solute + sd_k.solute
        
        sd_j.multi = multi_j_p
        sd_k.multi = multi_k_p
        
        update_rcubed!(sd_j, rcubed_new)
        update_rcubed!(sd_k, rcubed_new)
        
        sd_j.solute = solute_new
        sd_k.solute = solute_new
    end
end

"""
    collision_step!(droplets::Vector{Droplet}, model_config::ModelConfig) -> CollisionStepResult

Perform one collision-coalescence timestep.

This function implements the main collision algorithm:
1. Shuffle the droplet array
2. Pair droplets (first half with second half)
3. Compute collision probability for each pair
4. Perform coalescence when collision occurs

# Arguments
- `droplets`: Vector of droplets to process (modified in place)
- `model_config`: Model configuration

# Returns
- `CollisionStepResult` with statistics from this timestep
"""
function collision_step!(droplets::Vector{Droplet}, model_config::ModelConfig)::CollisionStepResult
    t_c = model_config.step_seconds
    delta_v = model_config.delta_v
    kernel = model_config.kernel
    
    # Shuffle the droplet list (Knuth shuffle is done by knuth_shuffle!)
    knuth_shuffle!(droplets)
    
    # Generate candidate pairs
    n_part = length(droplets)
    half_n_part = div(n_part, 2)
    scaling = (Float64(n_part) * (Float64(n_part) - 1.0) / 2.0) / Float64(half_n_part)
    
    counter = 0
    big_probs = 0
    max_prob = 0.0
    min_prob = 1.0
    
    # Process pairs
    for i in 1:half_n_part
        sd_j = droplets[i]
        sd_k = droplets[i + half_n_part]
        
        # Generate random number
        phi = rand()
        
        # Compute collision kernel
        k_ij = compute_kernel(kernel, sd_j, sd_k)
        
        max_xi = max(sd_j.multi, sd_k.multi)
        min_xi = min(sd_j.multi, sd_k.multi)
        
        if min_xi == 0
            continue
        end
        
        prob = scaling * Float64(max_xi) * (Float64(t_c) / delta_v) * k_ij
        
        # Update statistics
        max_prob = prob > max_prob ? prob : max_prob
        min_prob = prob < min_prob ? prob : min_prob
        if prob > 1.0
            big_probs += 1
        end
        
        # Check for collision and coalesce if necessary
        if (prob - floor(prob)) > phi
            gamma = floor(prob) + 1.0
            if sd_j.multi < sd_k.multi
                multi_coalesce!(sd_k, sd_j, gamma)
            else
                multi_coalesce!(sd_j, sd_k, gamma)
            end
            counter += 1
        end
    end
    
    # Compute final statistics
    total_xi = Float64(sum(d.multi for d in droplets))
    total_water_mass = total_water(droplets)
    
    if model_config.debug
        println("$counter collisions simulated")
        println("Max/min probabilities (count): $min_prob $max_prob $big_probs")
        println("Total ξ: $total_xi")
    end
    
    return CollisionStepResult(
        counter,
        big_probs,
        max_prob,
        min_prob,
        total_xi,
        total_water_mass
    )
end
