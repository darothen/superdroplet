__precompile__()

using Printf

const RHO_WATER = 1000.0
const RHO_AIR = 1.0
const THIRD = 1. / 3.
const MULTI_THRESH = 1e4

## droplet.f90

droplet_count = 0

mutable struct Droplet
    multi::Int
    rcubed::Float64
    solute::Float64
    density::Float64
    radius::Float64
    id::Int

    function Droplet(multi, rcubed)

        global droplet_count
        dc = droplet_count
        droplet_count += 1

        new(multi, rcubed, 0., 1000., rcubed^THIRD, dc)
    end
end

function calc_volume(droplet::Droplet)
    vol = droplet.rcubed*4.0*π/3.
    return vol
end

function calc_mass(droplet::Droplet)
    mass = droplet.density*calc_volume(droplet)
    return mass
end

function calc_terminal_v(droplet::Droplet)
    # Beard (1976)

    r = droplet.rcubed^THIRD
    d = 2. * r * 1e6 # diameter, m -> μm
    x = calc_mass(droplet) * 1e3 # mass, kg -> g

    if d <= 134.43
        alpha = 4.5795e5
        x_to_beta = x^(2. * THIRD)
    elseif 134.43 < d && d <= 1511.64
        alpha = 4962.0
        x_to_beta = x^THIRD
    elseif 1511.64 < d && d <= 3477.84
        alpha = 1732.0
        x_to_beta = x^(THIRD / 2.)
    else
        alpha = 917.0
        x_to_beta = 1.0
    end

    tv = 1e-2 * alpha * x_to_beta # cm/s -> m/s
    return tv
end

## math.f90

function exp_dist_moments(x::Float64, n0::Float64, x0::Float64, l::Int=0)::Float64
    x^l * (n0 / x0) * exp(-1. * x / x0)
end

function expon_pdf(x::Float64)::Float64
    scale::Float64 = x
    λ::Float64 = 1.0 / scale
    λ * exp(-1. * λ * x)
end

## util.f90

mutable struct Timestamp
    minutes::Number
    seconds::Number

    function Timestamp(seconds::Number)
        seconds_over = mod(seconds, 60)
        minutes = floor((seconds - seconds_over) / 60)
        new(minutes, seconds_over)
    end
end

function total_water(droplets::Array{Droplet})::Float64
    sum(map(d -> calc_mass(d)*d.multi / 1000., droplets))
end


function knuthshuffle!(v::AbstractVector)
    # https://www.rosettacode.org/wiki/Knuth_shuffle#Julia
    for i in length(v):-1:2
        j = rand(1:i)
        v[i], v[j] = v[j], v[i]
    end
    return v
end
# knuthshuffle!(v::AbstractVector) = knuthshuffle!(Base.Random.GLOBAL_RNG, v)
# v = collect(1:20)
# println("# v = $v\n   -> ", knuthshuffle!(v))

# collisions.f90

function calc_hydro_kernel(E_coal::Float64, E_coll::Float64,
                           r_sum::Float64, tv_diff::Float64)
    E_coal * E_coll * π * r_sum*r_sum * abs(tv_diff)
end

const GOLOVIN_B = 1500.
function golovin_kernel(sd_j::Droplet, sd_k::Droplet)
    GOLOVIN_B * (sd_j.rcubed + sd_k.rcubed) * 4. * π * THIRD
end

function hydro_kernel(sd_j::Droplet, sd_k::Droplet)
    r_j = sd_j.rcubed^THIRD
    r_k = sd_k.rcubed^THIRD
    tv_j = calc_terminal_v(sd_j)
    tv_k = calc_terminal_v(sd_k)

    tv_diff = tv_j - tv_k
    r_sum = r_j + r_k

    return calc_hydro_kernel(1., 1., r_sum, tv_diff)
end

function long_kernel(sd_j::Droplet, sd_k::Droplet)
    r_j = sd_j.rcubed^THIRD
    r_k = sd_k.rcubed^THIRD
    x_j = calc_mass(sd_j)
    x_k = calc_mass(sd_k)
    tv_j = calc_terminal_v(sd_j)
    tv_k = calc_terminal_v(sd_k)

    tv_diff = tv_j - tv_k
    r_sum = r_j + r_k

    r_small = min(r_j, r_k)*1e6 # convert to micron
    r_large = max(r_j, r_k)*1e6

    # Collection efficiency cut-off in limit of very
    # large drops
    if r_large >= 50
        E_coll = 1.
    else
        E_coll = 4.5e-4 * r_large*r_large * (1.0 - 3.0/(max(3.0, r_small) + 1e-2))
    end

    # Limit collection efficiency to 0 <= E_coll <= 1.0
    E_coll = min(E_coll, 1.)
    E_coll = max(0., E_coll)

    return calc_hydro_kernel(1., E_coll, r_sum, tv_diff)

end

function multi_coalesce!(sd_j::Droplet, sd_k::Droplet, γ::Float64)
    # γt = min(γ, round(Int, sd_j.multi ÷ sd_k.multi))
    γt = min(γ, floor(sd_j.multi / sd_k.multi))
    # excess = sd_j.multi - round(Int, γt * sd_k.multi)
    excess = sd_j.multi - floor(Int, γt * sd_k.multi)

    # println("mc ", excess, " ", sd_j.multi, " ", floor(γt * sd_k.multi))

    if excess > 0
        multi_j_p = excess
        multi_k_p = sd_k.multi

        rcubed_j_p = sd_j.rcubed
        rcubed_k_p = γt * rcubed_j_p + sd_k.rcubed

        solute_j_p = sd_j.solute
        solute_k_p = γt * solute_j_p + sd_k.solute
    else  # implies excess == 0
        # multi_j_p = round(Int, sd_k.multi / 2)
        multi_j_p = floor(Int, sd_k.multi / 2.0)
        multi_k_p = sd_k.multi - multi_j_p

        rcubed_j_p = γt * sd_j.rcubed + sd_k.rcubed
        rcubed_k_p = rcubed_j_p

        solute_j_p = γt * sd_j.solute + sd_k.solute
        solute_k_p = solute_j_p
    end

    sd_j.multi = multi_j_p
    sd_j.rcubed = rcubed_j_p
    sd_j.solute = solute_j_p

    sd_k.multi = multi_k_p
    sd_k.rcubed = rcubed_k_p
    sd_k.solute = solute_k_p
end

function collision_step!(droplets::Array{Droplet}, t_c::Float64,
                         delta_v::Float64, kernel::Symbol)
    println("PREP STEPS")

    # Make the random permutation of the droplet list
    println("    SHUFFLE LIST")
    knuthshuffle!(droplets)

    # Make the candidate pairs
    println("    GEN PAIRS")

    # Generate the uniform random numbers)
    println("PROBABILITY LOOP")
    n_part = length(droplets)
    scaling = (n_part * (n_part - 1) / 2.) / floor(n_part / 2.)

    println("PROB / COLLISION LOOP")
    counter = 0
    half_n_part = n_part ÷ 2

    big_probs = 0
    max_prob = 0.
    min_prob = 1.0

    for i in 1:half_n_part
        sd_j = droplets[i]
        sd_k = droplets[i + half_n_part]

        M₀ = calc_mass(sd_j) + calc_mass(sd_k)

        # Check for NAN in droplet sizes
        # TODO

        ϕ = rand()
        if kernel == :golovin
            K_ij = golovin_kernel(sd_j, sd_k)
        elseif kernel == :hydro
            K_ij = hydro_kernel(sd_j, sd_k)
        elseif kernel == :long
            K_ij = long_kernel(sd_j, sd_k)
        else  # Golovin default
            K_ij = golovin_kernel(sd_j, sd_k)
        end

        max_ξ = max(sd_j.multi, sd_k.multi)
        min_ξ = min(sd_j.multi, sd_k.multi)
        if min_ξ == 0
            continue
        end
        prob = scaling * max_ξ * (t_c / delta_v) * K_ij

        max_prob = prob > max_prob ? prob : max_prob
        min_prob = prob < min_prob ? prob : min_prob
        if prob > 1
            big_probs += 1
        end

        # Check for collision and coalesce if necessary
        if prob - floor(prob) > ϕ
            γ = floor(prob) + 1.

            if sd_j.multi < sd_k.multi
                multi_coalesce!(sd_k, sd_j, γ)
            else
                multi_coalesce!(sd_j, sd_k, γ)
            end

            # Need to copy back? I don't think so...

            counter += 1
        end

        M₁ = calc_mass(sd_j) + calc_mass(sd_k)

    end

    @printf("%5g collisions simulated\n", counter)
    @printf("Max/min probabilities (count): %8.2f %8.2f %5g\n",
            min_prob, max_prob, big_probs)
    @printf("Total ξ: %g\n", sum([d.multi for d in droplets]))
end
