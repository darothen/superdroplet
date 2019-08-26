
using Loess

const nRs = 250
const expons = range(0, stop=log10(5e3), length=nRs)
const Rs = 10 .^ expons
const RHO_WATER = 1e3

function gtilde(R, r_grid, xi, sigma, delta_V=1e6)
    nr = length(R)
    nrg = length(r_grid)

    ln_R = log.(R)
    ln_r_grid = log.(r_grid)

    A = 1. / sqrt(2*π) / sigma

    result = zeros(nr)
    for i in 1:nr
        s = 0.
        for j in 1:nrg
            Y = ln_R[i] - ln_r_grid[j]
            expon = (Y^2.) / 2. / (sigma^2.)
            W = A * exp(-expon)
            s += xi[j] * (r_grid[j]^3.) * W
        result[i] = s * RHO_WATER * 4 * π / 3. / delta_V
        end
    end

    return result
end

function kde_plot(droplets::Array{Droplet}, delta_V=1e6)
    r_grid = [d.rcubed for d in droplets]
    r_grid .^= (1. / 3.)
    xi = [d.multi for d in droplets]

    # Use the larger value because it generally smooths things
    # more nicely
    sigma₀ = 1.5 # 0.62
    sigma = sigma₀*(length(r_grid)^(-1. / 5.))

    xx = copy(Rs)
    yy = gtilde(Rs*1e-6, r_grid, xi, sigma, delta_V) * 1e3 # g m^-3

    return xx, yy
end

function bin_droplets(droplets::Array{Droplet})
    mids = 0.5*(Rs[1:end-1] + Rs[2:end])

    # Sort the droplet array... this should make the binning
    # easier
    sorted_droplets = sort(droplets, by=d->d.rcubed)
    n_droplets = length(sorted_droplets)

    bins = []
    bin_idx = 1
    droplet_idx = 1
    for bin in zip(Rs[1:end-1], Rs[2:end])

        # println("bin ", bin_idx, " ", bin)
        # println("droplet_idx is ", droplet_idx)
        _bin::Vector{Droplet} = []
        lo, hi = (bin .* 1e-6) .^ 3
        for idx in droplet_idx:n_droplets

            d = sorted_droplets[idx]
                # println("   ", idx, " ", d.rcubed, " ", hi)
            if d.rcubed < hi
                push!(_bin, d)
                droplet_idx += 1
            else
                break
            end

        end

        push!(bins, _bin)
        bin_idx += 1

    end

    total_waters = [total_water(dl) for dl in bins]
    counts = [length(bin) for bin in bins]
    probs = counts ./ n_droplets

    return mids, total_waters

end
