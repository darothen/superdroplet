include("sce_util.jl")
include("plot.jl")

using Distributions
using Printf
using Plots
using RollingFunctions

##################################################################
## SCRIPT PARAMETERS
##################################################################

const DEBUG = false
const PLOT = true

const plot_fn = :bin           # :kde or :bin
# const n_0 = 2.0^23
const n_0 = 27.0*(2^23)        # total number concentration
# const R_0 = 30.531e-6
const R_0 = 30.531e-6 / 3.     # total droplet radius
const X_0 = (4. * π / 3.) * (R_0^3)  # total droplet volume
const M_0 = X_0*RHO_WATER      # total droplet water mass
const delta_v = 1e6            # total parcel volume
const t_c = 1.0                # model timestep (seconds)
const kernel = :long          # Collision kernel

const n_part = 2^17   # Total number of superdroplets
const t_end = 3601.0   # Total simulation time (seconds)
const plot_dt = 1800.0/3 # Output interval time

##################################################################
## SIMULATION STARTS HERE
##################################################################

# Initialize the droplet array;
# First, pre-compute the superdroplet multiplicity
total_droplets = delta_v * n_0
xi_i = floor(Int, total_droplets / n_part)

# Compute droplet masses according to an exponential distribution
dist = Exponential(X_0)
x_grid = rand(dist, n_part)
sort!(x_grid)

# Compute the droplet radius grid
r_grid = (x_grid .* 3. ./ π ./ 4.).^THIRD

# Create the initial droplet array
droplets = [Droplet(xi_i, r^3) for r in r_grid]
m_grid = x_grid .* [d.density for d in droplets]

println("\nGRID SETUP")
@printf "   radii: %5.3g - %5.3g m\n" droplets[1].radius droplets[end].radius
vol_lo, vol_hi = [calc_volume(droplets[1]), calc_volume(droplets[end])]
@printf "  volume: %5.3g - %5.3g m^3\n" vol_lo vol_hi
mass_lo, mass_hi = [calc_mass(droplets[1]), calc_mass(droplets[end])]
@printf "    mass: %5.3g - %5.3g kg\n" mass_lo mass_hi

println("\nSD SETUP")
@printf "   N_s: %d\n" n_part
@printf "  xi_i: %d\n" xi_i
N_per_sd = total_droplets / xi_i / n_part
@printf " N per SD_xi: %d\n" N_per_sd
@printf "(Initialized %d droplets)\n" n_part

if plot_fn == :kde
    xx, yy = kde_plot(droplets)
else
    xx, yy = bin_droplets(droplets)
end
window = 9
yy = rollmedian(yy, window)
xx = xx[window÷2+1:end-window÷2]
if PLOT
    plt = plot(xx, yy, linewidth=2.5, xaxis=:log, label="Initial")
    display(plt)
end

println("\nBEGINNING MAIN ROUTINE")

wm0 = total_water(droplets)
@printf "Initial water mass = %g kg\n" wm0

## Main Loop
t = 0.
ti = 0
masses = [wm0, ]
while t < t_end
    global t, ti, t_c
    ts = Timestamp(t)

    @printf "\n STEP %5d ( %2g min %02g sec )\n" ti ts.minutes ts.seconds

    if PLOT && ((t % plot_dt == 0)) && ti > 0
        println("   PLOTTING ... ")
        s = @sprintf "t = %2g min %2g sec" ts.minutes ts.seconds
        println(s)

        if plot_fn == :kde
            xx, yy = kde_plot(droplets)
        else
            xx, yy = bin_droplets(droplets)
        end
        yy = rollmedian(yy, window)

        xx = xx[window÷2+1:end-window÷2]
        plot!(plt, xx, yy, linewidth=1.5, label=s)
        display(plt)
    end

    collision_step!(droplets, t_c, delta_v, kernel)
    push!(masses, total_water(droplets))

    ti += 1
    t = ti * t_c

end
