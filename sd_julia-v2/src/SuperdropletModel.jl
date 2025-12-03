"""
SuperdropletModel

A Julia implementation of the Shima et al. (2009) superdroplet method for
simulating stochastic collision-coalescence in clouds.

"""
module SuperdropletModel

# Export core types and functions
export Droplet, ModelConfig, Kernel
export Golovin, Hydro, Long
export run_simulation
export total_water, update_rcubed!
export collision_step!, CollisionStepResult
export Stopwatch, increment!, total_seconds
export bin_droplets, rolling_median

# Include files in dependency order:
# 1. constants.jl - no dependencies
# 2. droplet.jl - needs constants
# 3. kernels.jl - needs constants and Droplet, defines Kernel enum
# 4. config.jl - needs Kernel enum from kernels.jl
# 5. collision.jl - needs all of the above

include("core/constants.jl")
include("core/droplet.jl")
include("physics/kernels.jl")
include("core/config.jl")
include("physics/collision.jl")

# Include utility modules
include("utils/math.jl")
include("utils/time.jl")
include("utils/binning.jl")

# Include main simulation driver
include("sce.jl")

end # module SuperdropletModel
