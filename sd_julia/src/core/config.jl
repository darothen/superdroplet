"""
Configuration for the collision-coalescence model.
"""

"""
    ModelConfig

Configuration for the collision-coalescence model.

# Fields
- `step_seconds::Int`: Timestep duration in seconds
- `delta_v::Float64`: Total parcel volume (m³)
- `num_droplets::Int`: Number of superdroplets
- `kernel::Kernel`: Collision kernel to use
- `debug::Bool`: Enable debug output
"""
struct ModelConfig
    step_seconds::Int
    delta_v::Float64
    num_droplets::Int
    kernel::Kernel
    debug::Bool
end

"""
Custom display for ModelConfig.
"""
function Base.show(io::IO, config::ModelConfig)
    println(io, "\nMODEL SETUP")
    println(io, "ModelConfig(")
    println(io, "    step_seconds = $(config.step_seconds),")
    println(io, "    delta_v = $(config.delta_v),")
    println(io, "    num_droplets = $(config.num_droplets),")
    println(io, "    kernel = $(config.kernel),")
    println(io, "    debug = $(config.debug)")
    print(io, ")")
end
