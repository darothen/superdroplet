using Test
using SuperdropletModel

@testset "SuperdropletModel.jl" begin
    @testset "Constants" begin
        include("test_constants.jl")
    end
    
    @testset "Droplet" begin
        include("test_droplet.jl")
    end
    
    @testset "Kernels" begin
        include("test_kernels.jl")
    end
    
    @testset "Collision" begin
        include("test_collision.jl")
    end
    
    @testset "Time" begin
        include("test_time.jl")
    end
    
    @testset "Math" begin
        include("test_math.jl")
    end
    
    @testset "Binning" begin
        include("test_binning.jl")
    end
end
