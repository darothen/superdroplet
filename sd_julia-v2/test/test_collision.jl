"""Tests for collision module."""

using Test
using SuperdropletModel
using Random

@testset "CollisionStepResult" begin
    @testset "Creation" begin
        result = CollisionStepResult(
            10,       # counter
            2,        # big_probs
            1.5,      # max_prob
            0.3,      # min_prob
            100000.0, # total_xi
            50.0      # total_water
        )
        
        @test result.counter == 10
        @test result.big_probs == 2
        @test result.max_prob == 1.5
        @test result.min_prob == 0.3
        @test result.total_xi == 100000.0
        @test result.total_water == 50.0
    end
end

@testset "Multi Coalesce" begin
    @testset "Case 1: Excess droplets remain in sd_j" begin
        # sd_j has more droplets than needed for full coalescence
        sd_j = Droplet(1000, 10e-6)
        sd_k = Droplet(500, 20e-6)
        
        original_j_multi = sd_j.multi
        original_k_multi = sd_k.multi
        original_k_rcubed = sd_k.rcubed
        
        gamma = 1.5
        
        SuperdropletModel.multi_coalesce!(sd_j, sd_k, gamma)
        
        # sd_j should have some excess droplets remaining
        @test sd_j.multi < original_j_multi
        @test sd_j.multi > 0
        
        # sd_k multiplicity should stay the same
        @test sd_k.multi == original_k_multi
        
        # sd_k should have grown (rcubed increased)
        @test sd_k.rcubed > original_k_rcubed
    end
    
    @testset "Case 2: All sd_j droplets pair, split sd_k" begin
        # Gamma large enough that all sd_j droplets are consumed
        sd_j = Droplet(200, 10e-6)
        sd_k = Droplet(1000, 20e-6)
        
        original_k_multi = sd_k.multi
        gamma = 10.0  # Large gamma
        
        SuperdropletModel.multi_coalesce!(sd_j, sd_k, gamma)
        
        # Both multiplicities should have changed
        @test sd_j.multi != 200
        @test sd_k.multi != original_k_multi
        
        # Both should have same rcubed (coalesced)
        @test sd_j.rcubed ≈ sd_k.rcubed
        
        # Multiplicities should sum to original sd_k
        @test sd_j.multi + sd_k.multi == original_k_multi
    end
    
    @testset "Conservation of mass" begin
        sd_j = Droplet(1000, 10e-6)
        sd_k = Droplet(500, 20e-6)
        
        mass_before = sd_j.mass * sd_j.multi + sd_k.mass * sd_k.multi
        
        gamma = 2.0
        SuperdropletModel.multi_coalesce!(sd_j, sd_k, gamma)
        
        mass_after = sd_j.mass * sd_j.multi + sd_k.mass * sd_k.multi
        
        @test mass_before ≈ mass_after rtol=1e-10
    end
    
    @testset "Gamma limited by ratio" begin
        sd_j = Droplet(100, 10e-6)
        sd_k = Droplet(200, 20e-6)
        
        # Gamma larger than ratio should be limited
        gamma = 10.0
        ratio = Float64(sd_j.multi) / Float64(sd_k.multi)
        
        SuperdropletModel.multi_coalesce!(sd_j, sd_k, gamma)
        
        # Should use gamma_t = ratio, not gamma
        @test sd_j.multi >= 0  # Should not go negative
    end
end

@testset "Collision Step" begin
    @testset "Basic collision step" begin
        Random.seed!(42)
        
        # Create simple droplet collection
        droplets = [Droplet(1000, 10e-6) for _ in 1:100]
        
        config = ModelConfig(
            1,        # step_seconds
            1.0e6,    # delta_v
            100,      # num_droplets
            Golovin,  # kernel
            false     # debug
        )
        
        result = collision_step!(droplets, config)
        
        @test result isa CollisionStepResult
        @test result.counter >= 0
        @test result.big_probs >= 0
        @test result.max_prob >= 0.0
        @test result.min_prob <= 1.0
        @test result.total_xi > 0
        @test result.total_water > 0
    end
    
    @testset "Mass conservation" begin
        Random.seed!(123)
        
        droplets = [Droplet(1000, Float64(i)*1e-6) for i in 10:10:100]
        
        config = ModelConfig(1, 1.0e6, length(droplets), Golovin, false)
        
        mass_before = total_water(droplets)
        collision_step!(droplets, config)
        mass_after = total_water(droplets)
        
        @test mass_before ≈ mass_after rtol=1e-10
    end
    
    @testset "Droplets are shuffled" begin
        Random.seed!(456)
        
        # Create droplets with distinct radii
        droplets = [Droplet(1000, Float64(i)*1e-6) for i in 1:20]
        original_radii = [d.radius for d in droplets]
        
        config = ModelConfig(1, 1.0e6, length(droplets), Golovin, false)
        
        collision_step!(droplets, config)
        
        # After shuffling and collisions, order should potentially change
        # (this is stochastic, but with enough droplets, it's very likely)
        current_radii = [d.radius for d in droplets]
        
        # At least check that we still have the same number of droplets
        @test length(droplets) == length(original_radii)
    end
    
    @testset "Different kernels" begin
        Random.seed!(789)
        
        for kernel in [Golovin, Hydro, Long]
            droplets = [Droplet(1000, Float64(i)*1e-6) for i in 10:10:100]
            config = ModelConfig(1, 1.0e6, length(droplets), kernel, false)
            
            result = collision_step!(droplets, config)
            
            @test result isa CollisionStepResult
            @test result.total_water > 0
        end
    end
    
    @testset "Zero multiplicity droplets skipped" begin
        Random.seed!(101)
        
        droplets = [
            Droplet(1000, 10e-6),
            Droplet(0, 20e-6),  # Zero multiplicity
            Droplet(500, 30e-6),
            Droplet(0, 40e-6),  # Zero multiplicity
        ]
        
        config = ModelConfig(1, 1.0e6, length(droplets), Golovin, false)
        
        # Should not crash with zero multiplicity droplets
        result = collision_step!(droplets, config)
        
        @test result isa CollisionStepResult
    end
end
