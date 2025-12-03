"""Tests for droplet module."""

using Test
using SuperdropletModel

# Import constants for easier use
const PI = SuperdropletModel.PI
const FOUR_THIRD = SuperdropletModel.FOUR_THIRD
const THREE_FOURTH = SuperdropletModel.THREE_FOURTH
const RHO_WATER = SuperdropletModel.RHO_WATER

@testset "Droplet Creation" begin
    @testset "Basic properties" begin
        multi = 1000
        radius = 10e-6  # 10 microns
        
        droplet = Droplet(multi, radius)
        
        @test droplet.multi == multi
        @test droplet.radius == radius
        @test droplet.solute == 0.0
        @test droplet.density == RHO_WATER
    end
    
    @testset "Computed properties" begin
        radius = 10e-6
        droplet = Droplet(1, radius)
        
        expected_rcubed = radius^3
        expected_volume = FOUR_THIRD * PI * expected_rcubed
        expected_mass = expected_volume * RHO_WATER
        
        @test droplet.rcubed ≈ expected_rcubed
        @test droplet.volume ≈ expected_volume
        @test droplet.mass ≈ expected_mass
    end
    
    @testset "Terminal velocity" begin
        droplet = Droplet(1, 10e-6)
        @test droplet.terminal_velocity > 0
    end
end

@testset "Terminal Velocity" begin
    @testset "Small droplet (d <= 134.43 μm)" begin
        radius = 50e-6  # 50 microns, diameter = 100 μm
        volume = FOUR_THIRD * PI * radius^3
        mass = volume * RHO_WATER
        
        tv = SuperdropletModel.compute_terminal_velocity(radius, mass)
        
        @test tv > 0
        @test tv < 10.0  # Should be reasonable for small droplets
    end
    
    @testset "Medium droplet (134.43 < d <= 1511.64 μm)" begin
        radius = 500e-6  # 500 microns, diameter = 1000 μm
        volume = FOUR_THIRD * PI * radius^3
        mass = volume * RHO_WATER
        
        tv = SuperdropletModel.compute_terminal_velocity(radius, mass)
        
        @test tv > 0
        @test tv < 10.0
    end
    
    @testset "Large droplet (1511.64 < d <= 3477.84 μm)" begin
        radius = 1500e-6  # 1500 microns, diameter = 3000 μm
        volume = FOUR_THIRD * PI * radius^3
        mass = volume * RHO_WATER
        
        tv = SuperdropletModel.compute_terminal_velocity(radius, mass)
        
        @test tv > 0
        @test tv < 10.0
    end
    
    @testset "Very large droplet (d > 3477.84 μm)" begin
        radius = 2000e-6  # 2000 microns, diameter = 4000 μm
        volume = FOUR_THIRD * PI * radius^3
        mass = volume * RHO_WATER
        
        tv = SuperdropletModel.compute_terminal_velocity(radius, mass)
        
        @test tv > 0
        @test tv < 10.0
    end
    
    @testset "Increases with size" begin
        radii = [10e-6, 50e-6, 100e-6, 500e-6]
        terminal_velocities = Float64[]
        
        for radius in radii
            volume = FOUR_THIRD * PI * radius^3
            mass = volume * RHO_WATER
            tv = SuperdropletModel.compute_terminal_velocity(radius, mass)
            push!(terminal_velocities, tv)
        end
        
        # Check that terminal velocities generally increase
        for i in 1:length(terminal_velocities)-1
            @test terminal_velocities[i+1] >= terminal_velocities[i]
        end
    end
end

@testset "Update rcubed" begin
    @testset "Basic update" begin
        droplet = Droplet(1000, 10e-6)
        old_mass = droplet.mass
        old_radius = droplet.radius
        
        new_rcubed = (20e-6)^3
        update_rcubed!(droplet, new_rcubed)
        
        @test droplet.rcubed ≈ new_rcubed
        @test droplet.radius ≈ 20e-6 atol=1e-10
        @test droplet.mass > old_mass
    end
    
    @testset "Properties recalculated" begin
        droplet = Droplet(1000, 10e-6)
        new_rcubed = (15e-6)^3
        
        update_rcubed!(droplet, new_rcubed)
        
        expected_volume = new_rcubed * FOUR_THIRD * PI
        expected_mass = expected_volume * RHO_WATER
        
        @test droplet.volume ≈ expected_volume
        @test droplet.mass ≈ expected_mass
    end
end

@testset "Total Water" begin
    @testset "Empty vector" begin
        total = total_water(Droplet[])
        @test total == 0.0
    end
    
    @testset "Single droplet" begin
        droplet = Droplet(1000, 10e-6)
        total = total_water([droplet])
        expected = droplet.mass * droplet.multi
        @test total ≈ expected
    end
    
    @testset "Multiple droplets" begin
        droplets = [
            Droplet(1000, 10e-6),
            Droplet(500, 50e-6),
            Droplet(100, 100e-6)
        ]
        total = total_water(droplets)
        expected = sum(d.mass * d.multi for d in droplets)
        @test total ≈ expected
    end
    
    @testset "High multiplicity" begin
        droplets = [
            Droplet(10000, 5e-6),
            Droplet(5000, 20e-6),
            Droplet(1000, 100e-6)
        ]
        total = total_water(droplets)
        expected = sum(d.mass * d.multi for d in droplets)
        @test total ≈ expected
    end
end
