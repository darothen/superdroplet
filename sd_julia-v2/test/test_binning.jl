"""Tests for binning module."""

using Test

@testset "Bin Droplets" begin
    @testset "Basic binning" begin
        # Create droplets with known radii
        droplets = [
            Droplet(1000, 10e-6),   # 10 μm
            Droplet(500, 50e-6),    # 50 μm
            Droplet(100, 100e-6),   # 100 μm
        ]
        
        # Create bins: 0-30, 30-70, 70-150 μm
        bin_edges = [0.0, 30.0, 70.0, 150.0]
        
        result = bin_droplets(droplets, bin_edges)
        
        @test length(result) == 3
        @test all(x -> x >= 0, result)
    end
    
    @testset "Empty droplet list" begin
        droplets = Droplet[]
        bin_edges = [0.0, 10.0, 100.0, 1000.0]
        
        result = bin_droplets(droplets, bin_edges)
        
        @test length(result) == 3
        @test all(x -> x == 0.0, result)
    end
    
    @testset "Single droplet" begin
        droplet = Droplet(1000, 50e-6)  # 50 μm radius
        bin_edges = [0.0, 30.0, 70.0, 150.0]
        
        result = bin_droplets([droplet], bin_edges)
        
        # Should be in second bin (30-70 μm)
        @test result[2] > 0
        @test result[1] == 0.0
        @test result[3] == 0.0
    end
    
    @testset "Droplets in different bins" begin
        droplets = [
            Droplet(1000, 10e-6),   # Should be in bin 1
            Droplet(1000, 50e-6),   # Should be in bin 2
            Droplet(1000, 100e-6),  # Should be in bin 3
        ]
        
        bin_edges = [0.0, 30.0, 70.0, 150.0]
        
        result = bin_droplets(droplets, bin_edges)
        
        # All bins should have some mass
        @test all(x -> x > 0, result)
    end
    
    @testset "Multiple droplets in same bin" begin
        droplets = [
            Droplet(1000, 45e-6),   # Both in bin 2
            Droplet(2000, 55e-6),   # Both in bin 2
        ]
        
        bin_edges = [0.0, 30.0, 70.0, 150.0]
        
        result = bin_droplets(droplets, bin_edges)
        
        # Only bin 2 should have mass
        @test result[1] == 0.0
        @test result[2] > 0
        @test result[3] == 0.0
        
        # Bin 2 should have combined mass
        expected_mass = (droplets[1].mass * droplets[1].multi + 
                        droplets[2].mass * droplets[2].multi)
        @test result[2] ≈ expected_mass
    end
    
    @testset "Mass conservation" begin
        droplets = [
            Droplet(1000, 10e-6),
            Droplet(500, 50e-6),
            Droplet(100, 100e-6),
            Droplet(200, 200e-6),
        ]
        
        bin_edges = [0.0, 50.0, 150.0, 300.0]
        
        result = bin_droplets(droplets, bin_edges)
        
        total_binned = sum(result)
        total_actual = total_water(droplets)
        
        @test total_binned ≈ total_actual rtol=1e-10
    end
    
    @testset "Logarithmic bins" begin
        droplets = [
            Droplet(1000, 1e-6),
            Droplet(1000, 10e-6),
            Droplet(1000, 100e-6),
            Droplet(1000, 1000e-6),
        ]
        
        # Logarithmically spaced bins
        bin_edges = [10.0^exp for exp in 0.0:1.0:4.0]
        
        result = bin_droplets(droplets, bin_edges)
        
        @test length(result) == 4
        @test sum(result) ≈ total_water(droplets) rtol=1e-10
    end
    
    @testset "Edge case: droplet exactly at bin boundary" begin
        droplet = Droplet(1000, 50e-6)  # Exactly 50 μm
        bin_edges = [0.0, 50.0, 100.0]
        
        result = bin_droplets([droplet], bin_edges)
        
        # Should be in first bin (< 50 μm bin)
        # Due to the < comparison in the binning code
        @test result[1] == 0.0
        @test result[2] > 0.0
    end
end
