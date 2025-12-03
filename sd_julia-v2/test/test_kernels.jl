"""Tests for collision kernels module."""

using Test

@testset "Kernel Enum" begin
    @test Golovin isa Kernel
    @test Hydro isa Kernel
    @test Long isa Kernel
end

@testset "Golovin Kernel" begin
    @testset "Constant value" begin
        expected = 1500.0 * FOUR_THIRD * PI
        @test SuperdropletModel.GOLOVIN_CONSTANT ≈ expected
    end
    
    @testset "Basic computation" begin
        small = Droplet(1000, 10e-6)
        medium = Droplet(500, 50e-6)
        
        result = SuperdropletModel.golovin_kernel(small, medium)
        expected = SuperdropletModel.GOLOVIN_CONSTANT * (small.rcubed + medium.rcubed)
        
        @test result ≈ expected
    end
    
    @testset "Symmetry" begin
        small = Droplet(1000, 10e-6)
        large = Droplet(100, 100e-6)
        
        result1 = SuperdropletModel.golovin_kernel(small, large)
        result2 = SuperdropletModel.golovin_kernel(large, small)
        
        @test result1 ≈ result2
    end
    
    @testset "Positive values" begin
        small = Droplet(1000, 10e-6)
        medium = Droplet(500, 50e-6)
        
        result = SuperdropletModel.golovin_kernel(small, medium)
        @test result > 0
    end
    
    @testset "Increases with size" begin
        d1 = Droplet(1, 10e-6)
        d2 = Droplet(1, 50e-6)
        d3 = Droplet(1, 100e-6)
        
        k12 = SuperdropletModel.golovin_kernel(d1, d2)
        k23 = SuperdropletModel.golovin_kernel(d2, d3)
        
        @test k23 > k12
    end
end

@testset "Calc Hydro Kernel" begin
    @testset "Basic calculation" begin
        e_coal = 1.0
        e_coll = 1.0
        r_sum = 100e-6  # 100 microns
        tv_diff = 1.0  # 1 m/s
        
        result = SuperdropletModel.calc_hydro_kernel(e_coal, e_coll, r_sum, tv_diff)
        expected = e_coal * e_coll * PI * r_sum * r_sum * abs(tv_diff)
        
        @test result ≈ expected
    end
    
    @testset "Negative velocity difference" begin
        e_coal = 1.0
        e_coll = 1.0
        r_sum = 100e-6
        tv_diff = -1.0
        
        result = SuperdropletModel.calc_hydro_kernel(e_coal, e_coll, r_sum, tv_diff)
        @test result > 0  # Should be positive due to abs()
    end
    
    @testset "Zero velocity difference" begin
        e_coal = 1.0
        e_coll = 1.0
        r_sum = 100e-6
        tv_diff = 0.0
        
        result = SuperdropletModel.calc_hydro_kernel(e_coal, e_coll, r_sum, tv_diff)
        @test result ≈ 0.0
    end
    
    @testset "Efficiency factors" begin
        r_sum = 100e-6
        tv_diff = 1.0
        
        result_full = SuperdropletModel.calc_hydro_kernel(1.0, 1.0, r_sum, tv_diff)
        result_half = SuperdropletModel.calc_hydro_kernel(0.5, 1.0, r_sum, tv_diff)
        
        @test result_half ≈ result_full * 0.5
    end
end

@testset "Hydro Kernel" begin
    @testset "Basic computation" begin
        small = Droplet(1000, 10e-6)
        large = Droplet(100, 100e-6)
        
        result = SuperdropletModel.hydro_kernel(small, large)
        @test result > 0
    end
    
    @testset "Symmetry" begin
        small = Droplet(1000, 10e-6)
        large = Droplet(100, 100e-6)
        
        result1 = SuperdropletModel.hydro_kernel(small, large)
        result2 = SuperdropletModel.hydro_kernel(large, small)
        
        @test result1 ≈ result2
    end
    
    @testset "Same size droplets" begin
        d1 = Droplet(1, 50e-6)
        d2 = Droplet(1, 50e-6)
        
        result = SuperdropletModel.hydro_kernel(d1, d2)
        # Same size droplets should have same terminal velocity, so kernel ~0
        @test result ≈ 0.0 atol=1e-10
    end
    
    @testset "Increases with size difference" begin
        d_small = Droplet(1, 10e-6)
        d_medium = Droplet(1, 50e-6)
        d_large = Droplet(1, 200e-6)
        
        k_small_medium = SuperdropletModel.hydro_kernel(d_small, d_medium)
        k_small_large = SuperdropletModel.hydro_kernel(d_small, d_large)
        
        @test k_small_large > k_small_medium
    end
end

@testset "Long Kernel" begin
    @testset "Basic computation" begin
        small = Droplet(1000, 10e-6)
        large = Droplet(100, 100e-6)
        
        result = SuperdropletModel.long_kernel(small, large)
        @test result > 0
    end
    
    @testset "Symmetry" begin
        small = Droplet(1000, 10e-6)
        large = Droplet(100, 100e-6)
        
        result1 = SuperdropletModel.long_kernel(small, large)
        result2 = SuperdropletModel.long_kernel(large, small)
        
        @test result1 ≈ result2
    end
    
    @testset "Same size droplets" begin
        d1 = Droplet(1, 50e-6)
        d2 = Droplet(1, 50e-6)
        
        result = SuperdropletModel.long_kernel(d1, d2)
        @test result ≈ 0.0 atol=1e-10
    end
    
    @testset "Collection efficiency for large droplets" begin
        d_small = Droplet(1, 10e-6)
        d_large = Droplet(1, 60e-6)  # > 50 μm
        
        result = SuperdropletModel.long_kernel(d_small, d_large)
        @test result > 0
        
        # Verify it's less than or equal to hydro kernel
        hydro_result = SuperdropletModel.hydro_kernel(d_small, d_large)
        @test result <= hydro_result
    end
    
    @testset "Collection efficiency bounds" begin
        radii = [5e-6, 10e-6, 20e-6, 30e-6, 40e-6]
        
        for r1 in radii
            for r2 in radii
                if r1 == r2
                    continue
                end
                d1 = Droplet(1, r1)
                d2 = Droplet(1, r2)
                
                result = SuperdropletModel.long_kernel(d1, d2)
                @test result >= 0
            end
        end
    end
    
    @testset "Long kernel vs hydro kernel" begin
        radii = [10e-6, 20e-6, 30e-6, 40e-6, 60e-6, 100e-6]
        
        for r1 in radii
            for r2 in radii
                if r1 == r2
                    continue
                end
                d1 = Droplet(1, r1)
                d2 = Droplet(1, r2)
                
                long_result = SuperdropletModel.long_kernel(d1, d2)
                hydro_result = SuperdropletModel.hydro_kernel(d1, d2)
                
                # Long kernel should be <= hydro kernel
                @test long_result <= hydro_result + 1e-15
            end
        end
    end
end

@testset "Compute Kernel Dispatch" begin
    @testset "Golovin dispatch" begin
        d1 = Droplet(1, 10e-6)
        d2 = Droplet(1, 50e-6)
        
        result = SuperdropletModel.compute_kernel(Golovin, d1, d2)
        expected = SuperdropletModel.golovin_kernel(d1, d2)
        
        @test result ≈ expected
    end
    
    @testset "Hydro dispatch" begin
        d1 = Droplet(1, 10e-6)
        d2 = Droplet(1, 50e-6)
        
        result = SuperdropletModel.compute_kernel(Hydro, d1, d2)
        expected = SuperdropletModel.hydro_kernel(d1, d2)
        
        @test result ≈ expected
    end
    
    @testset "Long dispatch" begin
        d1 = Droplet(1, 10e-6)
        d2 = Droplet(1, 50e-6)
        
        result = SuperdropletModel.compute_kernel(Long, d1, d2)
        expected = SuperdropletModel.long_kernel(d1, d2)
        
        @test result ≈ expected
    end
end
