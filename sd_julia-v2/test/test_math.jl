"""Tests for math utilities module."""

using Test
using Random

@testset "Knuth Shuffle" begin
    @testset "Basic shuffle" begin
        Random.seed!(42)
        v = collect(1:10)
        original = copy(v)
        
        SuperdropletModel.knuth_shuffle!(v)
        
        # Should have same elements
        @test sort(v) == sort(original)
        # But different order (with high probability)
        @test v != original  # This is probabilistic but very likely
    end
    
    @testset "Single element" begin
        v = [1]
        SuperdropletModel.knuth_shuffle!(v)
        @test v == [1]
    end
    
    @testset "Two elements" begin
        Random.seed!(123)
        v = [1, 2]
        SuperdropletModel.knuth_shuffle!(v)
        # Should have same elements
        @test sort(v) == [1, 2]
    end
    
    @testset "Preserves all elements" begin
        Random.seed!(456)
        v = collect(1:100)
        original = copy(v)
        
        SuperdropletModel.knuth_shuffle!(v)
        
        @test sort(v) == sort(original)
        @test length(v) == length(original)
    end
end

@testset "Median" begin
    @testset "Single value" begin
        @test SuperdropletModel.median([5.0]) == 5.0
    end
    
    @testset "Two values" begin
        @test SuperdropletModel.median([3.0, 7.0]) == 5.0
    end
    
    @testset "Odd length" begin
        @test SuperdropletModel.median([1.0, 2.0, 3.0, 4.0, 5.0]) == 3.0
    end
    
    @testset "Even length" begin
        @test SuperdropletModel.median([1.0, 2.0, 3.0, 4.0]) == 2.5
    end
    
    @testset "Unsorted input" begin
        @test SuperdropletModel.median([5.0, 1.0, 3.0, 2.0, 4.0]) == 3.0
    end
    
    @testset "With duplicates" begin
        @test SuperdropletModel.median([1.0, 2.0, 2.0, 2.0, 3.0]) == 2.0
    end
    
    @testset "Negative values" begin
        @test SuperdropletModel.median([-5.0, -2.0, 0.0, 2.0, 5.0]) == 0.0
    end
end

@testset "Rolling Median" begin
    @testset "Basic rolling median" begin
        values = [1.0, 2.0, 3.0, 4.0, 5.0]
        window = 3
        
        result = rolling_median(values, window)
        
        @test length(result) == length(values)
        @test all(x -> !isnan(x), result)
    end
    
    @testset "Window size 1" begin
        values = [1.0, 2.0, 3.0, 4.0, 5.0]
        result = rolling_median(values, 1)
        
        # Window size 1 should be close to original (with small edge effects)
        @test length(result) == length(values)
    end
    
    @testset "Large window" begin
        values = collect(1.0:10.0)
        window = 9
        
        result = rolling_median(values, window)
        
        @test length(result) == length(values)
    end
    
    @testset "Constant values" begin
        values = fill(5.0, 10)
        window = 3
        
        result = rolling_median(values, window)
        
        @test all(x -> x ≈ 5.0, result)
    end
    
    @testset "Single value" begin
        values = [5.0]
        result = rolling_median(values, 3)
        
        @test length(result) == 1
        @test result[1] ≈ 5.0
    end
end
