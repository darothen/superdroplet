"""Tests for constants module."""

using Test

@testset "Mathematical Constants" begin
    @test PI ≈ π
    @test THIRD ≈ 1.0/3.0
    @test THREE_FOURTH ≈ 3.0/4.0
    @test FOUR_THIRD ≈ 4.0/3.0
end

@testset "Physical Constants" begin
    @test RHO_WATER == 1.0e3
    @test RHO_WATER > 0
end
