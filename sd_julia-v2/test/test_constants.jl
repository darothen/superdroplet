"""Tests for constants module."""

using Test
using SuperdropletModel

@testset "Mathematical Constants" begin
    @test SuperdropletModel.PI ≈ π
    @test SuperdropletModel.THIRD ≈ 1.0/3.0
    @test SuperdropletModel.THREE_FOURTH ≈ 3.0/4.0
    @test SuperdropletModel.FOUR_THIRD ≈ 4.0/3.0
end

@testset "Physical Constants" begin
    @test SuperdropletModel.RHO_WATER == 1.0e3
    @test SuperdropletModel.RHO_WATER > 0
end
