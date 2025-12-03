"""Tests for time module."""

using Test

@testset "Stopwatch Creation" begin
    @testset "Zero seconds" begin
        sw = Stopwatch(0)
        @test sw.minutes == 0
        @test sw.seconds == 0
    end
    
    @testset "Less than minute" begin
        sw = Stopwatch(45)
        @test sw.minutes == 0
        @test sw.seconds == 45
    end
    
    @testset "Exactly one minute" begin
        sw = Stopwatch(60)
        @test sw.minutes == 1
        @test sw.seconds == 0
    end
    
    @testset "Minutes and seconds" begin
        sw = Stopwatch(90)
        @test sw.minutes == 1
        @test sw.seconds == 30
    end
    
    @testset "Multiple minutes" begin
        sw = Stopwatch(185)  # 3 minutes and 5 seconds
        @test sw.minutes == 3
        @test sw.seconds == 5
    end
    
    @testset "Large value" begin
        sw = Stopwatch(3661)  # 61 minutes and 1 second
        @test sw.minutes == 61
        @test sw.seconds == 1
    end
end

@testset "Total Seconds" begin
    @testset "Zero time" begin
        sw = Stopwatch(0)
        @test total_seconds(sw) == 0
    end
    
    @testset "Seconds only" begin
        sw = Stopwatch(45)
        @test total_seconds(sw) == 45
    end
    
    @testset "Minutes and seconds" begin
        sw = Stopwatch(90)
        @test total_seconds(sw) == 90
    end
    
    @testset "Large time" begin
        sw = Stopwatch(3661)
        @test total_seconds(sw) == 3661
    end
end

@testset "Stopwatch Increment" begin
    @testset "Basic increment" begin
        sw = Stopwatch(0)
        increment!(sw, 10)
        @test sw.seconds == 10
        @test sw.minutes == 0
    end
    
    @testset "No overflow" begin
        sw = Stopwatch(30)
        increment!(sw, 15)
        @test sw.seconds == 45
        @test sw.minutes == 0
    end
    
    @testset "With overflow" begin
        sw = Stopwatch(50)
        increment!(sw, 20)
        @test sw.seconds == 10
        @test sw.minutes == 1
    end
    
    @testset "Exactly to minute" begin
        sw = Stopwatch(40)
        increment!(sw, 20)
        @test sw.seconds == 0
        @test sw.minutes == 1
    end
    
    @testset "Multiple minutes" begin
        sw = Stopwatch(30)
        increment!(sw, 150)  # 2 minutes and 30 seconds
        @test sw.seconds == 0
        @test sw.minutes == 3
    end
    
    @testset "From existing minutes" begin
        sw = Stopwatch(120)  # 2 minutes
        increment!(sw, 45)
        @test sw.seconds == 45
        @test sw.minutes == 2
    end
    
    @testset "Multiple increments" begin
        sw = Stopwatch(0)
        increment!(sw, 30)
        increment!(sw, 30)
        increment!(sw, 30)
        @test sw.seconds == 30
        @test sw.minutes == 1
    end
    
    @testset "Large increment" begin
        sw = Stopwatch(0)
        increment!(sw, 7325)  # 122 minutes and 5 seconds
        @test sw.seconds == 5
        @test sw.minutes == 122
    end
end

@testset "Stopwatch Round Trip" begin
    # Test that total_seconds is inverse of constructor
    for test_seconds in [0, 45, 60, 90, 185, 3600, 3661]
        sw = Stopwatch(test_seconds)
        @test total_seconds(sw) == test_seconds
    end
end
