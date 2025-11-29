!> Time tracking module
!!
!! This module provides a simple stopwatch type for tracking simulation time.
!! Uses object-oriented programming with type-bound procedures.
!!
!! @author Daniel Rothenberg
!! @date November 2025
module time

    implicit none
    private

    public :: stopwatch_t

    !> Stopwatch type for tracking elapsed time
    !!
    !! Tracks time in minutes and seconds with methods to increment
    !! and query total elapsed time.
    type :: stopwatch_t
        integer :: minutes   !< Elapsed minutes
        integer :: seconds   !< Elapsed seconds (0-59)
    contains
        procedure, public :: increment => stopwatch_increment       !< Increment by seconds
        procedure, public :: total_seconds => stopwatch_total_seconds !< Get total seconds
    end type stopwatch_t

    interface stopwatch_t
        module procedure init_stopwatch
    end interface

contains

    !> Constructor for stopwatch type
    !!
    !! Initializes a stopwatch with a given number of seconds,
    !! automatically converting to minutes and seconds format.
    !!
    !! @param[in] seconds - Initial elapsed time in seconds
    !! @return Initialized stopwatch_t instance
    type(stopwatch_t) function init_stopwatch(seconds)

        integer, intent(in) :: seconds
        integer :: seconds_over

        seconds_over = mod(seconds, 60)
        init_stopwatch%minutes = (seconds - seconds_over) / 60
        init_stopwatch%seconds = seconds_over

    end function init_stopwatch

    !> Increment stopwatch by seconds
    !!
    !! Adds delta_seconds to the stopwatch, automatically handling
    !! minute/second overflow.
    !!
    !! @param[inout] self - Stopwatch instance to increment
    !! @param[in] delta_seconds - Number of seconds to add
    subroutine stopwatch_increment(self, delta_seconds)
        class(stopwatch_t), intent(inout) :: self
        integer, intent(in) :: delta_seconds
        integer :: total

        total = self%minutes * 60 + self%seconds + delta_seconds
        self%minutes = total / 60
        self%seconds = mod(total, 60)
    end subroutine stopwatch_increment

    !> Get total elapsed time in seconds
    !!
    !! Returns the total elapsed time as a single integer value in seconds.
    !!
    !! @param[in] self - Stopwatch instance to query
    !! @return Total elapsed time in seconds
    integer pure function stopwatch_total_seconds(self)
        class(stopwatch_t), intent(in) :: self

        stopwatch_total_seconds = self%minutes * 60 + self%seconds
    end function stopwatch_total_seconds

end module