!> Physical and mathematical constants module
!!
!! This module defines precision parameters and commonly used physical
!! and mathematical constants for the superdroplet model.
!!
!! @author Daniel Rothenberg
!! @date November 2025
module constants
    use, intrinsic :: iso_fortran_env, only : real64

    implicit none

    private

    !> Double precision kind parameter (64-bit floating point)
    integer, parameter, public :: dp = real64

    !> Mathematical constant π
    !! Computed at compile time as 4*arctan(1)
    real(kind=dp), parameter, public :: PI = 4.0_dp * atan(1.0_dp)

    !> Density of water [kg/m³]
    real(kind=dp), parameter, public :: RHO_WATER = 1000.0_dp

    !> Mathematical constant 1/3
    real(kind=dp), parameter, public :: THIRD = 1.0_dp / 3.0_dp

    !> Mathematical constant 3/4
    real(kind=dp), parameter, public :: THREE_FOURTH = 3.0_dp / 4.0_dp

    !> Mathematical constant 4/3
    real(kind=dp), parameter, public :: FOUR_THIRD = 4.0_dp / 3.0_dp

    contains

end module constants