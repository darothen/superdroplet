module constants

    use iso_fortran_env
    implicit none

    integer, parameter :: rkind = real64
    integer, parameter :: ikind = int32
    integer, parameter :: lkind = int64

    real(kind=rkind), parameter :: PI = 4d0 * atan(1.0d0)
    real(kind=rkind), parameter :: RHO_WATER = 1e3 ! kg/m3
    real(kind=rkind), parameter :: RHO_AIR = 1.0
    real(kind=rkind), parameter :: THIRD = 1./3.
    real(kind=rkind), parameter :: MULTI_THRESH = 1e4

end module
