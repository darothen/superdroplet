module common

    integer, parameter :: rkind = kind(1.d0)
    integer, parameter :: ikind = selected_int_kind(8)
    integer, parameter :: lkind = selected_int_kind(16)

    real(kind=rkind), parameter :: PI = 4d0 * atan(1.0d0)
    real(kind=rkind), parameter :: RHO_WATER = 1e3 ! kg/m3
    real(kind=rkind), parameter :: RHO_AIR = 1.0
    real(kind=rkind), parameter :: THIRD = 1./3.
    real(kind=rkind), parameter :: MULTI_THRESH = 1e4

contains

    function ifloor(x)
        real(kind=rkind), intent(in) :: x
        integer(kind=ikind) :: ifloor

        ifloor = floor(x, ikind)

    end function ifloor

end module