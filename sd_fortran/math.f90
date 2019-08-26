module math

    use constants
    implicit none

contains

    real(kind=rkind) function exp_dist_moments(x, n0, x0, l)
        real(kind=rkind), intent(in) :: x, n0, x0
        integer(kind=ikind), optional :: l

        if ( .not. present(l) ) l = 0

        exp_dist_moments = (x**l)*(n0/x0)*exp(-x/x0)

    end function exp_dist_moments

    real(kind=rkind) function expon_pdf(x)
        real(kind=rkind), intent(in) :: x
        real(kind=rkind) :: lambda, scale

        scale = X
        lambda = 1.d0 / scale
        expon_pdf = lambda * exp(-lambda * x)

    end function expon_pdf

end module
