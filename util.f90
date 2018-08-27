module util

    use constants
    use droplet_class, only: droplet

    use fgsl

    implicit none

    type time_t
        integer :: minutes
        real :: seconds
    end type time_t

    integer, parameter :: file_unit = 10
    character (len=*), parameter :: suffix = "_output.txt"

    type(fgsl_rng) :: rng        ! FGSL random number generator
    type(fgsl_rng_type) :: rng_t ! Type of FGSL random number generaotr

contains

    function ifloor(x)
        real(kind=rkind), intent(in) :: x
        integer(kind=ikind) :: ifloor

        ifloor = floor(x, kind=ikind)

    end function ifloor

    subroutine shuffle_droplet_array(droplets)
        ! Implementation of Knuth Shuffle
        ! (http://rosettacode.org/wiki/Knuth_shuffle)
        type(droplet), dimension(:), intent(inout) :: droplets

        integer :: i, randpos
        type(droplet) :: temp
        real :: r

        do i = size(droplets), 2, -1
            call random_number(r)
            randpos = int(r*i) + 1
            temp = droplets(randpos)
            droplets(randpos) = droplets(i)
            droplets(i) = temp
        end do

    end subroutine shuffle_droplet_array

    real(kind=rkind) function total_water(droplets)
        type(droplet), dimension(:), intent(in) :: droplets

        integer(kind=ikind) :: n_droplets
        real(kind=rkind) :: droplet_water

        integer :: i

        total_water = 0.0
        n_droplets = size(droplets, kind=ikind)

        do i = 1, n_droplets
            droplet_water = droplets(i)%mass * droplets(i)%multi / 1d3
            total_water = total_water + droplet_water
        end do

    end function total_water

    type(time_t) function s_to_min_s(t) result(timestamp)
        real, intent(in) :: t ! time in seconds

        real :: seconds_over

        seconds_over = mod(t, 60.)
        timestamp%minutes = floor((t - seconds_over)/60.)
        timestamp%seconds = seconds_over

    end function s_to_min_s

    subroutine output_droplets(droplets, t)
        type(droplet), dimension(:), intent(in) :: droplets
        real, intent(in) :: t

        integer i
        integer(kind=ikind) :: n_droplets
        type(droplet) :: droplet_i
        character (len=1024) :: output_fn

        write (output_fn, "(F6.1,A)") t, suffix

        write (*,"(A,A,A)") "Writing output... (", trim(output_fn), ")"

        ! TODO: Safely open file (use iostatus)
        open (unit=file_unit, file=trim(adjustl(output_fn)))

        n_droplets = size(droplets, kind=ikind)
        do i = 1, n_droplets
            droplet_i = droplets(i)
            write (file_unit, "(ES16.5,',',I16)") &
                droplet_i%rcubed, droplet_i%multi
        end do

        close (file_unit)

        write (*,*) " done."

    end subroutine output_droplets

    real function urand()
        ! urand = fgsl_rng_uniform(rng)
        call random_number(urand)
    end function urand

end module
