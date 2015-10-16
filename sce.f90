
program sce

    use common
    use droplet_class, only: droplet, droplet_, droplet_count

    use fgsl

    implicit none

    real(kind=rkind), parameter ::    &
        n_0 = 2**23,                  &
        R_0 = 30.531d-6,              &
        X_0 = (4d0*PI/3d0)*(R_0**3.), &
        M_0 = X_0*RHO_WATER,          &
        delta_V = 1d6,                &
        t_c = 1.0

    integer(kind=ikind), parameter :: &
        n_part = 2**11,               &
        t_end = 3601,                 &
        plot_dt = 1200

    type(droplet), dimension(n_part) :: droplets
    real(kind=rkind), dimension(n_part) :: x_grid, r_grid, m_grid


    real(kind=rkind) :: total_droplets
    real(kind=rkind) :: rcubed_i ! Placeholder for cubed radius of a droplet
    type(droplet) :: droplet_i
    integer(kind=lkind) :: xi_i

    type(fgsl_rng) :: r      ! FGSL random number generator
    type(fgsl_rng_type) :: t ! Type of FGSL random number generaotr
    real(fgsl_double) :: phi ! Holder for a sample from the RNG

    integer(kind=lkind) :: N_per_SD

    integer :: i
! -----------------------------------------------------------------------------

    print *, "BEGINNING MAIN ROUTINE"

    ! Initialize the droplet array; note that we use a statically sized
    ! array since we know the number of droplets at compile time.

    ! Pre-compute the superdroplet multiplicity
    total_droplets = delta_V * n_0
    xi_i = floor(total_droplets / n_part, lkind)

    ! 1) Initialize the random number generator from FGSL

    t = fgsl_rng_env_setup()
    t = fgsl_rng_default
    r = fgsl_rng_alloc(t)

    ! Call the random number generator to compute masses.
    do i = 1, n_part
        ! The exponential dist is of the form
        !             1
        ! p(x) dx = ---- exp ( -x / mu ) dx
        !            mu
        ! where mu = 1/lambda from the Boost version, but is the same
        ! as the scipy form of the distribution!=
        phi = fgsl_ran_exponential(r, X_0) ! phi is a sample from the mass
                                           ! distribution
        x_grid(i) = phi
    end do

    ! Pre-sort the x_grid because it's a bit easier to sort doubles
    ! with FGSL than objects
    call fgsl_sort(x_grid, 1_fgsl_size_t, &
                           int(n_part, kind=fgsl_size_t))

    r_grid = (x_grid * 3. / PI / 4.)**THIRD

    ! Populate the initial droplet array
    do i = 1, n_part
        droplets(i) = droplet_(xi_i, r_grid(i)**3.)
        m_grid(i) = x_grid(i) * droplets(i)%density
    end do

    print *, "GRID SETUP"
    print *, "   radii: ", droplets(1)%radius, " - ", &
                           droplets(n_part-1)%radius, " m"
    print *, "  volume: ", droplets(1)%volume, " - ", &
                           droplets(n_part-1)%volume, " m^3"
    print *, "    mass: ", droplets(1)%mass, " - ", &
                           droplets(n_part-1)%mass, " kg"

    print *, "SD SETUP"
    print *, "   N_s: ", n_part
    print *, "  xi_i: ", xi_i
    N_per_SD = total_droplets / xi_i / n_part
    print *, " N per SD_xi":, N_per_SD
    print *, "(Initialized ", droplet_count, " droplets)"


    ! ) Free the memory from the random number generator
    call fgsl_rng_free(r)


end program sce