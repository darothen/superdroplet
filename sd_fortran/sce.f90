! ---------------------------------------------------------------
! Main stochastic collision/coalescence program, using the
! superdroplet method of Shima et al (2009)
! ---------------------------------------------------------------
program sce

    use collisions, only: collision_step, kernels, GOLOVIN, HYDRO, LONG, HALL
    use constants
    use droplet_class, only: droplet, droplet_, droplet_count
    use util

    use fgsl

    implicit none

    logical :: DEBUG = .false.

    !! Configuration parameters for the initial droplet size distribution
    real(kind=rkind), parameter ::    &
        ! Total number concentration
        ! n_0 = 2.**23,                 &
        n_0 = 27d0*(2d0**23),         &
        ! Total droplet radius
        ! R_0 = 30.531d-6,              &
        R_0 = 30.531e-6/3d0,              &
        ! Total droplet volume
        X_0 = (4.*PI/3.)*(R_0**3.),   &
        ! Total droplet water mass
        M_0 = X_0*RHO_WATER,          &
        ! Total parcel volume
        delta_V = 1e6,                &
        ! Model timestep (seconds)
        t_c = 1.0
        ! t_c = 0.1

    integer(kind=ikind), parameter :: &
        ! Total number of superdroplets
        n_part = 2**13,               &
        ! Total simulation time (seconds)
        t_end = 3601,                 &
        ! Output interval time (seconds)
        plot_dt = 1200

    integer(kind(kernels)), parameter :: &
        ! Collision kernel enumeration
        kern = GOLOVIN

    ! -- Nothing needs to be configured past here
    type(droplet), dimension(n_part) :: droplets
    real(kind=rkind), dimension(n_part) :: x_grid, r_grid, m_grid

    real(kind=rkind) :: total_droplets
    real(kind=rkind) :: rcubed_i ! Placeholder for cubed radius of a droplet
    type(droplet) :: droplet_i
    integer(kind=lkind) :: xi_i
    real(fgsl_double) :: phi     ! Holder for a sample from the RNG

    real(kind=rkind) :: N_per_SD, wm0

    integer :: i, ti
    real :: t
    type(time_t) :: timestamp

! -----------------------------------------------------------------------------

    ! Pre-initialize the random number generator, since it will be
    ! used both here and elsewhere in the program.
    rng_t = fgsl_rng_env_setup()
    rng_t = fgsl_rng_default
    rng   = fgsl_rng_alloc(rng_t)

    ! Initialize the droplet array; note that we use a statically sized
    ! array since we know the number of droplets at compile time.

    ! Pre-compute the superdroplet multiplicity
    total_droplets = delta_V * n_0
    xi_i = floor(total_droplets / n_part, lkind)

    ! Call the random number generator to compute masses.
    sample_loop : do i = 1, n_part
        ! The exponential dist is of the form
        !             1
        ! p(x) dx = ---- exp ( -x / mu ) dx
        !            mu
        ! where mu = 1/lambda from the Boost version, but is the same
        ! as the scipy form of the distribution!=
        phi = fgsl_ran_exponential(rng, X_0) ! phi is a sample from the mass
                                             ! distribution
        x_grid(i) = phi
    end do sample_loop

    ! Pre-sort the x_grid because it's a bit easier to sort doubles
    ! with FGSL than objects
    call fgsl_sort(x_grid, 1_fgsl_size_t, &
                           int(n_part, kind=fgsl_size_t))

    r_grid = (x_grid * 3. / PI / 4.)**THIRD

    ! Populate the initial droplet array
    populate_loop : do i = 1, n_part
        droplets(i) = droplet_(xi_i, r_grid(i)**3.)
        m_grid(i) = x_grid(i) * droplets(i)%density
    end do populate_loop

    write (*,*) "GRID SETUP"
    write (*,*) "   radii: ", droplets(1)%radius, " - ", &
                           droplets(n_part-1)%radius, " m"
    write (*,*) "  volume: ", droplets(1)%volume(), " - ", &
                           droplets(n_part-1)%volume(), " m^3"
    write (*,*) "    mass: ", droplets(1)%mass(), " - ", &
                           droplets(n_part-1)%mass(), " kg"

    write (*,*) "SD SETUP"
    write (*,*) "   N_s: ", n_part
    write (*,*) "  xi_i: ", xi_i
    N_per_SD = total_droplets / xi_i / n_part
    write (*,*) " N per SD_xi:", N_per_SD
    write (*,*) "(Initialized ", droplet_count, " droplets)"

    write (*,*) "BEGINNING MAIN ROUTINE"

    ! TODO: This is off by one order of magnitude, but it's not clear
    ! where the problem is coming from.
    wm0 = total_water(droplets)
    write (*,*) "Initial water mass = ", wm0, "kg"

    t = 0.
    ti = 0

    main_loop : do while (t < t_end)

        timestamp = s_to_min_s(t)

        write (*,'(/, A, I5," (", I2, " min")', advance='no') &
            "STEP", ti, timestamp%minutes
        if ( timestamp%seconds > 0 ) &
            write (*,'(1x,F5.2," sec")', advance='no') &
                timestamp%seconds
        write (*,'(")")')

        if (      ( floor( t / real(plot_dt) ) == ( t / plot_dt ) ) &
             .or. ( ti == 0 ) ) &
            call output_droplets(droplets, t)

        ! -------------------------------------------------------------
        call collision_step(droplets, t_c, delta_v, kern)
        ! -------------------------------------------------------------

        ti = ti + 1
        t = ti * t_c

        if ( DEBUG ) exit

    end do main_loop

    ! ) Free the memory from the random number generator
    call fgsl_rng_free(rng)

end program sce
