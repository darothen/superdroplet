!> @file main.f90
!! @brief Main driver for superdroplet collision-coalescence simulation
!!
!! This program implements the stochastic collision-coalescence algorithm
!! using the superdroplet method of Shima et al. (2009). The simulation
!! models cloud droplet collisions and coalescence in a well-mixed parcel.
!!
!!
!! @author Daniel Rothenberg
!! @date November 2025
!!
!! @see Shima et al. (2009), "The super-droplet method for the numerical
!!      simulation of clouds and precipitation"
program main

    ! Module imports
    use collisions, only : collision_step, collision_step_result_t, model_config_t, model_config_t_
    ! TODO: Switch to stdlib_kinds::dp instead of defining our own?
    use constants, only : dp, PI, THIRD, THREE_FOURTH
    use droplet_mod, only : droplet_t, droplet_t_, total_water
    use kernels, only : GOLOVIN, HYDRO, LONG
    use time, only : stopwatch_t
    use util, only : rolling_window_median

    use stdlib_kinds, only : stdlib_dp => dp
    use stdlib_math, only : logspace
    use stdlib_random, only : random_seed
    use stdlib_sorting, only : sort, sort_index
    use stdlib_stats_distribution_exponential, only : rvs_exp

    use csv_module, only : csv_file

    implicit none

    ! Configuration parameters
    logical, parameter :: DEBUG = .true.

    real(kind=dp), parameter :: &
        n_0 = 27d0 * (2d0**23),         & ! Total number concentration
        R_0 = 30.531e-6 / 3.0,          & ! Total droplet radius
        X_0 = (4. * PI / 3.0_dp) * (R_0**3.0_dp), & ! Total droplet volume
        delta_V = 1e6_dp                  ! Total parcel volume
    
    integer, parameter :: &
        n_part = 2**17,            & ! Total number of superdroplets
        t_c = 1,                   & ! Model timestep (seconds)
        t_end = 3601,              & ! Total simulation time (seconds)
        plot_dt = 600,             & ! Output interval time (seconds)
        kernel = LONG             ! Collision kernel enumeration    

    ! Workspace variables
    real(kind=dp) :: &
        x, wm0, wm_i, total_droplets, rcubed_max
    integer :: &
        i, j, droplet_idx, xi_i, seed_put, seed_get, step
    real(kind=dp), dimension(n_part) :: &
        r_grid, x_grid, droplet_rcubed
    integer, dimension(n_part) :: sorted_indices
    type(droplet_t), dimension(n_part) :: droplets  ! Droplet array
    type(stopwatch_t) :: stopwatch
    type(model_config_t) :: model_config
    type(collision_step_result_t):: collision_step_result

    ! Output configuration
    type(csv_file) :: csv_out
    logical :: status_ok
    integer, parameter :: NR = 250, smooth_window = 9
    real(kind=dp), parameter :: r_max_bin = 5d3
    real(kind=dp), dimension(NR+1) :: r_bin_edges ! Full radius grid (bin edges)
    real(kind=dp), dimension(NR) :: &
        r_bins, &  ! Radius bin centers
        bin_counts_raw, &  ! Raw (unsmoothed) bin counts
        bin_counts_smooth  ! Smoothed bin counts

    ! Pre-construct the grid for computing the droplet size distribution to output
    r_bin_edges = logspace(&
        real(0, kind=stdlib_dp), &
        real(log10(r_max_bin), kind=stdlib_dp), &
        NR+1, &
        base=10)
    r_bins = 0.5_dp * (r_bin_edges(1:NR) + r_bin_edges(2:NR+1))

    ! Initialize model configuration
    model_config = model_config_t_( &
        step_seconds = t_c, &
        num_droplets = n_part, &
        delta_v = delta_V, &
        kernel = kernel )

    ! -- Begin main program --
    seed_put = 1234567
    ! call random_seed(seed_put, seed_get)  ! Initialize RNG
    call random_seed()  ! Initialize RNG with default seeding

    ! Initialize CSV output
    ! call f%initialize(verbose = .true. )
    call csv_out%initialize(verbose = .true.)
    call csv_out%open('collision_output.csv', NR+1, status_ok)
    if (.not. status_ok) then
        print *, "ERROR: Could not open output CSV file."
        stop
    end if
    ! No header - write the radius grid centers as the first row, and prepend a -9999 for the 
    ! time column
    call csv_out%add(-9999)
    call csv_out%add(r_bins, real_fmt="(ES12.5)")
    call csv_out%next_row()

    ! Initialize the droplet array;
    ! first, pre-compute the superdroplet multiplicity
    total_droplets = delta_V * n_0
    xi_i = floor(total_droplets / real(n_part, kind=dp))

    ! Compute the droplet masses according to an exponential distribution
    x_grid = rvs_exp(1._dp / X_0, n_part)
    call sort(x_grid)
    r_grid = (x_grid * THREE_FOURTH / PI)**THIRD
    init_droplets_loop : do i = 1, n_part
        droplets(i) = droplet_t_(xi_i, r_grid(i))
    end do init_droplets_loop

    print '(/,A)', "GRID SETUP"
    print '(A,ES12.3,A,ES12.3,A)', "   radii: ", droplets(1)%radius, " ", droplets(n_part)%radius, " m"
    print '(A,ES12.3,A,ES12.3,A)', "  volume: ", droplets(1)%volume, " ", droplets(n_part)%volume, " m^3"
    print '(A,ES12.3,A,ES12.3,A)', "    mass: ", droplets(1)%mass, " ", droplets(n_part)%mass, " kg"

    print '(/,A)', "SD SETUP"
    print '(A,I15)', "          N_s:", n_part
    print '(A,I15)', "         xi_i:", xi_i
    print '(A,I15)', "  N per SD_xi:", INT(total_droplets / xi_i / n_part)

    wm0 = total_water(droplets)
    print '(/,A,F12.3,A)', "Initial total water mass = ", wm0, " kg"

    print '(/,A,/)', "BEGINNING MAIN SIMULATION LOOP"
    stopwatch = stopwatch_t(0)
    step = 0
    main_loop : do while (stopwatch%total_seconds() <= t_end)

        call collision_step(droplets, collision_step_result, model_config)

        if ( mod(step, plot_dt) == 0 ) then
            print '(A,I8,A,I8,A)', "Plotting. Stopwatch = ", stopwatch%minutes, " min ", stopwatch%seconds, " s"

            ! Compute droplet size distribution into bins
            bin_counts_raw = 0.0_dp
            droplet_idx = 1
            droplet_rcubed = (/ (droplets(i)%rcubed, i=1,n_part) /)
            call sort_index(droplet_rcubed, sorted_indices)
            ! Created a list of the indices of droplets sorted by radius
            ! Iterate over all the bins, sequentially
            outer_bin_loop : do i = 1, NR-1
                rcubed_max = (r_bin_edges(i+1)*1d-6)**3.0_dp
                ! Iterate over all the unprocessed droplets
                inner_droplet_loop : do j = droplet_idx, n_part
                    if ( droplets(sorted_indices(j))%rcubed < rcubed_max ) then
                        ! Add the total mass associated with this superdroplet to the bin
                        bin_counts_raw(i) = bin_counts_raw(i) + &
                            droplets(sorted_indices(j))%mass * real(droplets(sorted_indices(j))%multi, kind=dp)
                        droplet_idx = droplet_idx + 1
                    else
                        exit inner_droplet_loop
                    end if
                end do inner_droplet_loop
            end do outer_bin_loop

            ! Compute a rolling window median to smooth the binned counts.
            call rolling_window_median(bin_counts_raw, bin_counts_smooth, smooth_window)

            call csv_out%add(stopwatch%total_seconds())
            call csv_out%add(bin_counts_smooth, real_fmt="(ES12.5)")

            call csv_out%next_row()
        end if

        if ( DEBUG ) print '(A,I0,A,I0,A,I0,A,I0,A,ES12.3,A,ES12.3,A,I0,A,ES12.3,A)', &
            "STEP: ", step, " (", stopwatch%minutes, " min ", stopwatch%seconds, " s) | Collisions: ", &
            collision_step_result%counter, " | Probs: ", &
            collision_step_result%min_prob, " - ", collision_step_result%max_prob, " [", &
            collision_step_result%big_probs, "] | Total Water: ", total_water(droplets), " kg"

        step = step + 1
        call stopwatch%increment(t_c)

    end do main_loop

    print '(/,/,A)', "Simulation completed successfully."
    wm_i = total_water(droplets)
    print '(A,F12.3,A,F5.1,A)', "Remaining water mass: ", wm_i, " kg (", &
          (wm_i / wm0) * 100.d0, "% of initial)"

    call csv_out%close(status_ok)
    if (.not. status_ok) then
        print *, "ERROR: Could not close output CSV file."
        stop
    end if
    
end program main
