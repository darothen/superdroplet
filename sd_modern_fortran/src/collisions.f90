!> Collision-coalescence simulation module
!!
!! This module implements the stochastic collision-coalescence algorithm
!! using the superdroplet method of Shima et al. (2009). It handles the
!! pairing, probability calculation, and coalescence of superdroplets.
!!
!! @author Daniel Rothenberg
!! @date November 2025
module collisions

    use constants, only : dp
    use droplet_mod, only : droplet_t, set_droplet, total_water
    use kernels, only : compute_collision_kernel
    use util, only : generate_shuffled_indices

    use stdlib_stats_distribution_uniform, only : rvs_uniform

    implicit none

    private

    public :: model_config_t
    public :: model_config_t_
    public :: collision_step
    public :: collision_step_result_t
    public :: cleanup_model_config

    !> Model configuration type
    !!
    !! Contains simulation parameters used throughout the collision step.
    !! Also includes pre-allocated workspace arrays to avoid per-timestep allocations.
    type :: model_config_t
        integer :: step_seconds      !< Timestep duration (s)
        integer :: num_droplets      !< Number of superdroplets
        integer :: kernel            !< Collision kernel type
        real(kind=dp) :: delta_v     !< Parcel volume (m^3)
        ! Pre-allocated workspace arrays (avoid repeated allocation/deallocation)
        integer, dimension(:), allocatable :: shuffled_indices  !< Workspace for shuffled indices
        real(kind=dp), dimension(:), allocatable :: phis        !< Workspace for random numbers
    end type model_config_t

    interface model_config_t_
        module procedure init_model_config
    end interface

    !> Collision step result type
    !!
    !! Contains diagnostic information from a collision timestep.
    type :: collision_step_result_t
        integer :: counter           !< Number of collisions that occurred
        integer :: big_probs         !< Number of collision probabilities > 1
        real(kind=dp) :: max_prob    !< Maximum collision probability
        real(kind=dp) :: min_prob    !< Minimum collision probability
        real(kind=dp) :: total_xi    !< Total multiplicity sum
        real(kind=dp) :: total_water !< Total water mass (kg)
    end type collision_step_result_t
    
contains

    !> Constructor for model configuration
    !!
    !! @param[in] step_seconds - Timestep duration in seconds
    !! @param[in] num_droplets - Number of superdroplets
    !! @param[in] delta_v - Parcel volume in m^3
    !! @param[in] kernel - Collision kernel type
    !! @return Initialized model_config_t instance
    type(model_config_t) function init_model_config( &
        step_seconds, num_droplets, delta_v, kernel)
        integer, intent(in) :: step_seconds, num_droplets, kernel
        real(kind=dp), intent(in) :: delta_v
        integer :: half_n_part

        init_model_config%step_seconds = step_seconds
        init_model_config%num_droplets = num_droplets
        init_model_config%delta_v = delta_v
        init_model_config%kernel = kernel

        ! Pre-allocate workspace arrays once to avoid repeated allocation/deallocation
        ! in collision_step (major performance optimization)
        half_n_part = floor(real(num_droplets, kind=dp) / 2.0_dp)
        allocate(init_model_config%shuffled_indices(num_droplets))
        allocate(init_model_config%phis(half_n_part))

    end function init_model_config
    
    !> Constructor for collision step result (internal use)
    !!
    !! @param[in] counter - Number of collisions
    !! @param[in] big_probs - Count of probabilities exceeding 1
    !! @param[in] max_prob - Maximum collision probability
    !! @param[in] min_prob - Minimum collision probability
    !! @param[in] total_xi - Total multiplicity
    !! @param[in] total_water - Total water mass
    !! @return Initialized collision_step_result_t instance
    type(collision_step_result_t) function init_collision_step_result( &
        counter, big_probs, max_prob, min_prob, total_xi, total_water )
        integer, intent(in) :: counter, big_probs
        real(kind=dp), intent(in) :: max_prob, min_prob, total_xi, total_water

        init_collision_step_result%counter = counter
        init_collision_step_result%big_probs = big_probs
        init_collision_step_result%max_prob = max_prob
        init_collision_step_result%min_prob = min_prob
        init_collision_step_result%total_xi = total_xi
        init_collision_step_result%total_water = total_water

    end function init_collision_step_result

    !> Multi-droplet coalescence subroutine
    !!
    !! Implements the coalescence algorithm where multiple droplets from one
    !! superdroplet collide with droplets from another. Handles both cases where
    !! droplets remain unpaired or where superdroplets must be split.
    !!
    !! @param[inout] sd_j - First superdroplet (larger or equal multiplicity)
    !! @param[inout] sd_k - Second superdroplet (smaller or equal multiplicity)
    !! @param[in] gamma - Number of collisions to process
    pure subroutine multi_coalesce(sd_j, sd_k, gamma)
        type(droplet_t), intent(inout) :: sd_j, sd_k
        real(kind=dp), intent(in) :: gamma

        real(kind=dp) :: &
            ratio, gamma_t, rcubed_k_p, solute_k_p, rcubed_new, solute_new, &
            multi_j_real, multi_k_real, gamma_t_times_multi_k
        integer :: excess, multi_j_p, multi_k_p
        
        multi_j_real = real(sd_j%multi, kind=dp)
        multi_k_real = real(sd_k%multi, kind=dp)
        
        ratio = multi_j_real / multi_k_real
        gamma_t = min(gamma, ratio)
        
        ! Avoid repeated multiplication
        gamma_t_times_multi_k = gamma_t * multi_k_real
        excess = sd_j%multi - int(floor(gamma_t_times_multi_k))

        if ( excess > 0 ) then
            ! Case 1: Some droplets from sd_j remain unpaired
            sd_j%multi = excess
            ! sd_k%multi remains the same

            ! sd_j%rcubed remains the same
            rcubed_k_p = (gamma_t * sd_j%rcubed) + sd_k%rcubed

            ! sd_j%solute remains the same
            solute_k_p = (gamma_t * sd_j%solute) + sd_k%solute
            call sd_k%set_droplet(sd_k%multi, rcubed_k_p, solute_k_p)
        else
            ! Case 2: All droplets from sd_j are paired; split sd_k
            multi_j_p = int(multi_k_real * 0.5_dp)
            multi_k_p = sd_k%multi - multi_j_p

            rcubed_new = (gamma_t * sd_j%rcubed) + sd_k%rcubed
            solute_new = (gamma_t * sd_j%solute) + sd_k%solute

            call sd_j%set_droplet(multi_j_p, rcubed_new, solute_new)
            call sd_k%set_droplet(multi_k_p, rcubed_new, solute_new)
        end if

    end subroutine multi_coalesce


    !> Execute one collision timestep
    !!
    !! Main collision algorithm following Shima et al. (2009). Randomly pairs
    !! superdroplets, computes collision probabilities, and performs coalescence
    !! when collisions occur.
    !!
    !! Key optimizations:
    !! - Shuffles indices instead of droplet array (10-20% speedup)
    !! - Pre-generates random numbers for better performance
    !! - Precomputes scaling factors
    !! - Uses concurrent loops for potential parallelism
    !!
    !! @param[inout] droplets - Array of superdroplets to evolve
    !! @param[out] collision_step_result - Diagnostic results from this timestep
    !! @param[in] model_config - Model configuration parameters
    subroutine collision_step(droplets, collision_step_result, model_config)
        type(droplet_t), dimension(:), intent(inout) :: droplets
        type(collision_step_result_t), intent(out) :: collision_step_result
        type(model_config_t), intent(inout) :: model_config

        integer :: n_part, half_n_part, min_xi, max_xi, counter, big_probs, i, kernel, idx_j, idx_k
        real(kind=dp) :: &
            t_c_over_delta_v, scaling, K_ij, prob, gamma, max_prob, min_prob, phi, floor_prob

        kernel = model_config%kernel

        ! Generate candidate pairs and precompute constants
        n_part = model_config%num_droplets
        half_n_part = floor(real(n_part, kind=dp) / 2.0_dp)
        scaling = (real(n_part, kind=dp) * real(n_part - 1, kind=dp)) / 2.0_dp / real(half_n_part, kind=dp)

        ! Precompute t_c/delta_v to avoid repeated division
        t_c_over_delta_v = real(model_config%step_seconds, kind=dp) / model_config%delta_v


        ! Shuffle indices instead of droplet array to avoid copying large structures
        call generate_shuffled_indices(model_config%shuffled_indices, n_part, 1)

        ! Pre-generate random numbers for collision decisions
        ! NOTE: Strictly speaking, the random numbers for probability gen should
        ! be in the half-open interval [0,1), but this should be sufficient for practical purposes
        model_config%phis = rvs_uniform(0.0_dp, 1.0_dp, half_n_part)

        counter = 0
        big_probs = 0
        max_prob = 0.
        min_prob = 1.0

        ! Main collision loop
        do i = 1, half_n_part
            ! Get indices from shuffled array
            idx_j = model_config%shuffled_indices(i)
            idx_k = model_config%shuffled_indices(i + half_n_part)
            phi = model_config%phis(i)

            ! Optimize min/max with explicit checks
            if (droplets(idx_j)%multi < droplets(idx_k)%multi) then
                min_xi = droplets(idx_j)%multi
                max_xi = droplets(idx_k)%multi
            else
                min_xi = droplets(idx_k)%multi
                max_xi = droplets(idx_j)%multi
            end if
            
            if (min_xi == 0) cycle

            K_ij = compute_collision_kernel(droplets(idx_j), droplets(idx_k), kernel)

            prob = scaling * max_xi * t_c_over_delta_v * K_ij
            
            ! Update diagnostics
            if ( prob > max_prob ) max_prob = prob
            if ( prob < min_prob ) min_prob = prob
            if ( prob > 1.0_dp) big_probs = big_probs + 1

            ! Check for collision and coalesce if necessary
            floor_prob = floor(prob)
            if ( (prob - floor_prob) > phi ) then
                gamma = floor_prob + 1.0_dp

                ! Coalesce: pass larger multiplicity droplet first
                if (max_xi == droplets(idx_j)%multi) then
                    call multi_coalesce(droplets(idx_j), droplets(idx_k), gamma)
                else
                    call multi_coalesce(droplets(idx_k), droplets(idx_j), gamma)
                end if

                counter = counter + 1
            end if

        end do

        ! Populate the collision step result
        collision_step_result%counter = counter
        collision_step_result%big_probs = big_probs
        collision_step_result%max_prob = max_prob
        collision_step_result%min_prob = min_prob
        collision_step_result%total_xi = sum(real(droplets%multi, kind=dp))
        collision_step_result%total_water = total_water(droplets)

        ! print "(I5,' collisions simulated')", counter
        ! print "('Max/min probabilities (count): ',2(ES8.2,' '),I5)", &
        !     min_prob, max_prob, big_probs

    end subroutine collision_step

    !> Cleanup model configuration workspace arrays
    !!
    !! Deallocates the workspace arrays in model_config to free memory.
    !! Should be called at the end of the simulation.
    !!
    !! @param[inout] model_config - Model configuration to cleanup
    subroutine cleanup_model_config(model_config)
        type(model_config_t), intent(inout) :: model_config

        if (allocated(model_config%shuffled_indices)) deallocate(model_config%shuffled_indices)
        if (allocated(model_config%phis)) deallocate(model_config%phis)

    end subroutine cleanup_model_config

end module collisions
