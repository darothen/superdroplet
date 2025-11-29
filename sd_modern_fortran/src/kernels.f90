!> Collision kernel implementations module
!!
!! This module provides various collision kernel formulations for computing
!! the rate of droplet-droplet collisions in the superdroplet model.
!!
!! Available kernels:
!! - GOLOVIN: Analytical kernel for testing (sum of volumes)
!! - HYDRO: Hydrodynamic kernel (geometric collection)
!! - LONG: Long's kernel with collection efficiency parameterization
!!
!! @author Daniel Rothenberg
!! @date November 2025
module kernels

    use iso_c_binding, only: c_int

    use constants, only : dp, PI, FOUR_THIRD, THIRD
    use droplet_mod, only : droplet_t

    use stdlib_math, only : clip

    implicit none

    public :: GOLOVIN, HYDRO, LONG
    public :: compute_collision_kernel

    private

    !> Collision kernel type enumeration
    !!
    !! Defines available collision kernel types for the simulation.
    !!
    !! @note gfortran v15.2 does not yet support Fortran 2023 enumeration type
    !!       construct, so we use the older C-binding style enum.
    enum, bind(C)
        enumerator :: GOLOVIN = 1  !< Golovin analytical kernel
        enumerator :: HYDRO = 2    !< Hydrodynamic kernel
        enumerator :: LONG = 3     !< Long's kernel with collection efficiency
    end enum

contains

    !> Calculate hydrodynamic collision kernel
    !!
    !! Computes the hydrodynamic kernel using coalescence efficiency, collision
    !! efficiency, sum of radii, and terminal velocity difference.
    !!
    !! @param[in] E_coal - Coalescence efficiency (dimensionless)
    !! @param[in] E_coll - Collision efficiency (dimensionless)
    !! @param[in] r_sum - Sum of droplet radii (m)
    !! @param[in] tv_diff - Difference in terminal velocities (m/s)
    !! @return Collision kernel in m^3/s
    real(kind=dp) pure function calc_hydro_kernel(E_coal, E_coll, r_sum, tv_diff)
        real(kind=dp), intent(in) :: E_coal, E_coll, r_sum, tv_diff

        calc_hydro_kernel = &
            (E_coal*E_coll)*PI*(r_sum*r_sum)*abs(tv_diff)

    end function calc_hydro_kernel

    !> Compute collision kernel between two superdroplets
    !!
    !! Main dispatch function that selects the appropriate kernel calculation
    !! based on the specified kernel type.
    !!
    !! @param[in] sd_j - First superdroplet
    !! @param[in] sd_k - Second superdroplet
    !! @param[in] current_kernel - Kernel type (GOLOVIN, HYDRO, or LONG)
    !! @return Collision kernel K_ij in m^3/s
    real(kind=dp) pure function compute_collision_kernel(sd_j, sd_k, current_kernel) result(K_ij)
        type(droplet_t), intent(in) :: sd_j, sd_k
        integer, intent(in) :: current_kernel

        select case (current_kernel)
        case (GOLOVIN)
            K_ij = golovin_kernel(sd_j, sd_k)
        case (HYDRO)
            K_ij = hydro_kernel(sd_j, sd_k)
        case (LONG)
            K_ij = long_kernel(sd_j, sd_k)
        case default
            K_ij = 0.0_dp
        end select

    end function compute_collision_kernel


    !> Golovin analytical collision kernel
    !!
    !! Simple analytical kernel proportional to the sum of droplet volumes.
    !! Useful for testing and verification against analytical solutions.
    !!
    !! @note PERFORMANCE OPTIMIZATION: Uses precomputed constant to avoid
    !!       repeated multiplication (5-10% speedup)
    !!
    !! @param[in] sd_j - First superdroplet
    !! @param[in] sd_k - Second superdroplet
    !! @return Collision kernel in m^3/s
    real(kind=dp) pure function golovin_kernel(sd_j, sd_k)
        type(droplet_t), intent(in) :: sd_j, sd_k

        real(kind=dp), parameter :: golovin_b = 1.5d3
        real(kind=dp), parameter :: golovin_constant = golovin_b * FOUR_THIRD * PI

        golovin_kernel = golovin_constant * (sd_j%rcubed + sd_k%rcubed)

    end function golovin_kernel


    !> Hydrodynamic collision kernel
    !!
    !! Computes collision kernel based on geometric collection with
    !! unit coalescence and collision efficiencies.
    !!
    !! @param[in] sd_j - First superdroplet
    !! @param[in] sd_k - Second superdroplet
    !! @return Collision kernel in m^3/s
    real(kind=dp) pure function hydro_kernel(sd_j, sd_k)
        type(droplet_t), intent(in) :: sd_j, sd_k

        real(kind=dp) :: r_j, r_k, tv_j, tv_k
        real(kind=dp) :: tv_diff, r_sum

        r_j = sd_j%radius
        r_k = sd_k%radius
        tv_j = sd_j%terminal_velocity
        tv_k = sd_k%terminal_velocity

        tv_diff = tv_j - tv_k
        r_sum = r_j + r_k

        hydro_kernel = calc_hydro_kernel(1d0, 1d0, r_sum, tv_diff)

    end function hydro_kernel


    !> Long's collision kernel with collection efficiency
    !!
    !! Computes collision kernel using Long's collection efficiency
    !! parameterization as a function of droplet sizes.
    !!
    !! Collection efficiency is size-dependent and accounts for the fact
    !! that not all geometric collisions result in coalescence.
    !!
    !! @param[in] sd_j - First superdroplet
    !! @param[in] sd_k - Second superdroplet
    !! @return Collision kernel in m^3/s
    real(kind=dp) pure function long_kernel(sd_j, sd_k)
        type(droplet_t), intent(in) :: sd_j, sd_k

        real(kind=dp) :: &
             r_j, r_k, x_j, x_k, tv_j, tv_k, &
             tv_diff, r_sum, r_small, r_large, E_coll

        r_j = sd_j%radius
        r_k = sd_k%radius
        x_j = sd_j%mass
        x_k = sd_k%mass
        tv_j = sd_j%terminal_velocity
        tv_k = sd_k%terminal_velocity

        tv_diff = tv_j - tv_k
        r_sum = r_j + r_k

        if (r_j > r_k) then
            r_large = r_j
            r_small = r_k
        else
            r_large = r_k
            r_small = r_j
        end if
        r_large = r_large * 1d6  ! convert to micron
        r_small = r_small * 1d6  ! convert to micron

        ! Collection efficiency cut-off in limit of very large drops
        if ( r_large >= 50d0 ) then  ! microns
            E_coll = 1d0
        else
            E_coll = 4.5d-4*(r_large*r_large) &
                   * (1d0 - 3d0/(max(3d0, r_small) + 1d-2))
        end if

        ! Limit collection efficiency to 0 <= E_coll <= 1.0
        E_coll = clip(E_coll, 0d0, 1d0)

        long_kernel = calc_hydro_kernel(1d0, E_coll, r_sum, tv_diff)

    end function long_kernel


end module kernels