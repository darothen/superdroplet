!> Droplet data structure and operations module
!!
!! This module defines the core droplet_t type and associated operations
!! for the superdroplet method. Each superdroplet represents a collection
!! of identical real droplets.
!!
!! @author Daniel Rothenberg
!! @date November 2025
module droplet_mod

    use constants, only: dp, FOUR_THIRD, PI, RHO_WATER, THIRD
    implicit none
    private

    public :: droplet_t  ! droplet type
    public :: droplet_t_ ! droplet constructor bound to init_droplet
    public :: total_water
    public :: set_droplet

    !> Superdroplet data structure
    !!
    !! Represents a collection of identical real droplets in the superdroplet method.
    !! Contains both physical properties and derived quantities.
    !!
    !! @note Terminal velocity is cached and only recomputed when rcubed changes
    !!       for performance optimization.
    type :: droplet_t
        integer :: multi                !< Multiplicity (number of real droplets represented)
        real(kind=dp) :: &
            rcubed,                     & !< Cube of droplet radius (m^3)
            solute,                     & !< Solute mass (kg)
            density,                    & !< Droplet density (kg/m^3)
            radius,                     & !< Droplet radius (m)
            volume,                     & !< Droplet volume (m^3)
            mass,                       & !< Single droplet mass (kg)
            terminal_velocity             !< Terminal fall velocity (m/s)
    contains
        procedure, public :: set_droplet  !< Update droplet properties
    end type droplet_t

    interface droplet_t_
        module procedure init_droplet
    end interface

contains

    !> Constructor for droplet type
    !!
    !! Initializes a droplet with given multiplicity and radius.
    !! Computes all derived properties including volume, mass, and terminal velocity.
    !!
    !! @param[in] multi - Droplet multiplicity (number of real droplets represented)
    !! @param[in] radius - Droplet radius in meters
    !! @return Initialized droplet_t instance
    type(droplet_t) function init_droplet(multi, radius)

        integer, intent(in) :: multi
        real(kind=dp), intent(in) :: radius
        
        init_droplet%multi = multi
        init_droplet%radius = radius
        init_droplet%rcubed = radius**3.0_dp
        init_droplet%solute = 0.0_dp
        init_droplet%density = RHO_WATER
        init_droplet%volume = (4.0_dp / 3.0_dp) * PI * radius**3.0_dp
        init_droplet%mass = init_droplet%volume * init_droplet%density
        init_droplet%terminal_velocity = compute_terminal_velocity(radius, init_droplet%mass)

    end function init_droplet

    !> Mutator to update droplet properties
    !!
    !! Updates droplet multiplicity and rcubed value, recomputing derived
    !! properties only when rcubed has actually changed.
    !!
    !! @note Terminal velocity is only recomputed
    !!       when rcubed changes to avoid expensive calculations.
    !!
    !! @param[inout] self - Droplet instance to update
    !! @param[in] multi - New multiplicity value
    !! @param[in] rcubed - New radius cubed value (m^3)
    !! @param[in] solute - Optional new solute mass (kg)
    pure subroutine set_droplet(self, multi, rcubed, solute)

        class(droplet_t), intent(inout) :: self
        integer, intent(in) :: multi
        real(kind=dp), intent(in) :: rcubed
        real(kind=dp), intent(in), optional :: solute
        logical :: rcubed_changed

        ! Check if rcubed actually changed (avoid expensive terminal velocity recomputation)
        rcubed_changed = (abs(self%rcubed - rcubed) > 1e-20_dp)

        self%multi = multi
        self%rcubed = rcubed

        ! Only recompute derived properties if rcubed changed
        if (rcubed_changed) then
            self%radius = rcubed**(THIRD)
            self%volume = FOUR_THIRD * PI * rcubed  ! Optimized: use rcubed directly
            self%mass = self%volume * self%density
            self%terminal_velocity = compute_terminal_velocity(self%radius, self%mass)
        end if

        if ( present(solute) ) self%solute = solute

    end subroutine set_droplet

    !> Compute droplet terminal fall velocity
    !!
    !! Calculates terminal velocity using size-dependent power-law relationships
    !! based on Beard (1976) parameterization.
    !!
    !! @param[in] radius - Droplet radius (m)
    !! @param[in] mass - Droplet mass (kg)
    !! @return Terminal velocity in m/s
    pure function compute_terminal_velocity(radius, mass) result(tv)
        real(kind=dp), intent(in) :: radius, mass
        real(kind=dp) :: tv

        real(kind=dp) :: d, cbrt_x, alpha, x_to_beta

        d = 2.0_dp * radius * 1e6_dp ! diameter, m -> μm
        cbrt_x = (mass * 1e3_dp) ** THIRD ! cbr(mass), [kg -> g]^(1/3)

        if (d <= 134.43_dp) then 
            alpha = 4.5795e5_dp
            x_to_beta = cbrt_x * cbrt_x
        else if ((134.34_dp < d) .and. (d <= 1511.64_dp)) then
            alpha = 4962.0_dp
            x_to_beta = cbrt_x
        else if ((1511.64_dp < d) .and. (d <= 3477.84)) then
            alpha = 1732.0_dp
            x_to_beta = sqrt(cbrt_x)
        else
            alpha = 917.0_dp
            x_to_beta = 1.0_dp
        end if

        tv = 1e-2_dp * alpha * x_to_beta ! cm/s -> m/s

    end function compute_terminal_velocity

    !> Calculate total water mass across all superdroplets
    !!
    !! Sums the total water mass by multiplying each droplet's mass by its
    !! multiplicity and summing over all superdroplets.
    !!
    !! @param[in] droplets - Array of superdroplets
    !! @return Total water mass in kg
    pure function total_water(droplets) result(total_mass)
        type(droplet_t), dimension(:), intent(in) :: droplets
        real(kind=dp) :: total_mass
        integer :: i

        total_mass = 0.0_dp
        do i = 1, size(droplets)
            total_mass = total_mass + droplets(i)%mass * real(droplets(i)%multi, kind=dp)
        end do

    end function total_water

end module droplet_mod