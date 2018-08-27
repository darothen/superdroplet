module droplet_class

    use constants
    implicit none
    private

    public :: droplet  ! droplet type
    public :: droplet_ ! droplet constructor bound to init_droplet

    integer(kind=ikind), public :: droplet_count = 0

    type :: droplet
        integer(kind=lkind) :: multi
        real(kind=rkind)    :: rcubed
        real(kind=rkind)    :: solute
        real(kind=rkind)    :: density

        real(kind=rkind)    :: radius, volume, mass
        integer(kind=lkind) :: id
    contains
        procedure, public :: set_droplet
        procedure, public :: terminal_v=>calc_terminal_v
    end type droplet

    interface droplet_
        module procedure init_droplet
    end interface

! ----------------------------------------------------------------------------
contains

    ! Constructor
    type(droplet) function init_droplet(multi, rcubed, solute, density)

        integer(kind=lkind), intent(in) :: multi
        real(kind=rkind), intent(in) :: rcubed
        real(kind=rkind), optional, intent(in) :: solute, density

        init_droplet%multi = multi
        init_droplet%rcubed = rcubed

        if ( present(solute) ) then
            init_droplet%solute = solute
        else
            init_droplet%solute = 0.
        endif

        if ( present(density) ) then
            init_droplet%density = density
        else
            init_droplet%density = RHO_WATER
        endif

        init_droplet%radius = rcubed**(THIRD)
        init_droplet%volume = rcubed * 4. / PI / 3.
        init_droplet%mass   = init_droplet%volume * init_droplet%density

        droplet_count = droplet_count + 1_ikind
        init_droplet%id = droplet_count

    end function init_droplet

    ! Mutator - all attributes
    subroutine set_droplet(self, multi, rcubed, solute)

        class(droplet), intent(inout)  :: self
        integer(kind=lkind), intent(in) :: multi
        real(kind=rkind), intent(in)    :: rcubed, solute

        self%multi = multi
        self%rcubed = rcubed
        self%solute = solute

        self%radius = rcubed**THIRD

    end subroutine set_droplet

    function calc_terminal_v(self) result(terminal_v)
        class(droplet), intent(in) :: self
        real(kind=rkind) :: terminal_v

        real(kind=rkind) :: alpha, r, d, x, x_to_beta
        r = self%radius
        d = 2.*r*1e6 ! diameter, m-> micron
        x = self%mass * 1e3 ! mass, kg -> g

        if ( d <= 134.43 ) then
            alpha = 4.5795e5
            x_to_beta = x**(2*THIRD)
        else if ( (134.43 < d) .and. (d <= 1511.64) ) then
            alpha = 4962.0
            x_to_beta = x**(THIRD)
        else if ( (1511.64 < d) .and. (d <= 3477.84) ) then
            alpha = 1732.0
            x_to_beta = x**(THIRD/2)
        else
            alpha = 917.0
            x_to_beta = 1.0
        end if

        terminal_v = 1e-2 * alpha * x_to_beta ! cm/s -> m/s

    end function calc_terminal_v

end module
