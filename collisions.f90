module collisions

    use constants
    use droplet_class, only: droplet
    use util, only: urand, shuffle_droplet_array

    use fgsl

    implicit none

contains

    real(kind=rkind) function kernel(sd_j, sd_k)
        ! TODO: implement alternative collision kernels
        type(droplet), intent(in) :: sd_j, sd_k
        real(kind=rkind), parameter :: golovin_b = 1.5d3

        kernel = golovin_b*(sd_j%rcubed + sd_k%rcubed)*4.*PI/3.

    end function kernel

    subroutine multi_coalesce(sd_j, sd_k, gamma)
        type(droplet), intent(inout) :: sd_j, sd_k
        real(kind=rkind), intent(in) :: gamma

        real(kind=rkind) :: gamma_tilde, rcubed_j_p, rcubed_k_p, &
                            solute_j_p, solute_k_p, mjp, mkp
        integer(kind=lkind) :: excess, multi_j_p, multi_k_p

        mjp = real(sd_j%multi, kind=rkind)
        mkp = real(sd_k%multi, kind=rkind)
        gamma_tilde = min(gamma, real(floor(mjp/mkp), kind=rkind))
        excess = sd_j%multi - floor(gamma_tilde*sd_k%multi, &
                                    kind=lkind)

        if ( excess > 0 ) then
            multi_j_p = excess
            multi_k_p = sd_k%multi

            rcubed_j_p = sd_j%rcubed
            rcubed_k_p = gamma_tilde*rcubed_j_p + sd_k%rcubed

            solute_j_p = sd_j%solute
            solute_k_p = gamma_tilde*solute_j_p + sd_k%solute

            call sd_j%set_droplet(multi_j_p, rcubed_j_p, solute_j_p)
            call sd_k%set_droplet(multi_k_p, rcubed_k_p, solute_k_p)
        else
            multi_j_p = floor(sd_k%multi / 2.)
            multi_k_p = sd_k%multi - multi_j_p

            rcubed_j_p = gamma_tilde*sd_j%rcubed + sd_k%rcubed
            rcubed_k_p = rcubed_j_p

            solute_j_p = gamma_tilde*sd_j%solute + sd_k%solute
            solute_k_p = solute_j_p

            call sd_j%set_droplet(multi_k_p, rcubed_j_p, solute_j_p)
            call sd_k%set_droplet(multi_j_p, rcubed_k_p, solute_k_p)

        end if

    end subroutine multi_coalesce

    subroutine collision_step(droplets, t_c, delta_v)
        type(droplet), dimension(:), intent(inout) :: droplets
        real(kind=rkind), intent(in) :: t_c, delta_v

        integer(kind=lkind) :: n_part, half_n_part, max_xi
        real(kind=rkind) :: scaling, K_ij, prob, gamma

        real :: max_prob, min_prob, phi
        integer :: counter, big_probs, i
        type(droplet) :: sd_j, sd_k

        print *, "PREP STEPS"

        ! 1) Make the random permutation of the droplet list
        print *, "   SHUFFLE LIST"
        call shuffle_droplet_array(droplets)

        ! 2) Make the candidate pairs
        print *, "   GEN PAIRS"

        ! 3) Generate the uniform random numbers
        print *, "PROBABILITY LOOP"
        n_part = size(droplets)
        scaling = (n_part*(n_part-1)/2.)/floor(n_part/2., kind=lkind)

        print *, "PROB / COLLISION LOOP"
        counter = 0
        half_n_part = n_part/2

        big_probs = 0
        max_prob = 0.
        min_prob = 1.0

        do i = 1, half_n_part

            sd_j = droplets(i)
            sd_k = droplets(i + half_n_part)

            phi = urand()

            K_ij = kernel(sd_j, sd_k)
            max_xi = max(sd_j%multi, sd_k%multi)

            prob = scaling*max_xi*(t_c/delta_V)*K_ij

            if ( prob > max_prob ) max_prob = prob
            if ( prob < min_prob ) min_prob = prob
            if ( prob > 1) big_probs = big_probs + 1

            ! Check for collision and coalesce if necessary
            if ( (prob - floor(prob)) >= phi ) then
                gamma = floor(prob) + 1.

                ! print *, ">>>", i, sd_j%rcubed
                if (sd_j%multi < sd_k%multi) then
                    call multi_coalesce(sd_k, sd_j, gamma)
                else
                    call multi_coalesce(sd_j, sd_k, gamma)
                end if
                ! print *, "<<<", i, sd_j%rcubed

                ! Need to copy back into droplet array for
                ! updating
                droplets(i) = sd_j
                droplets(i + half_n_part) = sd_k

                counter = counter + 1
            end if

        end do

        print "(I5,' collisions simulated')", counter
        print "('Max/min probabilities (count): ',2(ES8.2,' '),I5)", &
            min_prob, max_prob, big_probs

    end subroutine collision_step

end module collisions
