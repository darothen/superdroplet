!> Utility functions module
!!
!! This module provides various utility functions for array manipulation,
!! shuffling, and statistical operations used throughout the simulation.
!!
!! @author Daniel Rothenberg
!! @date November 2025
module util

    use droplet_mod, only : droplet_t

    use stdlib_kinds, only : stdlib_dp => dp
    use stdlib_stats, only : median

    implicit none

    private

    public :: shuffle_droplet_array, generate_shuffled_indices, rolling_window_median

contains

    !> Shuffle droplet array in place (deprecated)
    !!
    !! Shuffles an array of droplets using the Fisher-Yates (Knuth) shuffle.
    !!
    !! @warning DEPRECATED: This function is inefficient as it copies large
    !!          droplet structures (~80 bytes each). Use generate_shuffled_indices
    !!          instead for 10-20% performance improvement.
    !!
    !! @param[inout] droplets - Array of droplets to shuffle
    subroutine shuffle_droplet_array(droplets)
        type(droplet_t), dimension(:), intent(inout) :: droplets

        integer :: i, randpos
        type(droplet_t) :: temp
        real :: r

        do i = size(droplets), 2, -1
            ! Use the built-in Fortran random number generator, which should return
            ! a value in [0.0, 1.0) generating using the xoshiro256** PRNG.
            call random_number(r)
            randpos = int(r*i) + 1
            temp = droplets(randpos)
            droplets(randpos) = droplets(i)
            droplets(i) = temp
        end do

    end subroutine shuffle_droplet_array

    !> Generate shuffled index array
    !!
    !! Creates a shuffled array of indices using the Fisher-Yates algorithm.
    !! This is much more efficient than shuffling the actual droplet array
    !! as it only shuffles small integers rather than large structures.
    !!
    !! @note PERFORMANCE OPTIMIZATION: Provides 10-20% speedup compared to
    !!       shuffling the droplet array directly.
    !!
    !! @param[out] indices - Output array of shuffled indices
    !! @param[in] n - Number of indices to generate
    !! @param[in] offset - Starting index (0 for 0-based, 1 for 1-based arrays)
    subroutine generate_shuffled_indices(indices, n, offset)
        integer, dimension(:), intent(out) :: indices
        integer, intent(in) :: n, offset

        integer :: i, j, temp
        real :: r

        ! Initialize indices sequentially
        do i = 1, n
            indices(i) = offset + (i - 1)
        end do

        ! Fisher-Yates shuffle on indices
        do i = n, 2, -1
            call random_number(r)
            j = int(r * i) + 1
            ! Swap indices(i) and indices(j)
            temp = indices(j)
            indices(j) = indices(i)
            indices(i) = temp
        end do

    end subroutine generate_shuffled_indices

    !> Compute rolling window median
    !!
    !! Applies a rolling window median filter to smooth the input array.
    !! Used for smoothing binned droplet size distribution output.
    !!
    !! @param[in] input_array - Input array to smooth
    !! @param[out] output_array - Smoothed output array (same size as input)
    !! @param[in] window_size - Size of rolling window (should be odd)
    subroutine rolling_window_median(input_array, output_array, window_size)
        real(kind=stdlib_dp), dimension(:), intent(in) :: input_array
        real(kind=stdlib_dp), dimension(:), intent(out) :: output_array
        integer, intent(in) :: window_size

        integer :: n, half_window, i, j, start_idx, end_idx, window_len
        real(kind=stdlib_dp), dimension(:), allocatable :: window_values
        real(kind=stdlib_dp) :: median_value

        n = size(input_array)
        half_window = window_size / 2

        allocate(window_values(window_size))

        do i = 1, n
            ! Determine the start and end indices for the current window
            start_idx = max(1, i - half_window)
            end_idx = min(n, i + half_window)

            ! Collect values in the current window
            window_len = 0
            do j = start_idx, end_idx
                window_len = window_len + 1
                window_values(window_len) = input_array(j)
            end do

            ! Compute median over the collected values
            median_value = median(window_values(1:window_len))

            output_array(i) = median_value
        end do

        deallocate(window_values)

    end subroutine rolling_window_median

end module