!> @file test_core.f90
!! @brief Small deterministic unit tests for core helper routines.
!> @brief Exercises basic derived-type and checkpoint helper behavior.
program test_core
    use model_types_mod, only: monte_carlo_study_config_t
    use model_types_mod, only: spatial_grid_t
    use model_types_mod, only: get_state_dimension
    use monte_carlo_study_mod, only: compute_monte_carlo_case_count
    use parameter_ml_estimation_mod, only: build_uniform_checkpoints
    implicit none

    call test_state_dimension()
    call test_uniform_checkpoints_even_spacing()
    call test_uniform_checkpoints_clamped_request()
    call test_monte_carlo_case_count()

    write (*, '(a)') "Unit tests passed."

contains

    subroutine assert_true(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (*, '(2a)') "Assertion failed: ", trim(message)
            error stop
        end if
    end subroutine assert_true

    subroutine test_state_dimension()
        type(spatial_grid_t) :: grid

        grid%nx = 4
        grid%ny = 3

        call assert_true(&
            get_state_dimension(grid) == 11, &
            "state dimension should be nx * ny - 1" &
        )
    end subroutine test_state_dimension

    subroutine test_uniform_checkpoints_even_spacing()
        integer, allocatable :: checkpoints(:)

        call build_uniform_checkpoints(10, 3, 4, checkpoints)

        call assert_true(size(checkpoints) == 3, "expected three checkpoints")
        call assert_true(checkpoints(1) == 4, "first checkpoint should be 4")
        call assert_true(checkpoints(2) == 7, "second checkpoint should be 7")
        call assert_true(checkpoints(3) == 10, "third checkpoint should be 10")
    end subroutine test_uniform_checkpoints_even_spacing

    subroutine test_uniform_checkpoints_clamped_request()
        integer, allocatable :: checkpoints(:)

        call build_uniform_checkpoints(6, 10, 4, checkpoints)

        call assert_true(size(checkpoints) == 3, "request should clamp to 3")
        call assert_true(all(checkpoints == [4, 5, 6]), &
            "clamped checkpoints should cover all available points")
    end subroutine test_uniform_checkpoints_clamped_request

    subroutine test_monte_carlo_case_count()
        type(monte_carlo_study_config_t) :: study_config

        allocate (study_config%basis_levels(2))
        allocate (study_config%observation_counts(3))
        allocate (study_config%time_steps(4))

        study_config%basis_levels = [10, 20]
        study_config%observation_counts = [1000, 5000, 20000]
        study_config%time_steps = [1.0e-6, 1.0e-5, 1.0e-4, 1.0e-3]

        call assert_true(&
            compute_monte_carlo_case_count(study_config) == 24, &
            "case count should equal the cartesian-product size" &
        )
    end subroutine test_monte_carlo_case_count

end program test_core
