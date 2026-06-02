!> @file test_core.f90
!! @brief Small deterministic unit tests for core helper routines.
!> @brief Exercises basic derived-type and checkpoint helper behavior.
program test_core
    use model_types_mod, only: dp
    use model_types_mod, only: monte_carlo_study_config_t
    use model_types_mod, only: spatial_grid_t
    use model_types_mod, only: get_state_dimension
    use monte_carlo_study_mod, only: compute_monte_carlo_case_count
    use parameter_ml_estimation_mod, only: build_uniform_checkpoints
    use solution_reconstruction_mod, only: build_endpoint_plot_coordinates
    use solution_reconstruction_mod, only: evaluate_cosine_mode
    use solution_reconstruction_mod, only: reconstruct_solution_field
    use solution_reconstruction_mod, only: reconstruct_solution_value
    implicit none

    call test_state_dimension()
    call test_uniform_checkpoints_even_spacing()
    call test_uniform_checkpoints_clamped_request()
    call test_monte_carlo_case_count()
    call test_endpoint_plot_coordinates()
    call test_cosine_mode_normalization()
    call test_modal_solution_reconstruction()

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

    subroutine assert_close(actual, expected, tolerance, message)
        real(dp), intent(in) :: actual
        real(dp), intent(in) :: expected
        real(dp), intent(in) :: tolerance
        character(len=*), intent(in) :: message

        if (abs(actual - expected) > tolerance) then
            write (*, '(2a)') "Assertion failed: ", trim(message)
            write (*, '(a,es18.8e3)') "  actual:   ", actual
            write (*, '(a,es18.8e3)') "  expected: ", expected
            error stop
        end if
    end subroutine assert_close

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

    subroutine test_endpoint_plot_coordinates()
        type(spatial_grid_t) :: grid
        real(dp), allocatable :: x_coordinates(:)
        real(dp), allocatable :: y_coordinates(:)

        grid%length_x = 5.0_dp
        grid%length_y = 2.0_dp

        call build_endpoint_plot_coordinates(&
            grid, 3, 5, x_coordinates, y_coordinates &
        )

        call assert_true(size(x_coordinates) == 3, "expected three x points")
        call assert_true(size(y_coordinates) == 5, "expected five y points")
        call assert_close(x_coordinates(1), 0.0_dp, 1.0e-12_dp, &
            "x grid should include left endpoint")
        call assert_close(x_coordinates(2), 2.5_dp, 1.0e-12_dp, &
            "x grid should include midpoint")
        call assert_close(x_coordinates(3), 5.0_dp, 1.0e-12_dp, &
            "x grid should include right endpoint")
        call assert_close(y_coordinates(5), 2.0_dp, 1.0e-12_dp, &
            "y grid should include top endpoint")
    end subroutine test_endpoint_plot_coordinates

    subroutine test_cosine_mode_normalization()
        type(spatial_grid_t) :: grid
        real(dp) :: mode_value

        grid%length_x = 5.0_dp
        grid%length_y = 5.0_dp

        mode_value = evaluate_cosine_mode(1, 0, 0.0_dp, 3.0_dp, grid)
        call assert_close(mode_value, sqrt(2.0_dp), 1.0e-12_dp, &
            "h_10 should carry sqrt(2) normalization")

        mode_value = evaluate_cosine_mode(1, 1, 0.0_dp, 0.0_dp, grid)
        call assert_close(mode_value, 2.0_dp, 1.0e-12_dp, &
            "h_11 should carry product normalization")
    end subroutine test_cosine_mode_normalization

    subroutine test_modal_solution_reconstruction()
        type(spatial_grid_t) :: grid
        integer :: mode_pairs(2, 2)
        real(dp) :: modal_coefficients(2)
        real(dp) :: value
        real(dp), allocatable :: field(:, :)
        real(dp), allocatable :: x_coordinates(:)
        real(dp), allocatable :: y_coordinates(:)

        grid%length_x = 5.0_dp
        grid%length_y = 5.0_dp
        mode_pairs = reshape([1, 0, 0, 1], shape(mode_pairs))
        modal_coefficients = [2.0_dp, -0.5_dp]

        call reconstruct_solution_value(&
            grid, mode_pairs, modal_coefficients, 0.0_dp, 0.0_dp, value &
        )
        call assert_close(value, 1.5_dp * sqrt(2.0_dp), 1.0e-12_dp, &
            "point reconstruction should sum both modal contributions")

        call build_endpoint_plot_coordinates(&
            grid, 2, 2, x_coordinates, y_coordinates &
        )
        call reconstruct_solution_field(&
            grid, mode_pairs, modal_coefficients, x_coordinates, y_coordinates, &
            field &
        )

        call assert_true(size(field, 1) == 2, "field x dimension should match")
        call assert_true(size(field, 2) == 2, "field y dimension should match")
        call assert_close(field(1, 1), value, 1.0e-12_dp, &
            "field reconstruction should match point reconstruction")
    end subroutine test_modal_solution_reconstruction

end program test_core
