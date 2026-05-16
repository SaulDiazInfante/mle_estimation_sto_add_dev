!> @file workflow_mod.f90
!! @brief Shared single-path workflow used by the application entry points.
!> @brief Runs operator assembly, simulation, validation, and final MLE.
module workflow_mod
    use driver_support_mod, only: read_wall_time_seconds
    use model_types_mod, only: dp
    use model_types_mod, only: parameter_estimates_t
    use model_types_mod, only: sde_parameters_t
    use model_types_mod, only: spatial_grid_t
    use model_types_mod, only: spectral_operator_set_t
    use parameter_ml_estimation_mod, only: estimate_model_parameters
    use sde_simulation_mod, only: set_random_seed
    use sde_simulation_mod, only: simulate_state_history
    use sde_simulation_mod, only: simulate_state_snapshots
    use spectral_operators_mod, only: assemble_problem_operators
    use spectral_operators_mod, only: build_cell_center_coordinates
    use spectral_operators_mod, only: reconstruct_field_from_modes
    use validation_mod, only: ensure_finite
    implicit none
    private

    public :: run_snapshot_comparison
    public :: run_simulation_and_estimation

contains

    !> Executes the maintained single-trajectory simulation and MLE workflow.
    subroutine run_simulation_and_estimation(&
        grid, parameters, seed_value, operators, state_history, estimates, &
        report_progress &
    )
        type(spatial_grid_t), intent(in) :: grid
        type(sde_parameters_t), intent(in) :: parameters
        integer, intent(in) :: seed_value
        type(spectral_operator_set_t), intent(out) :: operators
        real(dp), allocatable, intent(out) :: state_history(:, :)
        type(parameter_estimates_t), intent(out) :: estimates
        logical, intent(in), optional :: report_progress

        real(dp) :: finish_time
        logical :: report_mle_progress
        real(dp) :: start_time

        report_mle_progress = .false.
        if (present(report_progress)) report_mle_progress = report_progress

        start_time = read_wall_time_seconds()
        call assemble_problem_operators(grid, operators)
        call ensure_finite("initial_state", operators%initial_state)
        call ensure_finite("eigenvalues", operators%eigenvalues)
        call ensure_finite("diffusion_diagonal", operators%diffusion_diagonal)
        call ensure_finite("interaction_matrix", operators%interaction_matrix)

        call set_random_seed(seed_value)
        call simulate_state_history(operators, parameters, state_history)
        call ensure_finite("state_history", state_history)
        finish_time = read_wall_time_seconds()
        estimates%setup_time = finish_time - start_time

        start_time = read_wall_time_seconds()
        call estimate_model_parameters(&
            state_history, parameters%time_step, operators%interaction_matrix, &
            operators%eigenvalues, grid%gamma, estimates%sigma_hat, &
            estimates%beta_hat, estimates%theta_hat, &
            report_progress=report_mle_progress &
        )
        finish_time = read_wall_time_seconds()
        estimates%estimation_time = finish_time - start_time
    end subroutine run_simulation_and_estimation

    !> Simulates stochastic and deterministic snapshot paths and reconstructs
    !! them on the physical grid.
    subroutine run_snapshot_comparison(&
        grid, parameters, seed_value, snapshot_times, x_coordinates, &
        y_coordinates, deterministic_fields, stochastic_fields, &
        report_progress &
    )
        type(spatial_grid_t), intent(in) :: grid
        type(sde_parameters_t), intent(in) :: parameters
        integer, intent(in) :: seed_value
        real(dp), intent(in) :: snapshot_times(:)
        real(dp), allocatable, intent(out) :: x_coordinates(:)
        real(dp), allocatable, intent(out) :: y_coordinates(:)
        real(dp), allocatable, intent(out) :: deterministic_fields(:, :, :)
        real(dp), allocatable, intent(out) :: stochastic_fields(:, :, :)
        logical, intent(in), optional :: report_progress

        integer, allocatable :: checkpoints(:)
        type(sde_parameters_t) :: deterministic_parameters
        real(dp), allocatable :: reconstructed_field(:, :)
        real(dp), allocatable :: deterministic_snapshots(:, :)
        real(dp), allocatable :: stochastic_snapshots(:, :)
        type(spectral_operator_set_t) :: operators
        logical :: progress_enabled
        integer :: snapshot_index

        progress_enabled = .false.
        if (present(report_progress)) progress_enabled = report_progress

        call assemble_problem_operators(grid, operators)
        call ensure_finite("initial_state", operators%initial_state)
        call ensure_finite("eigenvalues", operators%eigenvalues)
        call ensure_finite("diffusion_diagonal", operators%diffusion_diagonal)
        call ensure_finite("interaction_matrix", operators%interaction_matrix)

        call build_snapshot_checkpoints(&
            parameters%time_step, parameters%n_observations, snapshot_times, &
            checkpoints &
        )

        call set_random_seed(seed_value)
        call simulate_state_snapshots(&
            operators, parameters, checkpoints, stochastic_snapshots, &
            report_progress=progress_enabled, &
            progress_label="Stochastic snapshot path" &
        )
        call ensure_finite("stochastic_snapshots", stochastic_snapshots)

        deterministic_parameters = parameters
        deterministic_parameters%sigma = 0.0_dp
        call simulate_state_snapshots(&
            operators, deterministic_parameters, checkpoints, &
            deterministic_snapshots, report_progress=progress_enabled, &
            progress_label="Deterministic snapshot path" &
        )
        call ensure_finite("deterministic_snapshots", deterministic_snapshots)

        call build_cell_center_coordinates(grid, x_coordinates, y_coordinates)
        allocate (deterministic_fields(grid%nx, grid%ny, size(snapshot_times)))
        allocate (stochastic_fields(grid%nx, grid%ny, size(snapshot_times)))

        do snapshot_index = 1, size(snapshot_times)
            call reconstruct_field_from_modes(&
                grid, operators%mode_pairs, &
                deterministic_snapshots(snapshot_index, :), &
                reconstructed_field &
            )
            deterministic_fields(:, :, snapshot_index) = reconstructed_field
            call ensure_finite("deterministic_field", reconstructed_field)

            call reconstruct_field_from_modes(&
                grid, operators%mode_pairs, &
                stochastic_snapshots(snapshot_index, :), reconstructed_field &
            )
            stochastic_fields(:, :, snapshot_index) = reconstructed_field
            call ensure_finite("stochastic_field", reconstructed_field)
        end do
    end subroutine run_snapshot_comparison

    !> Converts requested physical times into stored observation checkpoints.
    subroutine build_snapshot_checkpoints(&
        time_step, n_observations, snapshot_times, checkpoints &
    )
        integer, intent(in) :: n_observations
        real(dp), intent(in) :: time_step
        real(dp), intent(in) :: snapshot_times(:)
        integer, allocatable, intent(out) :: checkpoints(:)

        integer :: checkpoint_index
        integer :: nearest_step
        integer :: observation_index
        real(dp) :: aligned_time
        real(dp) :: tolerance

        if (time_step <= 0.0_dp) then
            write (*, '(a)') "time_step must be positive for snapshots."
            error stop
        end if

        if (n_observations < 1) then
            write (*, '(a)') "n_observations must be at least 1 for snapshots."
            error stop
        end if

        allocate (checkpoints(size(snapshot_times)))

        do checkpoint_index = 1, size(snapshot_times)
            nearest_step = nint(snapshot_times(checkpoint_index) / time_step)
            aligned_time = real(nearest_step, dp) * time_step
            tolerance = max(&
                1000.0_dp * epsilon(1.0_dp) * &
                max(1.0_dp, abs(snapshot_times(checkpoint_index))), &
                1.0e-12_dp &
            )

            if (abs(aligned_time - snapshot_times(checkpoint_index)) > tolerance) then
                write (*, '(a,1x,es12.5e3,a,1x,es12.5e3,a)') &
                    "Snapshot time", snapshot_times(checkpoint_index), &
                    "is not aligned with time_step", time_step, "."
                error stop
            end if

            observation_index = nearest_step + 1
            if (observation_index > n_observations) then
                write (*, '(a,1x,es12.5e3,a,1x,es12.5e3,a)') &
                    "Snapshot time", snapshot_times(checkpoint_index), &
                    "exceeds the simulated time horizon for time_step", &
                    time_step, "."
                error stop
            end if

            checkpoints(checkpoint_index) = observation_index
        end do
    end subroutine build_snapshot_checkpoints

end module workflow_mod
