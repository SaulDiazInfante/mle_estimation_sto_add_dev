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
    use spectral_operators_mod, only: assemble_problem_operators
    use validation_mod, only: ensure_finite
    implicit none
    private

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

end module workflow_mod
