!> @file main.f90
!! @brief Command-line entry point for the maintained simulation and estimation 
!! workflow.
!> @brief Runs operator assembly, simulation, estimation, and CSV export.
program max_likelihood_driver
    use csv_output_mod, only: write_estimator_history_csv
    use csv_output_mod, only: write_state_history_csv
    use driver_support_mod, only: assign_default_output_path
    use driver_support_mod, only: build_output_timestamp
    use driver_support_mod, only: default_estimator_history_name
    use driver_support_mod, only: default_minimum_trajectory_observations
    use driver_support_mod, only: default_requested_trajectory_points
    use driver_support_mod, only: default_seed_value
    use driver_support_mod, only: default_state_history_name
    use driver_support_mod, only: load_runtime_configuration
    use driver_support_mod, only: normalize_output_timestamp
    use driver_support_mod, only: read_wall_time_seconds
    use model_types_mod, only: dp
    use model_types_mod, only: parameter_estimates_t
    use model_types_mod, only: sde_parameters_t
    use model_types_mod, only: spatial_grid_t
    use model_types_mod, only: spectral_operator_set_t
    use parameter_ml_estimation_mod, only: build_uniform_checkpoints
    use parameter_ml_estimation_mod, only: estimate_model_parameters
    use parameter_ml_estimation_mod, only: estimate_parameter_history
    use parameter_ml_estimation_mod, only: print_estimation_report
    use sde_simulation_mod, only: set_random_seed
    use sde_simulation_mod, only: simulate_state_history
    use spectral_operators_mod, only: assemble_problem_operators
    use validation_mod, only: ensure_finite
    implicit none

    integer, allocatable :: checkpoints(:)
    integer :: minimum_trajectory_observations
    integer :: requested_trajectory_points
    integer :: seed_value
    logical :: write_state_history
    type(spatial_grid_t) :: grid
    type(sde_parameters_t) :: sde_parameters
    type(spectral_operator_set_t) :: operators
    type(parameter_estimates_t) :: estimates
    real(dp), allocatable :: beta_history(:)
    real(dp), allocatable :: sigma_history(:)
    real(dp), allocatable :: state_history(:, :)
    real(dp), allocatable :: theta_history(:)
    real(dp), allocatable :: times(:)
    real(dp) :: finish_time
    real(dp) :: start_time
    character(len=:), allocatable :: estimator_history_file
    character(len=:), allocatable :: output_timestamp
    character(len=:), allocatable :: state_history_file

    minimum_trajectory_observations = default_minimum_trajectory_observations
    requested_trajectory_points = default_requested_trajectory_points
    seed_value = default_seed_value
    write_state_history = .true.
    output_timestamp = build_output_timestamp()

    call load_runtime_configuration(&
        sde_parameters, minimum_trajectory_observations, &
        requested_trajectory_points, seed_value, write_state_history, &
        output_timestamp, state_history_file, estimator_history_file &
    )
    call normalize_output_timestamp(output_timestamp)
    call assign_default_output_path(&
        output_timestamp, default_state_history_name, state_history_file &
    )
    call assign_default_output_path(&
        output_timestamp, default_estimator_history_name, estimator_history_file &
    )

    start_time = read_wall_time_seconds()
    call assemble_problem_operators(grid, operators)
    call ensure_finite("initial_state", operators%initial_state)
    call ensure_finite("eigenvalues", operators%eigenvalues)
    call ensure_finite("diffusion_diagonal", operators%diffusion_diagonal)
    call ensure_finite("interaction_matrix", operators%interaction_matrix)

    call set_random_seed(seed_value)
    call simulate_state_history(operators, sde_parameters, state_history)
    call ensure_finite("state_history", state_history)
    if (write_state_history) then
        call write_state_history_csv(state_history_file, state_history)
        write (*, '(2a)') "Wrote state history to ", trim(state_history_file)
    end if
    finish_time = read_wall_time_seconds()
    estimates%setup_time = finish_time - start_time

    start_time = read_wall_time_seconds()
    call estimate_model_parameters(&
        state_history, sde_parameters%time_step, &
        operators%interaction_matrix, operators%eigenvalues, &
        grid%gamma, estimates%sigma_hat, estimates%beta_hat, &
        estimates%theta_hat, report_progress=.true. &
    )
    finish_time = read_wall_time_seconds()
    estimates%estimation_time = finish_time - start_time

    call build_uniform_checkpoints(&
        sde_parameters%n_observations, requested_trajectory_points, &
        minimum_trajectory_observations, checkpoints &
    )
    call estimate_parameter_history(&
        state_history, sde_parameters%time_step, &
        operators%interaction_matrix, operators%eigenvalues, &
        grid%gamma, checkpoints, times, sigma_history, beta_history, &
        theta_history, report_progress=.true. &
    )
    call write_estimator_history_csv(&
        estimator_history_file, checkpoints, times, sigma_history, &
        beta_history, theta_history, sde_parameters%sigma, &
        sde_parameters%beta, sde_parameters%theta &
    )
    write (*, '(2a)') "Wrote estimator history to ", trim(estimator_history_file)

    call print_estimation_report(&
        sde_parameters%sigma, sde_parameters%beta, &
        sde_parameters%theta, estimates &
    )

end program max_likelihood_driver
