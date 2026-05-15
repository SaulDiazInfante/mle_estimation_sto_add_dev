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
    use model_types_mod, only: dp
    use model_types_mod, only: parameter_estimates_t
    use model_types_mod, only: sde_parameters_t
    use model_types_mod, only: spatial_grid_t
    use model_types_mod, only: spectral_operator_set_t
    use parameter_ml_estimation_mod, only: build_uniform_checkpoints
    use parameter_ml_estimation_mod, only: estimate_parameter_history
    use parameter_ml_estimation_mod, only: print_estimation_report
    use workflow_mod, only: run_simulation_and_estimation
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
    character(len=:), allocatable :: estimator_history_file
    character(len=:), allocatable :: output_timestamp
    character(len=:), allocatable :: state_history_file

    minimum_trajectory_observations = default_minimum_trajectory_observations
    requested_trajectory_points = default_requested_trajectory_points
    seed_value = default_seed_value
    write_state_history = .true.
    output_timestamp = build_output_timestamp()

    call load_runtime_configuration(&
        grid, sde_parameters, minimum_trajectory_observations, &
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

    call run_simulation_and_estimation(&
        grid, sde_parameters, seed_value, operators, state_history, estimates, &
        report_progress=.true. &
    )
    if (write_state_history) then
        call write_state_history_csv(state_history_file, state_history)
        write (*, '(2a)') "Wrote state history to ", trim(state_history_file)
    end if

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
