!> @file monte_carlo_study.f90
!! @brief Entry point for Monte Carlo MLE validation studies.
!> @brief Sweeps basis size, step size, and observation count combinations.
program monte_carlo_study_driver
    use csv_output_mod, only: write_monte_carlo_replicates_csv
    use csv_output_mod, only: write_monte_carlo_summary_csv
    use driver_support_mod, only: assign_default_output_path
    use driver_support_mod, only: build_output_timestamp
    use driver_support_mod, only: default_monte_carlo_replicates_name
    use driver_support_mod, only: default_monte_carlo_summary_name
    use driver_support_mod, only: default_seed_value
    use driver_support_mod, only: load_core_runtime_configuration
    use driver_support_mod, only: load_monte_carlo_configuration
    use driver_support_mod, only: normalize_output_timestamp
    use model_types_mod, only: monte_carlo_case_summary_t
    use model_types_mod, only: monte_carlo_replicate_result_t
    use model_types_mod, only: monte_carlo_study_config_t
    use model_types_mod, only: sde_parameters_t
    use model_types_mod, only: spatial_grid_t
    use monte_carlo_study_mod, only: run_monte_carlo_study
    implicit none

    type(spatial_grid_t) :: grid
    type(sde_parameters_t) :: sde_parameters
    type(monte_carlo_study_config_t) :: study_config
    type(monte_carlo_replicate_result_t), allocatable :: replicate_results(:)
    type(monte_carlo_case_summary_t), allocatable :: case_summaries(:)
    integer :: seed_value
    character(len=:), allocatable :: output_timestamp
    character(len=:), allocatable :: replicate_file
    character(len=:), allocatable :: summary_file

    seed_value = default_seed_value
    output_timestamp = build_output_timestamp()

    call load_core_runtime_configuration(&
        grid, sde_parameters, seed_value, output_timestamp &
    )
    call load_monte_carlo_configuration(&
        grid, sde_parameters, study_config, summary_file, replicate_file &
    )
    call normalize_output_timestamp(output_timestamp)
    call assign_default_output_path(&
        output_timestamp, default_monte_carlo_summary_name, summary_file &
    )
    call assign_default_output_path(&
        output_timestamp, default_monte_carlo_replicates_name, replicate_file &
    )

    call run_monte_carlo_study(&
        grid, sde_parameters, study_config, seed_value, replicate_results, &
        case_summaries, report_progress=.true. &
    )
    call write_monte_carlo_summary_csv(summary_file, case_summaries)
    call write_monte_carlo_replicates_csv(replicate_file, replicate_results)

    write (*, '(2a)') "Wrote Monte Carlo summary to ", trim(summary_file)
    write (*, '(2a)') "Wrote Monte Carlo replicates to ", trim(replicate_file)
end program monte_carlo_study_driver
