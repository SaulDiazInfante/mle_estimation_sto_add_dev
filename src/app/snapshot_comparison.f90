!> @file snapshot_comparison.f90
!! @brief Entry point for reconstructing deterministic and stochastic solution
!! snapshots on the physical grid.
!> @brief Writes a long-form CSV suitable for article figures and comparisons.
program snapshot_comparison_driver
    use csv_output_mod, only: write_solution_snapshot_comparison_csv
    use driver_support_mod, only: assign_default_output_path
    use driver_support_mod, only: build_output_timestamp
    use driver_support_mod, only: default_seed_value
    use driver_support_mod, only: default_snapshot_comparison_name
    use driver_support_mod, only: load_core_runtime_configuration
    use driver_support_mod, only: load_snapshot_configuration
    use driver_support_mod, only: load_snapshot_plot_grid_configuration
    use driver_support_mod, only: normalize_output_timestamp
    use model_types_mod, only: dp
    use model_types_mod, only: sde_parameters_t
    use model_types_mod, only: spatial_grid_t
    use workflow_mod, only: run_snapshot_comparison
    implicit none

    type(spatial_grid_t) :: grid
    type(sde_parameters_t) :: sde_parameters
    integer :: seed_value
    integer :: plot_nx
    integer :: plot_ny
    real(dp), allocatable :: deterministic_fields(:, :, :)
    real(dp), allocatable :: snapshot_times(:)
    real(dp), allocatable :: stochastic_fields(:, :, :)
    real(dp), allocatable :: x_coordinates(:)
    real(dp), allocatable :: y_coordinates(:)
    character(len=:), allocatable :: output_timestamp
    character(len=:), allocatable :: snapshot_file

    seed_value = default_seed_value
    output_timestamp = build_output_timestamp()

    call load_core_runtime_configuration(&
        grid, sde_parameters, seed_value, output_timestamp &
    )
    call load_snapshot_configuration(&
        sde_parameters%time_step, sde_parameters%n_observations, snapshot_times, &
        snapshot_file &
    )
    call load_snapshot_plot_grid_configuration(grid, plot_nx, plot_ny)
    call normalize_output_timestamp(output_timestamp)
    call assign_default_output_path(&
        output_timestamp, default_snapshot_comparison_name, snapshot_file &
    )

    call run_snapshot_comparison(&
        grid, sde_parameters, seed_value, snapshot_times, x_coordinates, &
        y_coordinates, deterministic_fields, stochastic_fields, plot_nx, plot_ny, &
        report_progress=.true. &
    )
    call write_solution_snapshot_comparison_csv(&
        snapshot_file, grid, x_coordinates, y_coordinates, snapshot_times, &
        deterministic_fields, stochastic_fields &
    )

    write (*, '(2a)') "Wrote solution snapshot comparison to ", trim(snapshot_file)
end program snapshot_comparison_driver
