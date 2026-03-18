!> @file main.f90
!! @brief Command-line entry point for the maintained simulation and estimation workflow.
!> @brief Runs operator assembly, simulation, estimation, and CSV export.
program max_likelihood_driver
    use csv_output_mod, only: write_estimator_history_csv
    use csv_output_mod, only: write_state_history_csv
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

    character(len=*), parameter :: default_output_directory = "data/output"
    character(len=*), parameter :: default_estimator_history_name = &
        "estimator_trajectory.csv"
    character(len=*), parameter :: default_state_history_name = &
        "state_history.csv"
    integer, parameter :: default_minimum_trajectory_observations = 10000
    integer, parameter :: default_requested_trajectory_points = 100
    integer, parameter :: default_seed_value = 42

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

    call cpu_time(start_time)
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
    call cpu_time(finish_time)
    estimates%setup_time = finish_time - start_time

    call cpu_time(start_time)
    call estimate_model_parameters(&
        state_history, sde_parameters%time_step, &
        operators%interaction_matrix, operators%eigenvalues, &
        grid%gamma, estimates%sigma_hat, estimates%beta_hat, &
        estimates%theta_hat &
    )
    call cpu_time(finish_time)
    estimates%estimation_time = finish_time - start_time

    call build_uniform_checkpoints(&
        sde_parameters%n_observations, requested_trajectory_points, &
        minimum_trajectory_observations, checkpoints &
    )
    call estimate_parameter_history(&
        state_history, sde_parameters%time_step, &
        operators%interaction_matrix, operators%eigenvalues, &
        grid%gamma, checkpoints, times, sigma_history, beta_history, &
        theta_history &
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

contains

    subroutine load_runtime_configuration(&
        parameters, minimum_points, requested_points, seed, &
        write_state_history, output_timestamp, state_history_path, &
        estimator_history_path &
    )
        type(sde_parameters_t), intent(inout) :: parameters
        integer, intent(inout) :: minimum_points
        integer, intent(inout) :: requested_points
        integer, intent(inout) :: seed
        logical, intent(inout) :: write_state_history
        character(len=:), allocatable, intent(inout) :: output_timestamp
        character(len=:), allocatable, intent(inout) :: state_history_path
        character(len=:), allocatable, intent(inout) :: estimator_history_path

        call read_integer_env("SARGAZO_N_OBSERVATIONS", parameters%n_observations)
        call read_real_env("SARGAZO_TIME_STEP", parameters%time_step)
        call read_real_env("SARGAZO_BETA", parameters%beta)
        call read_real_env("SARGAZO_THETA", parameters%theta)
        call read_real_env("SARGAZO_SIGMA", parameters%sigma)
        call read_integer_env("SARGAZO_REQUESTED_TRAJECTORY_POINTS", requested_points)
        call read_integer_env(&
            "SARGAZO_MINIMUM_TRAJECTORY_OBSERVATIONS", minimum_points &
        )
        call read_integer_env("SARGAZO_SEED", seed)
        call read_logical_env("SARGAZO_WRITE_STATE_HISTORY", write_state_history)
        call read_string_env("SARGAZO_OUTPUT_TIMESTAMP", output_timestamp)
        call read_string_env("SARGAZO_STATE_HISTORY_FILE", state_history_path)
        call read_string_env("SARGAZO_ESTIMATOR_HISTORY_FILE", estimator_history_path)

        if (minimum_points > parameters%n_observations) then
            minimum_points = max(2, parameters%n_observations)
        end if
    end subroutine load_runtime_configuration

    subroutine assign_default_output_path(timestamp, base_name, output_path)
        character(len=*), intent(in) :: timestamp
        character(len=*), intent(in) :: base_name
        character(len=:), allocatable, intent(inout) :: output_path

        if (allocated(output_path)) return

        output_path = default_output_directory//"/"//trim(timestamp)//"_"// &
            trim(base_name)
    end subroutine assign_default_output_path

    subroutine read_integer_env(name, value)
        character(len=*), intent(in) :: name
        integer, intent(inout) :: value

        character(len=256) :: buffer
        integer :: ios
        integer :: length
        integer :: parsed_value
        integer :: status

        call get_environment_variable(name, buffer, length=length, status=status)
        if (status /= 0 .or. length == 0) return

        read (buffer(1:length), *, iostat=ios) parsed_value
        if (ios /= 0) then
            write (*, '(3a)') "Invalid integer value for ", trim(name), "."
            error stop
        end if

        value = parsed_value
    end subroutine read_integer_env

    subroutine read_real_env(name, value)
        character(len=*), intent(in) :: name
        real(dp), intent(inout) :: value

        character(len=256) :: buffer
        integer :: ios
        integer :: length
        real(dp) :: parsed_value
        integer :: status

        call get_environment_variable(name, buffer, length=length, status=status)
        if (status /= 0 .or. length == 0) return

        read (buffer(1:length), *, iostat=ios) parsed_value
        if (ios /= 0) then
            write (*, '(3a)') "Invalid real value for ", trim(name), "."
            error stop
        end if

        value = parsed_value
    end subroutine read_real_env

    subroutine read_logical_env(name, value)
        character(len=*), intent(in) :: name
        logical, intent(inout) :: value

        character(len=256) :: buffer
        integer :: length
        character(len=:), allocatable :: normalized
        integer :: status

        call get_environment_variable(name, buffer, length=length, status=status)
        if (status /= 0 .or. length == 0) return

        normalized = to_lower(buffer(1:length))

        select case (trim(normalized))
        case ("1", "true", "yes", "on")
            value = .true.
        case ("0", "false", "no", "off")
            value = .false.
        case default
            write (*, '(3a)') "Invalid logical value for ", trim(name), "."
            error stop
        end select
    end subroutine read_logical_env

    subroutine read_string_env(name, value)
        character(len=*), intent(in) :: name
        character(len=:), allocatable, intent(inout) :: value

        character(len=512) :: buffer
        integer :: length
        integer :: status

        call get_environment_variable(name, buffer, length=length, status=status)
        if (status /= 0 .or. length == 0) return

        value = trim(buffer(1:length))
    end subroutine read_string_env

    function build_output_timestamp() result(timestamp)
        character(len=:), allocatable :: timestamp

        integer :: values(8)

        call date_and_time(values=values)
        allocate (character(len=15) :: timestamp)
        write (timestamp, '(i4.4,i2.2,i2.2,"T",i2.2,i2.2,i2.2)') &
            values(1), values(2), values(3), values(5), values(6), values(7)
    end function build_output_timestamp

    subroutine normalize_output_timestamp(timestamp)
        character(len=:), allocatable, intent(inout) :: timestamp

        character(len=:), allocatable :: trimmed

        if (.not. allocated(timestamp)) return

        trimmed = trim(adjustl(timestamp))
        if (len(trimmed) == 0) then
            timestamp = build_output_timestamp()
            return
        end if

        if (len(trimmed) == 19 .and. trimmed(5:5) == "-" .and. &
            trimmed(8:8) == "-" .and. &
            (trimmed(11:11) == "T" .or. trimmed(11:11) == "t") .and. &
            trimmed(14:14) == ":" .and. trimmed(17:17) == ":") then
            timestamp = trimmed(1:4)//trimmed(6:7)//trimmed(9:10)//"T"// &
                trimmed(12:13)//trimmed(15:16)//trimmed(18:19)
            return
        end if

        if (len(trimmed) == 15 .and. &
            (trimmed(9:9) == "T" .or. trimmed(9:9) == "t")) then
            timestamp = trimmed
            timestamp(9:9) = "T"
            return
        end if

        timestamp = sanitize_filename_component(trimmed)
    end subroutine normalize_output_timestamp

    pure function sanitize_filename_component(text) result(sanitized)
        character(len=*), intent(in) :: text
        character(len=len(text)) :: sanitized

        integer :: index

        do index = 1, len(text)
            select case (text(index:index))
            case ('"', "*", "/", ":", "<", ">", "?", "\", "|", " ")
                sanitized(index:index) = "_"
            case default
                sanitized(index:index) = text(index:index)
            end select
        end do
    end function sanitize_filename_component

    pure function to_lower(text) result(lowered)
        character(len=*), intent(in) :: text
        character(len=len(text)) :: lowered

        integer :: ascii_code
        integer :: index

        do index = 1, len(text)
            ascii_code = iachar(text(index:index))
            if (ascii_code >= iachar('A') .and. ascii_code <= iachar('Z')) then
                lowered(index:index) = achar(ascii_code + 32)
            else
                lowered(index:index) = text(index:index)
            end if
        end do
    end function to_lower

end program max_likelihood_driver
