!> @file driver_support_mod.f90
!! @brief Helper routines and defaults used by the main application driver.
!> @brief Centralizes runtime configuration, timestamps, and output-path helpers.
module driver_support_mod
    use iso_fortran_env, only: int64
    use model_types_mod, only: dp
    use model_types_mod, only: monte_carlo_study_config_t
    use model_types_mod, only: sde_parameters_t
    use model_types_mod, only: spatial_grid_t
    implicit none
    private

    character(len=*), parameter, public :: default_estimator_history_name = &
        "estimator_trajectory.csv"
    character(len=*), parameter, public :: default_monte_carlo_replicates_name = &
        "monte_carlo_replicates.csv"
    integer, parameter, public :: default_monte_carlo_replicates = 1000
    character(len=*), parameter, public :: default_monte_carlo_summary_name = &
        "monte_carlo_summary.csv"
    integer, parameter, public :: default_minimum_trajectory_observations = 10000
    character(len=*), parameter, public :: default_output_directory = "data/output"
    integer, parameter, public :: default_requested_trajectory_points = 100
    integer, parameter, public :: default_seed_value = 42
    character(len=*), parameter, public :: default_state_history_name = &
        "state_history.csv"

    public :: assign_default_output_path
    public :: build_output_timestamp
    public :: load_core_runtime_configuration
    public :: load_monte_carlo_configuration
    public :: load_runtime_configuration
    public :: normalize_output_timestamp
    public :: read_wall_time_seconds

contains

    subroutine load_core_runtime_configuration(&
        grid, parameters, seed, output_timestamp &
    )
        type(spatial_grid_t), intent(inout) :: grid
        type(sde_parameters_t), intent(inout) :: parameters
        integer, intent(inout) :: seed
        character(len=:), allocatable, intent(inout) :: output_timestamp

        call read_integer_env("SARGAZO_GRID_NX", grid%nx)
        call read_integer_env("SARGAZO_GRID_NY", grid%ny)
        call read_integer_env("SARGAZO_VELOCITY_MODE_X", grid%velocity_mode_x)
        call read_integer_env("SARGAZO_VELOCITY_MODE_Y", grid%velocity_mode_y)
        call read_real_env("SARGAZO_LENGTH_X", grid%length_x)
        call read_real_env("SARGAZO_LENGTH_Y", grid%length_y)
        call read_real_env("SARGAZO_GAMMA", grid%gamma)
        call read_integer_env("SARGAZO_N_OBSERVATIONS", parameters%n_observations)
        call read_real_env("SARGAZO_TIME_STEP", parameters%time_step)
        call read_real_env("SARGAZO_BETA", parameters%beta)
        call read_real_env("SARGAZO_THETA", parameters%theta)
        call read_real_env("SARGAZO_SIGMA", parameters%sigma)
        call read_integer_env("SARGAZO_SEED", seed)
        call read_string_env("SARGAZO_OUTPUT_TIMESTAMP", output_timestamp)

        call validate_grid_configuration(grid)
        call validate_sde_configuration(parameters)
    end subroutine load_core_runtime_configuration

    subroutine load_runtime_configuration(&
        grid, parameters, minimum_points, requested_points, seed, &
        write_state_history, output_timestamp, state_history_path, &
        estimator_history_path &
    )
        type(spatial_grid_t), intent(inout) :: grid
        type(sde_parameters_t), intent(inout) :: parameters
        integer, intent(inout) :: minimum_points
        integer, intent(inout) :: requested_points
        integer, intent(inout) :: seed
        logical, intent(inout) :: write_state_history
        character(len=:), allocatable, intent(inout) :: output_timestamp
        character(len=:), allocatable, intent(inout) :: state_history_path
        character(len=:), allocatable, intent(inout) :: estimator_history_path

        call load_core_runtime_configuration(grid, parameters, seed, output_timestamp)
        call read_integer_env("SARGAZO_REQUESTED_TRAJECTORY_POINTS", requested_points)
        call read_integer_env(&
            "SARGAZO_MINIMUM_TRAJECTORY_OBSERVATIONS", minimum_points &
        )
        call read_logical_env("SARGAZO_WRITE_STATE_HISTORY", write_state_history)
        call read_string_env("SARGAZO_STATE_HISTORY_FILE", state_history_path)
        call read_string_env("SARGAZO_ESTIMATOR_HISTORY_FILE", estimator_history_path)

        if (minimum_points > parameters%n_observations) then
            minimum_points = max(2, parameters%n_observations)
        end if

        call validate_trajectory_configuration(&
            minimum_points, requested_points, parameters%n_observations &
        )
    end subroutine load_runtime_configuration

    subroutine load_monte_carlo_configuration(&
        base_grid, base_parameters, study_config, summary_path, replicate_path &
    )
        type(spatial_grid_t), intent(in) :: base_grid
        type(sde_parameters_t), intent(in) :: base_parameters
        type(monte_carlo_study_config_t), intent(inout) :: study_config
        character(len=:), allocatable, intent(inout) :: summary_path
        character(len=:), allocatable, intent(inout) :: replicate_path

        study_config%n_replicates = default_monte_carlo_replicates
        call set_default_integer_vector( &
            study_config%basis_levels, [base_grid%nx] &
        )
        call set_default_integer_vector( &
            study_config%observation_counts, [base_parameters%n_observations] &
        )
        call set_default_real_vector( &
            study_config%time_steps, [base_parameters%time_step] &
        )

        call read_integer_env("SARGAZO_MC_REPLICATES", study_config%n_replicates)
        call read_integer_vector_env(&
            "SARGAZO_MC_BASIS_LEVELS", study_config%basis_levels &
        )
        call read_integer_vector_env(&
            "SARGAZO_MC_N_OBSERVATIONS", study_config%observation_counts &
        )
        call read_real_vector_env( &
            "SARGAZO_MC_TIME_STEPS", study_config%time_steps &
        )
        call read_string_env("SARGAZO_MC_SUMMARY_FILE", summary_path)
        call read_string_env("SARGAZO_MC_REPLICATE_FILE", replicate_path)

        if (study_config%n_replicates < 1) then
            write (*, '(a)') "SARGAZO_MC_REPLICATES must be positive."
            error stop
        end if

        if (.not. allocated(study_config%basis_levels) .or. &
            size(study_config%basis_levels) < 1 .or. &
            any(study_config%basis_levels < 1)) then
            write (*, '(a)') "SARGAZO_MC_BASIS_LEVELS must contain positive integers."
            error stop
        end if

        if (.not. allocated(study_config%observation_counts) .or. &
            size(study_config%observation_counts) < 1 .or. &
            any(study_config%observation_counts < 2)) then
            write (*, '(a)') "SARGAZO_MC_N_OBSERVATIONS must be at least 2."
            error stop
        end if

        if (.not. allocated(study_config%time_steps) .or. &
            size(study_config%time_steps) < 1 .or. &
            any(study_config%time_steps <= 0.0_dp)) then
            write (*, '(a)') "SARGAZO_MC_TIME_STEPS must contain positive values."
            error stop
        end if
    end subroutine load_monte_carlo_configuration

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

    real(dp) function read_wall_time_seconds()
        integer(int64) :: clock_count
        integer(int64) :: clock_rate

        call system_clock(clock_count, clock_rate)
        if (clock_rate <= 0_int64) then
            read_wall_time_seconds = 0.0_dp
            return
        end if

        read_wall_time_seconds = real(clock_count, dp) / real(clock_rate, dp)
    end function read_wall_time_seconds

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

    subroutine read_integer_vector_env(name, values)
        character(len=*), intent(in) :: name
        integer, allocatable, intent(inout) :: values(:)

        character(len=2048) :: buffer
        integer :: length
        integer :: status

        call get_environment_variable(name, buffer, length=length, status=status)
        if (status /= 0 .or. length == 0) return

        call parse_integer_list(name, trim(buffer(1:length)), values)
    end subroutine read_integer_vector_env

    subroutine read_real_vector_env(name, values)
        character(len=*), intent(in) :: name
        real(dp), allocatable, intent(inout) :: values(:)

        character(len=2048) :: buffer
        integer :: length
        integer :: status

        call get_environment_variable(name, buffer, length=length, status=status)
        if (status /= 0 .or. length == 0) return

        call parse_real_list(name, trim(buffer(1:length)), values)
    end subroutine read_real_vector_env

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

    subroutine validate_grid_configuration(grid)
        type(spatial_grid_t), intent(in) :: grid

        if (grid%nx < 1 .or. grid%ny < 1) then
            write (*, '(a)') "Grid dimensions must be positive."
            error stop
        end if

        if (grid%length_x <= 0.0_dp .or. grid%length_y <= 0.0_dp) then
            write (*, '(a)') "Domain lengths must be positive."
            error stop
        end if

        if (grid%gamma < 0.0_dp) then
            write (*, '(a)') "gamma must be non-negative."
            error stop
        end if
    end subroutine validate_grid_configuration

    subroutine validate_sde_configuration(parameters)
        type(sde_parameters_t), intent(in) :: parameters

        if (parameters%n_observations < 1) then
            write (*, '(a)') "SARGAZO_N_OBSERVATIONS must be at least 1."
            error stop
        end if

        if (parameters%time_step <= 0.0_dp) then
            write (*, '(a)') "SARGAZO_TIME_STEP must be positive."
            error stop
        end if

        if (parameters%sigma < 0.0_dp) then
            write (*, '(a)') "SARGAZO_SIGMA must be non-negative."
            error stop
        end if
    end subroutine validate_sde_configuration

    subroutine validate_trajectory_configuration(&
        minimum_points, requested_points, total_points &
    )
        integer, intent(in) :: minimum_points
        integer, intent(in) :: requested_points
        integer, intent(in) :: total_points

        if (requested_points < 1) then
            write (*, '(a)') "SARGAZO_REQUESTED_TRAJECTORY_POINTS must be positive."
            error stop
        end if

        if (minimum_points < 2 .or. minimum_points > total_points) then
            write (*, '(a)') &
                "SARGAZO_MINIMUM_TRAJECTORY_OBSERVATIONS must belong to [2, n_observations]."
            error stop
        end if
    end subroutine validate_trajectory_configuration

    subroutine set_default_integer_vector(values, defaults)
        integer, allocatable, intent(inout) :: values(:)
        integer, intent(in) :: defaults(:)

        if (allocated(values)) deallocate (values)
        allocate (values(size(defaults)))
        values = defaults
    end subroutine set_default_integer_vector

    subroutine set_default_real_vector(values, defaults)
        real(dp), allocatable, intent(inout) :: values(:)
        real(dp), intent(in) :: defaults(:)

        if (allocated(values)) deallocate (values)
        allocate (values(size(defaults)))
        values = defaults
    end subroutine set_default_real_vector

    subroutine parse_integer_list(name, raw_text, values)
        character(len=*), intent(in) :: name
        character(len=*), intent(in) :: raw_text
        integer, allocatable, intent(inout) :: values(:)

        character(len=:), allocatable :: token
        integer :: ios
        integer :: parsed_value
        integer :: token_count
        integer :: token_index

        token_count = count_csv_tokens(raw_text)
        if (token_count < 1) then
            write (*, '(3a)') "No values were provided for ", trim(name), "."
            error stop
        end if

        if (allocated(values)) deallocate (values)
        allocate (values(token_count))

        do token_index = 1, token_count
            token = extract_csv_token(raw_text, token_index)
            read (token, *, iostat=ios) parsed_value
            if (ios /= 0) then
                write (*, '(3a)') "Invalid integer list for ", trim(name), "."
                error stop
            end if

            values(token_index) = parsed_value
        end do
    end subroutine parse_integer_list

    subroutine parse_real_list(name, raw_text, values)
        character(len=*), intent(in) :: name
        character(len=*), intent(in) :: raw_text
        real(dp), allocatable, intent(inout) :: values(:)

        character(len=:), allocatable :: token
        integer :: ios
        real(dp) :: parsed_value
        integer :: token_count
        integer :: token_index

        token_count = count_csv_tokens(raw_text)
        if (token_count < 1) then
            write (*, '(3a)') "No values were provided for ", trim(name), "."
            error stop
        end if

        if (allocated(values)) deallocate (values)
        allocate (values(token_count))

        do token_index = 1, token_count
            token = extract_csv_token(raw_text, token_index)
            read (token, *, iostat=ios) parsed_value
            if (ios /= 0) then
                write (*, '(3a)') "Invalid real list for ", trim(name), "."
                error stop
            end if

            values(token_index) = parsed_value
        end do
    end subroutine parse_real_list

    pure integer function count_csv_tokens(raw_text) result(token_count)
        character(len=*), intent(in) :: raw_text

        integer :: delimiter_index
        integer :: start_index
        character(len=:), allocatable :: trimmed_text

        trimmed_text = trim(raw_text)
        if (len(trimmed_text) == 0) then
            token_count = 0
            return
        end if

        start_index = 1
        token_count = 0
        do
            call skip_csv_separators(trimmed_text, start_index)
            if (start_index > len(trimmed_text)) exit

            delimiter_index = index(trimmed_text(start_index:), ",")
            token_count = token_count + 1
            if (delimiter_index == 0) exit

            start_index = start_index + delimiter_index
        end do
    end function count_csv_tokens

    pure function extract_csv_token(raw_text, token_index) result(token)
        character(len=*), intent(in) :: raw_text
        integer, intent(in) :: token_index
        character(len=:), allocatable :: token

        integer :: current_index
        integer :: delimiter_index
        integer :: start_index
        integer :: token_counter
        character(len=:), allocatable :: trimmed_text

        trimmed_text = trim(raw_text)
        start_index = 1
        token_counter = 0

        do
            call skip_csv_separators(trimmed_text, start_index)
            if (start_index > len(trimmed_text)) exit

            delimiter_index = index(trimmed_text(start_index:), ",")
            token_counter = token_counter + 1

            if (delimiter_index == 0) then
                if (token_counter == token_index) then
                    token = trim(adjustl(trimmed_text(start_index:)))
                    return
                end if
                exit
            end if

            current_index = start_index + delimiter_index - 2
            if (token_counter == token_index) then
                token = trim(adjustl(trimmed_text(start_index:current_index)))
                return
            end if

            start_index = start_index + delimiter_index
        end do

        token = ""
    end function extract_csv_token

    pure subroutine skip_csv_separators(text, start_index)
        character(len=*), intent(in) :: text
        integer, intent(inout) :: start_index

        do while (start_index <= len(text))
            if (text(start_index:start_index) /= "," .and. &
                text(start_index:start_index) /= " ") then
                exit
            end if

            start_index = start_index + 1
        end do
    end subroutine skip_csv_separators

end module driver_support_mod
