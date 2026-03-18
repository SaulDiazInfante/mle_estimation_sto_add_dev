!> @file driver_support_mod.f90
!! @brief Helper routines and defaults used by the main application driver.
!> @brief Centralizes runtime configuration, timestamps, and output-path helpers.
module driver_support_mod
    use iso_fortran_env, only: int64
    use model_types_mod, only: dp
    use model_types_mod, only: sde_parameters_t
    implicit none
    private

    character(len=*), parameter, public :: default_estimator_history_name = &
        "estimator_trajectory.csv"
    integer, parameter, public :: default_minimum_trajectory_observations = 10000
    character(len=*), parameter, public :: default_output_directory = "data/output"
    integer, parameter, public :: default_requested_trajectory_points = 100
    integer, parameter, public :: default_seed_value = 42
    character(len=*), parameter, public :: default_state_history_name = &
        "state_history.csv"

    public :: assign_default_output_path
    public :: build_output_timestamp
    public :: load_runtime_configuration
    public :: normalize_output_timestamp
    public :: read_wall_time_seconds

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

end module driver_support_mod
