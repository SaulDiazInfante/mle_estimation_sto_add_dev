!> @file progress_reporting_mod.f90
!! @brief Console progress reporting utilities for long-running loops.
!> @brief Maintains and redraws a single-line progress indicator with ETA.
module progress_reporting_mod
    use iso_fortran_env, only: int64
    use iso_fortran_env, only: output_unit
    use iso_fortran_env, only: real64
    implicit none
    private

    integer, parameter :: progress_bar_width = 24

    !> Mutable state used to report progress for a single loop or task.
    type, public :: progress_tracker_t
        logical :: active_line = .false.              !!< Tracks whether an in-place line is active.
        character(len=:), allocatable :: detail       !!< Optional case-specific details appended to the line.
        logical :: enabled = .false.                  !!< Enables or disables all output.
        integer :: last_completed_work = 0            !!< Most recent completed-work value seen by the tracker.
        integer :: last_line_length = 0               !!< Width of the most recent rendered line.
        integer :: next_report = 0                    !!< Next completed-work threshold to report.
        integer :: report_interval = 0                !!< Spacing between progress updates.
        integer :: total_work = 0                     !!< Total number of work units in the task.
        real(real64) :: start_time = 0.0_real64       !!< Wall time recorded at tracker initialization.
        character(len=:), allocatable :: label        !!< Human-readable task label.
    end type progress_tracker_t

    public :: finalize_progress_tracker
    public :: initialize_progress_tracker
    public :: set_progress_tracker_detail
    public :: update_progress_tracker

contains

    !> Initializes a tracker and optionally prints an initial 0% progress line.
    subroutine initialize_progress_tracker(&
        tracker, label, total_work, report_count, enabled, detail &
    )
        type(progress_tracker_t), intent(out) :: tracker
        character(len=*), intent(in) :: label
        integer, intent(in) :: total_work
        integer, intent(in) :: report_count
        logical, intent(in) :: enabled
        character(len=*), intent(in), optional :: detail

        tracker%enabled = enabled .and. total_work > 0
        if (present(detail)) tracker%detail = trim(detail)
        tracker%label = trim(label)
        tracker%total_work = max(0, total_work)

        if (.not. tracker%enabled) return

        tracker%report_interval = max( &
            1, ceiling_division(tracker%total_work, max(1, report_count)) &
        )
        tracker%last_completed_work = 0
        tracker%next_report = tracker%report_interval

        tracker%start_time = read_wall_time_seconds()
        call write_progress_line(tracker, 0)
    end subroutine initialize_progress_tracker

    !> Updates the tracker when the number of completed work units increases.
    subroutine update_progress_tracker(tracker, completed_work)
        type(progress_tracker_t), intent(inout) :: tracker
        integer, intent(in) :: completed_work

        integer :: bounded_work

        if (.not. tracker%enabled) return

        bounded_work = min(completed_work, tracker%total_work)
        tracker%last_completed_work = bounded_work
        if (bounded_work < tracker%next_report .and. &
            bounded_work < tracker%total_work) then
            return
        end if

        call write_progress_line(tracker, bounded_work)

        do while (tracker%next_report <= bounded_work .and. &
            tracker%next_report < tracker%total_work)
            tracker%next_report = tracker%next_report + tracker%report_interval
        end do

        if (bounded_work >= tracker%total_work) then
            tracker%next_report = tracker%total_work + 1
        end if
    end subroutine update_progress_tracker

    !> Finalizes the tracker and terminates the single-line progress output.
    subroutine finalize_progress_tracker(tracker)
        type(progress_tracker_t), intent(inout) :: tracker

        if (.not. tracker%enabled) return
        if (tracker%next_report <= tracker%total_work) then
            call update_progress_tracker(tracker, tracker%total_work)
        end if

        if (tracker%active_line) then
            write (output_unit, '()')
            flush (output_unit)
            tracker%active_line = .false.
            tracker%last_line_length = 0
        end if
    end subroutine finalize_progress_tracker

    !> Sets optional task details and optionally refreshes the current line.
    subroutine set_progress_tracker_detail(tracker, detail, refresh)
        type(progress_tracker_t), intent(inout) :: tracker
        character(len=*), intent(in) :: detail
        logical, intent(in), optional :: refresh

        logical :: refresh_line

        tracker%detail = trim(detail)
        if (.not. tracker%enabled) return

        refresh_line = .false.
        if (present(refresh)) refresh_line = refresh
        if (refresh_line) then
            call write_progress_line(tracker, tracker%last_completed_work)
        end if
    end subroutine set_progress_tracker_detail

    pure integer function ceiling_division(numerator, denominator) result(value)
        integer, intent(in) :: denominator
        integer, intent(in) :: numerator

        value = (numerator + denominator - 1) / denominator
    end function ceiling_division

    real(real64) function read_wall_time_seconds() result(seconds)
        integer(int64) :: clock_count
        integer(int64) :: clock_rate

        call system_clock(clock_count, clock_rate)
        if (clock_rate <= 0_int64) then
            seconds = 0.0_real64
            return
        end if

        seconds = real(clock_count, real64) / real(clock_rate, real64)
    end function read_wall_time_seconds

    subroutine write_progress_line(tracker, completed_work)
        type(progress_tracker_t), intent(inout) :: tracker
        integer, intent(in) :: completed_work

        character(len=320) :: base_line
        character(len=640) :: progress_line
        character(len=progress_bar_width) :: progress_bar
        character(len=32) :: eta_display
        real(real64) :: completion_fraction
        real(real64) :: current_time
        real(real64) :: elapsed_seconds
        real(real64) :: eta_seconds
        integer :: current_line_length
        integer :: filled_segments
        integer :: padding_length
        integer :: eta_hours, eta_minutes, eta_secs, eta_total

        current_time = read_wall_time_seconds()
        elapsed_seconds = current_time - tracker%start_time

        completion_fraction = 0.0_real64
        eta_seconds = 0.0_real64

        if (tracker%total_work > 0 .and. completed_work > 0) then
            completion_fraction = real(completed_work, real64) / &
                real(tracker%total_work, real64)

            if (completion_fraction > 0.0_real64 .and. &
                completion_fraction < 1.0_real64) then
                eta_seconds = elapsed_seconds * &
                    (1.0_real64 - completion_fraction) / completion_fraction
            end if
        end if

        ! Format ETA with hours, minutes, seconds
        eta_total = int(eta_seconds)
        eta_hours = eta_total / 3600
        eta_minutes = (eta_total - eta_hours * 3600) / 60
        eta_secs = eta_total - eta_hours * 3600 - eta_minutes * 60
        
        if (eta_hours > 0) then
            write (eta_display, '(i0,"h ",i0,"m ",i0,"s")') eta_hours, eta_minutes, eta_secs
        else if (eta_minutes > 0) then
            write (eta_display, '(i0,"m ",i0,"s")') eta_minutes, eta_secs
        else
            write (eta_display, '(i0,"s")') eta_secs
        end if

        filled_segments = int(completion_fraction * real(progress_bar_width, real64))
        filled_segments = max(0, min(progress_bar_width, filled_segments))
        progress_bar = repeat('#', filled_segments) // &
            repeat('-', progress_bar_width - filled_segments)

        write (base_line, '(a,1x,"[",a,"]",1x,f6.2,a,2x,'// &
            '"(",i0,"/",i0,")",2x,"eta ",a)') &
            trim(tracker%label), progress_bar, &
            100.0_real64 * completion_fraction, "%", completed_work, &
            tracker%total_work, trim(adjustl(eta_display))

        progress_line = trim(base_line)
        if (allocated(tracker%detail)) then
            if (len_trim(tracker%detail) > 0) then
                progress_line = trim(progress_line)//"  |  "// &
                    trim(tracker%detail)
            end if
        end if

        current_line_length = len_trim(progress_line)
        padding_length = max(0, tracker%last_line_length - current_line_length)

        write (output_unit, '(3a)', advance='no') achar(13), &
            progress_line(1:current_line_length), repeat(' ', padding_length)
        flush (output_unit)
        tracker%active_line = .true.
        tracker%last_line_length = current_line_length
    end subroutine write_progress_line

end module progress_reporting_mod
