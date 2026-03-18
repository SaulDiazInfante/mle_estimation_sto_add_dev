module csv_output_mod
    use model_types_mod, only: dp
    implicit none
    private

    public :: write_estimator_history_csv
    public :: write_state_history_csv

contains

    subroutine write_state_history_csv(file_name, state_history)
        character(len=*), intent(in) :: file_name
        real(dp), intent(in) :: state_history(:, :)

        integer :: unit_id
        integer :: column_index
        integer :: row_index

        open (newunit=unit_id, file=file_name, status="replace", action="write")

        write (unit_id, '(a)', advance='no') "observation"
        do column_index = 1, size(state_history, 2)
            write (unit_id, '(",",a,i0)', advance='no') &
                "state_", column_index
        end do
        write (unit_id, *)

        do row_index = 1, size(state_history, 1)
            write (unit_id, '(i0)', advance='no') row_index
            do column_index = 1, size(state_history, 2)
                write (unit_id, '(",",es26.17e3)', advance='no') &
                    state_history(row_index, column_index)
            end do
            write (unit_id, *)
        end do

        close (unit_id)
    end subroutine write_state_history_csv

    subroutine write_estimator_history_csv(&
        file_name, observation_counts, times, sigma_history, beta_history, &
        theta_history, sigma_true, beta_true, theta_true &
    )
        character(len=*), intent(in) :: file_name
        integer, intent(in) :: observation_counts(:)
        real(dp), intent(in) :: times(:)
        real(dp), intent(in) :: sigma_history(:)
        real(dp), intent(in) :: beta_history(:)
        real(dp), intent(in) :: theta_history(:)
        real(dp), intent(in) :: sigma_true
        real(dp), intent(in) :: beta_true
        real(dp), intent(in) :: theta_true

        integer :: point_index
        integer :: unit_id

        if (size(observation_counts) /= size(times)) then
            write (*, '(a)') "Observation counts and times must match"
            error stop
        end if

        if (size(times) /= size(sigma_history)) then
            write (*, '(a)') "Sigma history length does not match times"
            error stop
        end if

        if (size(beta_history) /= size(theta_history)) then
            write (*, '(a)') "Beta and theta histories must match"
            error stop
        end if

        if (size(times) /= size(beta_history)) then
            write (*, '(a)') "Estimator histories must have the same length"
            error stop
        end if

        open (newunit=unit_id, file=file_name, status="replace", action="write")
        write (unit_id, '(a)') &
            "n_obs,time,sigma_hat,sigma_true,beta_hat,beta_true,"// &
            "theta_hat,theta_true"

        do point_index = 1, size(observation_counts)
            write (unit_id, '(i0,",",es26.17e3,",",es26.17e3,",",'// &
                'es26.17e3,",",es26.17e3,",",es26.17e3,",",'// &
                'es26.17e3,",",es26.17e3)') &
                observation_counts(point_index), &
                times(point_index), &
                sigma_history(point_index), sigma_true, &
                beta_history(point_index), beta_true, &
                theta_history(point_index), theta_true
        end do

        close (unit_id)
    end subroutine write_estimator_history_csv

end module csv_output_mod
