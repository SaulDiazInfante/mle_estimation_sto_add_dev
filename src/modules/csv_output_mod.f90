!> @file csv_output_mod.f90
!! @brief CSV writers for generated trajectories and parameter estimates.
!> @brief Encapsulates the text output format used by downstream plotting scripts.
module csv_output_mod
    use model_types_mod, only: dp
    use model_types_mod, only: monte_carlo_case_summary_t
    use model_types_mod, only: monte_carlo_replicate_result_t
    implicit none
    private

    public :: write_estimator_history_csv
    public :: write_monte_carlo_replicates_csv
    public :: write_monte_carlo_summary_csv
    public :: write_state_history_csv

contains

    !> Writes the simulated modal state history to a CSV file.
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

    !> Writes time-indexed parameter estimates and their reference values to CSV.
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

    !> Writes replicate-level Monte Carlo estimates for downstream histograms.
    subroutine write_monte_carlo_replicates_csv(file_name, replicate_results)
        character(len=*), intent(in) :: file_name
        type(monte_carlo_replicate_result_t), intent(in) :: replicate_results(:)

        integer :: replicate_index
        integer :: unit_id

        open (newunit=unit_id, file=file_name, status="replace", action="write")
        write (unit_id, '(a)') &
            "basis_level,nx,ny,state_dimension,n_obs,time_step,replicate,"// &
            "sigma_hat,sigma_true,beta_hat,beta_true,theta_hat,theta_true,"// &
            "setup_time,estimation_time,total_time"

        do replicate_index = 1, size(replicate_results)
            write (unit_id, '(i0,",",i0,",",i0,",",i0,",",i0,",",'// &
                'es26.17e3,",",i0,",",es26.17e3,",",es26.17e3,",",'// &
                'es26.17e3,",",es26.17e3,",",es26.17e3,",",es26.17e3,'// &
                '",",es26.17e3,",",es26.17e3,",",es26.17e3)') &
                replicate_results(replicate_index)%basis_level, &
                replicate_results(replicate_index)%nx, &
                replicate_results(replicate_index)%ny, &
                replicate_results(replicate_index)%state_dimension, &
                replicate_results(replicate_index)%n_observations, &
                replicate_results(replicate_index)%time_step, &
                replicate_results(replicate_index)%replicate, &
                replicate_results(replicate_index)%sigma_hat, &
                replicate_results(replicate_index)%sigma_true, &
                replicate_results(replicate_index)%beta_hat, &
                replicate_results(replicate_index)%beta_true, &
                replicate_results(replicate_index)%theta_hat, &
                replicate_results(replicate_index)%theta_true, &
                replicate_results(replicate_index)%setup_time, &
                replicate_results(replicate_index)%estimation_time, &
                replicate_results(replicate_index)%total_time
        end do

        close (unit_id)
    end subroutine write_monte_carlo_replicates_csv

    !> Writes combination-level Monte Carlo summaries with mean and spread.
    subroutine write_monte_carlo_summary_csv(file_name, case_summaries)
        character(len=*), intent(in) :: file_name
        type(monte_carlo_case_summary_t), intent(in) :: case_summaries(:)

        integer :: case_index
        integer :: unit_id

        open (newunit=unit_id, file=file_name, status="replace", action="write")
        write (unit_id, '(a)') &
            "basis_level,nx,ny,state_dimension,n_obs,time_step,n_replicates,"// &
            "sigma_mean,sigma_sd,sigma_true,beta_mean,beta_sd,beta_true,"// &
            "theta_mean,theta_sd,theta_true,average_setup_time,"// &
            "average_estimation_time,average_total_time"

        do case_index = 1, size(case_summaries)
            write (unit_id, '(i0,",",i0,",",i0,",",i0,",",i0,",",'// &
                'es26.17e3,",",i0,",",es26.17e3,",",es26.17e3,",",'// &
                'es26.17e3,",",es26.17e3,",",es26.17e3,",",es26.17e3,'// &
                '",",es26.17e3,",",es26.17e3,",",es26.17e3,",",'// &
                'es26.17e3,",",es26.17e3,",",es26.17e3)') &
                case_summaries(case_index)%basis_level, &
                case_summaries(case_index)%nx, &
                case_summaries(case_index)%ny, &
                case_summaries(case_index)%state_dimension, &
                case_summaries(case_index)%n_observations, &
                case_summaries(case_index)%time_step, &
                case_summaries(case_index)%n_replicates, &
                case_summaries(case_index)%sigma_mean, &
                case_summaries(case_index)%sigma_sd, &
                case_summaries(case_index)%sigma_true, &
                case_summaries(case_index)%beta_mean, &
                case_summaries(case_index)%beta_sd, &
                case_summaries(case_index)%beta_true, &
                case_summaries(case_index)%theta_mean, &
                case_summaries(case_index)%theta_sd, &
                case_summaries(case_index)%theta_true, &
                case_summaries(case_index)%average_setup_time, &
                case_summaries(case_index)%average_estimation_time, &
                case_summaries(case_index)%average_total_time
        end do

        close (unit_id)
    end subroutine write_monte_carlo_summary_csv

end module csv_output_mod
