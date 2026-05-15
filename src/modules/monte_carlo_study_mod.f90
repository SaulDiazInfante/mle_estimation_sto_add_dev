!> @file monte_carlo_study_mod.f90
!! @brief Monte Carlo harness for repeated simulation-and-estimation runs.
!> @brief Sweeps basis sizes, step sizes, and observation counts.
module monte_carlo_study_mod
    use model_types_mod, only: dp
    use model_types_mod, only: monte_carlo_case_summary_t
    use model_types_mod, only: monte_carlo_replicate_result_t
    use model_types_mod, only: monte_carlo_study_config_t
    use model_types_mod, only: parameter_estimates_t
    use model_types_mod, only: sde_parameters_t
    use model_types_mod, only: spatial_grid_t
    use model_types_mod, only: spectral_operator_set_t
    use progress_reporting_mod, only: finalize_progress_tracker
    use progress_reporting_mod, only: initialize_progress_tracker
    use progress_reporting_mod, only: progress_tracker_t
    use progress_reporting_mod, only: update_progress_tracker
    use workflow_mod, only: run_simulation_and_estimation
    implicit none
    private

    integer, parameter :: default_study_progress_reports = 20

    public :: compute_monte_carlo_case_count
    public :: run_monte_carlo_study

contains

    !> Returns the number of parameter combinations in the study sweep.
    pure integer function compute_monte_carlo_case_count(study_config) &
        result(case_count)
        type(monte_carlo_study_config_t), intent(in) :: study_config

        if (.not. allocated(study_config%basis_levels) .or. &
            .not. allocated(study_config%observation_counts) .or. &
            .not. allocated(study_config%time_steps)) then
            case_count = 0
            return
        end if

        case_count = size(study_config%basis_levels) * &
            size(study_config%observation_counts) * &
            size(study_config%time_steps)
    end function compute_monte_carlo_case_count

    !> Runs all requested Monte Carlo combinations and aggregates summaries.
    subroutine run_monte_carlo_study(&
        base_grid, base_parameters, study_config, base_seed, replicate_results, &
        case_summaries, report_progress &
    )
        type(spatial_grid_t), intent(in) :: base_grid
        type(sde_parameters_t), intent(in) :: base_parameters
        type(monte_carlo_study_config_t), intent(in) :: study_config
        integer, intent(in) :: base_seed
        type(monte_carlo_replicate_result_t), allocatable, intent(out) :: &
            replicate_results(:)
        type(monte_carlo_case_summary_t), allocatable, intent(out) :: &
            case_summaries(:)
        logical, intent(in), optional :: report_progress

        integer :: basis_index
        integer :: case_count
        integer :: case_index
        integer :: observation_index
        type(spectral_operator_set_t) :: operators
        type(parameter_estimates_t) :: estimates
        type(progress_tracker_t) :: progress_tracker
        logical :: progress_enabled
        integer :: replicate
        integer :: replicate_index
        integer :: step_index
        type(spatial_grid_t) :: study_grid
        type(sde_parameters_t) :: study_parameters
        real(dp), allocatable :: state_history(:, :)
        real(dp) :: beta_mean
        real(dp) :: beta_sum
        real(dp) :: beta_sumsq
        real(dp) :: setup_sum
        real(dp) :: sigma_mean
        real(dp) :: sigma_sum
        real(dp) :: sigma_sumsq
        real(dp) :: theta_mean
        real(dp) :: theta_sum
        real(dp) :: theta_sumsq
        real(dp) :: total_sum
        real(dp) :: estimation_sum
        integer :: total_replicates

        case_count = compute_monte_carlo_case_count(study_config)
        if (case_count < 1) then
            write (*, '(a)') "The Monte Carlo study requires at least one case."
            error stop
        end if

        if (study_config%n_replicates < 1) then
            write (*, '(a)') "The Monte Carlo study requires at least one replicate."
            error stop
        end if

        total_replicates = case_count * study_config%n_replicates
        allocate (replicate_results(total_replicates))
        allocate (case_summaries(case_count))

        progress_enabled = .false.
        if (present(report_progress)) progress_enabled = report_progress
        call initialize_progress_tracker(&
            progress_tracker, "Monte Carlo study progress", total_replicates, &
            default_study_progress_reports, progress_enabled &
        )

        case_index = 0
        replicate_index = 0
        do basis_index = 1, size(study_config%basis_levels)
            do observation_index = 1, size(study_config%observation_counts)
                do step_index = 1, size(study_config%time_steps)
                    case_index = case_index + 1

                    study_grid = base_grid
                    study_grid%nx = study_config%basis_levels(basis_index)
                    study_grid%ny = study_config%basis_levels(basis_index)
                    study_parameters = base_parameters
                    study_parameters%n_observations = &
                        study_config%observation_counts(observation_index)
                    study_parameters%time_step = &
                        study_config%time_steps(step_index)

                    sigma_sum = 0.0_dp
                    sigma_sumsq = 0.0_dp
                    beta_sum = 0.0_dp
                    beta_sumsq = 0.0_dp
                    theta_sum = 0.0_dp
                    theta_sumsq = 0.0_dp
                    setup_sum = 0.0_dp
                    estimation_sum = 0.0_dp
                    total_sum = 0.0_dp

                    do replicate = 1, study_config%n_replicates
                        replicate_index = replicate_index + 1

                        call run_simulation_and_estimation(&
                            study_grid, study_parameters, &
                            derive_replicate_seed(base_seed, case_index, replicate), &
                            operators, state_history, estimates &
                        )

                        call fill_replicate_result(&
                            replicate_results(replicate_index), study_grid, &
                            study_parameters, operators%state_dimension, &
                            study_config%basis_levels(basis_index), replicate, &
                            estimates &
                        )

                        sigma_sum = sigma_sum + estimates%sigma_hat
                        sigma_sumsq = sigma_sumsq + estimates%sigma_hat**2
                        beta_sum = beta_sum + estimates%beta_hat
                        beta_sumsq = beta_sumsq + estimates%beta_hat**2
                        theta_sum = theta_sum + estimates%theta_hat
                        theta_sumsq = theta_sumsq + estimates%theta_hat**2
                        setup_sum = setup_sum + estimates%setup_time
                        estimation_sum = estimation_sum + estimates%estimation_time
                        total_sum = total_sum + estimates%setup_time + &
                            estimates%estimation_time

                        call update_progress_tracker(progress_tracker, replicate_index)
                    end do

                    sigma_mean = sigma_sum / real(study_config%n_replicates, dp)
                    beta_mean = beta_sum / real(study_config%n_replicates, dp)
                    theta_mean = theta_sum / real(study_config%n_replicates, dp)

                    call fill_case_summary(&
                        case_summaries(case_index), study_grid, study_parameters, &
                        study_config%basis_levels(basis_index), &
                        operators%state_dimension, study_config%n_replicates, &
                        sigma_mean, &
                        compute_sample_standard_deviation( &
                            sigma_sum, sigma_sumsq, study_config%n_replicates &
                        ), &
                        beta_mean, &
                        compute_sample_standard_deviation( &
                            beta_sum, beta_sumsq, study_config%n_replicates &
                        ), &
                        theta_mean, &
                        compute_sample_standard_deviation( &
                            theta_sum, theta_sumsq, study_config%n_replicates &
                        ), &
                        setup_sum / real(study_config%n_replicates, dp), &
                        estimation_sum / real(study_config%n_replicates, dp), &
                        total_sum / real(study_config%n_replicates, dp) &
                    )
                end do
            end do
        end do

        call finalize_progress_tracker(progress_tracker)
    end subroutine run_monte_carlo_study

    pure integer function derive_replicate_seed(&
        base_seed, case_index, replicate_index &
    ) result(seed_value)
        use iso_fortran_env, only: int64

        integer, intent(in) :: base_seed
        integer, intent(in) :: case_index
        integer, intent(in) :: replicate_index

        integer(int64) :: raw_seed

        raw_seed = int(base_seed, int64) + 104729_int64 * &
            int(case_index - 1, int64) + 13007_int64 * &
            int(replicate_index - 1, int64)
        seed_value = int(modulo(raw_seed, int(huge(1) - 1, int64))) + 1
    end function derive_replicate_seed

    pure real(dp) function compute_sample_standard_deviation(&
        sample_sum, sample_sumsq, n_samples &
    ) result(sample_sd)
        real(dp), intent(in) :: sample_sum
        real(dp), intent(in) :: sample_sumsq
        integer, intent(in) :: n_samples

        real(dp) :: mean_value
        real(dp) :: variance

        if (n_samples < 2) then
            sample_sd = 0.0_dp
            return
        end if

        mean_value = sample_sum / real(n_samples, dp)
        variance = (sample_sumsq - real(n_samples, dp) * mean_value**2) / &
            real(n_samples - 1, dp)
        sample_sd = sqrt(max(variance, 0.0_dp))
    end function compute_sample_standard_deviation

    pure subroutine fill_replicate_result(&
        replicate_result, grid, parameters, state_dimension, basis_level, &
        replicate, estimates &
    )
        type(monte_carlo_replicate_result_t), intent(out) :: replicate_result
        type(spatial_grid_t), intent(in) :: grid
        type(sde_parameters_t), intent(in) :: parameters
        integer, intent(in) :: state_dimension
        integer, intent(in) :: basis_level
        integer, intent(in) :: replicate
        type(parameter_estimates_t), intent(in) :: estimates

        replicate_result%basis_level = basis_level
        replicate_result%nx = grid%nx
        replicate_result%ny = grid%ny
        replicate_result%state_dimension = state_dimension
        replicate_result%n_observations = parameters%n_observations
        replicate_result%replicate = replicate
        replicate_result%time_step = parameters%time_step
        replicate_result%sigma_hat = estimates%sigma_hat
        replicate_result%sigma_true = parameters%sigma
        replicate_result%beta_hat = estimates%beta_hat
        replicate_result%beta_true = parameters%beta
        replicate_result%theta_hat = estimates%theta_hat
        replicate_result%theta_true = parameters%theta
        replicate_result%setup_time = estimates%setup_time
        replicate_result%estimation_time = estimates%estimation_time
        replicate_result%total_time = estimates%setup_time + &
            estimates%estimation_time
    end subroutine fill_replicate_result

    pure subroutine fill_case_summary(&
        case_summary, grid, parameters, basis_level, state_dimension, &
        n_replicates, sigma_mean, sigma_sd, beta_mean, beta_sd, theta_mean, &
        theta_sd, average_setup_time, average_estimation_time, &
        average_total_time &
    )
        type(monte_carlo_case_summary_t), intent(out) :: case_summary
        type(spatial_grid_t), intent(in) :: grid
        type(sde_parameters_t), intent(in) :: parameters
        integer, intent(in) :: basis_level
        integer, intent(in) :: state_dimension
        integer, intent(in) :: n_replicates
        real(dp), intent(in) :: sigma_mean
        real(dp), intent(in) :: sigma_sd
        real(dp), intent(in) :: beta_mean
        real(dp), intent(in) :: beta_sd
        real(dp), intent(in) :: theta_mean
        real(dp), intent(in) :: theta_sd
        real(dp), intent(in) :: average_setup_time
        real(dp), intent(in) :: average_estimation_time
        real(dp), intent(in) :: average_total_time

        case_summary%basis_level = basis_level
        case_summary%nx = grid%nx
        case_summary%ny = grid%ny
        case_summary%state_dimension = state_dimension
        case_summary%n_observations = parameters%n_observations
        case_summary%n_replicates = n_replicates
        case_summary%time_step = parameters%time_step
        case_summary%sigma_mean = sigma_mean
        case_summary%sigma_sd = sigma_sd
        case_summary%sigma_true = parameters%sigma
        case_summary%beta_mean = beta_mean
        case_summary%beta_sd = beta_sd
        case_summary%beta_true = parameters%beta
        case_summary%theta_mean = theta_mean
        case_summary%theta_sd = theta_sd
        case_summary%theta_true = parameters%theta
        case_summary%average_setup_time = average_setup_time
        case_summary%average_estimation_time = average_estimation_time
        case_summary%average_total_time = average_total_time
    end subroutine fill_case_summary

end module monte_carlo_study_mod
