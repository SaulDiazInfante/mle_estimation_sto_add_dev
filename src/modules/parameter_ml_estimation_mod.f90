!> @file parameter_ml_estimation_mod.f90
!! @brief Maximum-likelihood estimation routines for the spectral SDE model.
!> @brief Estimates drift and diffusion parameters from simulated state histories.
module parameter_ml_estimation_mod
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan
    use, intrinsic :: ieee_arithmetic, only: ieee_value
    use model_types_mod, only: dp
    use model_types_mod, only: parameter_estimates_t
    implicit none
    private

    public :: build_uniform_checkpoints
    public :: estimate_model_parameters
    public :: estimate_parameter_history
    public :: print_estimation_report

contains

    !> Estimates all model parameters from a complete state history.
    subroutine estimate_model_parameters(&
        state_history, time_step, interaction_matrix, eigenvalues, gamma, &
        sigma_hat, beta_hat, theta_hat &
    )
        real(dp), intent(in) :: state_history(:, :)
        real(dp), intent(in) :: time_step
        real(dp), intent(in) :: interaction_matrix(:, :)
        real(dp), intent(in) :: eigenvalues(:)
        real(dp), intent(in) :: gamma
        real(dp), intent(out) :: sigma_hat
        real(dp), intent(out) :: beta_hat
        real(dp), intent(out) :: theta_hat

        call estimate_sigma_from_quadratic_variation(&
            &state_history, &
            &time_step, &
            &eigenvalues, &
            &gamma, &
            &sigma_hat &
        &)
        call estimate_joint_drift_parameters(&
            &state_history, &
            &time_step, &
            &interaction_matrix, &
            &eigenvalues, gamma, &
            beta_hat, theta_hat &
        &)
    end subroutine estimate_model_parameters

    subroutine estimate_sigma_from_quadratic_variation(&
        &state_history, &
        &time_step, &
        &eigenvalues, &
        &gamma, &
        &sigma_hat &
    &)
        real(dp), intent(in) :: state_history(:, :)
        real(dp), intent(in) :: time_step
        real(dp), intent(in) :: eigenvalues(:)
        real(dp), intent(in) :: gamma
        real(dp), intent(out) :: sigma_hat

        integer :: mode_index
        integer :: n_observations
        integer :: n_state
        real(dp) :: increments_squared
        real(dp) :: sigma_sum
        real(dp) :: total_time

        n_observations = size(state_history, 1)
        n_state = size(state_history, 2)

        if (n_observations < 2) then
            write (*, '(a)') "At least two observations are required"
            error stop
        end if

        total_time = time_step * real(n_observations - 1, dp)
        sigma_sum = 0.0_dp

        do mode_index = 1, n_state
            increments_squared = sum(&
                (state_history(2:n_observations, mode_index) - &
                state_history(1:n_observations - 1, mode_index))**2 &
            )

            sigma_sum = sigma_sum + sqrt(&
                increments_squared * &
                eigenvalues(mode_index)**(2.0_dp * gamma) / total_time &
            )
        end do

        sigma_hat = sigma_sum / real(n_state, dp)
    end subroutine estimate_sigma_from_quadratic_variation

    subroutine estimate_joint_drift_parameters(&
        state_history, time_step, interaction_matrix, eigenvalues, gamma, &
        beta_hat, theta_hat &
    )
        real(dp), intent(in) :: state_history(:, :)
        real(dp), intent(in) :: time_step
        real(dp), intent(in) :: interaction_matrix(:, :)
        real(dp), intent(in) :: eigenvalues(:)
        real(dp), intent(in) :: gamma
        real(dp), intent(out) :: beta_hat
        real(dp), intent(out) :: theta_hat

        logical :: success

        call solve_joint_mle_system(&
            state_history, time_step, interaction_matrix, eigenvalues, gamma, &
            beta_hat, theta_hat, success &
        )

        if (.not. success) then
            write (*, '(a)') "The MLE normal equations are singular"
            error stop
        end if
    end subroutine estimate_joint_drift_parameters

    !> Builds evenly distributed observation checkpoints for trajectory summaries.
    subroutine build_uniform_checkpoints(&
        total_points, requested_points, minimum_points, checkpoints &
    )
        integer, intent(in) :: total_points
        integer, intent(in) :: requested_points
        integer, intent(in) :: minimum_points
        integer, allocatable, intent(out) :: checkpoints(:)

        integer :: available_points
        integer :: checkpoint_index
        integer :: n_checkpoints

        if (total_points < 2) then
            write (*, '(a)') "At least two observations are required"
            error stop
        end if

        if (requested_points < 1) then
            write (*, '(a)') "requested_points must be positive"
            error stop
        end if

        if (minimum_points < 2 .or. minimum_points > total_points) then
            write (*, '(a)') "minimum_points must belong to [2, total_points]"
            error stop
        end if

        available_points = total_points - minimum_points + 1
        n_checkpoints = min(requested_points, available_points)
        allocate (checkpoints(n_checkpoints))

        if (n_checkpoints == 1) then
            checkpoints(1) = total_points
            return
        end if

        do checkpoint_index = 1, n_checkpoints
            checkpoints(checkpoint_index) = minimum_points + int(&
                real(checkpoint_index - 1, dp) * &
                real(available_points - 1, dp) / &
                real(n_checkpoints - 1, dp) &
            )
        end do
    end subroutine build_uniform_checkpoints

    !> Computes the parameter estimates at several observation checkpoints.
    subroutine estimate_parameter_history(&
        state_history, time_step, interaction_matrix, eigenvalues, gamma, &
        checkpoints, times, sigma_history, beta_history, theta_history &
    )
        real(dp), intent(in) :: state_history(:, :)
        real(dp), intent(in) :: time_step
        real(dp), intent(in) :: interaction_matrix(:, :)
        real(dp), intent(in) :: eigenvalues(:)
        real(dp), intent(in) :: gamma
        integer, intent(in) :: checkpoints(:)
        real(dp), allocatable, intent(out) :: times(:)
        real(dp), allocatable, intent(out) :: sigma_history(:)
        real(dp), allocatable, intent(out) :: beta_history(:)
        real(dp), allocatable, intent(out) :: theta_history(:)

        integer :: checkpoint_index
        real(dp) :: not_a_number
        integer :: n_observations_at_checkpoint
        logical :: success

        allocate (times(size(checkpoints)))
        allocate (sigma_history(size(checkpoints)))
        allocate (beta_history(size(checkpoints)))
        allocate (theta_history(size(checkpoints)))

        not_a_number = ieee_value(0.0_dp, ieee_quiet_nan)

        do checkpoint_index = 1, size(checkpoints)
            n_observations_at_checkpoint = checkpoints(checkpoint_index)

            if (n_observations_at_checkpoint < 2 .or. &
                n_observations_at_checkpoint > size(state_history, 1)) then
                write (*, '(a)') "Checkpoint outside valid range"
                error stop
            end if

            times(checkpoint_index) = time_step * &
                real(n_observations_at_checkpoint - 1, dp)

            call estimate_sigma_from_quadratic_variation(&
                state_history(1:n_observations_at_checkpoint, :), time_step, &
                eigenvalues, gamma, sigma_history(checkpoint_index) &
            )

            beta_history(checkpoint_index) = not_a_number
            theta_history(checkpoint_index) = not_a_number

            call solve_joint_mle_system(&
                state_history(1:n_observations_at_checkpoint, :), time_step, &
                interaction_matrix, eigenvalues, gamma, &
                beta_history(checkpoint_index), &
                theta_history(checkpoint_index), success &
            )

            if (.not. success) then
                beta_history(checkpoint_index) = not_a_number
                theta_history(checkpoint_index) = not_a_number
            end if
        end do
    end subroutine estimate_parameter_history

    subroutine compute_mle_statistics(&
        state_history, time_step, interaction_matrix, eigenvalues, gamma, &
        statistics &
    )
        real(dp), intent(in) :: state_history(:, :)
        real(dp), intent(in) :: time_step
        real(dp), intent(in) :: interaction_matrix(:, :)
        real(dp), intent(in) :: eigenvalues(:)
        real(dp), intent(in) :: gamma
        real(dp), intent(out) :: statistics(5)

        real(dp), allocatable :: coupled_state_history(:, :)
        real(dp), allocatable :: statistic_1(:)
        real(dp), allocatable :: statistic_2(:)
        real(dp), allocatable :: statistic_3(:)
        real(dp), allocatable :: statistic_4(:)
        real(dp), allocatable :: statistic_5(:)
        real(dp), allocatable :: weight_1(:)
        real(dp), allocatable :: weight_2(:)
        real(dp), allocatable :: weight_3(:)
        integer :: mode_index
        integer :: n_observations
        integer :: n_state

        n_observations = size(state_history, 1)
        n_state = size(state_history, 2)

        allocate (coupled_state_history(n_observations, n_state))
        allocate (statistic_1(n_state), statistic_2(n_state))
        allocate (statistic_3(n_state), statistic_4(n_state))
        allocate (statistic_5(n_state))
        allocate (weight_1(n_state), weight_2(n_state), weight_3(n_state))

        coupled_state_history = matmul(state_history, transpose(interaction_matrix))

        !$omp parallel do default(none) schedule(static) &
        !$omp& shared(n_state, state_history, coupled_state_history, &
        !$omp& time_step, statistic_1, statistic_2, statistic_3, &
        !$omp& statistic_4, statistic_5) private(mode_index)
        do mode_index = 1, n_state
            call compute_ito_integral(&
                state_history(:, mode_index), &
                state_history(:, mode_index), &
                statistic_1(mode_index) &
            )
            call compute_ito_integral(&
                coupled_state_history(:, mode_index), &
                state_history(:, mode_index), &
                statistic_2(mode_index) &
            )
            call compute_trapezoidal_integral(&
                state_history(:, mode_index) * state_history(:, mode_index), &
                time_step, statistic_3(mode_index) &
            )
            call compute_trapezoidal_integral(&
                state_history(:, mode_index) * &
                coupled_state_history(:, mode_index), &
                time_step, statistic_4(mode_index) &
            )
            call compute_trapezoidal_integral(&
                coupled_state_history(:, mode_index) * &
                coupled_state_history(:, mode_index), &
                time_step, statistic_5(mode_index) &
            )
        end do
        !$omp end parallel do

        weight_1 = eigenvalues**(1.0_dp + 2.0_dp * gamma)
        weight_2 = eigenvalues**(2.0_dp * gamma)
        weight_3 = eigenvalues**(2.0_dp + 2.0_dp * gamma)

        statistics(1) = sum(statistic_1 * weight_1)
        statistics(2) = sum(statistic_2 * weight_2)
        statistics(3) = sum(statistic_3 * weight_3)
        statistics(4) = sum(statistic_4 * weight_1)
        statistics(5) = sum(statistic_5 * weight_2)
    end subroutine compute_mle_statistics

    subroutine solve_joint_mle_system(&
        state_history, time_step, interaction_matrix, eigenvalues, gamma, &
        beta_hat, theta_hat, success &
    )
        real(dp), intent(in) :: state_history(:, :)
        real(dp), intent(in) :: time_step
        real(dp), intent(in) :: interaction_matrix(:, :)
        real(dp), intent(in) :: eigenvalues(:)
        real(dp), intent(in) :: gamma
        real(dp), intent(out) :: beta_hat
        real(dp), intent(out) :: theta_hat
        logical, intent(out) :: success

        real(dp) :: denominator
        real(dp) :: scale
        real(dp) :: statistics(5)
        real(dp) :: tolerance

        call compute_mle_statistics(&
            state_history, time_step, interaction_matrix, eigenvalues, gamma, &
            statistics &
        )

        denominator = statistics(4)**2 - statistics(3) * statistics(5)
        scale = max(1.0_dp, abs(statistics(4)**2) + &
            abs(statistics(3) * statistics(5)))
        tolerance = sqrt(epsilon(1.0_dp)) * scale

        if (abs(denominator) <= tolerance) then
            success = .false.
            beta_hat = 0.0_dp
            theta_hat = 0.0_dp
            return
        end if

        success = .true.
        beta_hat = (statistics(1) * statistics(5) - &
            statistics(2) * statistics(4)) / denominator
        theta_hat = (statistics(2) * statistics(3) - &
            statistics(1) * statistics(4)) / denominator
    end subroutine solve_joint_mle_system

    pure subroutine compute_ito_integral(integrand, process, integral_value)
        real(dp), intent(in) :: integrand(:)
        real(dp), intent(in) :: process(:)
        real(dp), intent(out) :: integral_value

        integer :: n_observations

        n_observations = size(integrand)
        if (n_observations < 2) then
            integral_value = 0.0_dp
            return
        end if

        integral_value = dot_product(&
            integrand(1:n_observations - 1), &
            process(2:n_observations) - process(1:n_observations - 1) &
        )
    end subroutine compute_ito_integral

    pure subroutine compute_trapezoidal_integral(&
        path_values, time_step, integral_value &
    )
        real(dp), intent(in) :: path_values(:)
        real(dp), intent(in) :: time_step
        real(dp), intent(out) :: integral_value

        integer :: n_observations

        n_observations = size(path_values)
        if (n_observations < 2) then
            integral_value = 0.0_dp
            return
        end if

        integral_value = 0.5_dp * time_step * (&
            sum(path_values(1:n_observations - 1)) + &
            sum(path_values(2:n_observations)) &
        )
    end subroutine compute_trapezoidal_integral

    !> Prints a formatted report comparing true and estimated parameter values.
    subroutine print_estimation_report(&
        sigma_true, beta_true, theta_true, estimates &
    )
        real(dp), intent(in) :: sigma_true
        real(dp), intent(in) :: beta_true
        real(dp), intent(in) :: theta_true
        type(parameter_estimates_t), intent(in) :: estimates

        print '(a)', ''
        print '(a)', repeat('=', 86)
        print '(a)', 'MLE ESTIMATION REPORT'
        print '(a)', repeat('-', 86)
        write (*, '(a10,3x,a12,3x,a12,3x,a12,3x,a13)') &
            'Parameter', 'True Value', 'Estimate', 'Abs Error', &
            'Rel Error (%)'
        print '(a)', repeat('-', 86)

        call print_parameter_summary_line(&
            'sigma', sigma_true, estimates%sigma_hat &
        )
        call print_parameter_summary_line(&
            'beta', beta_true, estimates%beta_hat &
        )
        call print_parameter_summary_line(&
            'theta', theta_true, estimates%theta_hat &
        )

        print '(a)', repeat('-', 86)
        print '(a)', 'Runtime (seconds)'
        write (*, '(a20,2x,f12.6)') 'Set-up', estimates%setup_time
        write (*, '(a20,2x,f12.6)') 'MLE estimation', &
            estimates%estimation_time
        print '(a)', repeat('=', 86)
    end subroutine print_estimation_report

    subroutine print_parameter_summary_line(&
        parameter_name, true_value, estimate_value &
    )
        character(len=*), intent(in) :: parameter_name
        real(dp), intent(in) :: true_value
        real(dp), intent(in) :: estimate_value

        real(dp) :: absolute_error
        real(dp) :: relative_error

        absolute_error = abs(estimate_value - true_value)
        if (abs(true_value) > tiny(true_value)) then
            relative_error = 100.0_dp * absolute_error / abs(true_value)
        else
            relative_error = 0.0_dp
        end if

        write (*, '(a10,3x,f12.6,3x,f12.6,3x,es12.4,3x,f13.4)') &
            trim(parameter_name), true_value, estimate_value, &
            absolute_error, relative_error
    end subroutine print_parameter_summary_line

end module parameter_ml_estimation_mod
