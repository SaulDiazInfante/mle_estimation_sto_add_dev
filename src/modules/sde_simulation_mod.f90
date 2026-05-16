!> @file sde_simulation_mod.f90
!! @brief Time integration and random sampling for the spectral SDE model.
!> @brief Simulates state trajectories and random-number initialization.
module sde_simulation_mod
    use model_types_mod, only: dp
    use model_types_mod, only: sde_parameters_t
    use model_types_mod, only: spectral_operator_set_t
    use progress_reporting_mod, only: finalize_progress_tracker
    use progress_reporting_mod, only: initialize_progress_tracker
    use progress_reporting_mod, only: progress_tracker_t
    use progress_reporting_mod, only: update_progress_tracker
    implicit none
    private

    integer, parameter :: default_progress_reports = 20
    integer, parameter :: minimum_steps_for_progress = 10000
    real(dp), parameter :: pi_dp = acos(-1.0_dp)

    public :: set_random_seed
    public :: simulate_state_history
    public :: simulate_state_snapshots

contains

    !> Initializes the intrinsic Fortran RNG from a scalar seed.
    subroutine set_random_seed(seed_value)
        integer, intent(in) :: seed_value

        integer, allocatable :: seed_values(:)
        integer :: seed_index
        integer :: seed_size

        call random_seed(size=seed_size)
        allocate (seed_values(seed_size))

        do seed_index = 1, seed_size
            seed_values(seed_index) = modulo(&
                seed_value + 104729 * (seed_index - 1), huge(1) - 1 &
            ) + 1
        end do

        call random_seed(put=seed_values)
    end subroutine set_random_seed

    !> Simulates the full modal state history with Euler-Maruyama time stepping.
    subroutine simulate_state_history(&
        operators, parameters, state_history &
    )
        type(spectral_operator_set_t), intent(in) :: operators
        type(sde_parameters_t), intent(in) :: parameters
        real(dp), allocatable, intent(out) :: state_history(:, :)

        real(dp), allocatable :: diffusion_step(:)
        real(dp), allocatable :: drift(:)
        real(dp), allocatable :: gaussian_noise(:)
        real(dp), allocatable :: state(:)
        real(dp) :: sqrt_dt
        type(progress_tracker_t) :: simulation_progress
        integer :: step_index

        if (parameters%n_observations < 1) then
            write (*, '(a)') "n_observations must be at least 1"
            error stop
        end if

        allocate (state_history( &
            parameters%n_observations, operators%state_dimension &
        ))
        allocate (state(operators%state_dimension))
        allocate (drift(operators%state_dimension))
        allocate (gaussian_noise(operators%state_dimension))
        allocate (diffusion_step(operators%state_dimension))

        state = operators%initial_state
        state_history(1, :) = state

        if (parameters%n_observations == 1) return

        sqrt_dt = sqrt(parameters%time_step)
        diffusion_step = parameters%sigma * &
            operators%diffusion_diagonal * sqrt_dt
        call initialize_progress_tracker(&
            simulation_progress, "Simulation progress", &
            parameters%n_observations - 1, default_progress_reports, &
            parameters%n_observations - 1 >= minimum_steps_for_progress &
        )

        do step_index = 2, parameters%n_observations
            call evaluate_linear_drift(&
                parameters%beta, parameters%theta, &
                operators%eigenvalues, operators%interaction_matrix, &
                state, drift &
            )
            call generate_standard_normal_sample(&
                operators%state_dimension, gaussian_noise &
            )
            state = state + drift * parameters%time_step + &
                diffusion_step * gaussian_noise
            state_history(step_index, :) = state
            call update_progress_tracker(simulation_progress, step_index - 1)
        end do

        call finalize_progress_tracker(simulation_progress)
    end subroutine simulate_state_history

    !> Simulates only selected observation checkpoints to avoid storing full paths.
    subroutine simulate_state_snapshots(&
        operators, parameters, checkpoints, snapshot_states, report_progress, &
        progress_label &
    )
        type(spectral_operator_set_t), intent(in) :: operators
        type(sde_parameters_t), intent(in) :: parameters
        integer, intent(in) :: checkpoints(:)
        real(dp), allocatable, intent(out) :: snapshot_states(:, :)
        logical, intent(in), optional :: report_progress
        character(len=*), intent(in), optional :: progress_label

        real(dp), allocatable :: diffusion_step(:)
        real(dp), allocatable :: drift(:)
        real(dp), allocatable :: gaussian_noise(:)
        real(dp), allocatable :: state(:)
        logical :: progress_enabled
        character(len=:), allocatable :: progress_name
        type(progress_tracker_t) :: simulation_progress
        integer :: checkpoint_index
        integer :: final_checkpoint
        real(dp) :: sqrt_dt
        integer :: step_index

        if (parameters%n_observations < 1) then
            write (*, '(a)') "n_observations must be at least 1"
            error stop
        end if

        if (size(checkpoints) < 1) then
            write (*, '(a)') "At least one snapshot checkpoint is required."
            error stop
        end if

        if (any(checkpoints < 1) .or. any(checkpoints > parameters%n_observations)) then
            write (*, '(a)') "Snapshot checkpoints must lie within [1, n_observations]."
            error stop
        end if

        do checkpoint_index = 2, size(checkpoints)
            if (checkpoints(checkpoint_index) <= checkpoints(checkpoint_index - 1)) then
                write (*, '(a)') "Snapshot checkpoints must be strictly increasing."
                error stop
            end if
        end do

        allocate (snapshot_states(size(checkpoints), operators%state_dimension))
        allocate (state(operators%state_dimension))
        allocate (drift(operators%state_dimension))
        allocate (gaussian_noise(operators%state_dimension))
        allocate (diffusion_step(operators%state_dimension))

        state = operators%initial_state
        checkpoint_index = 1
        if (checkpoints(checkpoint_index) == 1) then
            snapshot_states(checkpoint_index, :) = state
            checkpoint_index = checkpoint_index + 1
        end if

        if (checkpoint_index > size(checkpoints)) return

        progress_enabled = .false.
        if (present(report_progress)) progress_enabled = report_progress
        if (present(progress_label)) then
            progress_name = trim(progress_label)
        else
            progress_name = "Snapshot simulation progress"
        end if

        final_checkpoint = checkpoints(size(checkpoints))
        sqrt_dt = sqrt(parameters%time_step)
        diffusion_step = parameters%sigma * operators%diffusion_diagonal * sqrt_dt
        call initialize_progress_tracker(&
            simulation_progress, progress_name, final_checkpoint - 1, &
            default_progress_reports, progress_enabled .and. &
            final_checkpoint - 1 >= minimum_steps_for_progress &
        )

        do step_index = 2, final_checkpoint
            call evaluate_linear_drift(&
                parameters%beta, parameters%theta, operators%eigenvalues, &
                operators%interaction_matrix, state, drift &
            )
            if (parameters%sigma /= 0.0_dp) then
                call generate_standard_normal_sample(&
                    operators%state_dimension, gaussian_noise &
                )
            else
                gaussian_noise = 0.0_dp
            end if

            state = state + drift * parameters%time_step + &
                diffusion_step * gaussian_noise

            if (step_index == checkpoints(checkpoint_index)) then
                snapshot_states(checkpoint_index, :) = state
                checkpoint_index = checkpoint_index + 1
                if (checkpoint_index > size(checkpoints)) exit
            end if

            call update_progress_tracker(simulation_progress, step_index - 1)
        end do

        call finalize_progress_tracker(simulation_progress)
    end subroutine simulate_state_snapshots

    pure subroutine evaluate_linear_drift(&
        beta, theta, eigenvalues, interaction_matrix, state, drift &
    )
        real(dp), intent(in) :: beta
        real(dp), intent(in) :: theta
        real(dp), intent(in) :: eigenvalues(:)
        real(dp), intent(in) :: interaction_matrix(:, :)
        real(dp), intent(in) :: state(:)
        real(dp), intent(out) :: drift(:)

        drift = -theta * matmul(interaction_matrix, state) - &
            beta * eigenvalues * state
    end subroutine evaluate_linear_drift

    subroutine generate_standard_normal_sample(sample_size, sample)
        integer, intent(in) :: sample_size
        real(dp), intent(out) :: sample(sample_size)

        real(dp) :: angle((sample_size + 1) / 2)
        real(dp) :: radius((sample_size + 1) / 2)
        real(dp) :: u1((sample_size + 1) / 2)
        real(dp) :: u2((sample_size + 1) / 2)

        call random_number(u1)
        call random_number(u2)

        u1 = max(u1, epsilon(1.0_dp))
        radius = sqrt(-2.0_dp * log(u1))
        angle = 2.0_dp * pi_dp * u2

        sample(1:sample_size:2) = radius * cos(angle)
        if (sample_size > 1) then
            sample(2:sample_size:2) = radius(1:sample_size / 2) * &
                sin(angle(1:sample_size / 2))
        end if
    end subroutine generate_standard_normal_sample

end module sde_simulation_mod
