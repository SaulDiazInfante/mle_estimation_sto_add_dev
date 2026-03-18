!===========================================================================
! File: aditive_v2_fortran_modular.f90
!
! Legacy consolidated reference file.
! The maintained implementation now lives in the split Fortran sources
! used by the Makefile.
!
! Build:
!   gfortran -std=f2008 -O3 aditive_v2_fortran_modular.f90 -o modular.out
!===========================================================================
module model_types
    use iso_fortran_env, only: real64
    implicit none
    private

    integer, parameter, public :: rk = real64

    type, public :: grid_config_t
        integer :: nx = 10
        integer :: ny = 10
        integer :: r_v = 1
        integer :: s_v = 1
        real(rk) :: lx = 5.0_rk
        real(rk) :: ly = 5.0_rk
        real(rk) :: gamma = 2.0_rk
    end type grid_config_t

    type, public :: sde_config_t
        integer :: nobs = 1000
        real(rk) :: delta = 1.0e-5_rk
        real(rk) :: beta = 0.1_rk
        real(rk) :: theta = 1.0_rk
        real(rk) :: sigma = 0.1_rk
    end type sde_config_t

    type, public :: operator_set_t
        integer :: dim = 0
        integer, allocatable :: mode_order(:, :)
        real(rk), allocatable :: initial_state(:)
        real(rk), allocatable :: lambdas(:)
        real(rk), allocatable :: diffusion_diag(:)
        real(rk), allocatable :: interaction(:, :)
    end type operator_set_t

    type, public :: estimation_result_t
        real(rk) :: sigma_hat = 0.0_rk
        real(rk) :: beta_hat = 0.0_rk
        real(rk) :: theta_hat = 0.0_rk
        real(rk) :: setup_time = 0.0_rk
        real(rk) :: mle_time = 0.0_rk
    end type estimation_result_t

    public :: active_dim

contains

    pure integer function active_dim(grid) result(dim)
        type(grid_config_t), intent(in) :: grid

        dim = grid%nx * grid%ny - 1
    end function active_dim

end module model_types

module diagnostics_mod
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use model_types, only: rk
    implicit none
    private

    interface assert_all_finite
        module procedure assert_all_finite_vector
        module procedure assert_all_finite_matrix
    end interface assert_all_finite

    public :: assert_all_finite, write_paths

contains

    subroutine assert_all_finite_vector(label, values)
        character(len=*), intent(in) :: label
        real(rk), intent(in) :: values(:)

        if (any(.not. ieee_is_finite(values))) then
            write (*, '(a)') trim(label)//" contains non-finite values"
            error stop
        end if
    end subroutine assert_all_finite_vector

    subroutine assert_all_finite_matrix(label, values)
        character(len=*), intent(in) :: label
        real(rk), intent(in) :: values(:, :)

        if (any(.not. ieee_is_finite(values))) then
            write (*, '(a)') trim(label)//" contains non-finite values"
            error stop
        end if
    end subroutine assert_all_finite_matrix

    subroutine write_paths(filename, paths)
        character(len=*), intent(in) :: filename
        real(rk), intent(in) :: paths(:, :)

        integer :: unit
        integer :: row_idx
        integer :: col_idx

        open (newunit=unit, file=filename, status="replace", action="write")
        write (unit, '(a)', advance='no') "obs_index"
        do col_idx = 1, size(paths, 2)
            write (unit, '(",",a,i0)', advance='no') "x_", col_idx
        end do
        write (unit, *)

        do row_idx = 1, size(paths, 1)
            write (unit, '(i0)', advance='no') row_idx
            do col_idx = 1, size(paths, 2)
                write (unit, '(",",es26.17e3)', advance='no') paths(row_idx, col_idx)
            end do
            write (unit, *)
        end do
        close (unit)
    end subroutine write_paths

end module diagnostics_mod

module spectral_model_mod
    use model_types, only: rk, grid_config_t, operator_set_t, active_dim
    implicit none
    private

    real(rk), parameter :: pi_rk = acos(-1.0_rk)

    public :: build_problem_data

contains

    subroutine build_problem_data(grid, operators)
        type(grid_config_t), intent(in) :: grid
        type(operator_set_t), intent(out) :: operators

        real(rk), allocatable :: u0_prime(:, :)

        operators%dim = active_dim(grid)

        call build_mode_table(grid, operators%mode_order, operators%lambdas)
        allocate (operators%initial_state(operators%dim))
        allocate (operators%interaction(operators%dim, operators%dim))
        allocate (operators%diffusion_diag(operators%dim))

        call build_fluctuation_field(grid, u0_prime)
        call project_field(grid, u0_prime, operators%mode_order, operators%initial_state)
        call build_operator_matrix(&
            & grid, &
            & operators%mode_order, &
            & operators%lambdas, &
            & operators%interaction, &
            & operators%diffusion_diag &
        &)
    end subroutine build_problem_data

    pure real(rk) function initial_bump(x, y, grid) result(u0)
        real(rk), intent(in) :: x, y
        type(grid_config_t), intent(in) :: grid

        real(rk), parameter :: support_fraction = 0.25_rk
        real(rk) :: ellipse_radius
        real(rk) :: x_center, y_center
        real(rk) :: x_normalized, y_normalized

        x_center = grid%lx / 4.0_rk
        y_center = 0.6_rk * grid%ly

        x_normalized = (x - x_center) / (support_fraction * grid%lx)
        y_normalized = (y - y_center) / (support_fraction * grid%ly)
        ellipse_radius = sqrt(x_normalized**2 + y_normalized**2)

        u0 = 0.0_rk
        if (ellipse_radius <= 1.0_rk) then
            u0 = 0.5_rk * (1.0_rk + cos(pi_rk * ellipse_radius))
        end if
    end function initial_bump

    subroutine build_fluctuation_field(grid, u0_prime)
        type(grid_config_t), intent(in) :: grid
        real(rk), allocatable, intent(out) :: u0_prime(:, :)

        real(rk), allocatable :: u0_grid(:, :)
        integer :: ix, iy
        real(rk) :: dx, dy, xi, yj
        real(rk) :: u0_bar

        allocate (u0_grid(grid%nx, grid%ny))
        allocate (u0_prime(grid%nx, grid%ny))

        dx = grid%lx / real(grid%nx, rk)
        dy = grid%ly / real(grid%ny, rk)

        do iy = 1, grid%ny
            yj = dy * (real(iy, rk) - 0.5_rk)
            do ix = 1, grid%nx
                xi = dx * (real(ix, rk) - 0.5_rk)
                u0_grid(ix, iy) = initial_bump(xi, yj, grid)
            end do
        end do

        u0_bar = sum(u0_grid) / real(grid%nx * grid%ny, rk)
        u0_prime = u0_grid - u0_bar
    end subroutine build_fluctuation_field

    subroutine build_mode_table(grid, mode_order, lambdas)
        type(grid_config_t), intent(in) :: grid
        integer, allocatable, intent(out) :: mode_order(:, :)
        real(rk), allocatable, intent(out) :: lambdas(:)

        integer :: dim
        integer :: mode_idx
        integer :: mode_i, mode_j

        dim = active_dim(grid)
        allocate (mode_order(dim, 2))
        allocate (lambdas(dim))

        mode_idx = 0
        do mode_i = 0, grid%nx - 1
            do mode_j = 0, grid%ny - 1
                if (mode_i + mode_j == 0) cycle
                mode_idx = mode_idx + 1
                mode_order(mode_idx, 1) = mode_i
                mode_order(mode_idx, 2) = mode_j
                lambdas(mode_idx) = eigenvalue(mode_i, mode_j, grid)
            end do
        end do

        call sort_modes_by_lambda(lambdas, mode_order)
    end subroutine build_mode_table

    pure real(rk) function eigenvalue(mode_i, mode_j, grid) result(lambda)
        integer, intent(in) :: mode_i, mode_j
        type(grid_config_t), intent(in) :: grid

        lambda = (pi_rk / grid%lx)**2 * real(mode_i, rk)**2 + &
                 (pi_rk / grid%ly)**2 * real(mode_j, rk)**2
    end function eigenvalue

    subroutine sort_modes_by_lambda(lambdas, mode_order)
        real(rk), intent(inout) :: lambdas(:)
        integer, intent(inout) :: mode_order(:, :)

        integer :: left_idx, right_idx, min_idx
        integer :: temp_mode(2)
        real(rk) :: temp_lambda

        do left_idx = 1, size(lambdas) - 1
            min_idx = left_idx
            do right_idx = left_idx + 1, size(lambdas)
                if (lambdas(right_idx) < lambdas(min_idx)) then
                    min_idx = right_idx
                end if
            end do

            if (min_idx /= left_idx) then
                temp_lambda = lambdas(left_idx)
                lambdas(left_idx) = lambdas(min_idx)
                lambdas(min_idx) = temp_lambda

                temp_mode = mode_order(left_idx, :)
                mode_order(left_idx, :) = mode_order(min_idx, :)
                mode_order(min_idx, :) = temp_mode
            end if
        end do
    end subroutine sort_modes_by_lambda

    subroutine project_field(grid, u0_prime, mode_order, coeffs)
        type(grid_config_t), intent(in) :: grid
        real(rk), intent(in) :: u0_prime(:, :)
        integer, intent(in) :: mode_order(:, :)
        real(rk), intent(out) :: coeffs(:)

        integer :: mode_idx

        do mode_idx = 1, size(coeffs)
            coeffs(mode_idx) = mode_integral(&
                & u0_prime, &
                & mode_order(mode_idx, 1), &
                & mode_order(mode_idx, 2), &
                & grid &
            &)
            coeffs(mode_idx) = coeffs(mode_idx) / (grid%lx * grid%ly)
        end do
    end subroutine project_field

    pure real(rk) function mode_integral(field, mode_i, mode_j, grid) result(intval)
        real(rk), intent(in) :: field(:, :)
        integer, intent(in) :: mode_i, mode_j
        type(grid_config_t), intent(in) :: grid

        integer :: ix, iy
        real(rk) :: dx, dy, xi, yj
        real(rk) :: basis_i, basis_j
        real(rk) :: wave_x, wave_y
        real(rk) :: cell_area

        dx = grid%lx / real(grid%nx, rk)
        dy = grid%ly / real(grid%ny, rk)
        cell_area = dx * dy
        wave_x = pi_rk * real(mode_i, rk) / grid%lx
        wave_y = pi_rk * real(mode_j, rk) / grid%ly

        intval = 0.0_rk
        do iy = 1, grid%ny
            yj = dy * (real(iy, rk) - 0.5_rk)
            basis_j = basis_normalization(mode_j) * cos(wave_y * yj)
            do ix = 1, grid%nx
                xi = dx * (real(ix, rk) - 0.5_rk)
                basis_i = basis_normalization(mode_i) * cos(wave_x * xi)
                intval = intval + cell_area * field(ix, iy) * basis_i * basis_j
            end do
        end do
    end function mode_integral

    pure real(rk) function basis_normalization(mode_number) result(norm_factor)
        integer, intent(in) :: mode_number

        if (mode_number == 0) then
            norm_factor = 1.0_rk
        else
            norm_factor = sqrt(2.0_rk)
        end if
    end function basis_normalization

    subroutine build_operator_matrix(grid, mode_order, lambdas, interaction, diffusion_diag)
        type(grid_config_t), intent(in) :: grid
        integer, intent(in) :: mode_order(:, :)
        real(rk), intent(in) :: lambdas(:)
        real(rk), intent(out) :: interaction(:, :)
        real(rk), intent(out) :: diffusion_diag(:)

        integer :: source_mode
        integer :: target_mode

        interaction = 0.0_rk

        do source_mode = 1, size(lambdas)
            diffusion_diag(source_mode) = lambdas(source_mode)**(-grid%gamma)
            do target_mode = 1, size(lambdas)
                interaction(target_mode, source_mode) = interaction_coefficient(&
                    & mode_order(source_mode, 1), &
                    & mode_order(source_mode, 2), &
                    & mode_order(target_mode, 1), &
                    & mode_order(target_mode, 2), &
                    & grid &
                &)
            end do
        end do
    end subroutine build_operator_matrix

    pure real(rk) function interaction_coefficient(mode_i, mode_j, mode_k, mode_l, grid) result(coeff)
        integer, intent(in) :: mode_i, mode_j, mode_k, mode_l
        type(grid_config_t), intent(in) :: grid

        real(rk) :: base_coefficient
        real(rk) :: interaction_term
        real(rk) :: sqrt2_factor
        real(rk) :: coefficient_i0, coefficient_j0
        real(rk) :: base_coeff_k0, base_coeff_l0
        logical :: i_minus_r_matches_k, i_plus_r_matches_k
        logical :: j_minus_s_matches_l, j_plus_s_matches_l

        coeff = 0.0_rk
        base_coefficient = pi_rk**2 / (2.0_rk * grid%lx * grid%ly)

        if (mode_i * mode_j * mode_k * mode_l > 0) then
            i_minus_r_matches_k = ((mode_i - grid%r_v)**2 == mode_k**2)
            i_plus_r_matches_k = ((mode_i + grid%r_v)**2 == mode_k**2)
            j_minus_s_matches_l = ((mode_j - grid%s_v)**2 == mode_l**2)
            j_plus_s_matches_l = ((mode_j + grid%s_v)**2 == mode_l**2)

            interaction_term = real(mode_j * grid%r_v - mode_i * grid%s_v, rk)

            if (i_minus_r_matches_k .and. j_minus_s_matches_l) then
                coeff = -base_coefficient * interaction_term
            end if

            if (i_plus_r_matches_k .and. j_minus_s_matches_l) then
                coeff = -base_coefficient * real(mode_j * grid%r_v + mode_i * grid%s_v, rk)
            end if

            if (i_plus_r_matches_k .and. j_plus_s_matches_l) then
                coeff = base_coefficient * interaction_term
            end if

            if (i_minus_r_matches_k .and. j_plus_s_matches_l) then
                coeff = base_coefficient * real(mode_j * grid%r_v + mode_i * grid%s_v, rk)
            end if
        end if

        if (mode_i == 0 .and. mode_j * mode_k * mode_l > 0) then
            sqrt2_factor = sqrt(2.0_rk) / 2.0_rk
            coefficient_i0 = sqrt2_factor * pi_rk**2 * real(mode_j * grid%r_v, rk) / (grid%lx * grid%ly)

            if (mode_k == grid%r_v .and. (mode_j - grid%s_v)**2 == mode_l**2) then
                coeff = -coefficient_i0
            end if

            if (mode_k == grid%r_v .and. (mode_j + grid%s_v)**2 == mode_l**2) then
                coeff = coefficient_i0
            end if
        end if

        if (mode_k == 0 .and. mode_i * mode_j * mode_l > 0) then
            sqrt2_factor = sqrt(2.0_rk) / 2.0_rk
            base_coeff_k0 = sqrt2_factor * pi_rk**2 / (grid%lx * grid%ly)

            if (mode_i == grid%r_v .and. (mode_j - grid%s_v)**2 == mode_l**2) then
                coeff = -base_coeff_k0 * real(mode_j * grid%r_v - mode_i * grid%s_v, rk)
            end if

            if (mode_i == grid%r_v .and. (mode_j + grid%s_v)**2 == mode_l**2) then
                coeff = base_coeff_k0 * real(mode_j * grid%r_v + mode_i * grid%s_v, rk)
            end if
        end if

        if (mode_j == 0 .and. mode_i * mode_k * mode_l > 0) then
            sqrt2_factor = sqrt(2.0_rk) / 2.0_rk
            coefficient_j0 = sqrt2_factor * pi_rk**2 * real(mode_i * grid%s_v, rk) / (grid%lx * grid%ly)

            if (mode_l == grid%s_v .and. (mode_i - grid%r_v)**2 == mode_k**2) then
                coeff = coefficient_j0
            end if

            if (mode_l == grid%s_v .and. (mode_i + grid%r_v)**2 == mode_k**2) then
                coeff = -coefficient_j0
            end if
        end if

        if (mode_l == 0 .and. mode_i * mode_j * mode_k > 0) then
            sqrt2_factor = sqrt(2.0_rk) / 2.0_rk
            base_coeff_l0 = sqrt2_factor * pi_rk**2 / (grid%lx * grid%ly)

            if (mode_j == grid%s_v .and. (mode_i - grid%r_v)**2 == mode_k**2) then
                coeff = -base_coeff_l0 * real(mode_j * grid%r_v - mode_i * grid%s_v, rk)
            end if

            if (mode_j == grid%s_v .and. (mode_i + grid%r_v)**2 == mode_k**2) then
                coeff = -base_coeff_l0 * real(mode_j * grid%r_v + mode_i * grid%s_v, rk)
            end if
        end if

        if (mode_i == 0 .and. mode_l == 0 .and. mode_j * mode_k > 0) then
            if (mode_k == grid%r_v .and. mode_j == grid%s_v) then
                coeff = -pi_rk**2 * real(mode_j * grid%r_v, rk) / (grid%lx * grid%ly)
            end if
        end if

        if (mode_j == 0 .and. mode_k == 0 .and. mode_i * mode_l > 0) then
            if (mode_i == grid%r_v .and. mode_l == grid%s_v) then
                coeff = pi_rk**2 * real(mode_i * grid%s_v, rk) / (grid%lx * grid%ly)
            end if
        end if
    end function interaction_coefficient

end module spectral_model_mod

module sde_solver_mod
    use model_types, only: rk, operator_set_t, sde_config_t
    implicit none
    private

    real(rk), parameter :: pi_rk = acos(-1.0_rk)

    public :: initialize_random_seed, simulate_paths

contains

    subroutine initialize_random_seed(base_seed)
        integer, intent(in) :: base_seed

        integer :: seed_size
        integer :: seed_idx
        integer, allocatable :: seed_values(:)

        call random_seed(size=seed_size)
        allocate (seed_values(seed_size))

        do seed_idx = 1, seed_size
            seed_values(seed_idx) = modulo(base_seed + 104729 * (seed_idx - 1), huge(1) - 1) + 1
        end do

        call random_seed(put=seed_values)
    end subroutine initialize_random_seed

    subroutine simulate_paths(operators, sde, paths)
        type(operator_set_t), intent(in) :: operators
        type(sde_config_t), intent(in) :: sde
        real(rk), allocatable, intent(out) :: paths(:, :)

        real(rk), allocatable :: state(:)
        real(rk), allocatable :: drift(:)
        real(rk), allocatable :: noise(:)
        real(rk), allocatable :: scaled_diffusion(:)
        real(rk) :: sqrt_delta
        integer :: step_idx

        if (sde%nobs < 1) then
            write (*, '(a)') "nobs must be at least 1"
            error stop
        end if

        allocate (paths(sde%nobs, operators%dim))
        allocate (state(operators%dim), drift(operators%dim), noise(operators%dim))
        allocate (scaled_diffusion(operators%dim))

        state = operators%initial_state
        paths(1, :) = state

        if (sde%nobs == 1) return

        sqrt_delta = sqrt(sde%delta)
        scaled_diffusion = sde%sigma * operators%diffusion_diag * sqrt_delta

        do step_idx = 2, sde%nobs
            call compute_linear_drift(&
                & sde%beta, &
                & sde%theta, &
                & operators%lambdas, &
                & operators%interaction, &
                & state, &
                & drift &
            &)
            call generate_normal_array(operators%dim, noise)
            state = state + drift * sde%delta + scaled_diffusion * noise
            paths(step_idx, :) = state
        end do
    end subroutine simulate_paths

    pure subroutine compute_linear_drift(beta, theta, lambdas, interaction, state, drift)
        real(rk), intent(in) :: beta, theta
        real(rk), intent(in) :: lambdas(:)
        real(rk), intent(in) :: interaction(:, :)
        real(rk), intent(in) :: state(:)
        real(rk), intent(out) :: drift(:)

        drift = -theta * matmul(interaction, state) - beta * lambdas * state
    end subroutine compute_linear_drift

    subroutine generate_normal_array(dim, noise)
        integer, intent(in) :: dim
        real(rk), intent(out) :: noise(dim)

        real(rk) :: u1((dim + 1) / 2)
        real(rk) :: u2((dim + 1) / 2)
        real(rk) :: radius((dim + 1) / 2)
        real(rk) :: angle((dim + 1) / 2)

        call random_number(u1)
        call random_number(u2)

        u1 = max(u1, epsilon(1.0_rk))
        radius = sqrt(-2.0_rk * log(u1))
        angle = 2.0_rk * pi_rk * u2

        noise(1:dim:2) = radius
        noise(1:dim:2) = noise(1:dim:2) * cos(angle)

        if (dim > 1) then
            noise(2:dim:2) = radius(1:dim/2)
            noise(2:dim:2) = noise(2:dim:2) * sin(angle(1:dim/2))
        end if
    end subroutine generate_normal_array

end module sde_solver_mod

module mle_estimation_mod
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use model_types, only: rk, estimation_result_t
    implicit none
    private

    public :: build_uniform_checkpoints, estimate_parameter_trajectory
    public :: estimate_parameters, print_estimation_report, write_estimator_trajectory

contains

    subroutine estimate_parameters(paths, delta, interaction, lambdas, gamma, sigma_hat, beta_hat, theta_hat)
        real(rk), intent(in) :: paths(:, :)
        real(rk), intent(in) :: delta, gamma
        real(rk), intent(in) :: interaction(:, :)
        real(rk), intent(in) :: lambdas(:)
        real(rk), intent(out) :: sigma_hat, beta_hat, theta_hat

        call estimate_sigma_from_qv(paths, delta, lambdas, gamma, sigma_hat)
        call estimate_joint_mle(paths, delta, interaction, lambdas, gamma, beta_hat, theta_hat)
    end subroutine estimate_parameters

    subroutine estimate_sigma_from_qv(paths, delta, lambdas, gamma, sigma_hat)
        real(rk), intent(in) :: paths(:, :)
        real(rk), intent(in) :: delta, gamma
        real(rk), intent(in) :: lambdas(:)
        real(rk), intent(out) :: sigma_hat

        integer :: dim
        integer :: mode_idx
        integer :: npoints
        real(rk) :: increments_squared
        real(rk) :: total_time
        real(rk) :: sigma_sum

        npoints = size(paths, 1)
        dim = size(paths, 2)

        if (npoints < 2) then
            write (*, '(a)') "At least two observations are required for sigma estimation"
            error stop
        end if

        total_time = delta * real(npoints - 1, rk)
        sigma_sum = 0.0_rk

        do mode_idx = 1, dim
            increments_squared = sum((paths(2:npoints, mode_idx) - paths(1:npoints - 1, mode_idx))**2)
            sigma_sum = sigma_sum + sqrt(increments_squared * lambdas(mode_idx)**(2.0_rk * gamma) / total_time)
        end do

        sigma_hat = sigma_sum / real(dim, rk)
    end subroutine estimate_sigma_from_qv

    subroutine estimate_joint_mle(paths, delta, interaction, lambdas, gamma, beta_hat, theta_hat)
        real(rk), intent(in) :: paths(:, :)
        real(rk), intent(in) :: delta, gamma
        real(rk), intent(in) :: interaction(:, :)
        real(rk), intent(in) :: lambdas(:)
        real(rk), intent(out) :: beta_hat, theta_hat

        logical :: success

        call try_estimate_joint_mle(paths, delta, interaction, lambdas, gamma, beta_hat, theta_hat, success)
        if (.not. success) then
            write (*, '(a)') "MLE normal equations are singular"
            error stop
        end if
    end subroutine estimate_joint_mle

    subroutine build_uniform_checkpoints(total_points, ncheckpoints, min_points, checkpoints)
        integer, intent(in) :: total_points, ncheckpoints, min_points
        integer, allocatable, intent(out) :: checkpoints(:)

        integer :: actual_points
        integer :: available_points
        integer :: point_idx

        if (total_points < 2) then
            write (*, '(a)') "At least two observations are required to build checkpoints"
            error stop
        end if

        if (ncheckpoints < 1) then
            write (*, '(a)') "ncheckpoints must be positive"
            error stop
        end if

        if (min_points < 2 .or. min_points > total_points) then
            write (*, '(a)') "min_points must belong to [2, total_points]"
            error stop
        end if

        available_points = total_points - min_points + 1
        actual_points = min(ncheckpoints, available_points)
        allocate (checkpoints(actual_points))

        if (actual_points == 1) then
            checkpoints(1) = total_points
            return
        end if

        do point_idx = 1, actual_points
            checkpoints(point_idx) = min_points + int(&
                & real(point_idx - 1, rk) * real(available_points - 1, rk) / &
                & real(actual_points - 1, rk) &
            )
        end do
    end subroutine build_uniform_checkpoints

    subroutine estimate_parameter_trajectory(&
        & paths, delta, interaction, lambdas, gamma, checkpoints, &
        & times, sigma_path, beta_path, theta_path &
    )
        real(rk), intent(in) :: paths(:, :)
        real(rk), intent(in) :: delta, gamma
        real(rk), intent(in) :: interaction(:, :)
        real(rk), intent(in) :: lambdas(:)
        integer, intent(in) :: checkpoints(:)
        real(rk), allocatable, intent(out) :: times(:)
        real(rk), allocatable, intent(out) :: sigma_path(:), beta_path(:), theta_path(:)

        integer :: point_idx
        integer :: nobs_at_point
        logical :: success
        real(rk) :: nan_value

        allocate (times(size(checkpoints)))
        allocate (sigma_path(size(checkpoints)))
        allocate (beta_path(size(checkpoints)))
        allocate (theta_path(size(checkpoints)))

        nan_value = ieee_value(0.0_rk, ieee_quiet_nan)

        do point_idx = 1, size(checkpoints)
            nobs_at_point = checkpoints(point_idx)

            if (nobs_at_point < 2 .or. nobs_at_point > size(paths, 1)) then
                write (*, '(a)') "Checkpoint index outside valid observation range"
                error stop
            end if

            times(point_idx) = delta * real(nobs_at_point - 1, rk)
            call estimate_sigma_from_qv(&
                & paths(1:nobs_at_point, :), delta, lambdas, gamma, sigma_path(point_idx) &
            )

            beta_path(point_idx) = nan_value
            theta_path(point_idx) = nan_value
            call try_estimate_joint_mle(&
                & paths(1:nobs_at_point, :), delta, interaction, lambdas, gamma, &
                & beta_path(point_idx), theta_path(point_idx), success &
            )
            if (.not. success) then
                beta_path(point_idx) = nan_value
                theta_path(point_idx) = nan_value
            end if
        end do
    end subroutine estimate_parameter_trajectory

    subroutine write_estimator_trajectory(&
        & filename, checkpoints, times, sigma_path, beta_path, theta_path, &
        & sigma_true, beta_true, theta_true &
    )
        character(len=*), intent(in) :: filename
        integer, intent(in) :: checkpoints(:)
        real(rk), intent(in) :: times(:), sigma_path(:), beta_path(:), theta_path(:)
        real(rk), intent(in) :: sigma_true, beta_true, theta_true

        integer :: point_idx
        integer :: unit

        if (size(checkpoints) /= size(times) .or. size(times) /= size(sigma_path) .or. &
            & size(sigma_path) /= size(beta_path) .or. size(beta_path) /= size(theta_path)) then
            write (*, '(a)') "Estimator trajectory arrays must have the same length"
            error stop
        end if

        open (newunit=unit, file=filename, status="replace", action="write")
        write (unit, '(a)') "n_obs,time,sigma_hat,sigma_true,beta_hat,beta_true,theta_hat,theta_true"

        do point_idx = 1, size(checkpoints)
            write (unit, '(i0,",",es26.17e3,",",es26.17e3,",",es26.17e3,",",es26.17e3,",",es26.17e3,",",es26.17e3,",",es26.17e3)') &
                & checkpoints(point_idx), &
                & times(point_idx), &
                & sigma_path(point_idx), sigma_true, &
                & beta_path(point_idx), beta_true, &
                & theta_path(point_idx), theta_true
        end do

        close (unit)
    end subroutine write_estimator_trajectory

    subroutine compute_sufficient_statistics(paths, delta, interaction, lambdas, gamma, stats)
        real(rk), intent(in) :: paths(:, :)
        real(rk), intent(in) :: delta, gamma
        real(rk), intent(in) :: interaction(:, :)
        real(rk), intent(in) :: lambdas(:)
        real(rk), intent(out) :: stats(5)

        real(rk), allocatable :: coupled_paths(:, :)
        real(rk), allocatable :: t1(:), t2(:), t3(:), t4(:), t5(:)
        real(rk), allocatable :: weight_1(:), weight_2(:), weight_3(:)
        integer :: dim
        integer :: mode_idx
        integer :: npoints

        npoints = size(paths, 1)
        dim = size(paths, 2)

        allocate (coupled_paths(npoints, dim))
        allocate (t1(dim), t2(dim), t3(dim), t4(dim), t5(dim))
        allocate (weight_1(dim), weight_2(dim), weight_3(dim))

        coupled_paths = matmul(paths, transpose(interaction))

        do mode_idx = 1, dim
            call integrate_ito(paths(:, mode_idx), paths(:, mode_idx), t1(mode_idx))
            call integrate_ito(coupled_paths(:, mode_idx), paths(:, mode_idx), t2(mode_idx))
            call integrate_trapezoidal(paths(:, mode_idx) * paths(:, mode_idx), delta, t3(mode_idx))
            call integrate_trapezoidal(paths(:, mode_idx) * coupled_paths(:, mode_idx), delta, t4(mode_idx))
            call integrate_trapezoidal(coupled_paths(:, mode_idx) * coupled_paths(:, mode_idx), delta, t5(mode_idx))
        end do

        weight_1 = lambdas**(1.0_rk + 2.0_rk * gamma)
        weight_2 = lambdas**(2.0_rk * gamma)
        weight_3 = lambdas**(2.0_rk + 2.0_rk * gamma)

        stats(1) = sum(t1 * weight_1)
        stats(2) = sum(t2 * weight_2)
        stats(3) = sum(t3 * weight_3)
        stats(4) = sum(t4 * weight_1)
        stats(5) = sum(t5 * weight_2)
    end subroutine compute_sufficient_statistics

    subroutine try_estimate_joint_mle(paths, delta, interaction, lambdas, gamma, beta_hat, theta_hat, success)
        real(rk), intent(in) :: paths(:, :)
        real(rk), intent(in) :: delta, gamma
        real(rk), intent(in) :: interaction(:, :)
        real(rk), intent(in) :: lambdas(:)
        real(rk), intent(out) :: beta_hat, theta_hat
        logical, intent(out) :: success

        real(rk) :: stats(5)
        real(rk) :: denominator
        real(rk) :: scale
        real(rk) :: tolerance

        call compute_sufficient_statistics(paths, delta, interaction, lambdas, gamma, stats)

        denominator = stats(4)**2 - stats(3) * stats(5)
        scale = max(1.0_rk, abs(stats(4)**2) + abs(stats(3) * stats(5)))
        tolerance = sqrt(epsilon(1.0_rk)) * scale

        if (abs(denominator) <= tolerance) then
            success = .false.
            beta_hat = 0.0_rk
            theta_hat = 0.0_rk
            return
        end if

        success = .true.
        beta_hat = (stats(1) * stats(5) - stats(2) * stats(4)) / denominator
        theta_hat = (stats(2) * stats(3) - stats(1) * stats(4)) / denominator
    end subroutine try_estimate_joint_mle

    pure subroutine integrate_ito(integrand, process, integral)
        real(rk), intent(in) :: integrand(:), process(:)
        real(rk), intent(out) :: integral

        integer :: npoints

        npoints = size(integrand)
        if (npoints < 2) then
            integral = 0.0_rk
            return
        end if

        integral = dot_product(integrand(1:npoints - 1), process(2:npoints) - process(1:npoints - 1))
    end subroutine integrate_ito

    pure subroutine integrate_trapezoidal(path, dt, integral)
        real(rk), intent(in) :: path(:)
        real(rk), intent(in) :: dt
        real(rk), intent(out) :: integral

        integer :: npoints

        npoints = size(path)
        if (npoints < 2) then
            integral = 0.0_rk
            return
        end if

        integral = 0.5_rk * dt * (sum(path(1:npoints - 1)) + sum(path(2:npoints)))
    end subroutine integrate_trapezoidal

    subroutine print_estimation_report(sigma_true, beta_true, theta_true, result)
        real(rk), intent(in) :: sigma_true, beta_true, theta_true
        type(estimation_result_t), intent(in) :: result

        print '(a)', ''
        print '(a)', repeat('=', 86)
        print '(a)', 'MLE ESTIMATION REPORT'
        print '(a)', repeat('-', 86)
        write (*, '(a10,3x,a12,3x,a12,3x,a12,3x,a13)') &
            & 'Parameter', 'True Value', 'Estimate', 'Abs Error', 'Rel Error (%)'
        print '(a)', repeat('-', 86)

        call print_parameter_line('sigma', sigma_true, result%sigma_hat)
        call print_parameter_line('beta', beta_true, result%beta_hat)
        call print_parameter_line('theta', theta_true, result%theta_hat)

        print '(a)', repeat('-', 86)
        print '(a)', 'Runtime (seconds)'
        write (*, '(a20,2x,f12.6)') 'Set-up', result%setup_time
        write (*, '(a20,2x,f12.6)') 'MLE estimation', result%mle_time
        print '(a)', repeat('=', 86)
    end subroutine print_estimation_report

    subroutine print_parameter_line(name, true_value, estimate)
        character(len=*), intent(in) :: name
        real(rk), intent(in) :: true_value, estimate

        real(rk) :: abs_err, rel_err

        abs_err = abs(estimate - true_value)
        if (abs(true_value) > tiny(true_value)) then
            rel_err = 100.0_rk * abs_err / abs(true_value)
        else
            rel_err = 0.0_rk
        end if

        write (*, '(a10,3x,f12.6,3x,f12.6,3x,es12.4,3x,f13.4)') &
            & trim(name), true_value, estimate, abs_err, rel_err
    end subroutine print_parameter_line

end module mle_estimation_mod
