!> @file spectral_operators_mod.f90
!! @brief Assembly of spectral operators, mode ordering, and initial conditions.
!> @brief Builds the modal representation needed by the solver and estimator.
module spectral_operators_mod
    use model_types_mod, only: dp
    use model_types_mod, only: get_state_dimension
    use model_types_mod, only: spatial_grid_t
    use model_types_mod, only: spectral_operator_set_t
    implicit none
    private

    real(dp), parameter :: pi_dp = acos(-1.0_dp)

    public :: assemble_problem_operators

contains

    !> Builds eigenvalues, diffusion terms, interactions, and the initial modal state.
    subroutine assemble_problem_operators(grid, operators)
        type(spatial_grid_t), intent(in) :: grid
        type(spectral_operator_set_t), intent(out) :: operators

        real(dp), allocatable :: fluctuation_field(:, :)

        operators%state_dimension = get_state_dimension(grid)

        call build_sorted_mode_pairs(&
            grid, operators%mode_pairs, operators%eigenvalues &
        )

        allocate (operators%initial_state(operators%state_dimension))
        allocate (operators%interaction_matrix( &
            operators%state_dimension, operators%state_dimension &
        ))
        allocate (operators%diffusion_diagonal(operators%state_dimension))

        call build_fluctuation_field(grid, fluctuation_field)
        call project_field_onto_modes(&
            grid, fluctuation_field, operators%mode_pairs, &
            operators%initial_state &
        )
        call assemble_interaction_and_diffusion(&
            grid, operators%mode_pairs, operators%eigenvalues, &
            operators%interaction_matrix, operators%diffusion_diagonal &
        )
    end subroutine assemble_problem_operators

    pure real(dp) function evaluate_initial_bump(x, y, grid) result(value)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: y
        type(spatial_grid_t), intent(in) :: grid

        real(dp), parameter :: support_fraction = 0.25_dp
        real(dp) :: ellipse_radius
        real(dp) :: x_center
        real(dp) :: x_scaled
        real(dp) :: y_center
        real(dp) :: y_scaled

        x_center = grid%length_x / 4.0_dp
        y_center = 0.6_dp * grid%length_y

        x_scaled = (x - x_center) / (support_fraction * grid%length_x)
        y_scaled = (y - y_center) / (support_fraction * grid%length_y)
        ellipse_radius = sqrt(x_scaled**2 + y_scaled**2)

        value = 0.0_dp
        if (ellipse_radius <= 1.0_dp) then
            value = 0.5_dp * (1.0_dp + cos(pi_dp * ellipse_radius))
        end if
    end function evaluate_initial_bump

    subroutine build_fluctuation_field(grid, fluctuation_field)
        type(spatial_grid_t), intent(in) :: grid
        real(dp), allocatable, intent(out) :: fluctuation_field(:, :)

        real(dp), allocatable :: initial_field(:, :)
        integer :: ix
        integer :: iy
        real(dp) :: dx
        real(dp) :: dy
        real(dp) :: field_mean
        real(dp) :: x_coord
        real(dp) :: y_coord

        allocate (initial_field(grid%nx, grid%ny))
        allocate (fluctuation_field(grid%nx, grid%ny))

        dx = grid%length_x / real(grid%nx, dp)
        dy = grid%length_y / real(grid%ny, dp)

        do iy = 1, grid%ny
            y_coord = dy * (real(iy, dp) - 0.5_dp)
            do ix = 1, grid%nx
                x_coord = dx * (real(ix, dp) - 0.5_dp)
                initial_field(ix, iy) = &
                    evaluate_initial_bump(x_coord, y_coord, grid)
            end do
        end do

        field_mean = sum(initial_field) / real(size(initial_field), dp)
        fluctuation_field = initial_field - field_mean
    end subroutine build_fluctuation_field

    subroutine build_sorted_mode_pairs(grid, mode_pairs, eigenvalues)
        type(spatial_grid_t), intent(in) :: grid
        integer, allocatable, intent(out) :: mode_pairs(:, :)
        real(dp), allocatable, intent(out) :: eigenvalues(:)

        integer :: mode_i
        integer :: mode_j
        integer :: mode_index
        integer :: n_state

        n_state = get_state_dimension(grid)
        allocate (mode_pairs(n_state, 2))
        allocate (eigenvalues(n_state))

        mode_index = 0
        do mode_i = 0, grid%nx - 1
            do mode_j = 0, grid%ny - 1
                if (mode_i + mode_j == 0) cycle

                mode_index = mode_index + 1
                mode_pairs(mode_index, 1) = mode_i
                mode_pairs(mode_index, 2) = mode_j
                eigenvalues(mode_index) = &
                    compute_eigenvalue(mode_i, mode_j, grid)
            end do
        end do

        call sort_modes_by_eigenvalue(eigenvalues, mode_pairs)
    end subroutine build_sorted_mode_pairs

    pure real(dp) function compute_eigenvalue(&
        mode_i, mode_j, grid &
    ) result(eigenvalue)
        integer, intent(in) :: mode_i
        integer, intent(in) :: mode_j
        type(spatial_grid_t), intent(in) :: grid

        eigenvalue = (pi_dp / grid%length_x)**2 * real(mode_i, dp)**2 + &
            (pi_dp / grid%length_y)**2 * real(mode_j, dp)**2
    end function compute_eigenvalue

    subroutine sort_modes_by_eigenvalue(eigenvalues, mode_pairs)
        real(dp), intent(inout) :: eigenvalues(:)
        integer, intent(inout) :: mode_pairs(:, :)

        integer :: candidate_index
        integer :: left_index
        integer :: minimum_index
        integer :: saved_pair(2)
        real(dp) :: saved_eigenvalue

        do left_index = 1, size(eigenvalues) - 1
            minimum_index = left_index

            do candidate_index = left_index + 1, size(eigenvalues)
                if (eigenvalues(candidate_index) < &
                    eigenvalues(minimum_index)) then
                    minimum_index = candidate_index
                end if
            end do

            if (minimum_index /= left_index) then
                saved_eigenvalue = eigenvalues(left_index)
                eigenvalues(left_index) = eigenvalues(minimum_index)
                eigenvalues(minimum_index) = saved_eigenvalue

                saved_pair = mode_pairs(left_index, :)
                mode_pairs(left_index, :) = mode_pairs(minimum_index, :)
                mode_pairs(minimum_index, :) = saved_pair
            end if
        end do
    end subroutine sort_modes_by_eigenvalue

    subroutine project_field_onto_modes(&
        grid, field, mode_pairs, modal_coefficients &
    )
        type(spatial_grid_t), intent(in) :: grid
        real(dp), intent(in) :: field(:, :)
        integer, intent(in) :: mode_pairs(:, :)
        real(dp), intent(out) :: modal_coefficients(:)

        integer :: mode_index

        do mode_index = 1, size(modal_coefficients)
            modal_coefficients(mode_index) = compute_mode_integral(&
                field, mode_pairs(mode_index, 1), mode_pairs(mode_index, 2), &
                grid &
            )
            modal_coefficients(mode_index) = &
                modal_coefficients(mode_index) / &
                (grid%length_x * grid%length_y)
        end do
    end subroutine project_field_onto_modes

    pure real(dp) function compute_mode_integral(&
        field, mode_i, mode_j, grid &
    ) result(integral_value)
        real(dp), intent(in) :: field(:, :)
        integer, intent(in) :: mode_i
        integer, intent(in) :: mode_j
        type(spatial_grid_t), intent(in) :: grid

        integer :: ix
        integer :: iy
        real(dp) :: basis_x
        real(dp) :: basis_y
        real(dp) :: cell_area
        real(dp) :: dx
        real(dp) :: dy
        real(dp) :: wave_number_x
        real(dp) :: wave_number_y
        real(dp) :: x_coord
        real(dp) :: y_coord

        dx = grid%length_x / real(grid%nx, dp)
        dy = grid%length_y / real(grid%ny, dp)
        cell_area = dx * dy

        wave_number_x = pi_dp * real(mode_i, dp) / grid%length_x
        wave_number_y = pi_dp * real(mode_j, dp) / grid%length_y

        integral_value = 0.0_dp
        do iy = 1, grid%ny
            y_coord = dy * (real(iy, dp) - 0.5_dp)
            basis_y = get_basis_normalization(mode_j) * &
                cos(wave_number_y * y_coord)

            do ix = 1, grid%nx
                x_coord = dx * (real(ix, dp) - 0.5_dp)
                basis_x = get_basis_normalization(mode_i) * &
                    cos(wave_number_x * x_coord)

                integral_value = integral_value + &
                    cell_area * field(ix, iy) * basis_x * basis_y
            end do
        end do
    end function compute_mode_integral

    pure real(dp) function get_basis_normalization(mode_number) &
        result(normalization)
        integer, intent(in) :: mode_number

        if (mode_number == 0) then
            normalization = 1.0_dp
        else
            normalization = sqrt(2.0_dp)
        end if
    end function get_basis_normalization

    subroutine assemble_interaction_and_diffusion(&
        grid, mode_pairs, eigenvalues, interaction_matrix, &
        diffusion_diagonal &
    )
        type(spatial_grid_t), intent(in) :: grid
        integer, intent(in) :: mode_pairs(:, :)
        real(dp), intent(in) :: eigenvalues(:)
        real(dp), intent(out) :: interaction_matrix(:, :)
        real(dp), intent(out) :: diffusion_diagonal(:)

        integer :: source_index
        integer :: target_index

        interaction_matrix = 0.0_dp

        do source_index = 1, size(eigenvalues)
            diffusion_diagonal(source_index) = &
                eigenvalues(source_index)**(-grid%gamma)

            do target_index = 1, size(eigenvalues)
                interaction_matrix(target_index, source_index) = &
                    compute_interaction_coefficient(&
                        mode_pairs(source_index, 1), &
                        mode_pairs(source_index, 2), &
                        mode_pairs(target_index, 1), &
                        mode_pairs(target_index, 2), &
                        grid &
                    )
            end do
        end do
    end subroutine assemble_interaction_and_diffusion

    pure real(dp) function compute_interaction_coefficient(&
        mode_i, mode_j, mode_k, mode_l, grid &
    ) result(coefficient)
        integer, intent(in) :: mode_i
        integer, intent(in) :: mode_j
        integer, intent(in) :: mode_k
        integer, intent(in) :: mode_l
        type(spatial_grid_t), intent(in) :: grid

        logical :: i_minus_r_matches_k
        logical :: i_plus_r_matches_k
        logical :: j_minus_s_matches_l
        logical :: j_plus_s_matches_l
        real(dp) :: base_coefficient
        real(dp) :: coefficient_i_zero
        real(dp) :: coefficient_j_zero
        real(dp) :: coefficient_k_zero
        real(dp) :: coefficient_l_zero
        real(dp) :: interaction_term
        real(dp) :: length_x
        real(dp) :: length_y
        integer :: r_mode
        integer :: s_mode
        real(dp) :: sqrt_half_two

        coefficient = 0.0_dp
        length_x = grid%length_x
        length_y = grid%length_y
        r_mode = grid%velocity_mode_x
        s_mode = grid%velocity_mode_y

        base_coefficient = pi_dp**2 / (2.0_dp * length_x * length_y)
        sqrt_half_two = sqrt(2.0_dp) / 2.0_dp

        if (mode_i * mode_j * mode_k * mode_l > 0) then
            i_minus_r_matches_k = ((mode_i - r_mode)**2 == mode_k**2)
            i_plus_r_matches_k = ((mode_i + r_mode)**2 == mode_k**2)
            j_minus_s_matches_l = ((mode_j - s_mode)**2 == mode_l**2)
            j_plus_s_matches_l = ((mode_j + s_mode)**2 == mode_l**2)

            interaction_term = real(mode_j * r_mode - mode_i * s_mode, dp)

            if (i_minus_r_matches_k .and. j_minus_s_matches_l) then
                coefficient = -base_coefficient * interaction_term
            end if

            if (i_plus_r_matches_k .and. j_minus_s_matches_l) then
                coefficient = -base_coefficient * &
                    real(mode_j * r_mode + mode_i * s_mode, dp)
            end if

            if (i_plus_r_matches_k .and. j_plus_s_matches_l) then
                coefficient = base_coefficient * interaction_term
            end if

            if (i_minus_r_matches_k .and. j_plus_s_matches_l) then
                coefficient = base_coefficient * &
                    real(mode_j * r_mode + mode_i * s_mode, dp)
            end if
        end if

        if (mode_i == 0 .and. mode_j * mode_k * mode_l > 0) then
            coefficient_i_zero = sqrt_half_two * pi_dp**2 * &
                real(mode_j * r_mode, dp) / (length_x * length_y)

            if (mode_k == r_mode .and. (mode_j - s_mode)**2 == mode_l**2) then
                coefficient = -coefficient_i_zero
            end if

            if (mode_k == r_mode .and. (mode_j + s_mode)**2 == mode_l**2) then
                coefficient = coefficient_i_zero
            end if
        end if

        if (mode_k == 0 .and. mode_i * mode_j * mode_l > 0) then
            coefficient_k_zero = sqrt_half_two * pi_dp**2 / &
                (length_x * length_y)

            if (mode_i == r_mode .and. (mode_j - s_mode)**2 == mode_l**2) then
                coefficient = -coefficient_k_zero * &
                    real(mode_j * r_mode - mode_i * s_mode, dp)
            end if

            if (mode_i == r_mode .and. (mode_j + s_mode)**2 == mode_l**2) then
                coefficient = coefficient_k_zero * &
                    real(mode_j * r_mode + mode_i * s_mode, dp)
            end if
        end if

        if (mode_j == 0 .and. mode_i * mode_k * mode_l > 0) then
            coefficient_j_zero = sqrt_half_two * pi_dp**2 * &
                real(mode_i * s_mode, dp) / (length_x * length_y)

            if (mode_l == s_mode .and. (mode_i - r_mode)**2 == mode_k**2) then
                coefficient = coefficient_j_zero
            end if

            if (mode_l == s_mode .and. (mode_i + r_mode)**2 == mode_k**2) then
                coefficient = -coefficient_j_zero
            end if
        end if

        if (mode_l == 0 .and. mode_i * mode_j * mode_k > 0) then
            coefficient_l_zero = sqrt_half_two * pi_dp**2 / &
                (length_x * length_y)

            if (mode_j == s_mode .and. (mode_i - r_mode)**2 == mode_k**2) then
                coefficient = -coefficient_l_zero * &
                    real(mode_j * r_mode - mode_i * s_mode, dp)
            end if

            if (mode_j == s_mode .and. (mode_i + r_mode)**2 == mode_k**2) then
                coefficient = -coefficient_l_zero * &
                    real(mode_j * r_mode + mode_i * s_mode, dp)
            end if
        end if

        if (mode_i == 0 .and. mode_l == 0 .and. mode_j * mode_k > 0) then
            if (mode_k == r_mode .and. mode_j == s_mode) then
                coefficient = -pi_dp**2 * real(mode_j * r_mode, dp) / &
                    (length_x * length_y)
            end if
        end if

        if (mode_j == 0 .and. mode_k == 0 .and. mode_i * mode_l > 0) then
            if (mode_i == r_mode .and. mode_l == s_mode) then
                coefficient = pi_dp**2 * real(mode_i * s_mode, dp) / &
                    (length_x * length_y)
            end if
        end if
    end function compute_interaction_coefficient

end module spectral_operators_mod
