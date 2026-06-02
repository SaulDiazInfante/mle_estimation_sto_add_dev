!> @file solution_reconstruction_mod.f90
!! @brief Physical-space reconstruction of modal SPDE solution paths.
!!
!! The solver evolves the Galerkin coefficients from equation (3) in the
!! manuscript. This module evaluates the expansion
!! u(t,x,y) = sum_k u_k(t) h_k(x,y) on plotting grids for snapshots and
!! animation frames.
module solution_reconstruction_mod
    use model_types_mod, only: dp
    use model_types_mod, only: spatial_grid_t
    use spectral_operators_mod, only: evaluate_basis_function
    implicit none
    private

    public :: build_endpoint_plot_coordinates
    public :: evaluate_cosine_mode
    public :: reconstruct_solution_field
    public :: reconstruct_solution_snapshots
    public :: reconstruct_solution_value

contains

    !> Builds endpoint-inclusive coordinates, matching the legacy field writer.
    subroutine build_endpoint_plot_coordinates(&
        grid, n_plot_x, n_plot_y, x_coordinates, y_coordinates &
    )
        type(spatial_grid_t), intent(in) :: grid
        integer, intent(in) :: n_plot_x
        integer, intent(in) :: n_plot_y
        real(dp), allocatable, intent(out) :: x_coordinates(:)
        real(dp), allocatable, intent(out) :: y_coordinates(:)

        integer :: ix
        integer :: iy
        real(dp) :: dx
        real(dp) :: dy

        if (n_plot_x < 2 .or. n_plot_y < 2) then
            write (*, '(a)') "Plot grids must have at least two points per axis."
            error stop
        end if

        allocate (x_coordinates(n_plot_x))
        allocate (y_coordinates(n_plot_y))

        dx = grid%length_x / real(n_plot_x - 1, dp)
        dy = grid%length_y / real(n_plot_y - 1, dp)

        do ix = 1, n_plot_x
            x_coordinates(ix) = real(ix - 1, dp) * dx
        end do

        do iy = 1, n_plot_y
            y_coordinates(iy) = real(iy - 1, dp) * dy
        end do
    end subroutine build_endpoint_plot_coordinates

    !> Evaluates the two-dimensional Neumann cosine basis h_ij(x,y).
    pure real(dp) function evaluate_cosine_mode(&
        mode_i, mode_j, x, y, grid &
    ) result(mode_value)
        integer, intent(in) :: mode_i
        integer, intent(in) :: mode_j
        real(dp), intent(in) :: x
        real(dp), intent(in) :: y
        type(spatial_grid_t), intent(in) :: grid

        mode_value = evaluate_basis_function(mode_i, x, grid%length_x) * &
            evaluate_basis_function(mode_j, y, grid%length_y)
    end function evaluate_cosine_mode

    !> Reconstructs u(t,x,y) from one vector of modal coefficients.
    subroutine reconstruct_solution_value(&
        grid, mode_pairs, modal_coefficients, x, y, value &
    )
        type(spatial_grid_t), intent(in) :: grid
        integer, intent(in) :: mode_pairs(:, :)
        real(dp), intent(in) :: modal_coefficients(:)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: y
        real(dp), intent(out) :: value

        integer :: mode_index

        call validate_reconstruction_inputs(mode_pairs, modal_coefficients)
        call validate_point_in_domain(grid, x, y)

        value = 0.0_dp
        do mode_index = 1, size(modal_coefficients)
            value = value + modal_coefficients(mode_index) * &
                evaluate_cosine_mode(&
                    mode_pairs(mode_index, 1), mode_pairs(mode_index, 2), &
                    x, y, grid &
                )
        end do
    end subroutine reconstruct_solution_value

    !> Reconstructs one modal state on an arbitrary tensor-product grid.
    subroutine reconstruct_solution_field(&
        grid, mode_pairs, modal_coefficients, x_coordinates, y_coordinates, &
        field &
    )
        type(spatial_grid_t), intent(in) :: grid
        integer, intent(in) :: mode_pairs(:, :)
        real(dp), intent(in) :: modal_coefficients(:)
        real(dp), intent(in) :: x_coordinates(:)
        real(dp), intent(in) :: y_coordinates(:)
        real(dp), allocatable, intent(out) :: field(:, :)

        real(dp), allocatable :: basis_x(:, :)
        real(dp), allocatable :: basis_y(:, :)

        call validate_reconstruction_inputs(mode_pairs, modal_coefficients)
        call validate_plot_coordinates(grid, x_coordinates, y_coordinates)
        call build_modal_basis_values(&
            grid, mode_pairs, x_coordinates, y_coordinates, basis_x, basis_y &
        )

        allocate (field(size(x_coordinates), size(y_coordinates)))
        call fill_reconstructed_field(modal_coefficients, basis_x, basis_y, field)
    end subroutine reconstruct_solution_field

    !> Reconstructs every saved modal snapshot on an arbitrary plotting grid.
    subroutine reconstruct_solution_snapshots(&
        grid, mode_pairs, snapshot_states, x_coordinates, y_coordinates, fields &
    )
        type(spatial_grid_t), intent(in) :: grid
        integer, intent(in) :: mode_pairs(:, :)
        real(dp), intent(in) :: snapshot_states(:, :)
        real(dp), intent(in) :: x_coordinates(:)
        real(dp), intent(in) :: y_coordinates(:)
        real(dp), allocatable, intent(out) :: fields(:, :, :)

        integer :: snapshot_index
        real(dp), allocatable :: basis_x(:, :)
        real(dp), allocatable :: basis_y(:, :)

        if (size(snapshot_states, 2) /= size(mode_pairs, 1)) then
            write (*, '(a)') &
                "Snapshot state width must match the number of mode pairs."
            error stop
        end if

        if (size(snapshot_states, 1) < 1) then
            write (*, '(a)') "At least one snapshot state is required."
            error stop
        end if

        call validate_reconstruction_inputs(mode_pairs, snapshot_states(1, :))
        call validate_plot_coordinates(grid, x_coordinates, y_coordinates)
        call build_modal_basis_values(&
            grid, mode_pairs, x_coordinates, y_coordinates, basis_x, basis_y &
        )

        allocate (&
            fields(size(x_coordinates), size(y_coordinates), size(snapshot_states, 1)) &
        )

        do snapshot_index = 1, size(snapshot_states, 1)
            call fill_reconstructed_field(&
                snapshot_states(snapshot_index, :), basis_x, basis_y, &
                fields(:, :, snapshot_index) &
            )
        end do
    end subroutine reconstruct_solution_snapshots

    subroutine build_modal_basis_values(&
        grid, mode_pairs, x_coordinates, y_coordinates, basis_x, basis_y &
    )
        type(spatial_grid_t), intent(in) :: grid
        integer, intent(in) :: mode_pairs(:, :)
        real(dp), intent(in) :: x_coordinates(:)
        real(dp), intent(in) :: y_coordinates(:)
        real(dp), allocatable, intent(out) :: basis_x(:, :)
        real(dp), allocatable, intent(out) :: basis_y(:, :)

        integer :: coordinate_index
        integer :: mode_index

        allocate (basis_x(size(x_coordinates), size(mode_pairs, 1)))
        allocate (basis_y(size(y_coordinates), size(mode_pairs, 1)))

        do mode_index = 1, size(mode_pairs, 1)
            do coordinate_index = 1, size(x_coordinates)
                basis_x(coordinate_index, mode_index) = evaluate_basis_function(&
                    mode_pairs(mode_index, 1), x_coordinates(coordinate_index), &
                    grid%length_x &
                )
            end do

            do coordinate_index = 1, size(y_coordinates)
                basis_y(coordinate_index, mode_index) = evaluate_basis_function(&
                    mode_pairs(mode_index, 2), y_coordinates(coordinate_index), &
                    grid%length_y &
                )
            end do
        end do
    end subroutine build_modal_basis_values

    subroutine fill_reconstructed_field(&
        modal_coefficients, basis_x, basis_y, field &
    )
        real(dp), intent(in) :: modal_coefficients(:)
        real(dp), intent(in) :: basis_x(:, :)
        real(dp), intent(in) :: basis_y(:, :)
        real(dp), intent(out) :: field(:, :)

        integer :: ix
        integer :: iy
        integer :: mode_index

        field = 0.0_dp
        do mode_index = 1, size(modal_coefficients)
            do iy = 1, size(basis_y, 1)
                do ix = 1, size(basis_x, 1)
                    field(ix, iy) = field(ix, iy) + &
                        modal_coefficients(mode_index) * &
                        basis_x(ix, mode_index) * basis_y(iy, mode_index)
                end do
            end do
        end do
    end subroutine fill_reconstructed_field

    subroutine validate_reconstruction_inputs(mode_pairs, modal_coefficients)
        integer, intent(in) :: mode_pairs(:, :)
        real(dp), intent(in) :: modal_coefficients(:)

        if (size(mode_pairs, 2) /= 2) then
            write (*, '(a)') "Mode pairs must have shape (:,2)."
            error stop
        end if

        if (size(mode_pairs, 1) /= size(modal_coefficients)) then
            write (*, '(a)') &
                "Modal coefficient count must match the number of mode pairs."
            error stop
        end if
    end subroutine validate_reconstruction_inputs

    subroutine validate_plot_coordinates(grid, x_coordinates, y_coordinates)
        type(spatial_grid_t), intent(in) :: grid
        real(dp), intent(in) :: x_coordinates(:)
        real(dp), intent(in) :: y_coordinates(:)

        integer :: coordinate_index

        if (size(x_coordinates) < 1 .or. size(y_coordinates) < 1) then
            write (*, '(a)') "Plot coordinate arrays must be non-empty."
            error stop
        end if

        do coordinate_index = 1, size(x_coordinates)
            call validate_coordinate(&
                "x", x_coordinates(coordinate_index), grid%length_x &
            )
        end do

        do coordinate_index = 1, size(y_coordinates)
            call validate_coordinate(&
                "y", y_coordinates(coordinate_index), grid%length_y &
            )
        end do
    end subroutine validate_plot_coordinates

    subroutine validate_point_in_domain(grid, x, y)
        type(spatial_grid_t), intent(in) :: grid
        real(dp), intent(in) :: x
        real(dp), intent(in) :: y

        call validate_coordinate("x", x, grid%length_x)
        call validate_coordinate("y", y, grid%length_y)
    end subroutine validate_point_in_domain

    subroutine validate_coordinate(axis_name, coordinate, domain_length)
        character(len=*), intent(in) :: axis_name
        real(dp), intent(in) :: coordinate
        real(dp), intent(in) :: domain_length

        real(dp) :: tolerance

        tolerance = max(1000.0_dp * epsilon(1.0_dp) * domain_length, 1.0e-12_dp)

        if (coordinate < -tolerance .or. coordinate > domain_length + tolerance) then
            write (*, '(3a)') &
                "Plot ", trim(axis_name), "-coordinate lies outside the domain."
            error stop
        end if
    end subroutine validate_coordinate

end module solution_reconstruction_mod
