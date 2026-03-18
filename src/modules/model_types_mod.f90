module model_types_mod
    use iso_fortran_env, only: real64
    implicit none
    private

    integer, parameter, public :: dp = real64

    type, public :: spatial_grid_t
        integer :: nx = 10
        integer :: ny = 10
        integer :: velocity_mode_x = 1
        integer :: velocity_mode_y = 1
        real(dp) :: length_x = 5.0_dp
        real(dp) :: length_y = 5.0_dp
        real(dp) :: gamma = 2.0_dp
    end type spatial_grid_t

    type, public :: sde_parameters_t
        integer :: n_observations = 1000000
        real(dp) :: time_step = 1.0e-5_dp
        real(dp) :: beta = 0.1_dp
        real(dp) :: theta = 1.0_dp
        real(dp) :: sigma = 0.1_dp
    end type sde_parameters_t

    type, public :: spectral_operator_set_t
        integer :: state_dimension = 0
        integer, allocatable :: mode_pairs(:, :)
        real(dp), allocatable :: initial_state(:)
        real(dp), allocatable :: eigenvalues(:)
        real(dp), allocatable :: diffusion_diagonal(:)
        real(dp), allocatable :: interaction_matrix(:, :)
    end type spectral_operator_set_t

    type, public :: parameter_estimates_t
        real(dp) :: sigma_hat = 0.0_dp
        real(dp) :: beta_hat = 0.0_dp
        real(dp) :: theta_hat = 0.0_dp
        real(dp) :: setup_time = 0.0_dp
        real(dp) :: estimation_time = 0.0_dp
    end type parameter_estimates_t

    public :: get_state_dimension

contains

    pure integer function get_state_dimension(grid) result(n_state)
        type(spatial_grid_t), intent(in) :: grid

        n_state = grid%nx * grid%ny - 1
    end function get_state_dimension

end module model_types_mod
