!> @file model_types_mod.f90
!! @brief Shared derived types and precision definitions for the maintained code.
!!
!! This file centralizes the public data structures used to configure the
!! spatial grid, define the SDE parameters, store assembled operators, and
!! return estimation results.
!> @brief Declares the core public types used across simulation and estimation.
module model_types_mod
    use iso_fortran_env, only: real64
    implicit none
    private

    !> Real kind used for floating-point work throughout the maintained solver.
    integer, parameter, public :: dp = real64

    !> Spatial and physical configuration for the rectangular computational 
    !! grid.
    type, public :: spatial_grid_t
        integer :: nx = 10                !!< Number of grid cells or modes in 
                                          !! x.
        integer :: ny = 10                !!< Number of grid cells or modes in 
                                          !! y.
        integer :: velocity_mode_x = 1    !!< Velocity forcing mode in x.
        integer :: velocity_mode_y = 1    !!< Velocity forcing mode in y.
        real(dp) :: length_x = 5.0_dp     !!< Domain length in x.
        real(dp) :: length_y = 5.0_dp     !!< Domain length in y.
        real(dp) :: gamma = 2.0_dp        !!< Fractional exponent used in 
                                          !! diffusion scaling.
    end type spatial_grid_t

    !> Parameters that control the SDE trajectory generation.
    type, public :: sde_parameters_t
        integer :: n_observations = 1000000 !!< Number of saved observations in the path.
        real(dp) :: time_step = 1.0e-5_dp   !!< Time increment between consecutive states.
        real(dp) :: beta = 0.1_dp           !!< Dissipative drift coefficient.
        real(dp) :: theta = 1.0_dp          !!< Coupling coefficient.
        real(dp) :: sigma = 0.1_dp          !!< Diffusion amplitude.
    end type sde_parameters_t

    !> Spectral operators and modal data assembled from the grid description.
    type, public :: spectral_operator_set_t
        integer :: state_dimension = 0                 !!< Number of active 
                                                       !! non-zero modes.
        integer, allocatable :: mode_pairs(:, :)       !!< Modal indices sorted 
                                                       !! by eigenvalue.
        real(dp), allocatable :: initial_state(:)      !!< Projected initial 
                                                       !! modal state.
        real(dp), allocatable :: eigenvalues(:)        !!< Eigenvalues of 
                                                       !! the spectral basis.
        real(dp), allocatable :: diffusion_diagonal(:) !!< Diagonal diffusion 
                                                       !! weights.
        real(dp), allocatable :: interaction_matrix(:, :) !!< Linear        
                                                       !! interaction operator.
    end type spectral_operator_set_t

    !> Aggregated parameter estimates and timing diagnostics.
    type, public :: parameter_estimates_t
        real(dp) :: sigma_hat = 0.0_dp       !!< Estimated diffusion amplitude.
        real(dp) :: beta_hat = 0.0_dp        !!< Estimated dissipative 
                                             !!  coefficient.
        real(dp) :: theta_hat = 0.0_dp       !!< Estimated coupling coefficient.
        real(dp) :: setup_time = 0.0_dp      !!< CPU time spent assembling and 
                                             !!  simulating.
        real(dp) :: estimation_time = 0.0_dp !!< CPU time spent in the 
                                             !!  estimation stage.
    end type parameter_estimates_t

    public :: get_state_dimension

    contains

    !> Computes the number of active modal coefficients for the chosen grid.
    pure integer function get_state_dimension(grid) result(n_state)
        type(spatial_grid_t), intent(in) :: grid

        n_state = grid%nx * grid%ny - 1
    end function get_state_dimension
end module model_types_mod
