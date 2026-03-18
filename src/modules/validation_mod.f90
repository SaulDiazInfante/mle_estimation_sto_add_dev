!> @file validation_mod.f90
!! @brief Finite-value validation helpers for vectors and matrices.
!> @brief Provides generic assertions that stop execution on non-finite values.
module validation_mod
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use model_types_mod, only: dp
    implicit none
    private

    !> Generic validation entry point for rank-1 and rank-2 real arrays.
    interface ensure_finite
        module procedure ensure_finite_vector
        module procedure ensure_finite_matrix
    end interface ensure_finite

    public :: ensure_finite

contains

    !> Validates that a vector contains only finite floating-point values.
    subroutine ensure_finite_vector(label, values)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: values(:)

        if (any(.not. ieee_is_finite(values))) then
            write (*, '(a)') trim(label)//" contains non-finite values"
            error stop
        end if
    end subroutine ensure_finite_vector

    !> Validates that a matrix contains only finite floating-point values.
    subroutine ensure_finite_matrix(label, values)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: values(:, :)

        if (any(.not. ieee_is_finite(values))) then
            write (*, '(a)') trim(label)//" contains non-finite values"
            error stop
        end if
    end subroutine ensure_finite_matrix

end module validation_mod
