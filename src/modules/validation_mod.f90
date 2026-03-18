module validation_mod
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use model_types_mod, only: dp
    implicit none
    private

    interface ensure_finite
        module procedure ensure_finite_vector
        module procedure ensure_finite_matrix
    end interface ensure_finite

    public :: ensure_finite

contains

    subroutine ensure_finite_vector(label, values)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: values(:)

        if (any(.not. ieee_is_finite(values))) then
            write (*, '(a)') trim(label)//" contains non-finite values"
            error stop
        end if
    end subroutine ensure_finite_vector

    subroutine ensure_finite_matrix(label, values)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: values(:, :)

        if (any(.not. ieee_is_finite(values))) then
            write (*, '(a)') trim(label)//" contains non-finite values"
            error stop
        end if
    end subroutine ensure_finite_matrix

end module validation_mod
