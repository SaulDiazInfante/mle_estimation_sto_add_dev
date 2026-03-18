program test_core
    use model_types_mod, only: spatial_grid_t
    use model_types_mod, only: get_state_dimension
    use parameter_ml_estimation_mod, only: build_uniform_checkpoints
    implicit none

    call test_state_dimension()
    call test_uniform_checkpoints_even_spacing()
    call test_uniform_checkpoints_clamped_request()

    write (*, '(a)') "Unit tests passed."

contains

    subroutine assert_true(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (*, '(2a)') "Assertion failed: ", trim(message)
            error stop
        end if
    end subroutine assert_true

    subroutine test_state_dimension()
        type(spatial_grid_t) :: grid

        grid%nx = 4
        grid%ny = 3

        call assert_true(&
            get_state_dimension(grid) == 11, &
            "state dimension should be nx * ny - 1" &
        )
    end subroutine test_state_dimension

    subroutine test_uniform_checkpoints_even_spacing()
        integer, allocatable :: checkpoints(:)

        call build_uniform_checkpoints(10, 3, 4, checkpoints)

        call assert_true(size(checkpoints) == 3, "expected three checkpoints")
        call assert_true(checkpoints(1) == 4, "first checkpoint should be 4")
        call assert_true(checkpoints(2) == 7, "second checkpoint should be 7")
        call assert_true(checkpoints(3) == 10, "third checkpoint should be 10")
    end subroutine test_uniform_checkpoints_even_spacing

    subroutine test_uniform_checkpoints_clamped_request()
        integer, allocatable :: checkpoints(:)

        call build_uniform_checkpoints(6, 10, 4, checkpoints)

        call assert_true(size(checkpoints) == 3, "request should clamp to 3")
        call assert_true(all(checkpoints == [4, 5, 6]), &
            "clamped checkpoints should cover all available points")
    end subroutine test_uniform_checkpoints_clamped_request

end program test_core
