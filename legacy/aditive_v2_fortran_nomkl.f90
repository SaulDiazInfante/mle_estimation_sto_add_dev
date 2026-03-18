!> @file aditive_v2_fortran_nomkl.f90
!! @brief Legacy single-file executable without MKL or external BLAS dependencies.
!===========================================================================
! File: aditive_v2_fortran_nomkl.f90
!
! Build (no MKL / no external BLAS required):
!   gfortran -O3 -fopenmp aditive_v2_fortran_nomkl.f90 -o b.out
!===========================================================================
!> @brief Legacy standalone executable retained for reference and comparison.
program multivariate
    use ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
    implicit none

    integer, parameter :: seed = 42
    integer, parameter :: nx = 10, ny = 10, nobs = 1000
    integer, parameter :: dim = nx * ny - 1
    integer, parameter :: r_v = 1, s_v = 1
    real, parameter :: delta = 1.0e-5, theta = 1.0, beta = 0.1
    real, parameter :: gamma = 2.0, sigma = 0.1
    real, parameter :: lx = 5.0, ly = 5.0
    real, parameter :: pi = 3.141592653589793

    integer :: i, unit, nseed
    integer, allocatable :: seed_values(:)
    real, allocatable :: a(:, :), b(:, :), g(:, :)
    real, allocatable :: g_diag(:), b_diag(:), lambdas(:)
    real, allocatable :: startx(:), path(:, :)
    real, allocatable :: u0_prime_grid(:, :)
    real :: sigma_hat, beta_hat, theta_hat
    real :: initial_run_time, final_run_time, run_time_set_up, run_time_ml

    !call random_seed(size=nseed)
    !allocate (seed_values(nseed))
    !seed_values = seed
    !call random_seed(put=seed_values)
    !deallocate (seed_values)
    allocate (a(dim, dim), b(dim, dim), g(dim, dim))
    allocate (g_diag(dim), b_diag(dim), lambdas(dim))
    allocate (startx(dim), path(nobs, dim))
    allocate (u0_prime_grid(nx, ny))

    run_time_set_up = 0.0
    run_time_ml = 0.0

    call cpu_time(initial_run_time)

    call build_velocity_fields(nx, ny, lx, ly, pi, r_v, s_v, u0_prime_grid)
    call project_u0_prime(u0_prime_grid, nx, ny, lx, ly, startx)
    call build_matrices(dim, nx, ny, lx, ly, r_v, s_v, gamma, a, b, g)

    if (any(ieee_is_nan(startx))) stop "startx contains NaN"
    if (any(ieee_is_nan(a))) stop "A contains NaN"
    if (any(ieee_is_nan(b))) stop "B contains NaN"
    if (any(ieee_is_nan(g))) stop "G contains NaN"
    if (any(.not. ieee_is_finite(startx))) then
        stop "startx contains non-finite values"
    end if
    if (any(.not. ieee_is_finite(a))) stop "A contains non-finite values"
    if (any(.not. ieee_is_finite(b))) stop "B contains non-finite values"
    if (any(.not. ieee_is_finite(g))) stop "G contains non-finite values"

    do i = 1, dim
        g_diag(i) = g(i, i)
        b_diag(i) = b(i, i)
    end do
    lambdas = g_diag

    call simulate_sde_paths(&
        &dim, nobs, &
        &delta, beta, theta, sigma, &
        &g_diag, a, &
        &b_diag, startx, pi, path&
    &)

    unit = 10
    open (unit, file="vector_output.txt", status="replace", action="write")
    do i = 1, nobs
        write (unit, *) path(i, :)
    end do
    close (unit)

    if (any(ieee_is_nan(path))) stop "path contains NaN"
    if (any(.not. ieee_is_finite(path))) stop "path contains non-finite values"
    call cpu_time(final_run_time)
    run_time_set_up = final_run_time - initial_run_time
    
    call cpu_time(initial_run_time)
    call estimate_sigma_from_qv(&
        &nobs, &
        &dim, &
        &path, &
        &delta, &
        &lambdas, &
        &gamma, &
        &sigma_hat &
    &)
    call estimate_mle_u(&
        &nobs, &
        &dim, &
        &path, &
        &delta, &
        &a, &
        &lambdas, &
        &gamma, &
        &theta, &
        &beta_hat, &
        &theta_hat &
    )

    call cpu_time(final_run_time)

    run_time_ml = final_run_time - initial_run_time
    call print_estimation_report(&
        &sigma, &
        &beta, &
        &theta, &
        &sigma_hat, &
        &beta_hat, &
        &theta_hat, &
        &run_time_set_up, &
        &run_time_ml &
    &)

    deallocate (a, b, g, g_diag, b_diag, lambdas, startx, path)
    deallocate (u0_prime_grid)

contains

    subroutine print_estimation_report(&
        &sigma_true, &
        &beta_true, &
        &theta_true, &
        &sigma_est, &
        &beta_est, &
        &theta_est, &
        &setup_time, &
        &mle_time &
    &)
        implicit none
        real, intent(in) :: sigma_true, beta_true, theta_true
        real, intent(in) :: sigma_est, beta_est, theta_est
        real, intent(in) :: setup_time, mle_time
        real :: abs_err, rel_err
        character(len=13) :: rel_err_txt

        print '(a)', ''
        print '(a)', repeat('=', 86)
        print '(a)', 'MLE ESTIMATION REPORT'
        print '(a)', repeat('-', 86)
        write (*, '(a10,3x,a12,3x,a12,3x,a12,3x,a13)') &
            & 'Parameter', 'True Value', 'Estimate',&
            & 'Abs Error', 'Rel Error (%)'
        print '(a)', repeat('-', 86)
!
        abs_err = abs(sigma_est - sigma_true)
        if (abs(sigma_true) > tiny(1.0)) then
            rel_err = 100.0 * abs_err / abs(sigma_true)
        else
            rel_err = 0.0
        end if
        write (rel_err_txt, '(f13.4)') rel_err
        rel_err_txt = adjustl(rel_err_txt)
        write (*, '(a10,3x,f12.6,3x,f12.6,3x,es12.4,3x,a13)') 'sigma', &
        & sigma_true, sigma_est, abs_err, rel_err_txt
!
        abs_err = abs(beta_est - beta_true)
        if (abs(beta_true) > tiny(1.0)) then
            rel_err = 100.0 * abs_err / abs(beta_true)
        else
            rel_err = 0.0
        end if
        write (rel_err_txt, '(f13.4)') rel_err
        rel_err_txt = adjustl(rel_err_txt)
        write (*, '(a10,3x,f12.6,3x,f12.6,3x,es12.4,3x,a13)') 'beta', &
            & beta_true, beta_est, abs_err, rel_err_txt

        abs_err = abs(theta_est - theta_true)
        if (abs(theta_true) > tiny(1.0)) then
            rel_err = 100.0 * abs_err / abs(theta_true)
        else
            rel_err = 0.0
        end if
        write (rel_err_txt, '(f13.4)') rel_err
        rel_err_txt = adjustl(rel_err_txt)
        write (*, '(a10,3x,f12.6,3x,f12.6,3x,es12.4,3x,a13)') 'theta', &
            & theta_true, theta_est, abs_err, rel_err_txt
        print '(a)', repeat('-', 86)
        print '(a)', 'Runtime (seconds)'
        write (*, '(a20,2x,f12.6)') 'Set-up', setup_time
        write (*, '(a20,2x,f12.6)') 'MLE estimation', mle_time
        print '(a)', repeat('=', 86)
    end subroutine print_estimation_report

!-----------------------------------------------------------------------
! Subroutine: compute_initial_bump
!
! Description:
!   Computes a compactly supported cosine bump function centered at
!   (Lx / 4, 0.6 * Ly) with elliptical support. The function is commonly
!   used as a smooth initial condition in PDE simulations.
!
! Input:
!   x, y  - spatial coordinates
!   Lx, Ly - domain lengths in x and y directions
!
! Output:
!   u0xy  - value of the initial condition at (x,y)
!
!-----------------------------------------------------------------------
    subroutine compute_initial_bump(x, y, cb_Lx, cb_Ly, cb_PI, u0xy)
        implicit none
        real, intent(in)  :: x, y, cb_Lx, cb_Ly, cb_PI
        real, intent(out) :: u0xy
        real :: frac, ellipse_radius, x_normalized, y_normalized
        real :: x_center, y_center, scale_factor

        frac = 0.25
        x_center = cb_Lx / 4.0
        y_center = 0.6 * cb_Ly
        scale_factor = frac

        ! Compute normalized distances in elliptical coordinates
        x_normalized = (x - x_center) / (scale_factor * cb_Lx)
        y_normalized = (y - y_center) / (scale_factor * cb_Ly)
        ellipse_radius = sqrt(x_normalized**2 + y_normalized**2)

        ! Cosine bump function with compact support
        u0xy = 0.0
        if (ellipse_radius <= 1.0) then
            u0xy = 0.5 * (1.0 + cos(cb_PI * ellipse_radius))
        end if
    end subroutine compute_initial_bump

!-----------------------------------------------------------------------
! Subroutine: compute_velocity_x
!
! Description:
!   Evaluates the analytical function
!       compute_velocity_x(x,y) = -2*(pi*s/Ly)*sin(pi*r*x/Lx)*cos(pi*s*y/Ly)
!   typically used in PDE test problems or manufactured solutions.
!
! Input:
!   x, y  - spatial coordinates
!   r, s  - mode parameters (typically integers but treated as real)
!   Lx, Ly - domain lengths in x and y directions
!
! Output:
!   v1xy  - function value at (x,y)
!
!-----------------------------------------------------------------------
    subroutine compute_velocity_x(x, y, r, s, vx_Lx, vx_Ly, v1xy)
        implicit none
        real, intent(in)  :: x, y, vx_Lx, vx_Ly
        integer, intent(in)  :: r, s
        real, intent(out) :: v1xy
        real, parameter :: vx_pi = 3.1415927
        real :: wave_x, wave_y, amplitude

        ! Compute wave numbers and amplitude
        wave_x = vx_pi * r * x / vx_Lx
        wave_y = vx_pi * s * y / vx_Ly
        amplitude = -2.0 * vx_pi * s / vx_Ly
        v1xy = amplitude * sin(wave_x) * cos(wave_y)

    end subroutine compute_velocity_x

!-----------------------------------------------------------------------
! Subroutine: compute_velocity_y
!
! Description:
!   Evaluates the analytical function
!       compute_velocity_y(x,y) = 2*(pi*r/Lx)*cos(pi*r*x/Lx)*sin(pi*s*y/Ly)
!   commonly used in manufactured solutions for PDE verification.
!
! Input:
!   x, y  - spatial coordinates
!   r, s  - mode parameters (typically integers but treated as real)
!   Lx, Ly - domain lengths in x and y directions
!
! Output:
!   v2xy  - function value at (x,y)
!
!-----------------------------------------------------------------------
    subroutine compute_velocity_y(x, y, r, s, vy_Lx, vy_Ly, v2xy)
        implicit none
        real, intent(in)  :: x, y, vy_Lx, vy_Ly
        integer, intent(in)  :: r, s
        real, intent(out) :: v2xy
        real, parameter :: vy_pi = 3.1415927
        real :: wave_x, wave_y, amplitude

        ! Compute wave numbers and amplitude
        wave_x = vy_pi * r * x / vy_Lx
        wave_y = vy_pi * s * y / vy_Ly
        amplitude = 2.0 * vy_pi * r / vy_Lx
        v2xy = amplitude * cos(wave_x) * sin(wave_y)
    end subroutine compute_velocity_y

!-----------------------------------------------------------------------
! Subroutine: compute_mode_integral
!
! Description:
!   Computes the weighted discrete integral of the field u0grid over a
!   uniform cell-centered mesh. The weights are cosine basis functions
!   determined by the mode indices (i,j). Typically used for spectral
!   or modal projections.
!
! Input:
!   u0grid(Nx,Ny) - field to be integrated
!   i, j          - mode indices
!   Nx, Ny        - number of grid cells
!   Lx, Ly        - domain lengths
!
! Output:
!   intval        - value of the weighted integral
!
!-----------------------------------------------------------------------
    subroutine compute_mode_integral(u0grid, mi, mj, mi_Nx, mi_Ny, mi_Lx, mi_Ly, intval)

        implicit none

        integer, intent(in) :: mi, mj, mi_Nx, mi_Ny
        real, intent(in) :: mi_Lx, mi_Ly
        real, intent(in) :: u0grid(mi_Nx, mi_Ny)
        real, intent(out) :: intval

        integer :: ix, iy
        real :: dx, dy, xi, yj
        real :: basis_i, basis_j, integrand, cell_area
        real :: wave_number_x, wave_number_y
        real :: normalization_i, normalization_j
        real, parameter :: mi_pi = 3.1415927

        ! Grid spacing
        dx = mi_Lx / real(mi_Nx)
        dy = mi_Ly / real(mi_Ny)
        cell_area = dx * dy

        ! Wave numbers
        wave_number_x = mi_pi * real(mi) / mi_Lx
        wave_number_y = mi_pi * real(mj) / mi_Ly

        ! Normalization factors for basis functions
        normalization_i = sqrt(1.0 + sign(1.0, real(mi)))
        normalization_j = sqrt(1.0 + sign(1.0, real(mj)))

        intval = 0.0

        do iy = 1, mi_Ny
            ! Cell center in y-direction
            yj = dy * (real(iy) - 0.5)
            basis_j = normalization_j * cos(wave_number_y * yj)
            do ix = 1, mi_Nx
                ! Cell center in x-direction
                xi = dx * (real(ix) - 0.5)
                basis_i = normalization_i * cos(wave_number_x * xi)

                ! Weighted integrand
                integrand = u0grid(ix, iy) * basis_i * basis_j
                intval = intval + cell_area * integrand
            end do
        end do
    end subroutine compute_mode_integral

!-----------------------------------------------------------------------
! Subroutine: build_velocity_fields
!
! Description:
!   Builds the spatial grids and evaluates the velocity components
!   (compute_velocity_x, compute_velocity_y) and the initial condition 
!   compute_initial_bump on a uniform cell-centered
!   mesh. It also computes the spatial mean of compute_initial_bump and the 
!        fluctuation
!   field compute_initial_bump'.
!
! Input:
!   Nx, Ny   - number of grid cells in x and y
!   Lx, Ly   - domain lengths
!   r_v, s_v - mode parameters for velocity field
!
! Output:
!   xgrid, ygrid        - cell-centered coordinates
!   v1grid, v2grid      - velocity fields
!   u0_grid             - initial condition
!   u0_prime_grid       - fluctuation field
!
! Note:
!   Requires external routines: compute_initial_bump, compute_velocity_x, compute_velocity_y, compute_mode_integral
!
!-----------------------------------------------------------------------
    subroutine build_velocity_fields(&
        & bvf_Nx, bvf_Ny, &
        & bvf_Lx, bvf_Ly, &
        & bvf_PI, bvf_rv, bvf_sv, &
        & bvf_u0_prime &
    &)

        implicit none

        integer, intent(in) :: bvf_Nx, bvf_Ny, bvf_rv, bvf_sv
        real, intent(in) :: bvf_Lx, bvf_Ly, bvf_PI

        real :: xgrid(bvf_Nx, bvf_Ny), ygrid(bvf_Nx, bvf_Ny)
        real :: v1grid(bvf_Nx, bvf_Ny), v2grid(bvf_Nx, bvf_Ny), u0_grid(bvf_Nx, bvf_Ny)
        real, intent(out) :: bvf_u0_prime(bvf_Nx, bvf_Ny)

        integer :: ix, iy
        real :: xi, yj, u0_bar, intval

        !----- Build grids and evaluate fields
        do iy = 1, bvf_Ny
            do ix = 1, bvf_Nx
                xgrid(ix, iy) = (bvf_Lx / real(bvf_Nx)) * (real(ix) - 0.5)
                ygrid(ix, iy) = (bvf_Ly / real(bvf_Ny)) * (real(iy) - 0.5)
                call compute_velocity_x(&
                    & xgrid(ix, iy), &
                    & ygrid(ix, iy), &
                    & bvf_rv, bvf_sv, bvf_Lx, bvf_Ly, &
                    & v1grid(ix, iy) &
                )

                call compute_velocity_y(&
                    & xgrid(ix, iy), &
                    & ygrid(ix, iy), &
                    & bvf_rv, bvf_sv, bvf_Lx, bvf_Ly, &
                    & v2grid(ix, iy) &
                )

                xi = xgrid(ix, iy)
                yj = ygrid(ix, iy)

                call compute_initial_bump(&
                    & xi, &
                    & yj, &
                    & bvf_Lx, &
                    & bvf_Ly, &
                    & bvf_PI, &
                    & u0_grid(ix, iy) &
                )
            end do
        end do

        !----- Compute spatial mean
        call compute_mode_integral(u0_grid, 0, 0, bvf_Nx, bvf_Ny, bvf_Lx, bvf_Ly, intval)
        u0_bar = (1.0 / (bvf_Lx * bvf_Ly)) * intval
        !----- Fluctuation field
        bvf_u0_prime = u0_grid - u0_bar
    end subroutine build_velocity_fields

!-----------------------------------------------------------------------
! Subroutine: build_mode_order
!
! Description:
!   Converts a linear index ordering into the corresponding (i,j)
!   grid indices. The mapping assumes lexicographic ordering and
!   accounts for the removal of the first eigenvalue.
!
! Input:
!   I(Nx*Ny-1)  - linear index vector (zero-based expected from MATLAB logic)
!   Nx, Ny      - grid dimensions
!
! Output:
!   Orden_ij(Nx*Ny-1,2) - array containing (i,j) index pairs
!
!-----------------------------------------------------------------------
    subroutine build_mode_order(bmo_Nx, bmo_Ny, bmo_Lx, bmo_Ly, Orden_ij)
        implicit none
        integer, intent(in) :: bmo_Nx, bmo_Ny
        real, intent(in):: bmo_Lx, bmo_Ly
        integer :: Iv(bmo_Nx * bmo_Ny - 1)
        integer, intent(out) :: Orden_ij(bmo_Nx * bmo_Ny - 1, 2)
        integer :: N2m1
        integer :: O, m, bmo_i, bmo_j
        N2m1 = bmo_Nx * bmo_Ny - 1

        call build_sorted_eigenvalues(bmo_Nx, bmo_Ny, bmo_Lx, bmo_Ly, N2m1, Iv)

        do O = 1, N2m1
            m = Iv(O)
            m = m + 1   ! first eigenvalue removed

            bmo_j = int((m - 1) / bmo_Nx)
            bmo_i = (m - 1) - bmo_j * bmo_Nx

            Orden_ij(O, 1) = bmo_i
            Orden_ij(O, 2) = bmo_j

        end do
    end subroutine build_mode_order

!-----------------------------------------------------------------------
! Subroutine: build_sorted_eigenvalues
!
! Description:
!   Computes the spectral eigenvalues lambda_ij for all modes with
!   i+j > 0 and sorts them in ascending order. Returns both the sorted
!   values and the permutation index vector.
!
! Input:
!   Nx, Ny   - grid dimensions
!   Lx, Ly   - domain lengths
!
! Output:
!   Eval_sorted(N2m1) - sorted eigenvalues
!
!-----------------------------------------------------------------------
    subroutine build_sorted_eigenvalues(&
        &bse_Nx, &
        &bse_Ny, &
        &bse_Lx, &
        bse_Ly, &
        &N2m1, &
        &mode_order&
    &)

        implicit none

        integer, intent(in) :: bse_Nx, bse_Ny, N2m1
        real, intent(in) :: bse_Lx, bse_Ly
        integer, intent(out) :: mode_order(N2m1)

        integer :: bse_i, bse_j, m
        integer :: mode_idx(N2m1)
        real :: Eval(N2m1)
        real, parameter :: bse_pi = 3.1415927

        !-------------------------------------------------
        m = 0
        do bse_i = 0, bse_Nx - 1
            do bse_j = 0, bse_Ny - 1

                if (bse_i + bse_j > 0) then
                    m = m + 1
                    Eval(m) = (bse_pi / bse_Lx)**2 * real(bse_i)**2 + &
                              (bse_pi / bse_Ly)**2 * real(bse_j)**2
                    mode_idx(m) = bse_i + bse_j * bse_Nx
                end if

            end do
        end do
        !-------------------------------------------------
        call sort_with_index(Eval, mode_idx, N2m1, mode_order)
    end subroutine build_sorted_eigenvalues

!-----------------------------------------------------------------------
! Subroutine: sort_with_index
!
! Description:
!   Sorts array in ascending order and returns permutation indices.
!   Simple O(n^2) version (sufficient for moderate N).
!
!-----------------------------------------------------------------------
    subroutine sort_with_index(swi_A, idx_in, n, idx_sorted)

        implicit none
        integer, intent(in) :: n
        real, intent(in) :: swi_A(n)
        integer, intent(in) :: idx_in(n)
        integer, intent(out) :: idx_sorted(n)

        integer :: swi_i, swi_j, min_idx
        integer :: temp_idx
        real :: temp
        real :: A_sorted(n)

        A_sorted = swi_A
        idx_sorted = idx_in

        do swi_i = 1, n - 1
            min_idx = swi_i
            do swi_j = swi_i + 1, n
                if (A_sorted(swi_j) < A_sorted(min_idx)) then
                    min_idx = swi_j
                end if
            end do

            if (min_idx /= swi_i) then
                temp = A_sorted(swi_i)
                A_sorted(swi_i) = A_sorted(min_idx)
                A_sorted(min_idx) = temp
                temp_idx = idx_sorted(swi_i)
                idx_sorted(swi_i) = idx_sorted(min_idx)
                idx_sorted(min_idx) = temp_idx
            end if
        end do
    end subroutine sort_with_index

!-----------------------------------------------------------------------
! Subroutine: project_u0_prime
!
! Description:
!   Computes the modal projection coefficients of the fluctuation field
!   u0_prime_grid onto the cosine basis indexed by Order_ij. Each
!   coefficient is obtained via the weighted integral normalized by
!   the domain area.
!
! Input:
!   u0_prime_grid(Nx,Ny) - fluctuation field
!   Order_ij(N2m1,2)     - mode index pairs (i,j)
!   Nx, Ny               - grid dimensions
!   Lx, Ly               - domain lengths
!
! Output:
!   u0_prime_proj_row(N2m1) - projection coefficients
!
! Requires:
!   subroutine compute_mode_integral
!
!-----------------------------------------------------------------------
    subroutine project_u0_prime(&
        & u0_prime_grid_in, &
        & pu0_Nx, pu0_Ny, &
        & pu0_Lx, pu0_Ly, &
        & u0_prime_proj_row&
    &)

        implicit none

        integer, intent(in) :: pu0_Nx, pu0_Ny
        integer :: Order_ij(pu0_Nx * pu0_Ny - 1, 2)
        real, intent(in) :: pu0_Lx, pu0_Ly
        real, intent(in) :: u0_prime_grid_in(pu0_Nx, pu0_Ny)
        real, intent(out) :: u0_prime_proj_row(pu0_Nx * pu0_Ny - 1)

        integer :: N2m1
        integer :: O, pu0_i, pu0_j
        real :: intval

        N2m1 = pu0_Nx * pu0_Ny - 1

        call build_mode_order(pu0_Nx, pu0_Ny, pu0_Lx, pu0_Ly, Order_ij)

        do O = 1, N2m1
            pu0_i = Order_ij(O, 1)
            pu0_j = Order_ij(O, 2)

            call compute_mode_integral(&
                &u0_prime_grid_in, &
                & pu0_i, pu0_j, &
                & pu0_Nx, pu0_Ny, &
                & pu0_Lx, pu0_Ly, &
                & intval &
            &)

            u0_prime_proj_row(O) = (1.0 / (pu0_Lx * pu0_Ly)) * intval
        end do

    end subroutine project_u0_prime

!-----------------------------------------------------------------------
! Subroutine: compute_interaction_coefficient
!
! Description:
!   Computes the analytical coefficient v_ijkl based on modal indices
!   (i,j,k,l) and parameters (r,s). The formula follows the piecewise
!   analytical expressions derived from the cosine basis interactions.
!
! Input:
!   i, j, k, l - mode indices
!   r, s       - velocity mode parameters
!
! Output:
!   v_ijkl     - coefficient value
!
!-----------------------------------------------------------------------
    subroutine compute_interaction_coefficient(&
        & mi, &
        & mj, &
        & mk, &
        & ml, &
        & mr, &
        & ms, &
        & Lx, &
        & Ly, &
        & v_ijkl &
    &)

        implicit none
        integer, intent(in) :: mi, mj, mk, ml, mr, ms
        real, intent(in) :: Lx, Ly
        real, intent(out) :: v_ijkl
        real, parameter :: pi = 3.1415927
        real :: base_coefficient, interaction_term
        real :: sqrt2_factor, coefficient_i0, base_coeff_k0
        real :: coefficient_j0, base_coeff_l0
        logical :: i_minus_r_matches_k, i_plus_r_matches_k
        logical :: j_minus_s_matches_l, j_plus_s_matches_l
        base_coefficient = pi**2 / (2.0 * Lx * Ly)
        v_ijkl = 0.0

        !-------------------------------------------------
        ! Case 1: All indices nonzero
        if (mi * mj * mk * ml > 0) then
            ! Precompute logical conditions
            i_minus_r_matches_k = ((mi - mr)**2 == mk**2)
            i_plus_r_matches_k = ((mi + mr)**2 == mk**2)
            j_minus_s_matches_l = ((mj - ms)**2 == ml**2)
            j_plus_s_matches_l = ((mj + ms)**2 == ml**2)

            interaction_term = mj * mr - mi * ms

            if (i_minus_r_matches_k .and. j_minus_s_matches_l) then
                v_ijkl = -base_coefficient * interaction_term
            end if

            if (i_plus_r_matches_k .and. j_minus_s_matches_l) then
                v_ijkl = -base_coefficient * (mj * mr + mi * ms)
            end if

            if (i_plus_r_matches_k .and. j_plus_s_matches_l) then
                v_ijkl = base_coefficient * interaction_term
            end if

            if (i_minus_r_matches_k .and. j_plus_s_matches_l) then
                v_ijkl = base_coefficient * (mj * mr + mi * ms)
            end if

        end if
        !-------------------------------------------------
        ! Case 2: mi = 0, others nonzero
        if (mi == 0 .and. mj * mk * ml > 0) then
            sqrt2_factor = sqrt(2.0) / 2.0
            coefficient_i0 = sqrt2_factor * pi**2 * mj * mr / (Lx * Ly)

            if (mk == mr .and. (mj - ms)**2 == ml**2) then
                v_ijkl = -coefficient_i0
            end if

            if (mk == mr .and. (mj + ms)**2 == ml**2) then
                v_ijkl = coefficient_i0
            end if
        end if
        !-------------------------------------------------
        ! Case 3: mk = 0, others nonzero
        if (mk == 0 .and. mi * mj * ml > 0) then
            sqrt2_factor = sqrt(2.0) / 2.0
            base_coeff_k0 = sqrt2_factor * pi**2 / (Lx * Ly)

            if (mi == mr .and. (mj - ms)**2 == ml**2) then
                v_ijkl = -base_coeff_k0 * (mj * mr - mi * ms)
            end if

            if (mi == mr .and. (mj + ms)**2 == ml**2) then
                v_ijkl = base_coeff_k0 * (mj * mr + mi * ms)
            end if

        end if
        !-------------------------------------------------
        ! Case 4: mj = 0, others nonzero
        if (mj == 0 .and. mi * mk * ml > 0) then
            sqrt2_factor = sqrt(2.0) / 2.0
            coefficient_j0 = sqrt2_factor * pi**2 * mi * ms / (Lx * Ly)

            if (ml == ms .and. (mi - mr)**2 == mk**2) then
                v_ijkl = coefficient_j0
            end if

            if (ml == ms .and. (mi + mr)**2 == mk**2) then
                v_ijkl = -coefficient_j0
            end if
        end if
        !-------------------------------------------------
        ! Case 5: ml = 0, others nonzero
        if (ml == 0 .and. mi * mj * mk > 0) then
            sqrt2_factor = sqrt(2.0) / 2.0
            base_coeff_l0 = sqrt2_factor * pi**2 / (Lx * Ly)

            if (mj == ms .and. (mi - mr)**2 == mk**2) then
                v_ijkl = -base_coeff_l0 * (mj * mr - mi * ms)
            end if

            if (mj == ms .and. (mi + mr)**2 == mk**2) then
                v_ijkl = -base_coeff_l0 * (mj * mr + mi * ms)
            end if
        end if
        !-------------------------------------------------
        ! Case 6: mi = 0 and ml = 0
        if (mi == 0 .and. ml == 0 .and. mj * mk > 0) then
            if (mk == mr .and. mj == ms) then
                v_ijkl = -pi**2 * mj * mr / (Lx * Ly)
            end if
        end if

        !-------------------------------------------------
        ! Case 7: mj = 0 and mk = 0
        if (mj == 0 .and. mk == 0 .and. mi * ml > 0) then
            if (mi == mr .and. ml == ms) then
                v_ijkl = pi**2 * mi * ms / (Lx * Ly)
            end if
        end if
    end subroutine compute_interaction_coefficient

!-----------------------------------------------------------------------
! Subroutine: build_matrices
!
! Description:
!   Builds matrices A, B and MLambda together with their flattened
!   row-wise storage vectors. Matrix A is filled using the analytical
!   coefficient routine, while B and MLambda are diagonal matrices
!   based on the spectral eigenvalues.
!
! Input:
!   N2m1                - number of active modes
!   Nx, Ny              - grid sizes
!   Lx, Ly              - domain lengths
!   r_v, s_v            - velocity parameters
!   gamma               - exponent for matrix B
!
! Output:
!   A, B, MLambda       - matrices (N2m1 x N2m1)
!
! Requires:
!   subroutine compute_interaction_coefficient
!
!-----------------------------------------------------------------------
    subroutine build_matrices(&
        & N2m1, bm_Nx, bm_Ny,&
        & bm_Lx, bm_Ly, bm_rv, bm_sv,&
        & bm_gamma, &
        & bm_A, bm_B, MLambda&
    &)

        implicit none

        integer, intent(in) :: N2m1, bm_Nx, bm_Ny
        integer, intent(in) :: bm_rv, bm_sv
        real, intent(in) :: bm_Lx, bm_Ly, bm_gamma
        integer :: Order_ij(N2m1, 2)
        real, intent(out) :: bm_A(N2m1, N2m1), bm_B(N2m1, N2m1)
        real, intent(out) :: MLambda(N2m1, N2m1)

        integer :: O1, O2, bi, bj, bk, bl
        integer :: r
        real :: lambda_ij
        real :: coeff
        real, parameter :: pi = 3.1415927

        ! Initialize
        bm_A = 0.0
        bm_B = 0.0
        MLambda = 0.0

        r = 0
        call build_mode_order(bm_Nx, bm_Ny, bm_Lx, bm_Ly, Order_ij)

        do O1 = 1, N2m1

            bi = Order_ij(O1, 1)
            bj = Order_ij(O1, 2)

            do O2 = 1, N2m1

                bk = Order_ij(O2, 1)
                bl = Order_ij(O2, 2)

                ! ---- Matrix A

                call compute_interaction_coefficient(&
                    & bi, bj, bk, bl, &
                    & bm_rv, bm_sv, &
                    & bm_Lx, bm_Ly, coeff&
                &)

                bm_A(O2, O1) = coeff
                r = r + 1
                ! ---- Matrices B and Lambda (diagonal)
                if (bi == bk .and. bj == bl) then
                    lambda_ij = (pi / bm_Lx)**2 * real(bi)**2 + &
                                (pi / bm_Ly)**2 * real(bj)**2
                    bm_B(O2, O1) = lambda_ij**(-bm_gamma)

                    MLambda(O2, O1) = lambda_ij
                end if
            end do
        end do

    end subroutine build_matrices

!! ============================================================
! SUBROUTINE: reshape_to_matrix
!
! Description:
!   Reshapes a vector of length DIM*DIM into a square
!   DIM x DIM matrix using column-major order (Fortran default).
!
! Arguments:
!   DIM   (integer, input)   : Matrix dimension
!   AM    (real*8, input)    : Vectorized matrix (length DIM*DIM)
!   A     (real*8, output)   : Reshaped DIM x DIM matrix
!
! ============================================================
    subroutine reshape_to_matrix(rtm_DIM, AM, rtm_A)
        implicit none
        integer, intent(in) :: rtm_DIM
        real, intent(in)  :: AM(rtm_DIM * rtm_DIM)
        real, intent(out) :: rtm_A(rtm_DIM, rtm_DIM)
        rtm_A = reshape(AM, [rtm_DIM, rtm_DIM])
    end subroutine reshape_to_matrix

! ============================================================
! SUBROUTINE: reshape_and_extract_diagonal
!
! Description:
!   Reshapes a vector into a DIM x DIM matrix and extracts
!   its diagonal elements into a separate vector.
!
! Arguments:
!   DIM      (integer, input)   : Matrix dimension
!   LM       (real*8, input)    : Vectorized matrix (length DIM*DIM)
!   LG       (real*8, output)   : Reshaped DIM x DIM matrix
!   lambdas  (real*8, output)   : Diagonal elements of LG
!
! ============================================================
    subroutine reshape_and_extract_diagonal(red_DIM, LM, LG, red_lambdas)
        implicit none

        integer, intent(in) :: red_DIM
        real, intent(in)  :: LM(red_DIM * red_DIM)
        real, intent(out) :: LG(red_DIM, red_DIM)
        real, intent(out) :: red_lambdas(red_DIM)

        integer :: red_i

        ! Reshape vector into matrix (Fortran column-major order)
        LG = reshape(LM, [red_DIM, red_DIM])
        do red_i = 1, red_DIM
            red_lambdas(red_i) = LG(red_i, red_i)
        end do

    end subroutine reshape_and_extract_diagonal

! ============================================================
! SUBROUTINE: reshape_matrix
!
! Description:
!   Reshapes a vectorized matrix into a DIM x DIM matrix.
!
! Arguments:
!   DIM   (integer, input)   : Matrix dimension
!   BM    (real*8, input)    : Vectorized matrix (length DIM*DIM)
!   B     (real*8, output)   : Reshaped DIM x DIM matrix
!
! ============================================================
    subroutine reshape_matrix(rm_DIM, BM, rm_B)
        implicit none
        integer, intent(in) :: rm_DIM
        real, intent(in)  :: BM(rm_DIM * rm_DIM)
        real, intent(out) :: rm_B(rm_DIM, rm_DIM)
        ! Reshape vector into matrix (column-major order)
        rm_B = reshape(BM, [rm_DIM, rm_DIM])
    end subroutine reshape_matrix
! ============================================================
! SUBROUTINE: simulate_sde_paths
!
! Description:
!   Simulates a multidimensional stochastic differential equation
!   using the Euler–Milstein scheme. The drift is linear with
!   diagonal and interaction terms, and the diffusion is diagonal.
!
! ============================================================
    subroutine simulate_sde_paths(&
        & sde_DIM, N, &
        & sde_delta, sde_beta, sde_theta, sde_sigma, &
        & sde_g_diag, sde_A, &
        & sde_b_diag, x0, sde_PI, paths &
    &)
        use ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
        implicit none

        integer, intent(in) :: sde_DIM, N
        real, intent(in)    :: sde_delta, sde_beta, sde_theta, sde_sigma, sde_PI
        real, intent(in)    :: sde_g_diag(sde_DIM), &
            &sde_A(sde_DIM, sde_DIM), sde_b_diag(sde_DIM)
        real, intent(in)    :: x0(sde_DIM)
        real, intent(out)   :: paths(N, sde_DIM)

        real :: x(sde_DIM)
        real :: alpha(sde_DIM), diff_diag(sde_DIM), W(sde_DIM)
        real :: std_dev
        integer :: k

        ! Initial condition
        x = x0
        paths(1, :) = x
        diff_diag = sde_sigma * sde_b_diag
        std_dev = sqrt(sde_delta)

        ! Time stepping
        do k = 2, N
            alpha = -sde_theta * matmul(sde_A, x) - sde_beta * sde_g_diag * x
            call generate_normal_array(sde_DIM, W)
            x = x + alpha * sde_delta + diff_diag * std_dev * W
            paths(k, :) = x
        end do
    end subroutine simulate_sde_paths

! ============================================================
! SUBROUTINE: euler_maruyama_step
!
! Description:
!   Performs one time step of a vector-valued SDE using
!   an Euler–Maruyama update (no derivative term included).
!
! Arguments:
!   DIM     (integer, input)  : Dimension of the state vector
!   delta   (real*8, input)   : Time step size
!   x       (real*8, input)   : Current state
!   alpha   (real*8, input)   : Drift evaluated at x
!   sigma   (real*8, input)   : Diffusion coefficient
!   W       (real*8, input)   : Brownian increment
!   xnew    (real*8, output)  : Updated state
!
! ============================================================
    subroutine euler_maruyama_step(DIM, delta, x, alpha, sigma, W, xnew)
        implicit none

        integer, intent(in) :: DIM
        real, intent(in)  :: delta
        real, intent(in)  :: x(DIM), alpha(DIM), sigma(DIM), W(DIM)
        real, intent(out) :: xnew(DIM)
        integer :: i

        !$omp parallel do default(shared) private(i)
        do i = 1, DIM
            xnew(i) = x(i) + alpha(i) * delta + sigma(i) * W(i)
        end do
        !$omp end parallel do
    end subroutine euler_maruyama_step

! ============================================================
! SUBROUTINE: generate_normal_array
!
! Description:
!   Generates an array of independent identically distributed
!   random numbers following a standard normal distribution 
!   N(0, 1) using the native Box-Muller transform.
!
!   By default, it produces mean 0 and standard deviation 1.
! ============================================================
    subroutine generate_normal_array(na_DIM, W)
        implicit none

        integer, intent(in) :: na_DIM
        real, intent(out)   :: W(na_DIM)

        real :: U1(na_DIM), U2(na_DIM)
        real, parameter :: epsilon = 1.0e-12
        real, parameter :: na_PI = 3.141592653589793

        ! Fortran's native uniform [0, 1) generator
        call random_number(U1)
        call random_number(U2)

        ! Avoid log(0) by clamping minimum value
        U1 = max(U1, epsilon)

        ! Native Array SIMD Box-Muller transform for N(0, 1)
        W = sqrt(-2.0 * log(U1)) * cos(2.0 * na_PI * U2)

    end subroutine generate_normal_array

! ============================================================
! SUBROUTINE: compute_linear_drift
!
! Description:
!   Computes the drift vector of a linear multidimensional SDE
!   using a diagonal mean-reversion term and a linear interaction
!   term expressed via matrix multiplication.
!
! ============================================================
    subroutine compute_linear_drift(&
        &ld_DIM, &
        &ld_beta, &
        &ld_theta, &
        &ld_g_diag, &
        &ld_A, &
        &x, &
        &alpha&
    &)
        implicit none
        integer, intent(in) :: ld_DIM
        real, intent(in)    :: ld_beta, ld_theta
        real, intent(in)    :: ld_g_diag(ld_DIM)
        real, intent(in)    :: ld_A(ld_DIM, ld_DIM)
        real, intent(in)    :: x(ld_DIM)
        real, intent(out)   :: alpha(ld_DIM)

        ! Drift: alpha = -beta*diag(g_diag)*x - theta*A*x
        alpha = -ld_theta * matmul(ld_A, x) - ld_beta * ld_g_diag * x

    end subroutine compute_linear_drift

! ============================================================
! SUBROUTINE: compute_diffusion_from_diagonal
!
! Description:
!   Computes a vector given by sigma times the diagonal of
!   matrix B. Used as a diffusion-related drift component.
!
! ============================================================
    subroutine compute_diffusion_from_diagonal(dfd_DIM, dfd_sigma, dfd_B, var)
        implicit none

        integer, intent(in) :: dfd_DIM
        real, intent(in)    :: dfd_sigma
        real, intent(in)    :: dfd_B(dfd_DIM, dfd_DIM)
        real, intent(out)   :: var(dfd_DIM)

        integer :: dfd_i

        do dfd_i = 1, dfd_DIM
            var(dfd_i) = dfd_sigma * dfd_B(dfd_i, dfd_i)
        end do

    end subroutine compute_diffusion_from_diagonal

! ============================================================
! SUBROUTINE: estimate_sigma_from_qv
!
! Description:
!   Computes the maximum likelihood estimator (MLE) of the
!   diffusion coefficient for each component of a multidimensional
!   process, and returns their average.
!
!   The estimator is based on quadratic variation of the
!   observed process Ys.
!
! Arguments:
!   npoints    (integer, input)  : Number of time increments
!   DIM        (integer, input)  : Dimension of the process
!   Ys         (real, input)     : Observations Ys(t_k) of the process
!                                  indexed as Ys(0:npoints, DIM)
!   delta      (real, input)     : Time step size
!   lambdas    (real, input)     : Scaling parameters per dimension
!   gamma      (real, input)     : Scaling exponent
!   sigma_hat  (real, output)    : Estimated diffusion coefficient
!
! ============================================================
    subroutine estimate_sigma_from_qv(&
        & npoints, &
        & esq_DIM, &
        & Ys, &
        & esq_delta, &
        & esq_lambdas, &
        & esq_gamma, &
        & esq_sigma_hat&
        &)
        use ieee_arithmetic, only: ieee_is_nan, IEEE_IS_FINITE
        implicit none
        integer, intent(in) :: npoints, esq_DIM
        real, intent(in)    :: esq_delta, esq_gamma
        real, intent(in)    :: Ys(npoints, esq_DIM), esq_lambdas(esq_DIM)
        real, intent(out)   :: esq_sigma_hat

        integer :: k
        real    :: del
        real    :: sigmas(esq_DIM)
        real    :: increments_squared, scaling_factor, quadratic_variation

        if (any(.not. ieee_is_finite(Ys))) then
            print *, "Ys with infinite"
            stop
        end if


        del = esq_delta * npoints

        ! Compute component-wise MLEs using quadratic variation
        do k = 1, esq_DIM
            ! Sum of squared increments
            increments_squared = sum( &
                (&
                    & Ys(2:npoints, k) &
                    & - Ys(1:(npoints - 1), k) &
                &)**2)

            scaling_factor = esq_lambdas(k)**(2.0 * esq_gamma)

            quadratic_variation = increments_squared * scaling_factor / del
            sigmas(k) = sqrt(quadratic_variation)
        end do

        if (any(ieee_is_nan(sigmas))) then
            print *, "Vector contains NaN"
        end if

        if (any(.not. ieee_is_finite(sigmas))) then
            print *, "sigmas with infinite"
            stop
        end if
        ! Average estimator across dimensions
        esq_sigma_hat = sum(sigmas) / esq_DIM
        return
    end subroutine estimate_sigma_from_qv

! ============================================================
! SUBROUTINE: integrate_ito
!
! Description:
!   Computes the Itô stochastic integral of a discrete-time
!   adapted process with respect to a Brownian motion using
!   a left-point (Itô) discretization.
!
! Arguments:
!   npoints   (integer, input)  : Number of discretization points
!   path      (real, input)     : Values of the integrand process
!                                 indexed as path(1:npoints)
!   W         (real, input)     : Brownian motion values
!                                 indexed as W(1:npoints)
!   integral  (real, output)    : Numerical approximation of the
!                                 Itô stochastic integral
!
! ============================================================
    subroutine integrate_ito(npoints, ito_path, W, integral)
        implicit none

        integer, intent(in) :: npoints
        real, intent(in) :: ito_path(npoints), W(npoints)
        real, intent(out):: integral

        integral = dot_product(&
            & ito_path(1:npoints - 1),&
            & W(2:npoints) - W(1:npoints - 1)&
        &)
        return
    end subroutine integrate_ito

! ============================================================
! SUBROUTINE: integrate_trapezoidal
!
! Description:
!   Computes the time integral of a discretely observed path
!   over the interval [0, npoints*delta] using the trapezoidal
!   quadrature rule.
!
! Arguments:
!   npoints   (integer, input)  : Number of discretization points
!   path      (real, input)     : Observations of the process
!                                 indexed as path(1:npoints)
!   delta     (real, input)     : Time step size
!   integral  (real, output)    : Numerical approximation of the
!                                 time integral of the path
!
! ============================================================
    subroutine integrate_trapezoidal(&
        & npoints,&
        & itr_path,&
        & itr_dt,&
        & integral&
    )
        implicit none

        integer, intent(in) :: npoints
        real, intent(in) :: itr_dt
        real, intent(in) :: itr_path(npoints)
        real, intent(out):: integral
        integral = 0.0
        integral = 0.5*itr_dt * &
            &(&
                & sum(itr_path(1:npoints - 1)) &
                & + sum(itr_path(2:npoints)) &
            &)
        return
    end subroutine integrate_trapezoidal

! ============================================================
! SUBROUTINE: compute_sufficient_statistics_aux
!
! Description:
!   Computes the reduced set of sufficient statistics used by the
!   auxiliary beta-only estimator (when theta is treated as known).
!
! Arguments:
!   npoints   (integer, input)  : Number of time increments
!   DIM       (integer, input)  : Dimension of the process
!   Us        (real, input)     : Observations of the process U,
!                                 indexed as Us(0:npoints, DIM)
!   delta     (real, input)     : Time step size
!   Vs        (real, input)     : Transformation matrix (DIM x DIM)
!   lambdas   (real, input)     : Scaling parameters per dimension
!   gamma     (real, input)     : Scaling exponent
!   SS        (real, output)    : Vector of sufficient statistics
!                                 SS(1:3)
!
! ============================================================
    subroutine compute_sufficient_statistics_u(&
        &npoints,&
        & csu_DIM,&
        & Us,&
        & csu_delta,&
        & Vs,&
        & csu_lambdas,&
        & csu_gamma,&
        & SS &
    &)
        implicit none

        integer, intent(in) :: npoints, csu_DIM
        real, intent(in) :: csu_delta, csu_gamma
        real, intent(in) :: Us(npoints, csu_DIM), &
            &Vs(csu_DIM, csu_DIM), csu_lambdas(csu_DIM)
        real, intent(out):: SS(5)

        real, allocatable :: Ik(:, :)

        integer :: k, milestone, next_report

        real :: T1(csu_DIM), T2(csu_DIM), &
            &T3(csu_DIM), T4(csu_DIM), T5(csu_DIM)
        real :: weight_1(csu_DIM), weight_2(csu_DIM), weight_3(csu_DIM)

        allocate (Ik(npoints, csu_DIM))
        Ik(:, :) = 0.0
        milestone = max(1, csu_DIM / 10)
        next_report = milestone

        !$omp parallel do default(shared) private(k) schedule(dynamic)
        do k = 1, csu_DIM
            Ik(:, k) = matmul(Us, Vs(k, :))
            ! T1 = ∫ U_k dU_k   (Itô)
            call integrate_ito(&
                & npoints,&
                & Us(1:npoints, k),&
                & Us(1:npoints, k),&
                & T1(k)&
            )
            call integrate_ito(&
                & npoints,&
                & Ik(1:npoints, k),&
                & Us(1:npoints, k),&
                & T2(k)&
            &)
            call integrate_trapezoidal(&
                & npoints,&
                & Us(1:npoints, k)**2,&
                & csu_delta,&
                & T3(k)&
            &)
            call integrate_trapezoidal(&
                & npoints,&
                & Us(1:npoints, k) * Ik(1:npoints, k),&
                & csu_delta,&
                & T4(k)&
            &)
            call integrate_trapezoidal(&
                & npoints,&
                & Ik(1:npoints, k)**2,&
                & csu_delta, &
                & T5(k)&
            &)

            if (k >= next_report .or. k == csu_DIM) then
                !$omp critical
                if (k >= next_report .or. k == csu_DIM) then
                    write(&
                        &*, &
                        &'(A, A, I4, A, I6, A, I6, A)',&
                         advance='no'&
                    &) &
                        &achar(13), &
                        &"[MLE] compute_sufficient_statistics_u: ", &
                        &int(100.0 * real(k) / real(csu_DIM)), &
                        &"% (", k, " /", csu_DIM, " )"
                    if (k == csu_DIM) print *
                    next_report = next_report + milestone
                end if
                !$omp end critical
            end if
        end do
        !$omp end parallel do

        weight_1 = csu_lambdas**(1.0 + 2.0 * csu_gamma)
        weight_2 = csu_lambdas**(2.0 * csu_gamma)
        weight_3 = csu_lambdas**(2.0 + 2.0 * csu_gamma)

        SS(1) = sum(T1 * weight_1)  ! ∫ U_k dU_k weighted
        SS(2) = sum(T2 * weight_2)  ! ∫ I_k dU_k weighted
        SS(3) = sum(T3 * weight_3)  ! ∫ U_k^2 dt weighted
        SS(4) = sum(T4 * weight_1)  ! ∫ I_k U_k dt weighted
        SS(5) = sum(T5 * weight_2)  ! ∫ I_k^2 dt weighted
        deallocate (Ik)
        return
    end subroutine compute_sufficient_statistics_u

! ============================================================
! SUBROUTINE: estimate_mle_aux
!
! Description:
!   Computes the auxiliary maximum likelihood estimator of beta
!   assuming theta is known.
!
! Arguments:
!   npoints    (integer, input)  : Number of time increments
!   DIM        (integer, input)  : Dimension of the process
!   Us         (real, input)     : Observations of the U-process,
!                                  indexed as Us(0:npoints, DIM)
!   delta      (real, input)     : Time step size
!   Vs         (real, input)     : Transformation matrix (DIM x DIM)
!   lambdas    (real, input)     : Scaling parameters per dimension
!   gamma      (real, input)     : Scaling exponent
!   theta      (real, input)     : Known drift interaction parameter
!   beta_hat   (real, output)    : MLE of parameter beta
!
! ============================================================
    subroutine estimate_mle_u(&
        &npoints, &
        &emu_DIM, &
        &Us, &
        &emu_delta, &
        &Vs, &
        &emu_lambdas, &
        &emu_gamma, &
        &emu_theta, &
        &emu_beta_hat, &
        &emu_theta_hat&
    &)
        implicit none

        integer, intent(in) :: npoints, emu_DIM
        real, intent(in) :: emu_delta, emu_gamma, emu_theta
        real, intent(in) :: Us(npoints, emu_DIM)
        real, intent(in) :: Vs(emu_DIM, emu_DIM), emu_lambdas(emu_DIM)
        real, intent(out):: emu_beta_hat, emu_theta_hat
        real :: SS(5)
        real :: dem, numerator_beta, numerator_theta
        call compute_sufficient_statistics_u(&
            & npoints,&
            & emu_DIM,&
            & Us,&
            & emu_delta,&
            & Vs,&
            & emu_lambdas,&
            & emu_gamma,&
            & SS&
        &)

        dem = SS(4)**2 - SS(3) * SS(5)
        numerator_beta = SS(1) * SS(5) - SS(2) * SS(4)
        numerator_theta = SS(2) * SS(3) - SS(1) * SS(4)
        emu_beta_hat = numerator_beta / dem
        emu_theta_hat = numerator_theta / dem
        return
    end subroutine estimate_mle_u
end program multivariate
