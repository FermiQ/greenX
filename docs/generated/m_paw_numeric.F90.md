# `m_paw_numeric.F90`

## Overview

The Fortran module `m_paw_numeric` provides a collection of numerical utility routines that are likely used in various parts of the `libPAW` library. These routines cover common numerical tasks such as spline interpolation, data smoothing, sorting, calculation of spherical Bessel functions, root finding for Bessel-related equations, and evaluation of the complementary error function.

## Key Components

### Public Subroutines and Functions

-   **`paw_spline(t, y, n, ybcbeg, ybcend, ypp)`**:
    -   **Purpose**: Computes the second derivatives (`ypp`) of a cubic spline that interpolates a given set of data points (`t`, `y`).
    -   **Inputs**:
        -   `t(n)`: `real(dp)`, array of knot values (abscissas), must be distinct and increasing.
        -   `y(n)`: `real(dp)`, array of data values at knots `t`.
        -   `n`: Integer, number of data points (must be \(\ge 2\)).
        -   `ybcbeg`, `ybcend`: `real(dp)`, boundary condition values for the first derivative at the beginning and end of the interval. If > 1.0e30, natural boundary conditions (second derivative = 0) are assumed.
    -   **Output**:
        -   `ypp(n)`: `real(dp)`, array of second derivatives of the spline at the knots.
    -   **Notes**: Adapted from a routine by John Burkardt, modified by Xavier Gonze. Solves a tridiagonal system for `ypp`.

-   **`paw_splint(nspline, xspline, yspline, ysplin2, nfit, xfit, yfit, ierr)`**:
    -   **Purpose**: Performs cubic spline interpolation. Given a tabulated function (`xspline`, `yspline`) and its pre-computed second derivatives (`ysplin2` from `paw_spline`), it evaluates the spline at new points `xfit`.
    -   **Inputs**:
        -   `nspline`: Integer, number of points in the tabulated function.
        -   `xspline(nspline)`: `real(dp)`, abscissas of the tabulated function.
        -   `yspline(nspline)`: `real(dp)`, values of the tabulated function.
        -   `ysplin2(nspline)`: `real(dp)`, second derivatives at `xspline`.
        -   `nfit`: Integer, number of points to interpolate at.
        -   `xfit(nfit)`: `real(dp)`, abscissas for interpolation.
    -   **Output**:
        -   `yfit(nfit)`: `real(dp)`, interpolated values at `xfit`.
        -   `ierr` (optional): Integer, error flag incremented if `xfit` points are outside `xspline` range.

-   **`paw_splint_der(nspline, xspline, yspline, ysplin2, nfit, xfit, dydxfit, ierr)`**:
    -   **Purpose**: Similar to `paw_splint`, but computes the first derivative of the spline at the points `xfit`.
    -   **Output**:
        -   `dydxfit(nfit)`: `real(dp)`, first derivatives of the spline at `xfit`.

-   **`paw_uniform_splfit(arg, derfun, fun, ider, newarg, newfun, numarg, numnew)`**:
    -   **Purpose**: Evaluates a cubic spline (function and/or derivatives) where the original data `fun` was defined on a uniformly spaced grid `arg`.
    -   **Inputs**:
        -   `arg(numarg)`: `real(dp)`, uniformly spaced abscissas of the original data.
        -   `fun(numarg,2)`: `real(dp)`, original function values `fun(:,1)` and their second derivatives `fun(:,2)`.
        -   `ider`: Integer, controls output: 0 for function only (`newfun`), 1 for function (`newfun`) and first derivative (`derfun`), 2 for second derivative only (`derfun`).
        -   `newarg(numnew)`: `real(dp)`, new abscissas for evaluation.
        -   `numarg`, `numnew`: Integers, sizes of `arg` and `newarg`.
    -   **Outputs**:
        -   `newfun(numnew)`: `real(dp)` (intent inout), function values at `newarg`.
        -   `derfun(numnew)`: `real(dp)`, derivative values at `newarg` (if `ider=1` or `ider=2`).
    -   **Notes**: If `newarg` points are outside the range of `arg`, extremal values are assigned.

-   **`paw_smooth(a, mesh, it)`**:
    -   **Purpose**: Smooths an array `a` of `mesh` equally spaced ordinates over `it` iterations.
    -   **Inputs**:
        -   `mesh`: Integer, size of the array `a`.
        -   `it`: Integer, number of smoothing iterations.
    -   **Input/Output**:
        -   `a(mesh)`: `real(dp)`, array to be smoothed.
    -   **Logic**: Uses a 10-point weighted average for interior points and simpler averages near boundaries. The weights seem to be fixed (e.g., `0.1*a(i) + 0.1*(a(i+1)+a(i-1)) + ...`).

-   **`paw_sort_dp(n, list, iperm, tol)`**:
    -   **Purpose**: Sorts a `real(dp)` array `list` into ascending order using the Heapsort algorithm. It also rearranges an integer array `iperm` (initially `iperm(i)=i`) to track the original indices.
    -   **Inputs**:
        -   `n`: Integer, dimension of `list` and `iperm`.
        -   `tol`: `real(dp)`, tolerance for considering two numbers equal. If numbers are within `tol`, their original relative order (from `iperm`) is maintained.
    -   **Input/Output**:
        -   `list(n)`: `real(dp)`, array to be sorted.
        -   `iperm(n)`: Integer, permutation array.

-   **`paw_jbessel(bes, besp, bespp, ll, order, xx)`**:
    -   **Purpose**: Computes the spherical Bessel function \(j_l(x)\) and optionally its first (`besp`) and second (`bespp`) derivatives.
    -   **Inputs**:
        -   `ll`: Integer, order \(l\) of the Bessel function.
        -   `order`: Integer, 0 for \(j_l\) only, 1 for \(j_l, j_l'\), 2 for \(j_l, j_l', j_l''\).
        -   `xx`: `real(dp)`, argument \(x\).
    -   **Outputs**: `bes`, `besp`, `bespp`.
    -   **Logic**: Uses series expansion for small \(x\) and standard recurrence relations or trigonometric forms for larger \(x\).

-   **`paw_solvbes(root, alpha, beta, ll, nq)`**:
    -   **Purpose**: Finds the first `nq` roots of the intrinsic equation \(\alpha j_l(Q) + \beta Q j_l'(Q) = 0\).
    -   **Inputs**:
        -   `alpha`, `beta`: `real(dp)`, coefficients in the equation.
        -   `ll`: Integer, order of the Bessel function.
        -   `nq`: Integer, number of roots to find.
    -   **Output**:
        -   `root(nq)`: `real(dp)`, array of found roots \(Q\).
    -   **Logic**: Uses a simple root-finding algorithm: steps along \(Q\) with `dh` until a sign change in the function is detected, then refines with bisection.

-   **`paw_jbessel_4spline(bes, besp, ll, order, xx, tol)`**:
    -   **Purpose**: Computes spherical Bessel functions \(j_l(x)\) and first derivatives, using polynomial approximations for \(x < tol\).
    -   **Inputs**:
        -   `ll`: Integer, order (0 to 3 handled with specific polynomials, \(\ge 4\) calls `paw_jbessel`).
        -   `order`: Integer, 0 for \(j_l\), 1 for \(j_l, j_l'\).
        -   `xx`: `real(dp)`, argument.
        -   `tol`: `real(dp)`, tolerance below which polynomial approximation is used.
    -   **Outputs**: `bes`, `besp`.
    -   **Notes**: Contains explicit polynomial expansions for \(j_0, j_1, j_2, j_3\) and their derivatives for small \(x\). The comment "Statement functions are obsolete. Sorry ..." suggests this might be older code or adapted from code that used statement functions.

-   **`paw_derfc(yy)`**: Elemental function.
    -   **Purpose**: Evaluates the complementary error function \(\text{erfc}(y)\).
    -   **Input**: `yy`: `real(dp)`.
    -   **Output**: `derfc_yy`: `real(dp)`.
    -   **Logic**: Uses different polynomial approximations based on the range of `abs(yy)` (0 to 0.477, 0.477 to 4.0, >4.0). Handles negative `yy` via \(\text{erfc}(-y) = 2 - \text{erfc}(y)\).

## Important Variables/Constants

-   The module uses constants like `zero`, `one`, `half`, `third`, `pi`, `tol*` from `m_libpaw_defs`.
-   Internal parameters in `paw_derfc` (`pp`, `qq`, `p1`, `q1`, `p2`, `q2`, `sqrpi`, `xbig`, etc.) are coefficients for the polynomial approximations of erfc.

## Usage Examples

These routines would be called by other `libPAW` modules for various numerical tasks.

```fortran
module m_example_numeric_usage
    use m_paw_numeric
    use m_libpaw_defs, only: dp
    implicit none

    subroutine test_spline_and_bessel
        real(dp), dimension(5) :: x_pts = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
        real(dp), dimension(5) :: y_pts = [0.0_dp, 0.8415_dp, 0.9093_dp, 0.1411_dp, -0.7568_dp] ! sin(x)
        real(dp), dimension(5) :: y_pp
        real(dp), dimension(3) :: x_interp = [0.5_dp, 1.5_dp, 2.5_dp]
        real(dp), dimension(3) :: y_interp
        integer :: ierr_splint

        real(dp) :: bessel_j1, bessel_j1_prime, bessel_j1_prime_prime

        ! Spline example
        call paw_spline(x_pts, y_pts, 5, 1.0_dp, -0.9899_dp, y_pp) ! Boundary derivatives: cos(0), cos(4)
        call paw_splint(5, x_pts, y_pts, y_pp, 3, x_interp, y_interp, ierr_splint)
        ! print *, "Interpolated sin(x) at ", x_interp, " are ", y_interp

        ! Bessel function example
        call paw_jbessel(bessel_j1, bessel_j1_prime, bessel_j1_prime_prime, 1, 2, 1.5_dp)
        ! print *, "j1(1.5) =", bessel_j1, "j1'(1.5) =", bessel_j1_prime

        ! Error function example
        ! print *, "erfc(1.0) =", paw_derfc(1.0_dp)

    end subroutine test_spline_and_bessel

end module m_example_numeric_usage
```

## Dependencies and Interactions

-   **`m_libpaw_defs` (`USE_DEFS`)**: For `dp` and various mathematical/physical constants.
-   **`m_libpaw_tools` (`USE_MSG_HANDLING`)**: For error reporting (`LIBPAW_ERROR`, `LIBPAW_BUG`).
-   **`m_libpaw_memory` (`USE_MEMORY_PROFILING`)**: For memory allocation macros (though not heavily used directly in this module, as most arrays are passed as arguments).

This module provides a set of general-purpose numerical tools that are likely fundamental for the accurate and efficient implementation of the PAW method and related physical calculations within `libPAW`.
