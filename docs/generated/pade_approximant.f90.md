# `pade_approximant.f90`

## Overview

This Fortran module, `pade_approximant`, provides routines for constructing and evaluating Padé approximants using Thiele's reciprocal differences algorithm. It supports both a standard Thiele interpolation and a modified "greedy" version that aims to optimize the selection of interpolation points. The module also incorporates functionality to enforce various symmetries on the function being approximated. All calculations are performed in double precision complex arithmetic.

## Key Components

- **Module `pade_approximant`**:
    - Uses `kinds`, only: `dp` for double precision type.
    - Publicly exposes:
        - `thiele_pade`: Subroutine to compute Padé approximant parameters.
        - `evaluate_thiele_pade`: Subroutine to evaluate a Padé approximant at a given point.
        - `evaluate_thiele_pade_tab`: Subroutine (likely an internal helper or variant for `evaluate_thiele_pade`) that evaluates the approximant using a tabulated approach (Wallis' method iteration).
        - `c_zero`, `c_one`: Public complex(dp) parameters for zero and one.
    - Contains private helper functions and internal logic.

- **`subroutine thiele_pade(...)`**:
    - **Purpose**: Computes the parameters `a_par` of the Thiele Padé approximant.
    - **Inputs**:
        - `n_par`: Order of the interpolant.
        - `x_ref`: Array of reference points (modified if `do_greedy` is true).
        - `y_ref`: Array of reference function values.
        - `do_greedy`: Logical flag to choose between the standard or greedy algorithm for point selection.
        - `enforce_symmetry`: Character string specifying the symmetry to enforce (e.g., "even", "odd", "none").
    - **Output**:
        - `a_par`: Array of the interpolant parameters.
    - **Logic**:
        - Applies symmetry to `y_ref` values using `apply_sym_y`.
        - If `do_greedy` is true:
            - Selects points iteratively to minimize `|P_n(x_{i+1}) - F(x_{i+1})|`.
            - The first point maximizes `|F|`.
            - Uses `thiele_pade_gcoeff` to compute generating function coefficients.
            - Uses `evaluate_thiele_pade_tab` to compute convergents during point selection.
            - Rescales Wallis coefficients (`acoef`, `bcoef`) to prevent overflow.
        - If `do_greedy` is false:
            - Directly interpolates using `thiele_pade_gcoeff` for each parameter.

- **`subroutine thiele_pade_gcoeff(x, y, g_func, n, symm)`**:
    - **Purpose**: Computes recurrence coefficients `g_func` for Thiele's continued fraction using tabulation.
    - **Inputs**:
        - `n`: Number of parameters/points considered so far.
        - `x`: Array of reference points.
        - `y`: Array of reference function values.
        - `symm`: Symmetry string.
    - **Input/Output**:
        - `g_func`: Recurrence matrix.
    - **Logic**:
        - `g_func(n, 1) = y(n)`.
        - Iteratively computes `g_func(n, idx)` based on previous values and `apply_sym_x` for point transformation.

- **`subroutine evaluate_thiele_pade_tab(n_par, x_ref, x, a_par, acoef, bcoef, symm)`**:
    - **Purpose**: Evaluates one step of Wallis' method for Padé approximant evaluation.
    - **Inputs**: `n_par`, `x_ref`, evaluation point `x`, parameters `a_par`, symmetry `symm`.
    - **Input/Output**: Wallis coefficients `acoef`, `bcoef`.
    - **Logic**: Updates `acoef(n_par)` and `bcoef(n_par)` based on `acoef(n_par-1)`, `bcoef(n_par-1)`, `a_par(n_par)`, and symmetrized `x` and `x_ref` values.

- **`subroutine evaluate_thiele_pade(n_par, x_ref, x, a_par, y, enforce_symmetry)`**:
    - **Purpose**: Evaluates the full Padé approximant at point `x`.
    - **Inputs**: `n_par`, `x_ref`, evaluation point `x`, parameters `a_par`, `enforce_symmetry`.
    - **Output**: `y` (value of the interpolant at `x`).
    - **Logic**:
        - Initializes Wallis coefficients `acoef`, `bcoef`.
        - Iteratively calls `evaluate_thiele_pade_tab` for `i_par = 1` to `n_par`.
        - Rescales `acoef` and `bcoef` at each step to prevent overflow.
        - Computes final `y = acoef(n_par) / bcoef(n_par)`.
        - Applies symmetry to the final result `y` using `apply_sym_y(x, y, enforce_symmetry)`.

- **`function apply_sym_x(x, symm) result(x_symm)`**:
    - **Purpose**: Applies symmetry transformation to a complex argument `x`.
    - **Inputs**: `x` (complex number), `symm` (symmetry string).
    - **Output**: `x_symm` (transformed complex number).
    - **Logic**: Uses a `select case` statement based on `symm` string ("mirror_real", "mirror_imag", "mirror_both", "even", "odd", "conjugate", "anti-conjugate", "none") to modify `x`. Prints an error and stops if symmetry is unknown.

- **`function apply_sym_y(x_original, y, symm) result(y_symm)`**:
    - **Purpose**: Applies symmetry transformation to a function value `y`, potentially based on whether its argument `x_original` was changed by `apply_sym_x`.
    - **Inputs**: `x_original`, `y` (function value), `symm` (symmetry string).
    - **Output**: `y_symm` (transformed function value).
    - **Logic**:
        - Calls `apply_sym_x(x_original, symm)` to get `x_projected`.
        - If `x_projected == x_original`, `y_symm = y`.
        - Otherwise (if `x` was projected), applies transformations to `y` based on `symm` for "odd", "conjugate", "anti-conjugate" cases. For other symmetries, `y_symm = y`.

## Important Variables/Constants

- `dp`: Parameter from `kinds` module, defining double precision.
- `c_zero`, `c_one`: Complex double precision constants for 0.0+0.0i and 1.0+0.0i.
- `tol`: Real(dp) parameter (e.g., `1.0E-6_dp`) used as a tolerance, typically for checking division by zero or near-zero in Wallis' method rescaling.
- `acoef`, `bcoef`: Arrays used in Wallis' method for evaluating continued fractions, representing coefficients of the numerator and denominator polynomials.
- `g_func`: 2D array storing recurrence coefficients in Thiele's algorithm.
- `a_par`: Array storing the computed parameters of the Padé approximant (the `a_i` in the continued fraction form).

## Usage Examples

```fortran
! module pade_approximant
!   use kinds, only: dp
!   implicit none
!   ! ... (declarations as in the file) ...
! contains
!   ! ... (all subroutines and functions) ...
! end module pade_approximant

program test_pade
    use kinds, only: dp
    use pade_approximant
    implicit none

    integer, parameter :: n_points = 5
    complex(dp) :: x_pts(n_points), y_pts(n_points)
    complex(dp) :: a_params(n_points)
    complex(dp) :: eval_x, result_y
    integer :: i

    ! Example data (e.g., for f(z) = 1/(z+1))
    x_pts(1) = cmplx(0.0_dp, 0.1_dp, kind=dp)
    x_pts(2) = cmplx(0.1_dp, 0.2_dp, kind=dp)
    x_pts(3) = cmplx(0.2_dp, 0.3_dp, kind=dp)
    x_pts(4) = cmplx(0.3_dp, 0.4_dp, kind=dp)
    x_pts(5) = cmplx(0.4_dp, 0.5_dp, kind=dp)

    do i = 1, n_points
        y_pts(i) = c_one / (x_pts(i) + c_one)
    end do

    ! Compute Padé parameters (greedy, no specific symmetry)
    call thiele_pade(n_points, x_pts, y_pts, a_params, .true., "none")

    ! Print parameters (example)
    ! do i = 1, n_points
    !    print *, "a_param(", i, ") = ", a_params(i)
    ! end do

    ! Evaluate the Padé approximant at a new point
    eval_x = cmplx(0.5_dp, 0.6_dp, kind=dp)
    call evaluate_thiele_pade(n_points, x_pts, eval_x, a_params, result_y, "none")

    ! print *, "Pade at", eval_x, " = ", result_y
    ! print *, "Actual at", eval_x, " = ", c_one / (eval_x + c_one)

end program test_pade
```

## Dependencies and Interactions

- **`kinds` module**: Used for defining `dp` (double precision).
- The module is self-contained in terms of Padé logic but provides functions `apply_sym_x` and `apply_sym_y` for symmetry. These are Fortran implementations similar to those in `Symmetry_pade.cpp/.hpp` but take a character string for symmetry type instead of an integer.
- The greedy algorithm in `thiele_pade` uses `evaluate_thiele_pade_tab` internally to make decisions about point selection.
- `evaluate_thiele_pade` uses `evaluate_thiele_pade_tab` to perform the core evaluation steps.
- This module is likely a core component of the `GX-AnalyticContinuation` library, providing the numerical engine for Padé approximation.

The symmetry functions `apply_sym_x` and `apply_sym_y` are defined locally within this Fortran module. Their behavior should be consistent with the C++ `Symmetry_pade` class for corresponding symmetry types.
The `kind=8` in `cmplx` calls within `apply_sym_x` assumes `dp` corresponds to `kind=8`. It's generally safer to use `kind=dp`. This is a minor point but good for consistency. (It seems `dp` is used correctly elsewhere).

The `enforce_symmetry` parameter is passed around as a `character(*)` or `character(len=*)` string (e.g., "none", "even", "odd").
The `stop` statement in `apply_sym_x` upon an unknown symmetry string will terminate the program. More robust error handling might be desired in a library context (e.g., returning an error code).
