# `gx_ac.F90`

## Overview

This Fortran module, `gx_ac`, serves as a high-level API for performing Padé analytic continuation. It allows users to create a Padé model from a set of function values at given complex points and then evaluate this model at new points. The module can utilize either standard double-precision Fortran routines or, if compiled with GMP support, leverage C++ multi-precision routines for increased accuracy. It handles parameter storage, choice of precision, symmetry enforcement, and memory management.

The core idea is to approximate a function \(\f$f(x)\f$\) using a rational Padé approximant \(\f$P_N(x)\f$\) constructed via Thiele's reciprocal-difference method.

## Key Components

- **Module `gx_ac`**:
    - Imports `kinds` for `dp` (double precision), `iso_c_binding` for C interoperability.
    - Imports `evaluate_thiele_pade`, `thiele_pade` from the local `pade_approximant` module.
    - Publicly exposes:
        - `params`: Derived type to store Padé model parameters and state.
        - `create_thiele_pade`: Function to initialize and compute the Padé model.
        - `evaluate_thiele_pade_at`: Function to evaluate the model at specified points.
        - `free_params`: Subroutine to deallocate memory associated with a `params` object.
        - `arbitrary_precision_available`: Logical parameter indicating if GMP support is compiled.

- **`type :: params`**:
    - **Purpose**: Stores all information related to a Padé model.
    - **Members**:
        - `initialized`: Logical flag, true if the model has been created.
        - `n_par`: Integer, number of Padé parameters.
        - `precision`: Integer, precision in bits (64 for double, >64 for multi-precision).
        - `multiprecision_used_internally`: Logical flag, true if GMP routines are used.
        - `enforced_symmetry`: `character(len=15)`, string describing the enforced symmetry (e.g., "even", "none").
        - `use_greedy`: Logical flag, true if the greedy algorithm was used for point selection.
        - `zero_everywhere`: Logical flag, true if all input function values `y` were zero.
        - `params_ptr`: `type(c_ptr)`, C pointer to the C++ `pade_model` struct if multi-precision is used.
        - `x_ref`: `complex(kind=dp), dimension(:), allocatable`, stores reference x-points for Fortran double precision.
        - `a_par`: `complex(kind=dp), dimension(:), allocatable`, stores Padé parameters for Fortran double precision.

- **`arbitrary_precision_available`**:
    - A `logical, parameter` that is `.true.` if the code is compiled with `GMPXX_FOUND` preprocessor macro defined, `.false.` otherwise.

- **Interface block (if `GMPXX_FOUND` is defined)**:
    - Declares Fortran interfaces to C++ functions (defined in `pade_mp.cpp` and exposed via `pade_mp.h`):
        - `thiele_pade_mp_aux` (maps to C++ `thiele_pade_mp`): Creates the multi-precision Padé model.
        - `evaluate_thiele_pade_mp_aux` (maps to C++ `evaluate_thiele_pade_mp`): Evaluates the multi-precision model.
        - `free_pade_model` (maps to C++ `free_pade_model`): Frees memory of the C++ model.

- **`function create_thiele_pade(...) result(par)`**:
    - **Purpose**: API function to construct a Padé model.
    - **Inputs**:
        - `n_par`: Integer, number of parameters.
        - `x`, `y`: `complex(dp)` arrays of reference points and function values.
        - `do_greedy` (optional): Logical, whether to use the greedy algorithm.
        - `precision` (optional): Integer, desired precision in bits.
        - `enforce_symmetry` (optional): Character string for symmetry type.
    - **Output**: `par` (a `type(params)` object).
    - **Logic**:
        1. Initializes the `par` object.
        2. Determines precision: uses provided `precision`, defaults to 128 (if GMP) or 64 (if not). Sets `multiprecision_used_internally`.
        3. Validates and sets `enforced_symmetry` and a corresponding integer `c_symmetry_label` for C++ calls. Stops if symmetry is unknown.
        4. Sets `use_greedy` and `c_do_greedy` for C++ calls.
        5. Checks for the edge case where all input `y` values are zero; if so, sets `par%zero_everywhere = .true.` and returns.
        6. If not using multi-precision:
            - Allocates `par%a_par` and `par%x_ref`.
            - Calls the Fortran `thiele_pade` subroutine.
        7. If using multi-precision (and GMP is available):
            - Calls `thiele_pade_mp_aux` to get a C pointer to the C++ model, stores it in `par%params_ptr`.

- **`function evaluate_thiele_pade_at(par, x) result(y)`**:
    - **Purpose**: API function to evaluate an existing Padé model.
    - **Inputs**:
        - `par`: `type(params)` object containing the model.
        - `x`: `complex(dp)` array of points to evaluate at.
    - **Output**: `y` (`complex(dp)` array of results).
    - **Logic**:
        1. Checks if `par` is initialized; if not, prints a warning and returns zeros.
        2. Checks if `par%zero_everywhere` is true; if so, returns zeros.
        3. Iterates through each point in `x`:
            - If not using multi-precision: Calls Fortran `evaluate_thiele_pade`.
            - If using multi-precision (and GMP is available): Calls `evaluate_thiele_pade_mp_aux` using `par%params_ptr`.

- **`subroutine free_params(par)`**:
    - **Purpose**: Deallocates memory held by a `params` object.
    - **Input/Output**: `par` (`type(params)` object).
    - **Logic**:
        1. If `par` is not initialized, returns.
        2. Deallocates `par%a_par` and `par%x_ref` if they were allocated (for Fortran double precision mode).
        3. If multi-precision was used and not `zero_everywhere` (and GMP is available): Calls `free_pade_model` with `par%params_ptr` to free the C++ allocated memory.

## Important Variables/Constants

- `dp`: Double precision kind parameter.
- `params`: The derived type itself is central to this API.
- `arbitrary_precision_available`: Compile-time constant indicating GMP availability.
- Preprocessor macro `GMPXX_FOUND`: Controls conditional compilation for GMP-related code sections.

## Usage Examples

```fortran
program demo_gx_ac
    use gx_ac
    use kinds, only: dp
    implicit none

    integer, parameter :: num_pts = 10
    complex(dp) :: x_data(num_pts), y_data(num_pts)
    complex(dp) :: eval_points(3), results(3)
    type(params) :: pade_model
    integer :: i

    ! 1. Prepare data
    do i = 1, num_pts
        x_data(i) = cmplx(real(i-1)*0.1_dp, 0.0_dp, kind=dp)
        y_data(i) = sin(x_data(i)%re) ! Example: sin(x)
    end do

    ! 2. Create Pade model
    ! Using default precision (128 if GMP available, else 64)
    ! Using default greedy algorithm, no symmetry
    pade_model = create_thiele_pade(num_pts, x_data, y_data)
    ! Or, with options:
    ! pade_model = create_thiele_pade(num_pts, x_data, y_data, &
    !                                 precision=256, enforce_symmetry="even", do_greedy=.true.)


    if (pade_model%initialized) then
        ! 3. Evaluate model
        eval_points(1) = cmplx(0.25_dp, 0.0_dp, kind=dp)
        eval_points(2) = cmplx(0.75_dp, 0.0_dp, kind=dp)
        eval_points(3) = cmplx(1.5_dp, 0.0_dp, kind=dp) ! Extrapolation

        results = evaluate_thiele_pade_at(pade_model, eval_points)

        do i = 1, size(eval_points)
            print *, "Pade at ", eval_points(i), " = ", results(i)
            print *, "Actual sin at ", eval_points(i)%re, " = ", sin(eval_points(i)%re)
        end do

        ! 4. Free model parameters
        call free_params(pade_model)
    else
        print *, "Failed to initialize Pade model."
    end if

end program demo_gx_ac
```

## Dependencies and Interactions

- **`kinds` module**: For `dp`.
- **`iso_c_binding`**: For C interoperability.
- **`pade_approximant` module**: For the double-precision Fortran implementations (`thiele_pade`, `evaluate_thiele_pade`).
- **C++ implementation (`pade_mp.cpp`, `pade_mp.h`, `ComplexGMP.cpp/.hpp`, `Symmetry_pade.cpp/.hpp`)**: Required if `GMPXX_FOUND` is defined and multi-precision capabilities are used. The Fortran module links to these C++ functions via `bind(C)`.
- **GMP library**: External dependency required for compiling the C++ multi-precision parts and for `arbitrary_precision_available` to be true.

This module acts as a comprehensive Fortran wrapper, abstracting away the direct calls to either local Fortran routines or external C++ routines based on user choice (precision) and compile-time availability (GMP). It provides a clean Fortranic interface for Padé analytic continuation.
The error handling for unknown symmetry strings in `create_thiele_pade` uses `stop`, which will terminate the program.
The `zero_everywhere` check is a practical way to handle cases where the Padé approximant is trivially zero, which might otherwise cause issues in the underlying algorithms.
