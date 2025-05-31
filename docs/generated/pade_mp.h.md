# `pade_mp.h`

## Overview

This header file defines the structure `pade_model` used to store parameters and results of a multi-precision Padé approximation. It also declares the C-style functions, intended for potential Fortran interface, that perform the Padé approximation (`thiele_pade_mp`), evaluate it (`evaluate_thiele_pade_mp`), and free the allocated memory for the model (`free_pade_model`).

## Key Components

- **`struct pade_model`**:
    - **Purpose**: A C-style struct to encapsulate all necessary information about a Padé model constructed using multi-precision arithmetic.
    - **Members**:
        - `int n_par`: Stores the number of Padé parameters (order of the approximant).
        - `int precision`: Stores the floating-point precision in bits used for the GMP calculations.
        - `Symmetry_pade symmetry`: An object of the `Symmetry_pade` class, used to handle and enforce symmetries on the Padé interpolant.
        - `ComplexGMP *a_par`: Pointer to a dynamically allocated array of `ComplexGMP` objects, storing the Padé coefficients.
        - `ComplexGMP *xref`: Pointer to a dynamically allocated array of `ComplexGMP` objects, storing the reference x-points used to construct the approximant. These points might be reordered if a greedy algorithm was used.

- **Function Declarations (within `extern "C"` block)**:
    - **`std::complex<double> evaluate_thiele_pade_mp(const std::complex<double> x, pade_model *params_ptr)`**:
        - **Purpose**: Declares the function (implemented in `pade_mp.cpp`) that evaluates a previously computed multi-precision Padé model at a given complex double-precision point `x`.
        - **Parameters**:
            - `x`: The `std::complex<double>` point at which to evaluate the model.
            - `params_ptr`: A pointer to the `pade_model` struct containing the approximant's data.
        - **Returns**: The evaluated result as a `std::complex<double>`.

    - **`pade_model *thiele_pade_mp(int n_par, const std::complex<double> *x_ref, const std::complex<double> *y_ref, int do_greedy, int precision, int symmetry)`**:
        - **Purpose**: Declares the function (implemented in `pade_mp.cpp`) that constructs a multi-precision Padé approximant.
        - **Parameters**:
            - `n_par`: Desired number of Padé parameters.
            - `x_ref`: Pointer to an array of `std::complex<double>` reference x-points.
            - `y_ref`: Pointer to an array of `std::complex<double>` reference y-values (function values at `x_ref`).
            - `do_greedy`: Integer flag; if non-zero (typically 1), a greedy algorithm for point selection is used.
            - `precision`: Desired precision in bits for GMP calculations.
            - `symmetry`: Integer label specifying the symmetry type to be enforced (corresponds to `Symmetry_pade` enum/constants).
        - **Returns**: A pointer to a newly allocated and populated `pade_model` struct. The caller is responsible for freeing this memory using `free_pade_model`.

    - **`void free_pade_model(pade_model *model)`**:
        - **Purpose**: Declares and implements an inline function to deallocate the memory associated with a `pade_model` struct.
        - **Parameters**:
            - `model`: A pointer to the `pade_model` struct to be freed.
        - **Logic**:
            - Checks if `model` is not null.
            - Deletes the dynamically allocated arrays `model->a_par` and `model->xref` using `delete[]`.
            - Deletes the `model` struct itself using `delete`.

## Important Variables/Constants

- The members of the `pade_model` struct are the primary "variables" defined by this header in terms of data structure.

## Usage Examples

```cpp
// Conceptual C++ usage (illustrating the interface)
// #include "pade_mp.h" // This file
// #include <vector>
// #include <complex>
// #include <iostream>

// int main() {
    // int num_points = 5;
    // std::vector<std::complex<double>> x_data(num_points);
    // std::vector<std::complex<double>> y_data(num_points);

    // // ... Populate x_data and y_data ...
    // for(int i=0; i<num_points; ++i) {
    //    x_data[i] = std::complex<double>(0.1 * i, 0.0);
    //    y_data[i] = std::sin(x_data[i].real()); // Example: sin(x)
    // }

    // int precision_bits = 128;
    // int use_greedy = 1;
    // int symmetry_type = 0; // e.g., sym_none from Symmetry_pade

    // // Construct the Pade model
    // pade_model *my_pade = thiele_pade_mp(num_points, x_data.data(), y_data.data(),
    //                                      use_greedy, precision_bits, symmetry_type);

    // if (my_pade) {
    //     std::complex<double> eval_point(0.25, 0.0);
    //     std::complex<double> result = evaluate_thiele_pade_mp(eval_point, my_pade);
    //     std::cout << "Pade approximation at " << eval_point << " is " << result << std::endl;

    //     // Free the model
    //     free_pade_model(my_pade);
    //     my_pade = nullptr; // Good practice
    // }

    // return 0;
// }
```

```fortran
! Conceptual Fortran usage (illustrating C interoperability)
! interface
!     function thiele_pade_mp_c(n_par, x_ref, y_ref, do_greedy, precision, symmetry) &
!                              bind(C, name='thiele_pade_mp')
!         use iso_c_binding
!         integer(c_int), value :: n_par, do_greedy, precision, symmetry
!         type(c_ptr) :: thiele_pade_mp_c
!         complex(c_double_complex), dimension(*) :: x_ref, y_ref
!     end function thiele_pade_mp_c

!     function evaluate_thiele_pade_mp_c(x, params_ptr) &
!                                      bind(C, name='evaluate_thiele_pade_mp')
!         use iso_c_binding
!         complex(c_double_complex) :: evaluate_thiele_pade_mp_c
!         complex(c_double_complex), value :: x
!         type(c_ptr), value :: params_ptr
!     end function evaluate_thiele_pade_mp_c

!     subroutine free_pade_model_c(model_ptr) bind(C, name='free_pade_model')
!         use iso_c_binding
!         type(c_ptr), value :: model_ptr
!     end subroutine free_pade_model_c
! end interface

! type(c_ptr) :: model_handle
! complex(c_double_complex) :: result_val
! complex(c_double_complex), target, allocatable :: x_f(:), y_f(:)
! integer :: n = 5
! ! ... allocate and fill x_f, y_f ...
! model_handle = thiele_pade_mp_c(n, x_f(1), y_f(1), 1, 128, 0)
! if (c_associated(model_handle)) then
!    result_val = evaluate_thiele_pade_mp_c(cmplx(0.25d0, 0.0d0, kind=c_double_complex), model_handle)
!    ! ... use result_val ...
!    call free_pade_model_c(model_handle)
! end if
```

## Dependencies and Interactions

- **`Symmetry_pade.hpp`**: Included for the `Symmetry_pade` class definition, which is a member of the `pade_model` struct.
- **`ComplexGMP.hpp`**: Included for the `ComplexGMP` class definition, used for the `a_par` and `xref` members of `pade_model`.
- **`<gmpxx.h>`**: Standard GMP C++ header, likely included by `ComplexGMP.hpp` as well, but explicitly included here.
- **`<complex>`**: Standard C++ header for `std::complex<double>`.
- **`pade_mp.cpp`**: This header file declares functions that are implemented in `pade_mp.cpp`.
- The `extern "C"` linkage specification indicates that these functions are intended to be callable from C and, by extension, Fortran (using `iso_c_binding`). This is crucial for the `GX-AnalyticContinuation` library's potential use in mixed-language environments.
- The `free_pade_model` function is critical for preventing memory leaks when using the dynamically allocated `pade_model` objects.

This header acts as the primary interface for the multi-precision Padé functionality when used from C++ or linked from Fortran.
