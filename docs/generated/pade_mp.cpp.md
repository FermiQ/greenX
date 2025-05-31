# `pade_mp.cpp`

## Overview

This C++ source file implements functions for computing and evaluating Padé approximants using Thiele's reciprocal differences algorithm with multi-precision arithmetic, leveraging the `ComplexGMP` class. It provides a "greedy" point selection strategy for constructing the approximant and includes mechanisms for enforcing symmetries. The results and parameters of the Padé model are stored in a `pade_model` struct.

## Key Components

- **`thiele_pade_gcoeff_mp(...)`**:
    - **Purpose**: Computes recurrence coefficients for Thiele's continued fraction using `ComplexGMP` multi-precision numbers.
    - **Parameters**:
        - `x`: Array of `ComplexGMP` reference points.
        - `y`: Vector of `ComplexGMP` reference function values.
        - `g_func`: 2D vector (matrix) to store `ComplexGMP` recurrence coefficients.
        - `n`: Current number of parameters/points being considered (0-indexed).
    - **Logic**: Populates the `n`-th row of `g_func` based on Thiele's recurrence relations.

- **`evaluate_thiele_pade_tab_mp<any_complex>(...)`**:
    - **Purpose**: Template function to perform one step of Wallis' method for evaluating the Padé approximant. It's templated to potentially work with types other than `ComplexGMP`, though used with `ComplexGMP` in this file.
    - **Parameters**:
        - `n_par`: Current parameter index in Wallis' method (0-indexed).
        - `x_ref`: Array of reference points of type `any_complex`.
        - `x`: The point at which to evaluate, of type `any_complex`.
        - `a_par`: Array of Padé parameters of type `any_complex`.
        - `acoef`, `bcoef`: Vectors to store Wallis' coefficients (numerator and denominator) of type `any_complex`.
    - **Logic**: Updates `acoef[n_par + 1]` and `bcoef[n_par + 1]`.

- **`evaluate_thiele_pade_mp(...)`**:
    - **Purpose**: Evaluates a pre-computed multi-precision Padé approximant at a given `std::complex<double>` point `x`.
    - **Parameters**:
        - `x`: The `std::complex<double>` point for evaluation.
        - `params`: Pointer to a `pade_model` struct containing the Padé model parameters (precision, coefficients, reference points, symmetry object).
    - **Returns**: `std::complex<double>` value of the interpolant at `x`.
    - **Logic**:
        1. Sets GMP default precision from `params->precision`.
        2. Applies symmetry to input `x` using `params->symmetry.apply_sym_x()`.
        3. Converts symmetrized `x` to `ComplexGMP`.
        4. Initializes Wallis' coefficients (`acoef`, `bcoef`) using `ComplexGMP`.
        5. Iteratively calls `evaluate_thiele_pade_tab_mp<ComplexGMP>` to compute the continued fraction.
        6. Rescales `acoef` and `bcoef` at each step if `bcoef[i_par + 1].abs()` is greater than a tolerance `tol` (which is set to "1" - this might be an issue, usually tol is small).
        7. Converts the `ComplexGMP` result to `std::complex<double>`.
        8. Applies symmetry to the output `y` using `params->symmetry.apply_sym_y()`.

- **`thiele_pade_mp(...)`**:
    - **Purpose**: Constructs the multi-precision Padé approximant.
    - **Parameters**:
        - `n_par`: Order of the interpolant.
        - `x_ref`: Array of `std::complex<double>` reference points.
        - `y_ref`: Array of `std::complex<double>` reference function values.
        - `do_greedy`: Integer flag (1 for greedy, 0 for naive).
        - `precision`: Integer for GMP floating-point precision in bits.
        - `symmetry`: Integer label for symmetry type.
    - **Returns**: Pointer to a `pade_model` struct populated with the computed model.
    - **Logic**:
        1. Sets GMP default precision.
        2. Creates and initializes a `pade_model` struct, including setting symmetry.
        3. Applies symmetry to input `x_ref` and `y_ref` using `params->symmetry.apply_sym_x` and `apply_sym_y`, converting them to `ComplexGMP` and storing in `x_projected`, `y_projected`. (Note: `x_ref` and `y_ref` are initially `std::complex<double>`, then converted to `ComplexGMP` for internal processing).
        4. Initializes `ComplexGMP` arrays/vectors for internal calculations (`a_par_mp`, `x_ref_mp`, `xtmp`, `ytmp`, `g_func`).
        5. **If `do_greedy == 1`**:
            - Implements a specific greedy algorithm:
                - Selects the first point `xtmp[0]` where `|y_ref_mp|` is maximum. Computes `a_par_mp[0]`.
                - Selects the second point `xtmp[1]` that maximizes `|x[i] - xtmp[0]|` among remaining points. Computes `a_par_mp[1]`.
                - Initializes Wallis coefficients.
                - Iteratively selects remaining points (`idx = 2` to `n_par-1`) by finding the point `x_in` that minimizes `|P_idx(x_in) - y_in|`, where `P_idx` is the current Padé approximant.
                - Uses `evaluate_thiele_pade_tab_mp` to calculate `P_idx`.
                - Updates reference points (`x_ref_mp`, `xtmp`, `ytmp`).
                - Rescales Wallis coefficients.
                - Calls `thiele_pade_gcoeff_mp` to update `g_func`.
                - Extracts `a_par_mp[idx]` from `g_func`.
        6. **Else (not greedy)**:
            - Iteratively calls `thiele_pade_gcoeff_mp` for each point to compute `a_par_mp`.
        7. Stores `a_par_mp` (Padé coefficients) and `x_ref_mp` (potentially reordered and symmetrized reference points) into the `params` struct.
        8. Cleans up temporary `xtmp`. Returns `params`.

## Important Variables/Constants

- `ComplexGMP`: Class used for multi-precision complex number arithmetic (defined in `ComplexGMP.hpp`).
- `pade_model`: Struct (defined in `pade_mp.h`) that holds all parameters and results of the Padé approximation, including precision, number of parameters, symmetry object, Padé coefficients (`a_par`), and reference points (`xref`).
- `precision`: Integer specifying the number of bits for GMP multi-precision floats.
- `tol`: An `mpf_class` constant used for tolerance checks. In `evaluate_thiele_pade_mp`, it's initialized as `mpf_class("1", params->precision)`. This is unusual as tolerance is typically a small value (e.g., 1e-16). A tolerance of 1 means rescaling happens only if the absolute value of `bcoef` is greater than 1. This might need review. In `thiele_pade_mp`, `tol` is `mpf_class("1e-6", precision)`.
- `g_func`: 2D vector storing recurrence coefficients (`ComplexGMP`).
- `a_par_mp`: Dynamically allocated array storing Padé coefficients (`ComplexGMP`).
- `x_ref_mp`: Dynamically allocated array storing reference x-values used for the model (`ComplexGMP`).

## Usage Examples

```cpp
// main.cpp (conceptual)
// #include "pade_mp.h"
// #include "ComplexGMP.hpp" // Assumed to be included by pade_mp.h or directly
// #include <vector>
// #include <complex>
// #include <iostream>

// int main() {
    // int n_points = 10;
    // int precision_bits = 256; // Example precision
    // std::vector<std::complex<double>> x_values(n_points);
    // std::vector<std::complex<double>> y_values(n_points);

    // // ... Fill x_values and y_values with data ...
    // for(int i=0; i<n_points; ++i) {
    //     x_values[i] = std::complex<double>(i*0.1, 0.05);
    //     // Example: y = 1.0 / (x + 1.0)
    //     y_values[i] = 1.0 / (x_values[i] + std::complex<double>(1.0,0.0));
    // }

    // // Construct Padé model (greedy, no symmetry example)
    // // Symmetry labels would come from Symmetry_pade.hpp (e.g., 0 for none)
    // pade_model* model = thiele_pade_mp(n_points, x_values.data(), y_values.data(),
    //                                    1, precision_bits, 0 /*sym_none*/);

    // if (model) {
    //     // Evaluate the model at a new point
    //     std::complex<double> test_x(0.5, 0.5);
    //     std::complex<double> result_y = evaluate_thiele_pade_mp(test_x, model);
    //     std::cout << "Pade at " << test_x << " = " << result_y << std::endl;

    //     // Remember to free memory allocated for the model
    //     delete[] model->a_par;
    //     delete[] model->xref;
    //     delete model;
    // }
    // return 0;
// }
```

## Dependencies and Interactions

- **`pade_mp.h`**: Header file containing the definition of `pade_model` struct and declarations for the functions in this `.cpp` file.
- **`ComplexGMP.hpp`**: Provides the `ComplexGMP` class for multi-precision arithmetic. This is a fundamental dependency.
- **`Symmetry_pade.hpp`**: (Assumed to be included via `pade_mp.h` where `pade_model::symmetry` is defined) Provides the `Symmetry_pade` class used within `pade_model` to handle symmetries.
- **GMP Library**: The GNU Multiple Precision Arithmetic Library is required by `ComplexGMP`. `mpf_set_default_prec` is a GMP function.
- **`<vector>`, `<complex>`, `<iostream>`, `<optional>`, `<algorithm>` (for `std::max_element`, `std::distance`, `std::swap`)**: Standard C++ libraries.

**Notes on Greedy Algorithm in `thiele_pade_mp`**:
- The greedy strategy seems to have specific choices for the first two points (`a_0`, `a_1` coefficients):
    - Point 0: Maximizes `|F(x)|`.
    - Point 1: Maximizes `|x - x_0|`.
- For subsequent points, it minimizes `|P_i(x_{i+1}) - F(x_{i+1})|`. This is a more common greedy criterion.
- The use of `std::vector<ComplexGMP> y_ref_mp;` and `y_ref_mp.push_back()` suggests `y_ref_mp` might not be pre-sized correctly if it's meant to be a direct copy of input `y_ref` after symmetrization. It's populated correctly during initialization.

**Memory Management**:
- `thiele_pade_mp` allocates memory for `pade_model`, `a_par_mp` (becomes `model->a_par`), and `x_ref_mp` (becomes `model->xref`). The caller is responsible for deleting these to prevent memory leaks (e.g., `delete[] model->a_par; delete[] model->xref; delete model;`).
- `xtmp` is allocated with `new` and `delete[]`d within `thiele_pade_mp`.

**Tolerance `tol` in `evaluate_thiele_pade_mp`**:
- The `tol` value set to `mpf_class("1", params->precision)` in `evaluate_thiele_pade_mp` is unusual for a typical floating-point comparison tolerance where small values are expected. If `bcoef[i_par+1].abs()` is, for example, 0.5, the rescaling `if (bcoef[i_par + 1].abs() > tol)` would be false, and division by a small `bcoef` might occur in the final step without this intermediate rescaling. This should be verified. The `thiele_pade_mp` function uses `mpf_class("1e-6", precision)` which is more conventional. This discrepancy is noteworthy.
