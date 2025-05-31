# `Symmetry_pade.cpp`

## Overview

This file implements the functionality for applying various symmetries to complex-valued functions, specifically in the context of Padé approximants. It defines how function arguments (`x`) and function values (`y`) are transformed based on a specified symmetry type.

## Key Components

- **`Symmetry_pade::Symmetry_pade(int symm)`**: Constructor that initializes the symmetry type.
    - Takes an integer `symm` (symmetry label) as input.
    - Validates if `symm` is within the defined range (0-7). Throws `std::invalid_argument` if out of range.
    - Stores the symmetry type in the `my_symmetry` member variable.

- **`Symmetry_pade::set_symmetry(int symm)`**: Sets or changes the symmetry type.
    - Takes an integer `symm` (symmetry label) as input.
    - Validates `symm` and updates `my_symmetry`. Throws `std::invalid_argument` if out of range.

- **`Symmetry_pade::apply_sym_x(std::complex<double> x)`**: Applies symmetry transformation to the function argument `x`.
    - Takes a `std::complex<double> x` as input.
    - Returns a `std::complex<double>` representing the transformed argument `x_symm`.
    - The transformation depends on the `my_symmetry` value:
        - `sym_none`: No change.
        - `sym_mirror_real`: `x_symm = (|Re(x)|, Im(x))`.
        - `sym_mirror_imag`: `x_symm = (Re(x), |Im(x)|)`.
        - `sym_mirror_both`: `x_symm = (|Re(x)|, |Im(x)|)`.
        - `sym_even`, `sym_odd`, `sym_conjugate`, `sym_anti_conjugate`: `x_symm = (|Re(x)|, copysign(1.0, Re(x)) * Im(x))`. (The logic seems identical for these cases in `apply_sym_x`).

- **`Symmetry_pade::apply_sym_y(std::complex<double> x_original, std::complex<double> y)`**: Applies symmetry transformation to the function value `y`, potentially based on the transformation of its argument `x`.
    - Takes `std::complex<double> x_original` (the argument before `apply_sym_x` was called) and `std::complex<double> y` (the function value at `x_original`).
    - Returns a `std::complex<double>` representing the transformed function value `y_symm`.
    - It first calls `apply_sym_x(x_original)` to get `x_projected`.
    - If `x_projected == x_original` (meaning `x` was not changed by `apply_sym_x`), then `y_symm = y`.
    - Otherwise, the transformation of `y` depends on `my_symmetry`:
        - `sym_odd`: `y_symm = -y`.
        - `sym_conjugate`: `y_symm = conj(y)`.
        - `sym_anti_conjugate`: `y_symm = -conj(y)`.
        - For other symmetry types (`sym_none`, `sym_mirror_real`, `sym_mirror_imag`, `sym_mirror_both`, `sym_even`), if `x_original` was changed, `y_symm = y` (no change to `y` itself, only `x` was projected).

## Important Variables/Constants

- `my_symmetry`: (Private member of `Symmetry_pade` class, defined in `.hpp`) Stores the integer label for the current symmetry type. The values correspond to enums like `sym_none`, `sym_mirror_real`, etc. (presumably defined in `Symmetry_pade.hpp`).

## Usage Examples

```cpp
// Assuming Symmetry_pade class and enums are defined in Symmetry_pade.hpp
// #include "Symmetry_pade.hpp"
// #include <iostream>
// #include <complex>

// int main() {
    // // Initialize with a symmetry type, e.g., odd symmetry
    // Symmetry_pade pade_sym(Symmetry_pade::sym_odd); // Assuming sym_odd is an enum/const

    // std::complex<double> x_val(-2.0, 3.0);
    // std::complex<double> y_val(1.0, -1.0);

    // // Apply symmetry to x
    // std::complex<double> x_transformed = pade_sym.apply_sym_x(x_val);
    // // For sym_odd and x_val(-2.0, 3.0), x_transformed would be (2.0, -3.0)
    // // std::cout << "Original x: " << x_val << ", Transformed x: " << x_transformed << std::endl;

    // // Apply symmetry to y based on original x and y
    // std::complex<double> y_transformed = pade_sym.apply_sym_y(x_val, y_val);
    // // For sym_odd, if x changed, y_transformed would be -y_val = (-1.0, 1.0)
    // // std::cout << "Original y: " << y_val << ", Transformed y: " << y_transformed << std::endl;

    // // Change symmetry
    // pade_sym.set_symmetry(Symmetry_pade::sym_conjugate); // Assuming sym_conjugate

    // x_transformed = pade_sym.apply_sym_x(x_val);
    // y_transformed = pade_sym.apply_sym_y(x_val, y_val);
    // // For sym_conjugate and x_val(-2.0, 3.0), x_transformed would be (2.0, -3.0)
    // // If x changed, y_transformed would be conj(y_val) = (1.0, 1.0)
    // // std::cout << "Original x: " << x_val << ", Transformed x (conjugate): " << x_transformed << std::endl;
    // // std::cout << "Original y: " << y_val << ", Transformed y (conjugate): " << y_transformed << std::endl;
// }
```

## Dependencies and Interactions

- **`Symmetry_pade.hpp`**: This source file includes `Symmetry_pade.hpp`, which contains the definition of the `Symmetry_pade` class, its member variable `my_symmetry`, and presumably the enum definitions for symmetry types (e.g., `sym_none`, `sym_odd`).
- **`<complex>`**: Uses `std::complex<double>` for representing complex numbers and functions like `std::abs`, `std::copysign`, `std::conj`, `real()`, `imag()`.
- **`<stdexcept>`**: Uses `std::invalid_argument` for error handling in the constructor and `set_symmetry` method.
- **`<cmath>`**: Implicitly used by `std::abs` and `std::copysign` for `double` types.

This module is crucial for ensuring that Padé approximants respect certain known symmetries of the function being approximated. By transforming the input points and/or the function values, it helps in constructing a more accurate and physically meaningful approximant.
The logic for `sym_even`, `sym_odd`, `sym_conjugate`, `sym_anti_conjugate` in `apply_sym_x` appears to be identical: `x_symm = std::complex<double>(std::abs(x.real()), std::copysign(1.0, x.real()) * x.imag());`. This might be intentional if the distinction for these symmetries is primarily handled in `apply_sym_y`.
