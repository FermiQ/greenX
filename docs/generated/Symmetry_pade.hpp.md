# `Symmetry_pade.hpp`

## Overview

This header file defines the `Symmetry_pade` class, which provides utilities to enforce specific symmetries on a Padé model. It declares the class structure, including private member variables for storing the symmetry type, static constants representing different symmetry types, and public methods for initializing, setting, and applying these symmetries.

## Key Components

- **`class Symmetry_pade`**:
    - **Private Members**:
        - `int my_symmetry`: Stores an integer label representing the active symmetry type for an instance of the class.
    - **Static Constants (Symmetry Types)**:
        - `static const int sym_none = 0;`: No symmetry.
        - `static const int sym_mirror_real = 1;`: Mirror symmetry at the imaginary axis (i.e., `f(z) = f(-Re(z) + i Im(z))`, argument becomes `|Re(z)| + i Im(z)` if restricted to positive real part).
        - `static const int sym_mirror_imag = 2;`: Mirror symmetry at the real axis (i.e., `f(z) = f(Re(z) - i Im(z))`, argument becomes `Re(z) + i |Im(z)|` if restricted to positive imag part).
        - `static const int sym_mirror_both = 3;`: Mirror at both real and imaginary axes (argument becomes `|Re(z)| + i |Im(z)|`).
        - `static const int sym_even = 4;`: Even function, `f(z) = f(-z)`.
        - `static const int sym_odd = 5;`: Odd function, `f(z) = -f(-z)`.
        - `static const int sym_conjugate = 6;`: Conjugate symmetry, often defined as `f(z) = conj(f(-z*))` or related forms. The comment states `f(z) = \overline{f(-z)}`.
        - `static const int sym_anti_conjugate = 7;`: Anti-conjugate symmetry, `f(z) = -conj(f(-z*))` or related. The comment states `f(z) = -\overline{f(-z)}`.
    - **Public Methods**:
        - `Symmetry_pade()`: Default constructor.
        - `Symmetry_pade(int symm)`: Constructor to initialize with a specific symmetry type. (Implemented in `Symmetry_pade.cpp`)
        - `void set_symmetry(int symm)`: Method to set or change the symmetry type. (Implemented in `Symmetry_pade.cpp`)
        - `std::complex<double> apply_sym_x(std::complex<double> x)`: Declares a method to apply symmetry transformation to a function argument `x`. (Implemented in `Symmetry_pade.cpp`)
        - `std::complex<double> apply_sym_y(std::complex<double> x_original, std::complex<double> y)`: Declares a method to apply symmetry transformation to a function value `y`, potentially based on the original argument `x_original`. (Implemented in `Symmetry_pade.cpp`)

## Important Variables/Constants

- `my_symmetry`: (Private member) An integer that holds the currently selected symmetry type for an object of `Symmetry_pade`.
- `sym_none`, `sym_mirror_real`, `sym_mirror_imag`, `sym_mirror_both`, `sym_even`, `sym_odd`, `sym_conjugate`, `sym_anti_conjugate`: Static constant integers used as labels or identifiers for the different supported symmetry types. These are used as arguments to the constructor and `set_symmetry` method.

## Usage Examples

```cpp
#include "Symmetry_pade.hpp"
#include <complex>
#include <iostream> // For demonstration

int main() {
    // Default constructor, symmetry might be uninitialized or sym_none by default
    Symmetry_pade sym_handler1;
    // It's good practice to set symmetry explicitly if not done by default constructor
    // sym_handler1.set_symmetry(Symmetry_pade::sym_none); // Assuming sym_none is accessible or use integer 0

    // Constructor with specific symmetry
    // Symmetry_pade sym_handler2(Symmetry_pade::sym_odd); // Or sym_handler2(5);
    Symmetry_pade sym_handler2(5); // Using the integer value for sym_odd

    std::complex<double> point_x(-1.0, 2.0);
    std::complex<double> value_y(3.0, 4.0);

    // Apply symmetry transformation to x
    // std::complex<double> transformed_x = sym_handler2.apply_sym_x(point_x);
    // std::cout << "Original x: " << point_x << ", Transformed x: " << transformed_x << std::endl;

    // Apply symmetry transformation to y
    // std::complex<double> transformed_y = sym_handler2.apply_sym_y(point_x, value_y);
    // std::cout << "Original y: " << value_y << ", Transformed y: " << transformed_y << std::endl;

    // Change symmetry
    // sym_handler2.set_symmetry(Symmetry_pade::sym_conjugate); // Or sym_handler2.set_symmetry(6);
    sym_handler2.set_symmetry(6); // Using the integer value for sym_conjugate

    // transformed_x = sym_handler2.apply_sym_x(point_x);
    // transformed_y = sym_handler2.apply_sym_y(point_x, value_y);
    // std::cout << "After changing to conjugate symmetry:" << std::endl;
    // std::cout << "Original x: " << point_x << ", Transformed x: " << transformed_x << std::endl;
    // std::cout << "Original y: " << value_y << ", Transformed y: " << transformed_y << std::endl;

    return 0;
}
```
*(Note: Accessing static const members like `Symmetry_pade::sym_odd` directly in `main` might require them to be defined outside the class or be public static const. If they are private static const as implied, their integer values would be used directly, or a public accessor/enum would be needed.)*
*Correction: The static const int members are private. This means their direct use like `Symmetry_pade::sym_odd` from outside the class or its friends is not allowed. The integer values (0-7) must be used when calling constructors or `set_symmetry` from external code, or public enums/aliases should be provided.*

## Dependencies and Interactions

- **`<complex>`**: Uses `std::complex<double>` for method arguments and return types.
- **`Symmetry_pade.cpp`**: This header file declares methods that are implemented in `Symmetry_pade.cpp`.
- This class is intended to be used by other parts of the `GX-AnalyticContinuation` system that deal with Padé approximants, allowing them to enforce specific symmetry properties on the functions being approximated.

The definitions of symmetries like `sym_conjugate` and `sym_anti_conjugate` as `f(z) = \overline{f(-z)}` and `f(z) = -\overline{f(-z)}` respectively, imply a transformation on both the argument (negation) and the value (conjugation/negated conjugation). The actual implementation in `Symmetry_pade.cpp` will clarify how these are precisely handled.
