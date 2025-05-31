# `ComplexGMP.hpp`

## Overview

This header file defines the `ComplexGMP` class, which represents complex numbers using arbitrary precision arithmetic provided by the GMP (GNU Multiple Precision Arithmetic Library). It includes constructors for various initialization scenarios and declares overloaded operators for basic arithmetic operations (addition, subtraction, multiplication, division) and a method to calculate the absolute value. It also provides a type casting operator to convert `ComplexGMP` objects to standard `std::complex<double>`.

## Key Components

- **`class ComplexGMP`**:
    - **Public Members**:
        - `mpf_class real_mp`: Stores the real part of the complex number using GMP's arbitrary precision float.
        - `mpf_class imag_mp`: Stores the imaginary part of the complex number using GMP's arbitrary precision float.
    - **Constructors**:
        - `ComplexGMP()`: Default constructor, initializes to `0 + 0i`.
        - `ComplexGMP(const mpf_class &real_mp, const mpf_class &imag_mp)`: Initializes from arbitrary precision real and imaginary parts.
        - `ComplexGMP(std::complex<double> z)`: Initializes from a standard `std::complex<double>`.
    - **Operator Overloads**:
        - `operator+`: Adds two `ComplexGMP` numbers. (Implemented inline)
        - `operator-`: Subtracts one `ComplexGMP` number from another. (Implemented inline)
        - `operator*`: Declares multiplication of two `ComplexGMP` numbers (implemented in `ComplexGMP.cpp`).
        - `operator/`: Declares division of two `ComplexGMP` numbers (implemented in `ComplexGMP.cpp`).
    - **Type Casting**:
        - `operator std::complex<double>() const`: Converts a `ComplexGMP` object to `std::complex<double>`. (Implemented inline)
    - **Methods**:
        - `mpf_class abs() const`: Declares the calculation of the absolute value (implemented in `ComplexGMP.cpp`).

## Important Variables/Constants

- `real_mp`: Member variable of `ComplexGMP` storing the real part.
- `imag_mp`: Member variable of `ComplexGMP` storing the imaginary part.

## Usage Examples

```cpp
#include "ComplexGMP.hpp"
#include <iostream> // For printing

int main() {
    // Default constructor
    ComplexGMP c1; // 0 + 0i

    // Constructor from mpf_class
    mpf_class r_mp("1.234567890123456789");
    mpf_class i_mp("9.876543210987654321");
    ComplexGMP c2(r_mp, i_mp);

    // Constructor from std::complex<double>
    std::complex<double> std_c(1.5, -2.5);
    ComplexGMP c3(std_c);

    // Addition
    ComplexGMP sum = c2 + c3;
    // std::cout << "Sum: " << static_cast<std::complex<double>>(sum) << std::endl;

    // Subtraction
    ComplexGMP diff = c2 - c3;
    // std::cout << "Difference: " << static_cast<std::complex<double>>(diff) << std::endl;

    // Multiplication (implementation in ComplexGMP.cpp)
    // ComplexGMP prod = c2 * c3;

    // Division (implementation in ComplexGMP.cpp)
    // ComplexGMP quot = c2 / c3;

    // Absolute value (implementation in ComplexGMP.cpp)
    // mpf_class abs_val = c2.abs();

    // Type casting to std::complex<double>
    std::complex<double> c2_std = static_cast<std::complex<double>>(c2);
    // std::cout << "c2 as std::complex<double>: " << c2_std << std::endl;

    return 0;
}
```

## Dependencies and Interactions

- **`<gmpxx.h>`**: Essential for `mpf_class`, providing the arbitrary precision floating-point capabilities.
- **`<complex>`**: Used for the standard `std::complex<double>` type, particularly in one of the constructors and the type casting operator.
- **`ComplexGMP.cpp`**: This header file declares methods (`operator*`, `operator/`, `abs`) that are implemented in `ComplexGMP.cpp`.
- The class is designed to be a high-precision alternative to `std::complex<double>` and would be used in parts of the `GX-AnalyticContinuation` system that require calculations with higher precision than standard double-precision floating-point numbers can offer.
