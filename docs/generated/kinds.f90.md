# `kinds.f90`

## Overview

The `kinds` Fortran module in the `GX-common` component serves a crucial role in defining standardized kind parameters for data types used throughout the GreenX library. This ensures consistency in numerical precision for real numbers and defines conventional lengths for character strings.

By centralizing these definitions, the library can easily be adapted to different compiler environments or precision requirements if necessary, although `dp` (double precision) is the standard for most scientific computations.

## Key Components

This module defines public integer parameters.

## Important Variables/Constants

- **`sp`**:
    - Type: `integer, parameter, public`
    - Value: `selected_real_kind(6, 30)`
    - Description: Defines the kind parameter for single-precision real numbers. It typically selects a real type that has at least 6 decimal digits of precision and can represent numbers up to approximately 10<sup>30</sup>.

- **`dp`**:
    - Type: `integer, parameter, public`
    - Value: `selected_real_kind(14, 200)`
    - Description: Defines the kind parameter for double-precision real numbers. This is the most commonly used precision for floating-point variables in the library, aiming for at least 14-15 decimal digits of precision and a range up to approximately 10<sup>200</sup>.

- **`medium_char`**:
    - Type: `integer, parameter, public`
    - Value: `100`
    - Description: Defines a standardized "medium" length for character strings. This can be used for declaring character variables that are expected to hold moderately sized strings.

- **`long_char`**:
    - Type: `integer, parameter, public`
    - Value: `200`
    - Description: Defines a standardized "long" length for character strings, suitable for longer informational messages or file paths.

## Usage Examples

```fortran
module example_using_kinds
  use kinds
  implicit none

  ! Declare variables with defined kinds
  real(kind=sp) :: single_precision_variable
  real(kind=dp) :: double_precision_variable
  character(len=medium_char) :: medium_string
  character(len=long_char) :: long_string

  double_precision_variable = 1.234567890123456_dp ! Using _dp suffix for literals
  single_precision_variable = 0.123456_sp         ! Using _sp suffix for literals

  medium_string = "This is a medium length string."
  long_string = "This is a longer string, perhaps for a detailed message or a file path that might exceed 100 characters."

  print *, "Double precision var: ", double_precision_variable
  print *, "Single precision var: ", single_precision_variable
  print *, "Medium string length: ", len(medium_string)
  print *, "Long string length: ", len(long_string)

end module example_using_kinds
```

## Dependencies and Interactions

- This module has no dependencies on other custom modules within the GreenX library. It uses the intrinsic Fortran function `selected_real_kind`.
- Numerous other modules throughout the GreenX library (e.g., `constants`, `minimax_utils`, `minimax_tau`, `minimax_omega`, etc.) depend on `kinds` to define the precision of their real variables (almost universally `dp`) and sometimes character variable lengths.

The `kinds` module is foundational for maintaining numerical stability and data type consistency across the entire codebase.
