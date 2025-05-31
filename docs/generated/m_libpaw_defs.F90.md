# `m_libpaw_defs.F90`

## Overview

The Fortran module `m_libpaw_defs` serves as a central repository for various fundamental definitions used within the `libPAW` library, particularly when `libPAW` is not being used with ABINIT (which would use its own `defs_basis` module). This module provides:

-   Definitions for double precision real and complex kinds.
-   A comprehensive list of named mathematical constants (e.g., `zero`, `one`, `pi`).
-   Various tolerance values.
-   Complex constants (`czero`, `cone`).
-   Character constants and string length parameters.
-   Standard UNIX I/O unit numbers.
-   Key physical constants relevant to atomic and electronic structure calculations.
-   Simple derived data types for creating "ragged" arrays (arrays where one or more dimensions are determined at runtime).
-   Small helper functions (`to_arrayN`) likely used by preprocessor macros for memory allocation, especially in the BigDFT context.

The module notes indicate that this file is derived from ABINIT's `defs_basis.F90` module, ensuring a consistent set of basic definitions when `libPAW` is used independently or with codes other than ABINIT.

## Key Components

### Named Constants

-   **Precision Kinds**:
    -   `dp`: Integer parameter, `kind(1.0d0)` for double precision real numbers.
    -   `dpc`: Integer parameter, `kind((1.0_dp,1.0_dp))` for double precision complex numbers.

-   **Real Constants (`real(dp), parameter`)**:
    -   Basic numbers: `zero`, `one`, `two`, `three`, `four`.
    -   Fractions: `half`, `third`, `quarter`, `eighth`.
    -   Square roots: `sqrt2`, `sqrt3`, `sqrthalf`.
    -   Pi: `pi`, `two_pi`, `four_pi`.
    -   Tolerances: `tol3` (0.001) down to `tol16` (1.0e-16).

-   **Complex Constants (`complex(dpc), parameter`)**:
    -   `czero`: Complex zero `(0.0_dp, 0.0_dp)`.
    -   `cone`: Complex one `(1.0_dp, 0.0_dp)`.

-   **Character Constants**:
    -   `ch10`: `char(10)` (carriage return/newline).
    -   `fnlen`: Integer parameter, `264` (max length for file names).
    -   `strlen`: Integer parameter, `2000000` (max length for input strings).

-   **UNIX Unit Numbers (`integer, parameter` or `save`)**:
    -   `ab_out`: `7` (output file, ABINIT convention).
    -   `std_out`: `6` (standard output).
    -   `std_err`: `0` (standard error).
    -   `tmp_unit`, `tmp_unit2`: `9`, `10` (temporary files).

-   **Physical Constants (`real(dp), parameter`)**:
    -   `Bohr_Ang`: `0.52917720859_dp` (Bohr radius in Angstroms).
    -   `Ha_eV`: `27.21138386_dp` (Hartree in electronvolts).
    -   `InvFineStruct`: `137.035999679_dp` (Inverse fine structure constant).
    -   `FineStructureConstant2`: `0.000053251354478_dp` (Square of fine structure constant).

### Derived Data Types for Ragged Arrays

These types are simple wrappers around allocatable arrays, allowing them to be components of other derived types, effectively creating ragged arrays or arrays of arrays.

-   **`coeffi1_type`**: Contains `integer, allocatable :: value(:)`.
-   **`coeff1_type`**: Contains `real(dp), allocatable :: value(:)`.
-   **`coeff2_type`**: Contains `real(dp), allocatable :: value(:,:)`.
-   **`coeff3_type`**: Contains `real(dp), allocatable :: value(:,:,:)`.

### Helper Functions (`to_arrayN`)

-   **`interface to_array`**: Generic interface for `to_array1` through `to_array6`.
-   **`function to_array1(i1) result(arr)`**: Takes 1 integer, returns `arr(1) = (/i1/)`.
-   **`function to_array2(i1,i2) result(arr)`**: Takes 2 integers, returns `arr(2) = (/i1,i2/)`.
-   ...and so on up to:
-   **`function to_array6(i1,i2,i3,i4,i5,i6) result(arr)`**: Takes 6 integers, returns `arr(6) = (/i1,i2,i3,i4,i5,i6/)`.
    These functions are likely used by the BigDFT memory allocation macros (`f_malloc(to_array SIZE)`) seen in `libpaw.h` to convert a list of dimension sizes into an array suitable for the `f_malloc` function.

## Important Variables/Constants

The module primarily *defines* constants rather than declaring variables. The parameters listed above (e.g., `dp`, `pi`, `Bohr_Ang`) are the key entities provided by this module for use throughout `libPAW`.

## Usage Examples

Other Fortran modules within `libPAW` would `use m_libpaw_defs` to access these definitions:

```fortran
module m_another_paw_module
    use m_libpaw_defs, only: dp, pi, Bohr_Ang, coeff1_type, zero
    implicit none

    real(dp) :: radius_bohr, radius_angstrom
    type(coeff1_type) :: some_coefficients

    subroutine calculate_something
        radius_bohr = 5.0_dp
        radius_angstrom = radius_bohr * Bohr_Ang

        if (abs(radius_angstrom) < tol8) then
            ! Do something if radius is very small
        end if

        allocate(some_coefficients%value(10))
        some_coefficients%value = zero ! Initialize with zero from m_libpaw_defs
        ! ...
        if (allocated(some_coefficients%value)) deallocate(some_coefficients%value)
    end subroutine calculate_something

end module m_another_paw_module
```

The `to_array` functions are used more implicitly through macros when `libPAW` is configured for BigDFT, for example:
`LIBPAW_ALLOCATE(ARR, (dim1, dim2))` might expand to `ARR = f_malloc(to_array2(dim1, dim2))` after C preprocessing.

## Dependencies and Interactions

-   This module has no dependencies on other `libPAW` modules, making it a foundational definitions module.
-   Many other modules in `libPAW` (those not specific to ABINIT's environment) will `use m_libpaw_defs` to ensure consistent precision and access to fundamental constants and types.
-   It plays a role in the conditional compilation set up by `libpaw.h`:
    -   If `HAVE_LIBPAW_ABINIT` is defined, ABINIT's `defs_basis` is used instead.
    -   If `HAVE_LIBPAW_BIGDFT` is defined or for the default configuration, `m_libpaw_defs` is specified by the `USE_DEFS` macro.

This module is essential for maintaining consistency and providing fundamental building blocks for the `libPAW` library, especially in diverse compilation environments.
