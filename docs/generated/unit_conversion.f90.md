# `unit_conversion.f90`

## Overview

The `unit_conversion` Fortran module, part of the `GX-common` component, is dedicated to providing standardized and accurate conversion factors between various units of measurement. The module stresses the importance of using values from CODATA 2018 recommendations to ensure the reliability of physical calculations within the GreenX library.

Guidelines for constants in this module include:
- Adherence to CODATA 2018 values.
- Definition using double precision (`dp`).
- Use of descriptive names for clarity.

## Key Components

This module directly defines parameters for unit conversion.

## Important Variables/Constants

- **`bohr_to_angstrom`**:
    - Type: `real(kind=dp), parameter, public`
    - Value: `0.529177210903_dp`
    - Description: The conversion factor to convert lengths from Bohr (atomic unit of length, a<sub>0</sub>) to Angstroms (Å). 1 Bohr = 0.529177210903 Å.

## Usage Examples

```fortran
module example_using_unit_conversion
  use unit_conversion
  use kinds, only: dp ! Assuming dp is from a module like 'kinds'
  implicit none

  real(kind=dp) :: length_bohr, length_angstrom

  length_bohr = 2.0_dp ! Example length in Bohr

  ! Convert Bohr to Angstroms
  length_angstrom = length_bohr * bohr_to_angstrom

  print *, length_bohr, " Bohr is equal to ", length_angstrom, " Angstroms."

end module example_using_unit_conversion
```

## Dependencies and Interactions

- **`kinds`**: The `unit_conversion` module imports the `dp` kind parameter from the `kinds` module to define its conversion factors in double precision.
- This module would be used by any part of the GreenX library that needs to convert values between Bohr and Angstroms, or other units if more conversion factors were added.

This module plays a role in maintaining consistency and accuracy in calculations involving physical quantities that might be expressed in different common units.
