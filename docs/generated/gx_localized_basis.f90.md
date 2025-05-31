# `gx_localized_basis.f90`

## Overview

This Fortran module, `gx_localized_basis`, serves as the public Application Programming Interface (API) for the localized basis set component of the GreenX library. Its primary purpose is to expose specific functionalities from internal modules to the end-user or other components of the library in a controlled manner. Currently, it re-exports the `gx_rirs_coefficients` subroutine from the `separable_ri` module.

## Key Components

- **Module `gx_localized_basis`**:
    - **Imports**:
        - `gx_rirs_coefficients` from the `separable_ri` module.
    - **Re-exports**:
        - `gx_rirs_coefficients`. By using the `only` keyword, it makes this specific subroutine available to any code that uses `gx_localized_basis`.

## Important Variables/Constants

This module itself does not define any variables or constants. It only acts as a pass-through for the `gx_rirs_coefficients` subroutine.

## Usage Examples

The module's own comment provides a usage example:

```fortran
program use_gx_localized_basis_api
    use gx_localized_basis
    use kinds, only: dp ! Assuming kinds is available for dp
    implicit none

    ! Parameters for gx_rirs_coefficients
    integer :: n_basis_pairs_val = 100 ! Example value
    integer :: n_loc_basbas_val = 200  ! Example value
    integer :: n_rk_points_val = 500   ! Example value

    real(dp), dimension(n_basis_pairs_val, n_rk_points_val) :: ovlp2fn_data
    real(dp), dimension(n_basis_pairs_val, n_loc_basbas_val) :: ovlp3fn_data
    real(dp) :: error_metric

    ! Initialize input data (example)
    ovlp2fn_data = 0.01_dp
    ovlp3fn_data = 0.001_dp

    ! Call the API function
    call gx_rirs_coefficients(n_basis_pairs_val, n_loc_basbas_val, n_rk_points_val, &
                              ovlp2fn_data, ovlp3fn_data, error_metric)

    print *, "gx_rirs_coefficients executed. Error: ", error_metric

end program use_gx_localized_basis_api
```

## Dependencies and Interactions

- **`separable_ri` module**: This is a direct dependency, as `gx_localized_basis` imports `gx_rirs_coefficients` from it.
- Any user code that `use gx_localized_basis` will gain access to `gx_rirs_coefficients` and will therefore indirectly depend on the `separable_ri` module and its own dependencies (like `kinds`, `lapack_interfaces`, `localized_basis_types`, `localized_basis_environments`).

This module simplifies the public interface of the `GX-LocalizedBasis` component by providing a single point of access for selected routines. This approach helps in managing complexity and allows internal module structures to be refactored without necessarily changing the public API, as long as the exposed subroutines maintain their signatures. Currently, the API surface is minimal, exposing only one subroutine.
