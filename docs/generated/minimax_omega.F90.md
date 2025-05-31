# `minimax_omega.F90`

## Overview

The `minimax_omega` Fortran module provides tabulated minimax coefficients for approximating the function `1/x`. This approximation is typically used in the context of frequency-dependent calculations, where `x` would represent frequency (`omega`). The specific form of the approximation used is `1/x ~ (1/pi) * sum_i [w_i * x^2 / (x^2 + a_i^2)^2]`, valid for `x` in a normalized range `[1:rc]`, where `rc` is the energy range factor.

This module stores pre-calculated sets of coefficients (`a_i`, referred to as `points` or `aw_erange_matrix` part 1 in code) and corresponding weights (`w_i`, `aw_erange_matrix` part 2) for various supported grid sizes and energy ranges. These are essential for the frequency integration schemes implemented in the GX-TimeFrequency module, particularly for methods like GW calculations.

The module relies on data published in:
- [https://doi.org/10.1021/ct5001268](https://doi.org/10.1021/ct5001268)
- [https://doi.org/10.1103/PhysRevB.94.165109](https://doi.org/10.1103/PhysRevB.94.165109)

## Key Components

### Public Subroutines:

- **`get_points_weights_omega(grid_size, e_range, ac_we, ierr)`**:
    - This is the primary public interface of the module.
    - It retrieves the appropriate minimax coefficients (`a_i`) and weights (`w_i`) for a given `grid_size` (number of frequency points) and `e_range` (ratio of max to min energy, `E_max/E_min`).
    - The retrieved coefficients and weights are scaled by `1/e_ratio` (where `e_ratio` is related to `e_min`) before being returned in the `ac_we` array.
    - **Inputs**:
        - `grid_size`: Integer, number of mesh points.
        - `e_range`: Real(dp), the energy range `E_max/E_min`.
    - **Outputs**:
        - `ac_we`: Real(dp) array of dimension `(2 * grid_size)`, storing coefficients (`a_i`) in the first half and weights (`w_i`) in the second half.
        - `ierr`: Integer, error status (0 for success, 1 if `grid_size` is not supported).

### Private Subroutines:

- **`set_aw_array_omega(grid_size, aw)`**:
    - Populates the `aw` derived type (of `er_aw_aux` type) with hardcoded, pre-tabulated minimax coefficients and weights.
    - It uses a large `select case(grid_size)` block to load the specific data for the requested `grid_size`.
    - For each `grid_size`, it defines:
        - `aw%energy_range(:)`: An array of discrete energy range values for which coefficients are tabulated.
        - `aw%aw_erange_matrix(:, k)`: A 2D array where the k-th column contains the `2*grid_size` values (coefficients `a_i` followed by weights `w_i`) corresponding to the k-th energy range in `aw%energy_range`.
    - **Inputs**:
        - `grid_size`: Integer, the requested grid size.
    - **In/Outputs**:
        - `aw`: `type(er_aw_aux)`, derived type instance to be filled with coefficients.

## Important Variables/Constants

- **`dp`**: Imported from `kinds`, defines the double precision kind for real numbers.
- **`ngrids`**: Integer parameter, defines the number of different grid sizes for which coefficients are stored.
- **`omega_npoints_supported(ngrids)`**: Public integer parameter array. Stores the list of supported grid sizes (e.g., 6, 8, 10, ..., 34).
- **`energy_ranges_grids(ngrids)`**: Integer parameter array. For each corresponding grid size in `omega_npoints_supported`, this array stores the number of discrete energy range values for which coefficients are tabulated.
- **`er_aw_aux`**: A derived type imported from `minimax_utils`. This type is used to store and manage the energy ranges, coefficients, and weights. It likely contains methods like `get_coeff_weight` to interpolate or select the correct set of `a_i` and `w_i` based on the input `e_range`.

## Usage Examples

The `get_points_weights_omega` subroutine is typically called by other routines (like those in `minimax_grids.F90`) that need these coefficients for constructing frequency grids or performing numerical integrations.

```fortran
module use_minimax_omega_example
  use minimax_omega
  use kinds, only: dp
  implicit none

  integer, parameter :: num_freq_points = 16
  real(dp) :: energy_ratio_val = 100.0_dp ! Example E_max/E_min
  real(dp), allocatable :: omega_coeffs_weights(:)
  integer :: err_status

  allocate(omega_coeffs_weights(2 * num_freq_points))

  call get_points_weights_omega(grid_size=num_freq_points, &
                                e_range=energy_ratio_val, &
                                ac_we=omega_coeffs_weights, &
                                ierr=err_status)

  if (err_status == 0) then
    print *, "Minimax omega coefficients and weights retrieved successfully."
    ! omega_coeffs_weights(1:num_freq_points) contains the 'a_i' coefficients
    ! omega_coeffs_weights(num_freq_points+1 : 2*num_freq_points) contains the 'w_i' weights
    ! These are scaled by 1/e_min (done via e_ratio internally)
  else
    print *, "Error retrieving minimax omega data: ", err_status
  end if

  deallocate(omega_coeffs_weights)

end module use_minimax_omega_example
```

## Dependencies and Interactions

- **`gx_common.h`**: Included for common definitions, potentially error handling macros like `_REGISTER_EXC`.
- **`kinds`**: Provides type definitions (e.g., `dp` for double precision).
- **`error_handling`**: Used for registering exceptions (e.g., `register_exc` macro call).
- **`minimax_utils`**: Provides the `er_aw_aux` derived type which is crucial for storing and accessing the tabulated coefficients and weights. The `er_aw_aux` type likely has a method (e.g., `get_coeff_weight`) that handles the logic of finding the closest tabulated `e_range` and returning the corresponding `a_i` and `w_i` values, possibly with some scaling or interpolation if the exact `e_range` is not tabulated.

This module serves as a data repository for the minimax frequency scheme. The actual mathematical details of how these coefficients are used would be found in modules that consume this data, such as `minimax_grids.F90`. The hardcoded nature of the coefficients in `set_aw_array_omega` means that to support new grid sizes or refine existing coefficients, this module's source code must be updated directly.
The module is designed to be extensible by adding new entries to `omega_npoints_supported`, `energy_ranges_grids`, and then adding a corresponding `case` in `set_aw_array_omega`.
