# `minimax_utils.F90`

## Overview

The `minimax_utils` Fortran module provides auxiliary data structures and procedures that support the main minimax routines (`minimax_tau.F90` and `minimax_omega.F90`) within the GX-TimeFrequency module. Its central component is the `er_aw_aux` derived type, designed to store and manage tabulated coefficients and weights used in minimax approximations. It also defines parameters for transformation types.

## Key Components

### Derived Type: `er_aw_aux`

- **Purpose**: To encapsulate the tabulated data for minimax approximations. This includes discrete energy ranges and the corresponding sets of coefficients (often denoted as `a_i` or points) and weights (`w_i`).
- **Members**:
    - `real(kind=dp), dimension(:), allocatable :: energy_range`: A 1D array storing sorted, discrete energy range values (typically `E_max/E_min`).
    - `real(kind=dp), dimension(:, :), allocatable :: aw_erange_matrix`: A 2D array. Each column `k` of this matrix holds the `2*grid_size` values (coefficients `a_i` followed by weights `w_i`) corresponding to the `k`-th energy range in the `energy_range` array.
- **Type-bound Procedure**:
    - `procedure :: get_coeff_weight => coeffs_and_weights`: This procedure is responsible for retrieving the appropriate set of coefficients and weights from the `aw_erange_matrix` based on a given target energy range.

### Public Parameters:

- **`cosine_tw = 1`**: Integer parameter, likely a flag for cosine transform from time to frequency.
- **`cosine_wt = 2`**: Integer parameter, likely a flag for cosine transform from frequency to time.
- **`sine_tw = 3`**: Integer parameter, likely a flag for sine transform from time to frequency.
  (These parameters are used in modules like `minimax_grids.F90`.)

### Subroutines/Functions:

- **`function find_erange(length, einter, eval) result(idx)`** (Private):
    - **Purpose**: Searches an array `einter` (representing `energy_range`) to find the index (`idx`) of the first element that is strictly greater than the input `eval` (the target `e_range`).
    - **Logic**: It iterates through `einter` to find the smallest value in `einter` that is still larger than `eval`.
    - **Return**: Returns the 1-based index `idx`. If `eval` is greater than or equal to all elements in `einter`, it returns `length + 1`.
    - *Note*: The comment mentions `einter` is unsorted and the algorithm is O(n). However, the `energy_range` member of `er_aw_aux` is described as "Sorted array". If sorted, a more efficient search like binary search could be used.

- **`subroutine coeffs_and_weights(this, grid_size, bup, e_range, ac_we, e_ratio)`** (Private, Type-bound to `er_aw_aux` as `get_coeff_weight`):
    - **Purpose**: Selects the appropriate column of coefficients and weights from `this%aw_erange_matrix` based on the input `e_range`. It also calculates a scaling factor `e_ratio`.
    - **Inputs**:
        - `this`: A `class(er_aw_aux)` instance containing the tabulated data.
        - `grid_size`: Integer, the number of grid points.
        - `bup`: Integer, the number of tabulated energy ranges (i.e., `size(this%energy_range)`).
        - `e_range`: Real(dp), the target energy range for which coefficients are desired.
    - **In/Outputs**:
        - `ac_we`: Real(dp) array, dimension `(2*grid_size)`. On output, it's filled with the selected coefficients and weights.
        - `e_ratio`: Real(dp). Initialized to `1.0_dp`. It's modified with a heuristic correction if `e_range` falls below the smallest tabulated range and `grid_size` is large (`>20`). This factor is then used by the calling routine to scale `ac_we`.
    - **Logic**:
        1. Calls `find_erange` to get the index `ien` for the column in `aw_erange_matrix`.
        2. If `ien` indicates that `e_range` is smaller than the smallest tabulated range (`ien == 1`) AND `grid_size > 20`, it computes an `e_ratio = this%energy_range(1) / e_range`. This `e_ratio` is further adjusted by dividing by `1.5_dp` if it's greater than `1.5_dp`.
        3. Copies the data from `this%aw_erange_matrix(:, ien)` to `ac_we`.

## Important Variables/Constants

- **`dp`**: Imported from `kinds`, defines double precision for real numbers.

## Usage Examples

The `er_aw_aux` type and its `get_coeff_weight` method are primarily used internally by `minimax_tau.F90` and `minimax_omega.F90`.

```fortran
! In minimax_tau.F90 or minimax_omega.F90
module some_module
  use kinds, only: dp
  use minimax_utils, only: er_aw_aux
  ! ...
  type(er_aw_aux) :: aw_data
  real(dp), allocatable :: coeffs_and_weights_array(:)
  real(dp) :: target_e_range, scaling_factor
  integer :: n_grid, n_energy_intervals, ierr
  ! ...
  ! (1) Allocate and fill aw_data%energy_range and aw_data%aw_erange_matrix
  !     (Done in set_aw_array_tau or set_aw_array_omega)
  ! ...
  n_grid = 16
  n_energy_intervals = size(aw_data%energy_range)
  target_e_range = 100.0_dp
  scaling_factor = 1.0_dp ! Will be modified by get_coeff_weight

  allocate(coeffs_and_weights_array(2 * n_grid))

  call aw_data%get_coeff_weight(grid_size=n_grid, bup=n_energy_intervals, &
                                e_range=target_e_range, &
                                ac_we=coeffs_and_weights_array, &
                                e_ratio=scaling_factor)
  ! ...
  ! coeffs_and_weights_array now contains data for target_e_range.
  ! It needs to be scaled by scaling_factor by the caller.
  ! ...
end module
```

## Dependencies and Interactions

- **`gx_common.h`**: Likely for common preprocessor definitions.
- **`kinds`**: For `dp` precision.
- **`minimax_tau.F90` / `minimax_omega.F90`**: These modules declare instances of `er_aw_aux` and call its `get_coeff_weight` method after populating the instance with hardcoded data.
- The parameters `cosine_tw`, `cosine_wt`, `sine_tw` are used by `minimax_grids.F90` to control behavior in transformation weight calculations.

This module provides essential utilities that prevent code duplication in `minimax_tau` and `minimax_omega` by centralizing the logic for storing and retrieving pre-tabulated coefficients based on energy ranges.
