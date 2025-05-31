# `minimax_tau.F90`

## Overview

The `minimax_tau` Fortran module is a data module containing tabulated minimax coefficients. These coefficients are used to approximate the function `1/x` with a sum of exponentials, typically in the form `1/x ~ sum_i [w_i * exp(-a_i * x)]`. This type of approximation is commonly employed in quantum chemistry and physics for calculations involving imaginary time (`tau`), where `x` would represent `tau` scaled by some energy.

This module stores pre-calculated sets of coefficients (`a_i`, also referred to as points) and corresponding weights (`w_i`) for various supported grid sizes and energy ranges. These are crucial for the imaginary time integration schemes within the GX-TimeFrequency module.

The underlying data and methodology are based on the works referenced in:
- [https://doi.org/10.1021/ct5001268](https://doi.org/10.1021/ct5001268)
- [https://doi.org/10.1103/PhysRevB.94.165109](https://doi.org/10.1103/PhysRevB.94.165109)

## Key Components

### Public Subroutines:

- **`get_points_weights_tau(grid_size, e_range, ac_we, ierr)`**:
    - This is the primary public routine for accessing the minimax coefficients and weights.
    - It retrieves the appropriate `a_i` (coefficients/points) and `w_i` (weights) for a specified `grid_size` (number of time points) and `e_range` (ratio `E_max/E_min`).
    - The retrieved coefficients and weights are scaled by `e_ratio` (where `e_ratio` is related to `e_min`) before being returned. Note the difference in scaling factor application compared to `minimax_omega` (multiplied here, divided in omega).
    - **Inputs**:
        - `grid_size`: Integer, the number of mesh points.
        - `e_range`: Real(dp), the energy range `E_max/E_min`.
    - **Outputs**:
        - `ac_we`: Real(dp) array of dimension `(2 * grid_size)`. The first half stores coefficients (`a_i`), and the second half stores weights (`w_i`).
        - `ierr`: Integer, error status (0 for success, 1 if `grid_size` is not supported).

### Private Subroutines:

- **`set_aw_array_tau(grid_size, aw)`**:
    - This subroutine populates the `aw` derived type (instance of `er_aw_aux`) with the hardcoded, pre-tabulated minimax coefficients and weights for the imaginary time approximation.
    - A large `select case(grid_size)` structure is used to load the specific data sets corresponding to the requested `grid_size`.
    - For each `grid_size`, it defines:
        - `aw%energy_range(:)`: An array of discrete energy range values for which coefficients are tabulated.
        - `aw%aw_erange_matrix(:, k)`: A 2D array. The k-th column contains `2*grid_size` values (coefficients `a_i` then weights `w_i`) for the k-th energy range specified in `aw%energy_range`.
    - **Inputs**:
        - `grid_size`: Integer, the grid size for which coefficients are needed.
    - **In/Outputs**:
        - `aw`: `type(er_aw_aux)`, the derived type instance to be filled.

## Important Variables/Constants

- **`dp`**: Imported from `kinds`, specifies the double precision kind for real numbers.
- **`ngrids`**: Integer parameter, indicating the number of different grid sizes supported.
- **`tau_npoints_supported(ngrids)`**: Public integer parameter array. Lists the supported grid sizes (e.g., 6, 8, ..., 34), identical to `omega_npoints_supported`.
- **`energy_ranges_grids(ngrids)`**: Integer parameter array. For each grid size, this stores the number of discrete `e_range` values for which coefficients have been tabulated.
- **`er_aw_aux`**: A derived type imported from `minimax_utils`. This type is essential for organizing and accessing the tabulated data (energy ranges, coefficients, and weights). It likely includes a method (e.g., `get_coeff_weight`) to select or interpolate coefficients based on the input `e_range`.

## Usage Examples

The `get_points_weights_tau` subroutine is primarily used internally by other parts of the GX-TimeFrequency module, for example, by `minimax_grids.F90` to set up imaginary time grids.

```fortran
module use_minimax_tau_example
  use minimax_tau
  use kinds, only: dp
  implicit none

  integer, parameter :: num_time_points = 16
  real(dp) :: energy_ratio_val = 100.0_dp ! Example E_max/E_min
  real(dp), allocatable :: tau_coeffs_weights(:)
  integer :: err_status

  allocate(tau_coeffs_weights(2 * num_time_points))

  call get_points_weights_tau(grid_size=num_time_points, &
                              e_range=energy_ratio_val, &
                              ac_we=tau_coeffs_weights, &
                              ierr=err_status)

  if (err_status == 0) then
    print *, "Minimax tau coefficients and weights retrieved successfully."
    ! tau_coeffs_weights(1:num_time_points) contains the 'a_i' coefficients
    ! tau_coeffs_weights(num_time_points+1 : 2*num_time_points) contains the 'w_i' weights
    ! These are scaled by e_min (done via e_ratio internally)
  else
    print *, "Error retrieving minimax tau data: ", err_status
  end if

  deallocate(tau_coeffs_weights)

end module use_minimax_tau_example
```

## Dependencies and Interactions

- **`gx_common.h`**: Included for common definitions, potentially error handling macros (`_REGISTER_EXC`).
- **`kinds`**: Provides the `dp` type definition.
- **`error_handling`**: Used for the `register_exc` macro for error reporting.
- **`minimax_utils`**: Supplies the `er_aw_aux` derived type, which is fundamental for managing the tabulated coefficient data. The `er_aw_aux%get_coeff_weight` method is used to extract the relevant data for a given `e_range`.

Similar to `minimax_omega.F90`, this module acts as a data store. The coefficients are hardcoded into `set_aw_array_tau`. Any updates or additions to supported grids/coefficients require direct modification of this source file. The module is designed for extensibility by updating `tau_npoints_supported`, `energy_ranges_grids`, and adding new `case` blocks in `set_aw_array_tau`.
The scaling of `ac_we` by `e_ratio` (effectively `e_min`) in `get_points_weights_tau` adapts the normalized coefficients to the specific energy scale of the problem.
