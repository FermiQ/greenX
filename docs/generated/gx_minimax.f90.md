# `gx_minimax.f90`

## Overview

The `gx_minimax` Fortran module serves as the primary public Application Programming Interface (API) for the minimax grid generation capabilities within the GX-TimeFrequency library. It consolidates and re-exports key routines from underlying modules, offering a simplified and unified interface for end-users.

The main purpose of this module is to provide functions for computing minimax grids suitable for RPA (Random Phase Approximation) energy calculations and GW calculations on imaginary time and frequency domains.

## Key Components

This module does not define new procedures or types itself but rather makes components from other modules publicly accessible.

### Re-exported Public Subroutines:

- **`gx_minimax_grid(...)`** (from `minimax_grids` module):
    - This is the main routine for generating a comprehensive set of minimax grids (imaginary time and frequency points and weights) and the transformation weights between these domains.
    - For detailed parameters, refer to the documentation for `minimax_grids.F90.md`.

- **`gx_minimax_grid_frequency(...)`** (from `minimax_grids` module):
    - A specialized routine to retrieve only the imaginary frequency grid points and weights.
    - For detailed parameters, refer to the documentation for `minimax_grids.F90.md`.

- **`gx_check_ntau(ntau, msg, ierr)`** (from `api_utilites` module):
    - Checks if a given number of imaginary time points (`ntau`) is supported by the library.
    - For detailed parameters, refer to the documentation for `api_utilities.f90.md`.

- **`gx_get_error_message(msg)`** (from `api_utilites` module):
    - Retrieves the last error message logged by the library if a previous API call resulted in an error.
    - For detailed parameters, refer to the documentation for `api_utilities.f90.md`.

### Re-exported Public Variables:

- **`tau_npoints_supported(:)`** (from `minimax_tau` module):
    - An integer array listing the supported number of grid points for imaginary time grids. Users can inspect this array to determine valid grid sizes.

## Important Variables/Constants

(None defined directly in this module. All relevant constants are part of the imported modules.)

## Usage Examples

The intended usage pattern, as suggested by the module's introductory comment:

```fortran
module my_calculation
  use gx_minimax
  use kinds, only: dp ! Assuming dp is needed for real type declarations
  implicit none

  integer :: num_points, status
  real(dp) :: min_energy, max_energy
  real(dp), allocatable :: tau_p(:), tau_w(:), omega_p(:), omega_w(:)
  real(dp), allocatable :: cosft_wt(:,:), cosft_tw(:,:), sinft_wt(:,:)
  real(dp) :: errs(3), dual_err
  character(len=256) :: error_msg

  ! Define parameters for the grid
  num_points = 16  ! Example: number of time/frequency points
  min_energy = 0.1_dp
  max_energy = 10.0_dp

  ! First, check if the desired number of points is supported
  call gx_check_ntau(num_points, error_msg, status)
  if (status /= 0) then
    print *, "Error: ntau = ", num_points, " is not supported."
    print *, "Details: ", trim(error_msg)
    ! Potentially print available points from tau_npoints_supported
    print *, "Supported ntau values are: ", tau_npoints_supported
    stop 1
  end if

  ! Allocate necessary arrays (gx_minimax_grid will allocate these if not already)
  allocate(tau_p(num_points), tau_w(num_points))
  allocate(omega_p(num_points), omega_w(num_points))
  allocate(cosft_wt(num_points, num_points))
  allocate(cosft_tw(num_points, num_points))
  allocate(sinft_wt(num_points, num_points))

  ! Call the main grid generation routine
  call gx_minimax_grid(num_points, min_energy, max_energy, &
       tau_p, tau_w, omega_p, omega_w, &
       cosft_wt, cosft_tw, sinft_wt, &
       errs, dual_err, status)

  ! Check for errors
  if (status /= 0) then
    call gx_get_error_message(error_msg)
    print *, "Error during gx_minimax_grid execution: ", trim(error_msg)
    ! Handle error appropriately
    stop 2
  end if

  print *, "Minimax grids and weights generated successfully."
  ! ... proceed with using the generated grids and weights ...

  deallocate(tau_p, tau_w, omega_p, omega_w, cosft_wt, cosft_tw, sinft_wt)

end module my_calculation
```

## Dependencies and Interactions

- **`minimax_grids`**: Provides the core grid generation subroutines (`gx_minimax_grid`, `gx_minimax_grid_frequency`).
- **`minimax_tau`**: Provides the `tau_npoints_supported` array, which lists valid grid sizes for imaginary time.
- **`api_utilites`**: Provides helper functions for input validation (`gx_check_ntau`) and error reporting (`gx_get_error_message`).

This module simplifies access to the minimax functionalities by bundling them under a single module. Users primarily interact with the GX-TimeFrequency library's minimax features through this API.
