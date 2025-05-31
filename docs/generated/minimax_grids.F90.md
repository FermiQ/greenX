# `minimax_grids.F90`

## Overview

The `minimax_grids` Fortran module is designed to calculate optimized frequency and imaginary time grids, along with corresponding integration points and weights. These grids are essential for numerical correlation methods, particularly in the context of Green's function approaches (likely GW calculations, given the naming convention). The module also computes weights for inhomogeneous cosine and sine transforms between the time and frequency domains. The methods employed are based on minimax principles to ensure accuracy, as referenced by scientific publications.

## Key Components

The module is structured around a few main public subroutines and several private helper routines.

### Public Subroutines:

- **`gx_minimax_grid(...)`**:
    - This is the main entry point for generating a comprehensive set of minimax grids and transformation weights.
    - **Inputs**: `num_points` (number of mesh points), `e_min` (minimum transition energy), `e_max` (maximum transition energy), optional `bare_cos_sin_weights`, optional `regularization`.
    - **Outputs**: `tau_points`, `tau_weights` (imaginary time grid and weights), `omega_points`, `omega_weights` (imaginary frequency grid and weights), `cosft_wt` (weights for tau -> omega cosine transform), `cosft_tw` (weights for omega -> tau cosine transform), `sinft_wt` (weights for tau -> omega sine transform), `max_errors` (max error for transforms), `cosft_duality_error` (error for `cosft_wt * cosft_tw - I`), `ierr` (exit status).
- **`gx_minimax_grid_frequency(...)`**:
    - A more specialized routine to retrieve only the imaginary frequency grid and weights.
    - **Inputs**: `num_points`, `e_min`, `e_max`.
    - **Outputs**: `omega_points`, `omega_weights`, `ierr`.

### Private Helper Subroutines:

- **`get_transformation_weights(...)`**:
    - Calculates the actual transformation weights for cosine (time to frequency, frequency to time) and sine (time to frequency) transforms using Singular Value Decomposition (SVD via LAPACK's `dgesdd`).
    - It involves constructing an auxiliary matrix `mat_A` and a vector `psi` based on the transformation type and then solving for the weights.
    - **Inputs**: `num_points`, `tau_points`, `omega_points`, `e_min`, `e_max`, `transformation_type`, `regularization`.
    - **Outputs**: `weights`, `max_error`, `ierr`.
- **`calculate_psi_and_mat_A(...)`**:
    - Constructs the `psi` vector and `mat_A` matrix required by `get_transformation_weights`. The exact forms depend on the `transformation_type`.
    - **Inputs**: `num_points`, `tau_points`, `omega_points`, `num_x_nodes`, `x_mu` (logarithmic grid for auxiliary calculations), `i_point`, `transformation_type`.
    - **Outputs**: `psi`, `mat_A`, `current_point`.
- **`calculate_max_error(...)`**:
    - Estimates the maximum error of the fitting procedure used to obtain the transformation weights by comparing the fitted function to the target function `psi` on a dense grid `x_mu`.
    - **Inputs**: `num_points`, `tau_points`, `omega_points`, `weights_work`, `num_x_nodes`, `x_mu`, `psi`, `current_point`, `transformation_type`.
    - **Outputs**: `max_error`.

## Important Variables/Constants

- **`dp`**: Imported from `kinds`, likely defines the double precision kind for real numbers.
- **`pi`**: Imported from `constants`.
- **`cos_t_to_cos_w`, `cos_w_to_cos_t`, `sin_t_to_sin_w`**: Integer parameters used as flags to specify the type of transformation in `get_transformation_weights` and related routines.
- **`e_min`, `e_max`**: Input parameters defining the energy range for grid generation. The ratio `e_max/e_min` is crucial.
- **`num_points`**: The number of desired grid points for both time and frequency meshes.
- **`regularization`**: An optional Tikhonov regularization parameter used in the SVD step of `get_transformation_weights` to improve numerical stability.
- **`bare_cos_sin_weights`**: An optional logical flag. If true, the output `cosft_wt`, `cosft_tw`, `sinft_wt` are "bare" transformation weights. If false (default), these weights are multiplied by the `cos(tau*omega)` or `sin(tau*omega)` terms respectively.

## Usage Examples

```fortran
module use_minimax_grids_example
  use minimax_grids
  use kinds, only: dp
  implicit none

  integer, parameter :: n_pts = 64
  real(dp) :: min_energy = 0.1_dp
  real(dp) :: max_energy = 10.0_dp
  real(dp), allocatable :: tau_p(:), tau_w(:), omega_p(:), omega_w(:)
  real(dp), allocatable :: c_wt(:,:), c_tw(:,:), s_wt(:,:)
  real(dp) :: errs(3), dual_err
  integer :: status

  ! Allocate arrays (simplified, gx_minimax_grid handles allocation for outputs)
  allocate(tau_p(n_pts), tau_w(n_pts), omega_p(n_pts), omega_w(n_pts))
  allocate(c_wt(n_pts, n_pts), c_tw(n_pts, n_pts), s_wt(n_pts, n_pts))

  call gx_minimax_grid(num_points=n_pts, e_min=min_energy, e_max=max_energy, &
       tau_points=tau_p, tau_weights=tau_w, &
       omega_points=omega_p, omega_weights=omega_w, &
       cosft_wt=c_wt, cosft_tw=c_tw, sinft_wt=s_wt, &
       max_errors=errs, cosft_duality_error=dual_err, ierr=status)

  if (status == 0) then
    print *, "Minimax grids and weights generated successfully."
    print *, "Max errors (cos_wt, cos_tw, sin_wt): ", errs
    print *, "Duality error for cosine transforms: ", dual_err
    ! Further use of tau_p, omega_p, c_wt etc.
  else
    print *, "Error generating minimax grids: ", status
  end if

  ! Example for frequency grid only
  call gx_minimax_grid_frequency(num_points=n_pts, e_min=min_energy, e_max=max_energy, &
                                 omega_points=omega_p, omega_weights=omega_w, ierr=status)
  if (status == 0) then
      print *, "Minimax frequency grid generated successfully."
  else
      print *, "Error generating frequency grid: ", status
  end if

  deallocate(tau_p, tau_w, omega_p, omega_w, c_wt, c_tw, s_wt)

end module use_minimax_grids_example
```

## Dependencies and Interactions

- **`gx_common.h`**: Included for common definitions, potentially error handling macros like `_REGISTER_EXC`.
- **`kinds`**: Provides type definitions (e.g., `dp` for double precision).
- **`error_handling`**: Used for registering exceptions (e.g., `register_exc`).
- **`constants`**: Provides physical or mathematical constants (e.g., `pi`).
- **`minimax_tau`**: Provides `get_points_weights_tau` for generating the initial time grid points and weights based on `e_range`.
- **`minimax_omega`**: Provides `get_points_weights_omega` for generating the initial frequency grid points and weights based on `e_range`.
- **`minimax_utils`**: (Assumed, based on typical structure, though not directly listed as `use`d for `cosine_wt`, `cosine_tw`, `sine_tw` which appear to be internal parameter constants in this file with names like `cos_t_to_cos_w`). The comment mentions `cosine_wt, cosine_tw, sine_tw` from `minimax_utils` but these are defined as local integer parameters. The module likely provides other utility functions used in the broader context.
- **`lapack_interfaces`**: Provides LAPACK routines, specifically `dgemm` (matrix multiplication) and `dgesdd` (singular value decomposition), which are critical for calculating transformation weights.

The module follows a pattern of first getting "raw" points and weights from `minimax_tau` and `minimax_omega` (which likely handle the core minimax optimization for the grid point distribution), then scaling them to the desired `[e_min, e_max]` range. Subsequently, it computes the transformation matrices between these grids using SVD for robustness. Error handling is done by setting an `ierr` flag and returning. Allocatable arrays are used for outputs, and they are expected to be deallocated by the caller or automatically when they go out of scope (F2008 standard).
The module relies on external libraries (LAPACK) for numerically intensive operations.
The references in the initial comment block point to the theoretical background for these grid generation and transformation techniques.
- [https://doi.org/10.1021/ct5001268](https://doi.org/10.1021/ct5001268)
- [https://doi.org/10.1103/PhysRevB.94.165109](https://doi.org/10.1103/PhysRevB.94.165109)
