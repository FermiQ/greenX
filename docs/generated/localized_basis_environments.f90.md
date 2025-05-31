# `localized_basis_environments.f90`

## Overview

This Fortran module, `localized_basis_environments`, provides subroutines for initializing, deallocating, and managing various data structures (environments) related to localized basis sets. These environments appear to be used in calculations involving Resolution of Identity (RI), polarizability, and an entity referred to as "w_engine". The module also includes utility functions for error calculation, fetching machine precision, and matrix operations like diagonalization and matrix power, primarily using LAPACK routines.

## Key Components

- **Module `localized_basis_environments`**:
    - Imports types (`separable_ri_types`, `polarizability_types`, `w_engine_types`) from `localized_basis_types`.
    - Imports LAPACK routines (`dgemm`, `dsyevx`, `dlaset`, `dlamch`).
    - Imports `kinds` for `dp` and `error_handling` for `register_exc`.

- **Initialization Subroutines**:
    - **`initialization(ri_rs)`**:
        - Initializes a `separable_ri_types` variable (`ri_rs`).
        - Sets hardcoded dimensions for basis sets and points.
        - Allocates arrays within `ri_rs` (e.g., `ovlp2fn`, `ovlp3fn`, `z_coeff`).
    - **`initialize_kohn_sham(pi_pq)`**:
        - Initializes Kohn-Sham related arrays within a `polarizability_types` variable (`pi_pq`).
        - Allocates `eigenvalues`, `eigenvectors`, and `wave` arrays in `pi_pq%ks`.
    - **`initialize_minimax_grids(pi_pq)`**:
        - Initializes minimax time-frequency grid arrays within `pi_pq%minimax`.
        - Sets a hardcoded number of points (`n_points = 6`).
        - Allocates and initializes `cos_tf`, `omega`, `tau`, `weights`.
    - **`initialize_polarizability(pi_pq)`**:
        - Initializes polarizability-specific arrays within `pi_pq`.
        - Allocates `tau`, `omega`, and `chi%matrix`. Initializes them to zero.
    - **`initialize_w_engine(we)`**:
        - Initializes arrays within a `w_engine_types` variable (`we`).
        - Allocates `omega` and `work`. Initializes `omega` to zero and `work` using `dlaset`.

- **Deallocation Subroutines**:
    - **`deallocations(ri_rs, keep_coeff)`**:
        - Deallocates arrays in `ri_rs`.
        - Optionally keeps `ri_rs%z_coeff` allocated based on `keep_coeff` logical.
    - **`deallocate_minimax_grids(pi_pq)`**: Deallocates minimax grid arrays in `pi_pq`.
    - **`deallocate_kohn_sham(pi_pq)`**: Deallocates Kohn-Sham arrays in `pi_pq`.
    - **`deallocate_polarizability(pi_pq, keep_pi)`**:
        - Deallocates arrays in `pi_pq`.
        - Calls `deallocate_minimax_grids` and `deallocate_kohn_sham`.
        - Optionally keeps `pi_pq%omega` allocated based on `keep_pi`.
    - **`deallocate_w_engine(we)`**: Deallocates arrays in `we`.

- **Utility Subroutines**:
    - **`calculate_error(ri_rs, n_basis_pairs, n_basbas, ovlp3fn, error)`**:
        - Calculates the maximum absolute difference between `ri_rs%ovlp3fn` (presumably RI-RS coefficients) and a provided `ovlp3fn` array (presumably RI-V coefficients).
    - **`get_machine_precision(safe_minimum)`**:
        - Retrieves the machine safe minimum (epsilon) using LAPACK's `dlamch('S')`.
    - **`power_genmat(matrix, n_dim, power, threshold)`**:
        - Computes `matrix = matrix ** power` for a symmetric matrix.
        - It first negates the matrix, diagonalizes it, negates eigenvalues back.
        - Checks for unphysical eigenvalues (if `n_nonsingular /= n_dim` after diagonalization of negated matrix).
        - Filters eigenvalues based on a `threshold`.
        - Reconstructs the matrix using `(V * (diag(evals)**(power/2))) * (diag(evals)**(power/2)) * V^T` effectively, but by modifying eigenvectors `work(j,i) = work(j,i) * sqrt(eigenvalues(i))**power` and then `matrix = work * work^T`.
    - **`diagonalize_genmat(matrix, n_dim, threshold, n_nonsingular, eigenvalues, work)`**:
        - Diagonalizes a symmetric matrix using LAPACK's `dsyevx`.
        - `threshold` is used in `dsyevx` for eigenvalue range selection (though here it seems to be an absolute threshold for values to be considered, rather than a range `vl, vu`).
        - `abs_tol` for `dsyevx` is set based on machine precision.
        - Includes error checking for `dsyevx` info flag.

## Important Variables/Constants

- **Derived Types (from `localized_basis_types`)**:
    - `separable_ri_types`: Holds data for separable RI, like overlap integrals (`ovlp2fn`, `ovlp3fn`) and coefficients (`z_coeff`).
    - `polarizability_types`: Holds data for polarizability calculations, including Kohn-Sham info (`ks`), minimax grids (`minimax`), and polarizability matrices (`tau`, `omega`, `chi`).
    - `w_engine_types`: Holds data for "w_engine", likely related to screened Coulomb interaction.
- **Hardcoded Dimensions**: Several dimensions are hardcoded in `initialization` (e.g., `n_basbas = 303`, `n_basis = 58`) and `initialize_minimax_grids` (`n_points = 6`). This suggests these might be specific to a particular system or a default configuration.
- `dp`: Double precision kind parameter from `kinds` module.

## Usage Examples

These subroutines are primarily for internal use by the `GX-LocalizedBasis` module to manage data structures.

```fortran
module my_calculation_module
    use kinds, only: dp
    use localized_basis_types
    use localized_basis_environments
    implicit none

    type(separable_ri_types) :: ri_data
    type(polarizability_types) :: polarizability_data
    type(w_engine_types) :: w_engine_data
    real(dp) :: matrix_to_power(10,10)
    real(dp) :: error_metric
    real(dp), dimension(3030,303) :: some_other_ovlp3fn ! Example dimensions

    ! 1. Initialize environments
    call initialization(ri_data)

    ! ... (Set values in ri_data.basis, ri_data.n_points if not done by default) ...
    ! ... (Initialize polarizability_data.ks, polarizability_data.ri_rs if needed) ...
    polarizability_data%ri_rs = ri_data ! shallow copy, be careful if ri_data is deallocated separately
    polarizability_data%ks%n_states = 50
    polarizability_data%ks%n_basis = ri_data%basis%n_basis

    call initialize_kohn_sham(polarizability_data)
    call initialize_minimax_grids(polarizability_data)
    call initialize_polarizability(polarizability_data)

    w_engine_data%pi_pq = polarizability_data ! shallow copy
    call initialize_w_engine(w_engine_data)

    ! 2. Perform some calculations ...
    ! For example, compute matrix^0.5
    ! call power_genmat(matrix_to_power, 10, 0.5_dp, 1e-8_dp)

    ! 3. Calculate error (assuming ri_data%ovlp3fn is populated)
    ! call calculate_error(ri_data, ri_data%basis%n_basis_pairs, ri_data%basis%n_basbas, &
    !                      some_other_ovlp3fn, error_metric)
    ! print *, "Error in ovlp3fn: ", error_metric

    ! 4. Deallocate environments when done
    call deallocations(ri_data)
    call deallocate_polarizability(polarizability_data)
    call deallocate_w_engine(w_engine_data)

contains
    ! ... other subroutines ...
end module my_calculation_module
```

## Dependencies and Interactions

- **`kinds`**: For `dp` precision.
- **`error_handling`**: For `register_exc` (though not explicitly used in the provided code snippet, it's imported).
- **`localized_basis_types`**: Crucial, as this module defines the derived types (`separable_ri_types`, `polarizability_types`, `w_engine_types`) that are initialized and managed here.
- **`lapack_interfaces`**: Provides LAPACK routines (`dgemm`, `dsyevx`, `dlaset`, `dlamch`) used for matrix operations and numerical utilities.
- This module is foundational for setting up the necessary data structures (environments) before actual physics/chemistry calculations involving localized basis sets can be performed by other modules in `GX-LocalizedBasis`.
- The hardcoded values suggest that modifying system size or basis characteristics might require changes in this module, unless these are overridden or further configured by higher-level routines.

**Notes**:
- The `power_genmat` routine has a specific way of calculating matrix powers involving negation and then re-negation of eigenvalues. The comment `matrix = -matrix` at the beginning and `eigenvalues = -eigenvalues` after diagonalization suggests it might be tailored for matrices that are negative definite or for a specific transformation. The subsequent check `if (eigenvalues(i_element) <= threshold) exit` implies it expects positive eigenvalues after the double negation, for `sqrt` to be valid.
- The threshold usage in `diagonalize_genmat` with `dsyevx('V','V','U', ...)` and `threshold` as `vl` (lower bound of interval for eigenvalues to be found) and `1.d5` as `vu` (upper bound) might be specific. If `threshold` is very small or negative, and `vu` is large, it effectively asks for all eigenvalues in that wide range.
- The `stop` statements in `power_genmat` and `diagonalize_genmat` will halt execution on certain errors, which is common in scientific computing but might be handled by error codes in a library setting.
- Shallow copies of types (e.g., `polarizability_data%ri_rs = ri_data`) mean that the underlying allocated arrays are shared. Care must be taken with deallocation order or by implementing proper deep copy / reference counting if independent lifetimes are needed. The current deallocation routines seem to deallocate components independently.
