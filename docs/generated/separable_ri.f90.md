# `separable_ri.f90`

## Overview

This Fortran module, `separable_ri`, implements routines related to the Separable Resolution of the Identity (RI) method, specifically for computing three-center overlap integrals and the corresponding RI fitting coefficients (`Z_Pk`, also referred to as `M_Pk` in comments). The module uses LAPACK for matrix operations and relies on types and utility functions defined in other modules of the `GX-LocalizedBasis` component.

## Key Components

- **Module `separable_ri`**:
    - Imports `dp` from `kinds`.
    - Imports `dgemm` from `lapack_interfaces`.
    - Imports `separable_ri_types` from `localized_basis_types`.
    - Imports `calculate_error`, `power_genmat`, `initialization`, `deallocations` from `localized_basis_environments`.

- **`subroutine gx_rirs_coefficients(...)`**:
    - **Purpose**: A high-level driver routine that computes RI-RS three-center overlap integrals (`ri_rs%ovlp3fn`), calculates the error against provided RI-V (`ovlp3fn`) integrals, and then deallocates temporary data.
    - **Inputs**:
        - `n_basis_pairs`, `n_loc_basbas`, `n_rk_points`: Dimensions.
        - `ovlp2fn`: Real array, product of two NAO basis functions \(\Phi_i(r_k) \Phi_j(r_k)\) evaluated on a real-space grid (dimension: `n_basis_pairs, n_rk_points`). Let's denote this as \(D_{ij}(r_k)\).
        - `ovlp3fn`: Real array, RI-V three-center overlap integrals \((ij|P)\) (dimension: `n_basis_pairs, n_loc_basbas`). This is the reference to compare against.
    - **Output**:
        - `error`: Real number, maximum absolute difference between the computed RI-RS `ri_rs%ovlp3fn` and the input `ovlp3fn`.
    - **Logic**:
        1. Initializes a local `separable_ri_types` variable `ri_rs` using `initialization`.
        2. Copies the input `ovlp2fn` into `ri_rs%ovlp2fn`.
        3. Calls `compute_ovlp3fn` to calculate the RI-RS three-center overlaps and store them in `ri_rs%ovlp3fn`. This internally first computes `ri_rs%z_coeff`.
        4. Calls `calculate_error` to compare `ri_rs%ovlp3fn` with the input `ovlp3fn` and compute the `error`.
        5. Deallocates `ri_rs` using `deallocations`.

- **`subroutine get_rirs_coefficients(...)`**:
    - **Purpose**: Computes the RI-RS fitting coefficients `ri_rs%z_coeff` and stores them in the provided `ri_rs` structure. It keeps these coefficients allocated.
    - **Inputs**:
        - `ri_rs`: A `separable_ri_types` variable (intent inout).
        - Dimensions: `n_basis_pairs`, `n_loc_basbas`, `n_rk_points`.
        - `ovlp2fn`: Input \(D_{ij}(r_k)\) array.
        - `ovlp3fn`: Input RI-V \((ij|P)\) array (used to compute `z_coeff`).
    - **Logic**:
        1. Calls `initialization(ri_rs)` (Note: This will re-initialize `ri_rs` including its dimensions based on hardcoded values in `localized_basis_environments`, which might conflict if `ri_rs` was already partially set up. It also allocates arrays in `ri_rs`).
        2. Copies input `ovlp2fn` to `ri_rs%ovlp2fn`.
        3. Calls `get_coeff_zrs` to compute `ri_rs%z_coeff` using `ovlp3fn` (RI-V) and `ri_rs%ovlp2fn`.
        4. Calls `deallocations(ri_rs, keep_coeff=.true.)` to deallocate parts of `ri_rs` but explicitly keeps `ri_rs%z_coeff`.

- **`subroutine compute_ovlp3fn(...)`**:
    - **Purpose**: Computes the RI-RS approximated three-center overlap integrals \((ij|P)_{RS} = \sum_k D_{ij}(r_k) Z_{Pk}\) and stores them in `ri_rs%ovlp3fn`.
    - **Inputs**:
        - `ri_rs`: A `separable_ri_types` variable (intent inout).
        - Dimensions: `n_basis_pairs`, `n_loc_basbas`.
        - `ovlp3fn`: Input RI-V \((ij|P)\) array (passed to `get_coeff_zrs`).
    - **Logic**:
        1. Calls `get_coeff_zrs` to compute `ri_rs%z_coeff` using the input `ovlp3fn` (RI-V reference) and `ri_rs%ovlp2fn`.
        2. Computes `ri_rs%ovlp3fn = ri_rs%ovlp2fn * (ri_rs%z_coeff)^T` using `dgemm`.
           `(ij|P)_{RS} = \sum_k D_{ij}(r_k) Z_{Pk}`.
           `dgemm('N', 'T', n_basis_pairs, n_loc_basbas, ri_rs%n_points, ... ri_rs%ovlp2fn, ri_rs%z_coeff, ... ri_rs%ovlp3fn)`.
           This matches `ovlp2fn(pair, k) * z_coeff(P, k)`.

- **`subroutine get_coeff_zrs(...)`**:
    - **Purpose**: Computes the least-squares RI fitting coefficients \(Z_{Pk}\) (stored in `ri_rs%z_coeff`).
    - **Formula**: \(Z = A B^{-1}\), where \(A_{Pk'} = \sum_{ij} (ij|P) D_{ij}(r_{k'})\) and \(B_{k k'} = \sum_{ij} D_{ij}(r_k) D_{ij}(r_{k'})\).
    - **Inputs**:
        - `ri_rs`: A `separable_ri_types` variable (intent inout, `z_coeff` is output, `ovlp2fn` is input).
        - Dimensions: `n_basis_pairs`, `n_loc_basbas`.
        - `ovlp3fn`: Input RI-V \((ij|P)\) array.
    - **Logic**:
        1. Allocates temporary matrices `aux_mata` and `aux_matb`.
        2. Computes \(A = (ovlp3fn)^T * (ri_rs\%ovlp2fn)\) using `dgemm`.
           `aux_mata(P,k') = \sum_{pair} ovlp3fn(pair,P) * ri_rs%ovlp2fn(pair,k')`. This corresponds to \(A_{Pk'}\).
           `dgemm('T','N',n_loc_basbas,ri_rs%n_points,n_basis_pairs, ... ovlp3fn, ri_rs%ovlp2fn, ... aux_mata)`.
        3. Computes \(B = (ri_rs\%ovlp2fn)^T * (ri_rs\%ovlp2fn)\) using `dgemm`.
           `aux_matb(k,k') = \sum_{pair} ri_rs%ovlp2fn(pair,k) * ri_rs%ovlp2fn(pair,k')`. This corresponds to \(B_{kk'}\).
           `dgemm('T', 'N',ri_rs%n_points,ri_rs%n_points,n_basis_pairs, ... ri_rs%ovlp2fn, ri_rs%ovlp2fn, ... aux_matb)`.
        4. Computes \(B^{-1}\) using `power_genmat(aux_matb, ri_rs%n_points, -1.0_dp, 1.0d-10)`.
        5. Computes \(Z = A * B^{-1}\) using `dgemm` and stores in `ri_rs%z_coeff`.
           `ri_rs%z_coeff(P,k) = \sum_{k'} aux_mata(P,k') * aux_matb(k',k)`. (Assuming `aux_matb` now stores B_inv(k',k)).
           `dgemm('N', 'N',n_loc_basbas,ri_rs%n_points,ri_rs%n_points, ... aux_mata, aux_matb, ... ri_rs%z_coeff)`.
        6. Deallocates `aux_mata`, `aux_matb`.

## Important Variables/Constants

- `separable_ri_types`: Derived type holding RI-related data, especially `z_coeff`, `ovlp2fn`, `ovlp3fn`.
- `ri_rs`: Instance of `separable_ri_types` used within the routines.
- `ovlp2fn`: Represents \(D_{ij}(r_k)\), the product of two basis functions on a real-space grid point \(r_k\).
- `ovlp3fn`: In `gx_rirs_coefficients` and `get_rirs_coefficients`, this is the input RI-V three-center overlap integrals \((ij|P)\). In `ri_rs%ovlp3fn`, it stores the computed RI-RS three-center overlaps.
- `z_coeff`: The RI fitting coefficients \(Z_{Pk}\), where \(P\) is an auxiliary basis index and \(k\) is a real-space grid point index.

## Usage Examples

```fortran
program test_separable_ri
    use kinds, only: dp
    use localized_basis_types, only: separable_ri_types
    use separable_ri
    ! Assuming localized_basis_environments is available for full setup if needed
    implicit none

    ! --- Dummy parameters ---
    integer, parameter :: n_basis_pairs_val = 1711
    integer, parameter :: n_loc_basbas_val = 303
    integer, parameter :: n_rk_points_val = 638

    real(dp) :: error_val
    real(dp), dimension(n_basis_pairs_val, n_rk_points_val) :: ovlp2fn_data
    real(dp), dimension(n_basis_pairs_val, n_loc_basbas_val) :: ovlp3fn_ref_data ! RI-V reference

    type(separable_ri_types) :: ri_environment_for_coeffs

    ! Initialize dummy arrays (replace with actual data)
    ovlp2fn_data = 0.01_dp
    ovlp3fn_ref_data = 0.001_dp

    ! Option 1: Compute RI-RS overlaps and error vs RI-V
    call gx_rirs_coefficients(n_basis_pairs_val, n_loc_basbas_val, n_rk_points_val, &
                              ovlp2fn_data, ovlp3fn_ref_data, error_val)
    print *, "Error between RI-RS and RI-V ovlp3fn: ", error_val

    ! Option 2: Get Z_Pk coefficients stored in ri_environment_for_coeffs
    ! Note: get_rirs_coefficients calls initialization internally, which uses hardcoded sizes.
    ! For this to work as intended with external ri_environment_for_coeffs,
    ! initialization logic might need adjustment or ri_environment_for_coeffs
    ! should be populated by localized_basis_environments routines first if custom
    ! dimensions are needed beyond the hardcoded ones.
    ! However, get_rirs_coefficients primarily *outputs* to ri_rs%z_coeff after re-init.

    ! To properly use get_rirs_coefficients and retain z_coeff:
    call initialization(ri_environment_for_coeffs) ! Initialize with default/hardcoded sizes
    ! Override necessary dimensions if they differ from defaults and are used by get_coeff_zrs logic
    ri_environment_for_coeffs%n_points = n_rk_points_val
    ! Ensure basis component of ri_environment_for_coeffs is consistent if used by get_coeff_zrs
    ! (though it seems n_loc_basbas is passed directly)

    call get_rirs_coefficients(ri_environment_for_coeffs, n_basis_pairs_val, &
                               n_loc_basbas_val, n_rk_points_val, &
                               ovlp2fn_data, ovlp3fn_ref_data)

    ! Now ri_environment_for_coeffs%z_coeff contains the Z_Pk coefficients.
    ! print *, "Z_Pk(1,1) = ", ri_environment_for_coeffs%z_coeff(1,1) ! Example access

    ! Clean up the externally managed type
    if (allocated(ri_environment_for_coeffs%z_coeff)) then
        print *, "z_coeff is allocated."
        ! ... use z_coeff ...
        deallocate(ri_environment_for_coeffs%z_coeff)
    end if
    ! Other components of ri_environment_for_coeffs might have been deallocated by
    ! deallocations(ri_rs, keep_coeff=.true.) inside get_rirs_coefficients.

end program test_separable_ri
```

## Dependencies and Interactions

- **`kinds`**: For `dp` precision.
- **`lapack_interfaces`**: `dgemm` is used for matrix multiplications.
- **`localized_basis_types`**: Defines `separable_ri_types`.
- **`localized_basis_environments`**:
    - `initialization` and `deallocations`: Used to manage the `ri_rs` type instance. The `initialization` subroutine uses hardcoded dimensions, which might be a limitation or a feature for default scenarios.
    - `calculate_error`: Used to compare RI-RS and RI-V three-center overlaps.
    - `power_genmat`: Used to compute the inverse of a matrix (\(B^{-1}\)) in `get_coeff_zrs`.
- The module provides a key step in the overall RI-based calculations, producing the fitting coefficients \(Z_{Pk}\) which are then used in other parts of the `GX-LocalizedBasis` (e.g., in the `polarizability` module to transform quantities from real-space grid to auxiliary basis representation).

**Important Considerations**:
- The `initialization` call within `get_rirs_coefficients` will reset `ri_rs` based on potentially hardcoded default dimensions from `localized_basis_environments`. If `ri_rs` is passed in with pre-allocated members or custom dimensions intended to be used by `get_coeff_zrs`, this could lead to unexpected behavior or errors. It seems `get_rirs_coefficients` is designed to populate a fresh or temporary `ri_rs` and output `z_coeff` (which is kept).
- The formula for \(Z_{Pk}\) coefficients involves matrix inversion (\(B^{-1}\)), which is handled by `power_genmat`. The stability and accuracy of this step can be sensitive to the conditioning of matrix \(B\).
- The RI-V `ovlp3fn` serves as the reference data for deriving the \(Z_{Pk}\) coefficients. The quality of these coefficients, and subsequently the RI-RS approximation, depends on this input and the choice of the real-space grid and auxiliary basis.
