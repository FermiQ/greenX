# `polarizability.f90`

## Overview

This Fortran module, `polarizability`, is dedicated to computing the irreducible polarizability \(\Pi_0(i\omega)\) in the auxiliary basis representation (`[PI(iw)]_PQ`). It primarily uses the Separable Resolution of the Identity (RI-RS) method. The calculation involves several steps: obtaining RI coefficients, processing Kohn-Sham wave functions, utilizing minimax grids for time and frequency transformations, evaluating the polarizability in real space (\(\chi_0\)), and then transforming it to the auxiliary basis representation.

## Key Components

- **Module `polarizability`**:
    - Imports types and routines from `kinds`, `error_handling`, `lapack_interfaces`, `localized_basis_types`, `localized_basis_environments`, `separable_ri`, and `gx_minimax`.

- **`subroutine gx_rirs_polarizability(...)`**:
    - **Purpose**: High-level routine to compute RI-RS polarizability and then immediately deallocate most working arrays.
    - **Inputs**: Basis/grid dimensions (`n_basis`, `n_basis_pairs`, etc.), Kohn-Sham eigenvalues and eigenvectors, 2-center and 3-center overlap integrals (`ovlp2fn`, `ovlp3fn`), and Kohn-Sham wave functions on a real-space grid (`wave`).
    - **Output**: `error` (though currently set to `0.0_dp` without calculation).
    - **Logic**:
        1. Calls `get_rirs_coefficients` to populate `pi_pq%ri_rs`.
        2. Calls `get_ks_wave` to process and store Kohn-Sham data into `pi_pq%ks`.
        3. Calls `get_minimax_grids` to set up `pi_pq%minimax`.
        4. Calls `initialize_polarizability` to allocate arrays in `pi_pq`.
        5. Calls `evaluate_polarizability` to compute `pi_pq%omega`.
        6. Calls `deallocate_polarizability` to free memory.

- **`subroutine get_rirs_polarizability(...)`**:
    - **Purpose**: Similar to `gx_rirs_polarizability` but intended to keep the computed polarizability (`pi_pq%omega`) by using `deallocate_polarizability(pi_pq, keep_pi=.true.)`.
    - **Inputs & Logic**: Almost identical to `gx_rirs_polarizability`, but passes `keep_pi=.true.` to the deallocation routine.

- **`subroutine get_ks_wave(...)`**:
    - **Purpose**: Initializes the Kohn-Sham part (`pi_pq%ks`) of the `polarizability_types` structure.
    - **Inputs**: `pi_pq` (polarizability environment), dimensions (`n_basis`, `n_states`, `n_rk_points`), Kohn-Sham eigenvalues, eigenvectors, and wave functions on the grid.
    - **Logic**:
        1. Sets dimensions in `pi_pq%ks` and `pi_pq%ri_rs`.
        2. Calls `initialize_kohn_sham` to allocate arrays within `pi_pq%ks`.
        3. Copies input KS data into `pi_pq%ks%eigenvalues`, `pi_pq%ks%eigenvectors`, `pi_pq%ks%wave`.

- **`subroutine get_minimax_grids(pi_pq)`**:
    - **Purpose**: Initializes the minimax grid data (`pi_pq%minimax`).
    - **Input**: `pi_pq`.
    - **Logic**:
        1. Calls `initialize_minimax_grids` to allocate arrays.
        2. Sets hardcoded `e_tran_max = 0.260_dp` and `e_tran_min = 31.725_dp`. (Note: `e_tran_min` is larger than `e_tran_max`, which might be an issue for `gx_minimax_grid` depending on its internal logic, or these variable names might be misleading).
        3. Checks if `e_tran_min <= 0.0_dp` (for metal systems, not supported).
        4. Calls `gx_minimax_grid` (from `gx_minimax` module) to get time/frequency points, weights, and transformation matrices.
        5. Stores results in `pi_pq%minimax`.

- **`subroutine evaluate_polarizability(pi_pq)`**:
    - **Purpose**: Computes the polarizability `pi_pq%omega` in the auxiliary basis, frequency domain.
    - **Input**: `pi_pq`.
    - **Logic**:
        1. Loops over imaginary time grid points (`i_tau` from `pi_pq%minimax%n_points`).
            a. Calls `evaluate_chi` to get \(\chi_0(i\tau_{i\_tau})\) in `pi_pq%chi%matrix`.
            b. Calls `transform_chi_to_pi` to convert `pi_pq%chi%matrix` to `pi_pq%tau` (auxiliary basis, time domain).
            c. Performs a cosine transform from time to frequency:
               `pi_pq%omega(:,:,i_omega) += pi_pq%tau(:,:) * pi_pq%minimax%cos_tf(i_omega, i_tau)` for each frequency `i_omega`.

- **`subroutine evaluate_chi(pi_pq, i_tau)`**:
    - **Purpose**: Computes the irreducible polarizability \(\chi_0(i\tau)\) on the real-space grid.
    - **Inputs**: `pi_pq`, current time index `i_tau`.
    - **Logic**:
        1. Gets `tau = pi_pq%minimax%tau(i_tau)`.
        2. Allocates `pi_pq%chi%matrix`.
        3. Loops `i_spin` from 1 to `pi_pq%ks%n_spin`:
            a. Calls `get_green_forward` to compute \(G(r, r', i\tau)\) in `pi_pq%chi%green_forward`.
            b. Calls `get_green_backward` to compute \(G(r, r', -i\tau)\) in `pi_pq%chi%green_backward`.
        4. Computes \(\chi_0(r, r', i\tau) = -G(r, r', i\tau) \times G(r', r, -i\tau)\) (Note: usually \(G(r,r',-i\tau)\) but here it is \(G(r',r, -i\tau)\) due to `dgemm` with `"n","t"` for G and G_backward, effectively G(r,r') and G(r',r) where the second G is backward in time). The product is element-wise if `green_forward` and `green_backward` are \(G(r_i,r_j,i\tau)\) and \(G(r_i,r_j,-i\tau)\). If it's a matrix product, the interpretation changes. Given the loop structure, it appears to be `pi_pq%chi%matrix(i_point,j_point) = -pi_pq%chi%green_forward(i_point,j_point) * pi_pq%chi%green_backward(i_point,j_point)`.
        5. Deallocates `green_forward` and `green_backward`.

- **`subroutine get_green_forward(pi_pq, i_spin, tau)`**:
    - **Purpose**: Computes the occupied part of the Green's function \(G_{occ}(r, r', i\tau)\) on the real-space grid.
    - **Inputs**: `pi_pq`, spin index `i_spin`, time `tau`.
    - **Logic**:
        1. Allocates `pi_pq%chi%green_forward` and temporary `wave_occ`.
        2. Scales occupied KS wave functions: `wave_occ(:, i_state) = pi_pq%ks%wave(:,i_state,i_spin) * exp(-0.5*tau*(pi_pq%ks%e_fermi - pi_pq%ks%eigenvalues(i_state,i_spin)))`.
        3. Computes \(G_{occ}(r_k, r_{k'}) = \sum_{m \in occ} \psi_m(r_k) \psi_m(r_{k'})\) (scaled) using `dgemm("n","t", ... wave_occ, wave_occ, ... pi_pq%chi%green_forward)`.
        4. Deallocates `wave_occ`.

- **`subroutine get_green_backward(pi_pq, i_spin, tau)`**:
    - **Purpose**: Computes the virtual part of the Green's function \(G_{virt}(r, r', -i\tau)\) on the real-space grid.
    - **Inputs**: `pi_pq`, spin index `i_spin`, time `tau`.
    - **Logic**:
        1. Allocates `pi_pq%chi%green_backward` and temporary `wave_virt`.
        2. Scales virtual KS wave functions: `wave_virt(:, i_state) = pi_pq%ks%wave(:,pi_pq%ks%n_occ + i_state,i_spin) * exp(-0.5*tau*(pi_pq%ks%eigenvalues(pi_pq%ks%n_occ + i_state,i_spin) - pi_pq%ks%e_fermi))`.
        3. Computes \(G_{virt}(r_k, r_{k'}) = \sum_{a \in virt} \psi_a(r_k) \psi_a(r_{k'})\) (scaled) using `dgemm("n","t", ... wave_virt, wave_virt, ... pi_pq%chi%green_backward)`.
        4. Deallocates `wave_virt`.

- **`subroutine transform_chi_to_pi(pi_pq)`**:
    - **Purpose**: Transforms \(\chi_0(r, r', i\tau)\) (stored in `pi_pq%chi%matrix`) to \(\Pi_0(P, Q, i\tau)\) (stored in `pi_pq%tau`) using RI coefficients.
    - **Input**: `pi_pq`.
    - **Logic**:
        1. Allocates temporary `mat_aux`.
        2. Right multiplication: `mat_aux(k,Q) = \sum_{k'} \chi_0(k,k') Z(Q,k')` via `dgemm("n","t", ..., pi_pq%chi%matrix, pi_pq%ri_rs%z_coeff, ..., mat_aux)`. Note: `z_coeff` is (aux_basis, rk_point), so `Z_Qk'` is `z_coeff(Q, k')`.
        3. Left multiplication: `Pi(P,Q) = \sum_{k} Z(P,k) mat_aux(k,Q)` via `dgemm("n","n", ..., pi_pq%ri_rs%z_coeff, mat_aux, ..., pi_pq%tau)`.
        4. Deallocates `mat_aux`.

## Important Variables/Constants

- `polarizability_types`: Derived type from `localized_basis_types` that holds all necessary data.
- `pi_pq`: Variable of `polarizability_types`, central to most subroutines.
- `e_tran_min`, `e_tran_max`: Hardcoded energy transition values for minimax grids. Their ordering (`e_tran_min > e_tran_max`) is unusual and might indicate specific requirements or a misunderstanding of the variable names' common interpretation.
- `dp`: Double precision kind parameter.

## Usage Examples

The main entry points are `gx_rirs_polarizability` or `get_rirs_polarizability`.

```fortran
program test_polarizability_calc
    use kinds, only: dp
    use localized_basis_types ! For polarizability_types if used directly
    use polarizability        ! The module being documented
    implicit none

    ! --- Dummy parameters for gx_rirs_polarizability ---
    integer :: n_basis_val = 58
    integer :: n_basis_pairs_val = 1711 ! n_basis * (n_basis+1)/2
    integer :: n_basbas_val = 303 ! Auxiliary basis size
    integer :: n_rk_points_val = 638 ! Real-space grid points for RI
    integer :: n_states_val = 100 ! Number of Kohn-Sham states

    real(dp) :: error_val
    real(dp), dimension(n_states_val) :: eigenvalues_val
    real(dp), dimension(n_basis_val, n_states_val) :: eigenvectors_val
    real(dp), dimension(n_basis_pairs_val, n_rk_points_val) :: ovlp2fn_val
    real(dp), dimension(n_basis_pairs_val, n_basbas_val) :: ovlp3fn_val
    real(dp), dimension(n_states_val, n_rk_points_val) :: wave_val
    ! --- End dummy parameters ---

    ! Initialize dummy arrays (replace with actual data loading/computation)
    eigenvalues_val = [(0.1_dp * real(i), i=1,n_states_val)]
    eigenvectors_val = 0.1_dp
    ovlp2fn_val = 0.01_dp
    ovlp3fn_val = 0.001_dp
    wave_val = 0.05_dp

    ! Call the main polarizability calculation routine
    call gx_rirs_polarizability(n_basis_val, n_basis_pairs_val, n_basbas_val, &
                                n_rk_points_val, n_states_val, eigenvalues_val, &
                                eigenvectors_val, ovlp2fn_val, ovlp3fn_val, &
                                wave_val, error_val)

    print *, "Polarizability calculation finished. Error metric (if implemented): ", error_val

    ! If get_rirs_polarizability were used to keep pi_pq%omega:
    ! type(polarizability_types) :: stored_polarizability
    ! call get_rirs_polarizability(stored_polarizability, n_basis_val, n_basis_pairs_val, &
    !                              n_basbas_val, n_rk_points_val, n_states_val, &
    !                              eigenvalues_val, eigenvectors_val, ovlp2fn_val, &
    !                              ovlp3fn_val, wave_val)
    ! ! ... use stored_polarizability%omega ...
    ! ! Then manually deallocate:
    ! ! call deallocate_polarizability(stored_polarizability, .true.) ! to keep omega
    ! ! if (allocated(stored_polarizability%omega)) deallocate(stored_polarizability%omega)
    ! ! call deallocate_polarizability(stored_polarizability) ! to clean everything else related to it.
    ! ! Proper deallocation would need careful handling of what was kept vs what was already freed.
    ! ! The `deallocate_polarizability` from `localized_basis_environments` already handles `keep_pi`.

end program test_polarizability_calc
```

## Dependencies and Interactions

- **`localized_basis_types`**: Defines `polarizability_types`.
- **`localized_basis_environments`**: Provides initialization and deallocation routines for components of `polarizability_types`.
- **`separable_ri`**: Provides `get_rirs_coefficients` to compute RI fitting coefficients (`z_coeff`).
- **`gx_minimax`**: Provides `gx_minimax_grid` for generating time and frequency points and transformation matrices.
- **`lapack_interfaces`**: Specifically `dgemm` is used for matrix multiplications in Green's function calculation and basis transformations.
- **`error_handling`**: `register_exc` is used to report errors (e.g., metal system detection, minimax grid errors).

This module orchestrates the complex calculation of polarizability by integrating various components like RI, Kohn-Sham data, and time-frequency transformations. The hardcoded energy transition values in `get_minimax_grids` are a point of attention for general applicability. The `error` output of `gx_rirs_polarizability` is currently not meaningfully computed. The distinction between `gx_rirs_polarizability` and `get_rirs_polarizability` lies in whether the computed `pi_pq%omega` is intended to be kept or immediately deallocated.
