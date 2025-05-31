# `w_engine.f90`

## Overview

This Fortran module, `w_engine`, is designed to compute the screened Coulomb interaction, often denoted as \(W(i\omega)\), within the framework of the Random Phase Approximation (RPA) or similar theories. It utilizes the irreducible polarizability \(\Pi_0(i\omega)\) computed elsewhere (likely by the `polarizability` module) and the Separable Resolution of the Identity (RI-RS) method. The calculation involves matrix operations, including inversion, to obtain \(W\) from \(\Pi_0\).

## Key Components

- **Module `w_engine`**:
    - Imports `dp` from `kinds`.
    - Imports `dgemm` from `lapack_interfaces` (though not directly used in the visible code, it might be used by `power_genmat`).
    - Imports `w_engine_types` from `localized_basis_types`.
    - Imports `initialize_w_engine`, `deallocate_w_engine`, `power_genmat` from `localized_basis_environments`.
    - Imports `get_rirs_polarizability` from the `polarizability` module.

- **`subroutine gx_w_engine(...)`**:
    - **Purpose**: This is the main public interface to compute the screened Coulomb interaction \(W\). It orchestrates the calculation of the polarizability and then uses it to compute \(W\) for each frequency point on a minimax grid.
    - **Inputs**:
        - Dimensions: `n_basis`, `n_basis_pairs`, `n_basbas`, `n_rk_points`, `n_states`.
        - Kohn-Sham data: `eigenvalues`, `eigenvectors`.
        - RI data: `ovlp2fn` (2-center overlaps on grid), `ovlp3fn` (3-center RI-V overlaps).
        - Kohn-Sham wave functions on grid: `wave`.
    - **Output**:
        - `error`: Real number, currently set to `0.0_dp` without any actual error calculation within this routine.
    - **Logic**:
        1. Initializes a local `w_engine_types` variable `we` using `initialize_w_engine`. This sets up `we%work` as an identity matrix and allocates `we%omega`.
        2. Calls `get_rirs_polarizability` to compute the irreducible polarizability \(\Pi_0(P,Q,i\omega)\) and store it in `we%pi_pq%omega`. This populates the `we%pi_pq` component of the `w_engine_types` structure.
        3. Loops over each imaginary frequency point `i_omega` from the minimax grid (`we%pi_pq%minimax%n_points`):
            a. Calls `get_screened_coulomb(we, i_omega)` to compute \(W(P,Q,i\omega)\) and store it in `we%omega(:,:,i_omega)`.
        4. Deallocates working arrays in `we` using `deallocate_w_engine`.

- **`subroutine get_screened_coulomb(we, i_omega)`**:
    - **Purpose**: Computes the screened Coulomb interaction \(W(P,Q,i\omega)\) for a single frequency point `i_omega`, given the irreducible polarizability \(\Pi_0(P,Q,i\omega)\).
    - **Formula (Conceptual)**: \(W = (1 - V \Pi_0)^{-1} V\). If working in an auxiliary basis where \(V\) is implicitly identity (or absorbed into \(\Pi_0\)), it might simplify. Given the steps, it seems to compute the correlation part of \(W_c = W - V\), as \(W_c = (1 - \Pi_0 V)^{-1} \Pi_0 V^2 - V\) or more directly \(W_c = \epsilon^{-1} V - V\), where \(\epsilon = 1 - V \Pi_0\).
    The code computes:
        1. `matrix_to_invert = I - pi_pq%omega(:,:,i_omega)` (Here `pi_pq%omega` is \(\Pi_0(i\omega)\) in auxiliary basis. `we%work` is an identity matrix \(I\)). Assuming \(V\) is implicitly handled or is the bare Coulomb interaction in this basis.
        2. `inverted_matrix = matrix_to_invert ^ -1` using `power_genmat`.
        3. `we%omega(:,:,i_omega) = inverted_matrix - I`. This suggests that `we%omega` stores the correlation part of the screened interaction, \(W_c = W - V\), if \(V\) is represented by \(I\) in this basis or context.
    - **Inputs**:
        - `we`: A `w_engine_types` variable containing \(\Pi_0\) in `we%pi_pq%omega` and workspace `we%work`.
        - `i_omega`: Integer index for the current frequency.
    - **Logic**:
        1. Computes \([I - \Pi_0(i\omega)]_{PQ}\) and stores it in `we%omega(:,:,i_omega)`.
           `we%omega(:,:,i_omega) = we%work(:,:) - we%pi_pq%omega(:,:,i_omega)`.
           (Note: `we%work` is initialized to identity in `initialize_w_engine`).
        2. Computes the inverse: \([I - \Pi_0(i\omega)]^{-1}_{PQ}\) using `power_genmat` with an exponent of -1.0. The result overwrites `we%omega(:,:,i_omega)`.
        3. Computes the final result for this frequency: `we%omega(:,:,i_omega) = we%omega(:,:,i_omega) - we%work(:,:)`. This step subtracts the identity matrix, yielding \([I - \Pi_0]^{-1} - I\).
        4. Deallocates a local `work` array (which was declared but not used in the snippet; this seems to be a leftover or mistake, as `we%work` is the one used).

## Important Variables/Constants

- `w_engine_types`: Derived type from `localized_basis_types` that holds the polarizability data (`pi_pq`), the computed screened Coulomb interaction (`omega`), and a workspace matrix (`work`).
- `we`: Instance of `w_engine_types`.
- `we%pi_pq%omega`: Input irreducible polarizability \(\Pi_0(P,Q,i\omega)\).
- `we%omega`: Output screened Coulomb interaction \(W(P,Q,i\omega)\) (specifically, the correlation part \(W_c = W-V\)).
- `we%work`: An identity matrix used in the calculation.
- `dp`: Double precision kind parameter.

## Usage Examples

The primary way to use this module is by calling `gx_w_engine`.

```fortran
program test_w_engine_calc
    use kinds, only: dp
    use w_engine        ! The module being documented
    implicit none

    ! --- Dummy parameters for gx_w_engine ---
    integer :: n_basis_val = 58
    integer :: n_basis_pairs_val = (n_basis_val * (n_basis_val + 1)) / 2
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

    ! Call the main W engine routine
    call gx_w_engine(n_basis_val, n_basis_pairs_val, n_basbas_val, &
                     n_rk_points_val, n_states_val, eigenvalues_val, &
                     eigenvectors_val, ovlp2fn_val, ovlp3fn_val, &
                     wave_val, error_val)

    print *, "W engine calculation finished. Error metric (if implemented): ", error_val

    ! After the call, 'we%omega' (within gx_w_engine) would have held the result,
    ! but it's deallocated internally by gx_w_engine.
    ! To use W, one would need a modified gx_w_engine that returns 'we' or
    ! passes it as an intent(out) and manages deallocation externally.

end program test_w_engine_calc
```

## Dependencies and Interactions

- **`kinds`**: For `dp` precision.
- **`lapack_interfaces`**: `dgemm` is likely used by `power_genmat`.
- **`localized_basis_types`**: Defines `w_engine_types` and its components like `polarizability_types`.
- **`localized_basis_environments`**:
    - `initialize_w_engine`: Sets up the `w_engine_types` instance, notably preparing `we%work` as an identity matrix.
    - `deallocate_w_engine`: Cleans up memory.
    - `power_genmat`: Used to compute the matrix inverse \([I - \Pi_0]^{-1}\).
- **`polarizability` module**: `get_rirs_polarizability` is called to compute the irreducible polarizability \(\Pi_0\), which is the essential input for the screening calculation.

This module performs a crucial step in many-body perturbation theory calculations (like GW approximation) by computing the screened Coulomb interaction. It builds directly upon the polarizability computed by other modules. The current `gx_w_engine` subroutine computes \(W\) but deallocates it internally, suggesting it might be part of a larger workflow where \(W\) is immediately consumed, or it needs modification to output the `we` structure for later use. The `error` output parameter is not currently used to reflect any computational error from this stage. The deallocation of a local `work` variable in `get_screened_coulomb` that was never allocated or used in that scope appears to be a minor code cleanup issue.
