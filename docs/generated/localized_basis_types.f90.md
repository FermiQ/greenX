# `localized_basis_types.f90`

## Overview

This Fortran module, `localized_basis_types`, defines a set of derived data types (structs) that are used throughout the `GX-LocalizedBasis` component of the GreenX library. These types encapsulate data related to basis sets, grids, atomic species, Kohn-Sham electronic structure, minimax approximations, polarizability, and Resolution of Identity (RI) techniques. The module serves as a centralized place for data structure definitions, promoting modularity and clear data organization.

## Key Components (Derived Types)

### Sub-Types (Building Blocks)

- **`basis_types`**:
    - **Purpose**: Stores counts related to basis functions.
    - **Members**:
        - `n_basbas`: Integer, likely number of auxiliary basis functions.
        - `n_basis`: Integer, number of primary basis functions.
        - `n_loc_basbas`: Integer, likely number of localized auxiliary basis functions.
        - `n_basis_pairs`: Integer, number of pairs of primary basis functions.

- **`grids_types`**:
    - **Purpose**: Defines real-space integration grids.
    - **Members**:
        - `m_points`: Integer, (purpose not immediately clear, perhaps max points or number of grid types).
        - `n_points`: `integer, dimension(:), allocatable`, number of points for each grid.
        - `r_points`: `real(kind=8), dimension(:,:,:), allocatable`, coordinates of grid points (e.g., grid_idx, point_idx, xyz).
        - `w_points`: `real(kind=8), dimension(:,:), allocatable`, weights of grid points.

- **`species_types`**:
    - **Purpose**: Stores information about atomic species and their coordinates.
    - **Members**:
        - `n_atoms`: Integer, total number of atoms.
        - `n_species`: Integer, number of unique atomic species.
        - `species`: `character, dimension(:), allocatable`, names/symbols of species.
        - `coords`: `real(kind=8), dimension(:,:), allocatable`, atomic coordinates (atom_idx, xyz).

- **`minimax_types`**:
    - **Purpose**: Holds data for minimax approximations, likely for frequency/time representations.
    - **Members**:
        - `n_points`: Integer, number of minimax points.
        - `cos_tf`: `real(kind=8), dimension(:,:), allocatable`, transformation matrix (e.g., cosine transformation coefficients).
        - `omega`: `real(kind=8), dimension(:), allocatable`, frequency points.
        - `tau`: `real(kind=8), dimension(:), allocatable`, time points.
        - `weights`: `real(kind=8), dimension(:), allocatable`, weights for RPA (Random Phase Approximation) or minimax quadrature.

- **`kohn_sham_types`**:
    - **Purpose**: Encapsulates data from Kohn-Sham DFT calculations.
    - **Members**:
        - `n_basis`: Integer, number of basis functions (consistent with `basis_types`).
        - `n_homo`: Integer, index of the Highest Occupied Molecular Orbital.
        - `n_lumo`: Integer, index of the Lowest Unoccupied Molecular Orbital.
        - `n_occ`: Integer, number of occupied states.
        - `n_states`: Integer, total number of Kohn-Sham states.
        - `n_spin`: Integer, number of spin channels.
        - `n_virt`: Integer, number of virtual (unoccupied) states.
        - `e_fermi`: Integer (likely should be `real(kind=8)`), Fermi energy.
        - `eigenvalues`: `real(kind=8), dimension(:,:), allocatable`, Kohn-Sham eigenvalues (state_idx, spin_idx).
        - `eigenvectors`: `real(kind=8), dimension(:,:,:), allocatable`, Kohn-Sham eigenvectors (basis_idx, state_idx, spin_idx).
        - `occupied`: `real(kind=8), dimension(:,:), allocatable`, occupied wave functions (basis_idx, state_idx).
        - `virtual`: `real(kind=8), dimension(:,:), allocatable`, virtual wave functions (basis_idx, state_idx).
        - `wave`: `real(kind=8), dimension(:,:,:), allocatable`, all Kohn-Sham wave functions (basis_idx, state_idx, spin_idx or point_idx, state_idx, spin_idx if on a grid).

- **`real_space_chi_types`**:
    - **Purpose**: Stores components of the polarizability \(\chi\) in real space or related representations.
    - **Members**:
        - `matrix`: `real(kind=8), dimension(:,:), allocatable`, Polarizability matrix, possibly on real-space grid points (`r_k`).
        - `green_forward`: `real(kind=8), dimension(:,:), allocatable`, Green's function in forward time/frequency.
        - `green_backward`: `real(kind=8), dimension(:,:), allocatable`, Green's function in backward time/frequency.

### Global (Main) Types

- **`separable_ri_types`**:
    - **Purpose**: Top-level type for Separable Resolution of Identity (RI) calculations.
    - **Members**:
        - `n_points`: Integer, number of RI integration points.
        - `error`: `real(kind=8)`, stores an error metric associated with the RI approximation.
        - `ovlp2fn`: `real(kind=8), dimension(:,:), allocatable`, 2-center overlap integrals (basis_pair_idx, ri_point_idx).
        - `ovlp3fn`: `real(kind=8), dimension(:,:), allocatable`, 3-center overlap integrals (basis_pair_idx, aux_basis_idx).
        - `z_coeff`: `real(kind=8), dimension(:,:), allocatable`, RI fitting coefficients (aux_basis_idx, ri_point_idx).
        - `basis`: `type(basis_types)`, nested basis information.
        - `grids`: `type(grids_types)`, nested grid information.
        - `species`: `type(species_types)`, nested atomic species information.

- **`polarizability_types`**:
    - **Purpose**: Top-level type for polarizability calculations.
    - **Members**:
        - `tau`: `real(kind=8), dimension(:,:), allocatable`, Polarizability or related quantity in time domain (aux_basis_idx_1, aux_basis_idx_2).
        - `omega`: `real(kind=8), dimension(:,:,:), allocatable`, Polarizability or related quantity in frequency domain (aux_basis_idx_1, aux_basis_idx_2, freq_idx).
        - `chi`: `type(real_space_chi_types)`, nested real-space polarizability components.
        - `ks`: `type(kohn_sham_types)`, nested Kohn-Sham data.
        - `minimax`: `type(minimax_types)`, nested minimax approximation data.
        - `ri_rs`: `type(separable_ri_types)`, nested separable RI data, indicating polarizability calculations might use RI.

- **`w_engine_types`**:
    - **Purpose**: Top-level type for "w_engine", likely related to calculating the screened Coulomb interaction W.
    - **Members**:
        - `pi_pq`: `type(polarizability_types)`, nested polarizability data, suggesting W is derived from \(\chi\).
        - `omega`: `real(kind=8), dimension(:,:,:), allocatable`, Screened Coulomb interaction W in frequency domain (aux_basis_idx_1, aux_basis_idx_2, freq_idx).
        - `work`: `real(kind=8), dimension(:,:), allocatable`, A workspace array.

## Important Variables/Constants

This module primarily defines types, not variables or constants with specific values, except for the implicit `real(kind=8)` which specifies double precision for real numbers. The structure and naming of type members (e.g., `n_basis`, `eigenvalues`) are themselves important for understanding data flow.

## Usage Examples

These types are intended to be used in variable declarations within other modules of the `GX-LocalizedBasis` library.

```fortran
module calculation_example
    use localized_basis_types
    use kinds, only: dp ! Assuming kinds module provides dp, and kind=8 is dp
    implicit none

    type(separable_ri_types) :: my_ri_data
    type(polarizability_types) :: my_polarizability
    type(w_engine_types) :: my_w_interaction

    subroutine setup_calculation
        ! Example: Initialize basis counts for RI
        my_ri_data%basis%n_basis = 100
        my_ri_data%basis%n_basbas = 500
        my_ri_data%basis%n_loc_basbas = 400
        my_ri_data%basis%n_basis_pairs = (my_ri_data%basis%n_basis * (my_ri_data%basis%n_basis + 1)) / 2
        my_ri_data%n_points = 2000

        ! Allocate arrays (dimensions based on counts above)
        allocate(my_ri_data%ovlp2fn(my_ri_data%basis%n_basis_pairs, my_ri_data%n_points))
        allocate(my_ri_data%z_coeff(my_ri_data%basis%n_basbas, my_ri_data%n_points))
        ! ... and so on for other allocatable members ...

        ! Populate Kohn-Sham data for polarizability calculation
        my_polarizability%ks%n_basis = my_ri_data%basis%n_basis
        my_polarizability%ks%n_states = 150
        my_polarizability%ks%n_occ = 50
        allocate(my_polarizability%ks%eigenvalues(my_polarizability%ks%n_states, 1))
        ! ... initialize eigenvalues ...

        ! Link RI data to polarizability calculation
        my_polarizability%ri_rs = my_ri_data ! This is a shallow copy

        ! Link polarizability data to W engine
        my_w_interaction%pi_pq = my_polarizability ! Shallow copy

        ! ... further allocations and initializations ...
    end subroutine setup_calculation

    ! ... other subroutines to perform calculations using these data structures ...

    subroutine cleanup_calculation
        if (allocated(my_ri_data%ovlp2fn)) deallocate(my_ri_data%ovlp2fn)
        if (allocated(my_ri_data%z_coeff)) deallocate(my_ri_data%z_coeff)
        ! ... and so on for other allocatable members ...
    end subroutine cleanup_calculation

end module calculation_example
```

## Dependencies and Interactions

- **`implicit none`**: Enforces explicit type declaration for all variables.
- **`real(kind=8)`**: Specifies that all real numbers are of kind 8 (typically double precision). This implies a dependency on a system or `kinds` module that defines this kind.
- This module is a foundational part of `GX-LocalizedBasis`. Other modules within this component will use these types to declare variables, pass data between subroutines, and organize the results of their computations.
- The structure of these types (e.g., `polarizability_types` containing `kohn_sham_types` and `separable_ri_types`) indicates logical dependencies and workflows within the larger calculation scheme (e.g., Kohn-Sham results and RI approximations are inputs or components for polarizability calculations).
- The use of `allocatable` members means that memory management (allocation and deallocation) for these arrays must be handled by the routines that use these types, typically in initialization and finalization steps (as seen in `localized_basis_environments.f90`).

This module provides the blueprint for the data containers used in the localized basis set calculations, ensuring consistency and a structured approach to managing complex data. The comment `!> \brief This module contains the types for the localized basis set component of the library` accurately describes its role. The note about `e_fermi` being an integer in `kohn_sham_types` is a potential point of attention, as Fermi energy is typically a real value.
