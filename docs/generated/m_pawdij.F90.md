# `m_pawdij.F90`

## Overview

The Fortran module `m_pawdij` is a core component of `libPAW` responsible for computing the PAW pseudopotential strength matrix elements, commonly denoted as \(D_{ij}\). These \(D_{ij}\) values define the strength of the non-local PAW operator, \(\hat{V}_{NL} = \sum_{ij} D_{ij} |\tilde{p}_i\rangle\langle\tilde{p}_j|\), where \(|\tilde{p}_i\rangle\) are the PAW projector functions.

The module provides a main driver routine `pawdij` that orchestrates the calculation of the total \(D_{ij}\) by summing various physical contributions. Each contribution is typically calculated by a dedicated public subroutine within this module (e.g., `pawdijhartree`, `pawdijxc`, `pawdijfock`). Additionally, the module includes routines for symmetrizing the computed \(D_{ij}\) matrices and for printing them for debugging or analysis.

## Key Components

### Main Driver Subroutine

-   **`pawdij(...)`**:
    -   **Purpose**: Computes the total \(D_{ij}\) matrix elements by summing various contributions. It can also handle first-order \(D_{ij}^{(1)}\) for Response Function (RF) calculations.
    -   **Inputs**: A comprehensive list of inputs including:
        -   Flags controlling calculation details (`cplex`, `ipert`, `pawxcdev`).
        -   System parameters (`natom`, `nspden`, `ntypat`, `ucvol`, `charge`).
        -   Grid information (`gprimd`, `nfft`, `nfftot`, `qphon`, `xred`).
        -   PAW data structures: `paw_an` (angular grid data), `paw_ij` (for storing \(D_{ij}\) results), `pawang` (angular mesh setup), `pawfgrtab` (fine grid data), `pawrad` (radial grids), `pawrhoij` (on-site density matrices), `pawtab` (PAW dataset info).
        -   Potentials: `vtrial` (total local potential), `vxc` (XC potential on grid).
        -   Optional arguments for MPI, DFT+U shifts, electron-positron DFT, hybrid functional parameters, nuclear dipole moments.
    -   **Output**: Updates the relevant `dij*` arrays within the `paw_ij` data structure (e.g., `paw_ij(iatom)%dij`, `paw_ij(iatom)%dijhartree`).
    -   **Logic**:
        1.  Performs consistency checks on input arguments.
        2.  Initializes MPI communication details if relevant.
        3.  Determines which \(D_{ij}\) components need to be computed based on `has_*` flags in `paw_ij` and input parameters.
        4.  Iterates over atoms assigned to the current process.
        5.  For each atom, it calls specialized subroutines (see below) to compute individual contributions to \(D_{ij}\) if they are needed and their prerequisites are met.
        6.  Sums these contributions into `paw_ij(iatom)%dij`.
        7.  Handles printing of computed \(D_{ij}\) matrices based on `pawprtvol`.

### Subroutines for Specific \(D_{ij}\) Contributions

-   **`pawdijhartree(dijhartree, qphase, nspden, pawrhoij, pawtab)`**:
    -   Computes the Hartree contribution \(D_{ij}^H = \sum_{kl} \langle ij | kl \rangle \rho_{kl}\), where \(\langle ij | kl \rangle\) are PAW four-center electron repulsion integrals (`pawtab%eijkl`) and \(\rho_{kl}\) is the on-site density matrix (`pawrhoij%rhoijp`).
-   **`pawdijfock(dijfock_vv, dijfock_cv, cplex_dij, qphase, hyb_mixing, hyb_mixing_sr, ndij, pawrhoij, pawtab)`**:
    -   Computes the Fock exact-exchange contribution \(D_{ij}^X\), separating valence-valence (`dijfock_vv`) and core-valence (`dijfock_cv`) parts. Uses `pawtab%eijkl` (or `eijkl_sr` for short-range) and `pawrhoij`.
-   **`pawdijxc(dijxc, ..., vxc1, vxct1, ..., vxctau1, vxcttau1)`**:
    -   Computes the exchange-correlation contribution \(D_{ij}^{XC}\) using densities and potentials represented on a real-space angular grid (`vxc1`, `vxct1`). Handles XC functionals up to meta-GGA (using `vxctau1`, `vxcttau1` if `usekden=1`).
    -   Formula: \(D_{ij}^{XC} = \langle \phi_i | V_{xc}^{AE} | \phi_j \rangle - \langle \tilde{\phi}_i | V_{xc}^{PS} | \tilde{\phi}_j \rangle - \int V_{xc}^{PS,loc} \hat{Q}_{ij} dr \).
-   **`pawdijxcm(dijxc, ..., vxc1, vxct1, ...)`**:
    -   Similar to `pawdijxc`, but for XC potentials represented as spherical harmonic moments (`vxc1`, `vxct1` on LM basis).
-   **`pawdijhat(dijhat, ..., Pot, ...)`**:
    -   Computes the "hat" contribution \(D_{ij}^{\hat{n}} = \int V_{eff}^{loc}(r) \hat{Q}_{ij}(r) dr\), where \(\hat{Q}_{ij}\) is the PAW compensation charge operator and \(V_{eff}^{loc}\) is the local effective potential on the real-space grid (`Pot`).
-   **`pawdijnd(dijnd, ..., nucdipmom, ...)`**:
    -   Computes the contribution from nuclear magnetic dipole moments: \(D_{ij}^{ND} = \alpha^2 \langle \phi_i | \frac{\vec{L} \cdot \vec{m}}{|r-R|^3} | \phi_j \rangle - \langle \tilde{\phi}_i | \dots | \tilde{\phi}_j \rangle\).
-   **`pawdijso(dijso, ..., vh1, vxc1)`**:
    -   Computes the spin-orbit coupling contribution \(D_{ij}^{SO} = \langle \phi_i | \frac{\alpha^2}{2} \frac{1}{r} \frac{dV}{dr} \vec{L} \cdot \vec{S} | \phi_j \rangle - \langle \tilde{\phi}_i | \dots | \tilde{\phi}_j \rangle\). Uses spherical components of \(V_H\) and \(V_{xc}\) from `paw_an_type` (`vh1`, `vxc1`).
-   **`pawdiju(dijpawu, ..., vpawu, ...)`**:
    -   Computes the DFT+U contribution \(D_{ij}^U\) using an on-site PAW+U potential `vpawu` (moments \(V_{m_1m_2}^{l\sigma}\)).
    -   \(D_{ij}^U = \sum_{m_1m_2} V_{m_1m_2}^{l\sigma} \langle \phi_i^l | Y_{lm_1}Y_{lm_2}^* | \phi_j^l \rangle\).
-   **`pawdiju_euijkl(diju, ..., pawrhoij, pawtab)`**:
    -   Alternative way to compute DFT+U contribution using \(D_{ij}^U = \sum_{kl} \rho_{kl} E_{ijkl}^U\), where \(E_{ijkl}^U\) are effective U-interaction parameters (`pawtab%euijkl`).
-   **`pawdijexxc(dijexxc, ..., vpawx, vxc_ex)`**:
    -   Computes contribution for local exact exchange (EXX) within PAW: \(D_{ij}^{EXXC} = \alpha ( \langle \phi_i | V_{Fock}^{corr} | \phi_j \rangle - \langle \phi_i | V_{xc}^{corr} | \phi_j \rangle ) \). Uses `vpawx` (EXX potential moments) and `vxc_ex` (on-site XC for correlated electrons).
-   **`pawdijfr(...)`**:
    -   Computes the "frozen" part of the first-order \(D_{ij}^{(1)}\) for response function calculations. This includes terms from the first-order change in the local potential \(V_{loc}^{(1)}\) and the first-order change in the compensation charge \(\hat{Q}_{ij}^{(1)}\).

### Utility Subroutines

-   **`pawpupot(..., noccmmp, nocctot, ..., vpawu)`**:
    -   Computes the on-site DFT+U potential moments `vpawu` from the occupation matrix `noccmmp` and total occupation `nocctot`.
-   **`pawxpot(..., pawrhoij, pawtab, vpawx)`**:
    -   Computes the on-site local exact-exchange potential moments `vpawx`.
-   **`symdij(...)`**: Symmetrizes a specific component of the \(D_{ij}\) matrix (selected by `option_dij`) using crystal symmetries.
-   **`symdij_all(...)`**: Calls `symdij` for all available \(D_{ij}\) components that have been computed.
-   **`pawdij_gather(...)`**: Gathers distributed `coeff2_type` arrays (used for storing \(D_{ij}\)-like data before it's in `paw_ij_type`) across MPI processes.
-   **`pawdij_print_dij(...)`**: Prints a \(D_{ij}\) matrix (passed as a raw array) with formatting, used by `paw_ij_print`.

## Important Variables/Constants

-   `paw_ij_type`: Defined in `m_paw_ij`, this is the primary data structure holding the results. Its `has_*` flags determine which components are computed by `pawdij`.
-   `pawtab_type`: From `m_pawtab`, provides essential PAW dataset parameters (projectors, partial waves, \(D_{ij}^0\), \(E_{ijkl}\), \(E_{ijkl}^U\), etc.).
-   `paw_an_type`: From `m_paw_an`, provides on-site potentials (Hartree, XC) on angular grids or as LM-moments, needed for some \(D_{ij}\) contributions like spin-orbit.
-   `pawrhoij_type`: From `m_pawrhoij`, provides on-site occupation matrices \(\rho_{ij}\), crucial for Hartree, Fock, and DFT+U contributions.
-   `pawfgrtab_type`, `pawang_type`: Provide grid and angular integration data.

## Usage Examples

The main routine `pawdij` is typically called by higher-level DFT or response function codes after necessary prerequisites (like on-site density matrices \(\rho_{ij}\) and on-site potentials \(V_{xc}^{on-site}\)) have been computed.

```fortran
! Conceptual call within a larger program
! Assuming all input PAW data structures (paw_an_data, paw_ij_data, etc.)
! and potentials (eff_potential, xc_potential_grid) are already prepared.

call pawdij(cplex=current_cplex, enunit=output_energy_unit, gprimd=reciprocal_vectors, &
            ipert=current_perturbation_index, my_natom=num_local_atoms, natom=total_atoms, &
            nfft=num_fft_points_local, nfftot=num_fft_points_total, nspden=num_spin_density_components, &
            ntypat=num_atom_types, paw_an=paw_an_data, paw_ij=paw_ij_data, &
            pawang=angular_grid_setup, pawfgrtab=fine_grid_data_atomic, &
            pawprtvol=print_verbosity, pawrad=radial_grids_atomic, &
            pawrhoij=onsite_density_matrices, pawspnorb=spin_orbit_flag, &
            pawtab=paw_datasets, pawxcdev=xc_representation_flag, &
            qphon=phonon_q_vector, spnorbscl=so_scaling_factor, ucvol=unit_cell_volume, &
            charge=total_charge, vtrial=eff_potential, vxc=xc_potential_grid, &
            xred=atom_positions_reduced, comm_atom=mpi_atom_communicator)

! After this call, paw_ij_data will be updated with the computed D_ij components.
! Symmetrization might be called next:
! call symdij_all(...)
```

## Dependencies and Interactions

-   **`m_libpaw_defs`, `m_libpaw_tools`, `m_libpaw_mpi`, `m_libpaw_memory`**: For basic definitions, utilities, MPI, and memory.
-   **`m_paral_atom`**: For atom distribution logic.
-   **`m_paw_io`**: `pawio_print_ij` is used by `pawdij_print_dij`.
-   **`m_pawang`, `m_pawrad`, `m_pawtab`**: Define the core PAW data structures used as input.
-   **`m_paw_an`, `m_paw_ij`, `m_pawfgrtab`, `m_pawrhoij`**: Define other PAW data structures that are either input to or output from this module.
-   **`m_paw_finegrid`**: Provides `pawgylm` (for \(g_l Y_{lm}\) on fine grid) and `pawexpiqr` (for \(e^{i\vec{q}\cdot\vec{r}}\) factors), used in `pawdijhat` and potentially `pawdijfr`.
-   **`m_paw_sphharm`**: Provides `slxyzs` (for \(\langle S_{l'm'} | L_z | S_{lm} \rangle\), though not directly visible in `pawdijso`'s direct calls, it's used by underlying routines for SO coupling).

This module is central to the PAW method implementation, as it assembles the non-local operator strengths \(D_{ij}\) from their various physical origins. The modular design where each contribution is handled by a specific subroutine aids in clarity and maintainability.
