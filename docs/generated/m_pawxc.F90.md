# `m_pawxc.F90`

## Overview

The Fortran module `m_pawxc` is a crucial part of `libPAW` that handles Exchange-Correlation (XC) functional calculations within the PAW atomic augmentation spheres. It provides interfaces to compute XC energy densities, XC potentials, and XC kernels (second derivatives of the XC energy). The module is designed to work with densities represented either on a real-space angular grid (r, \(\theta\), \(\phi\)) or as spherical harmonic (LM) moments.

Key functionalities include:
-   Calculating XC quantities for ground-state DFT.
-   Handling electron-positron correlation functionals.
-   Calculating first-order changes in XC potentials and second-order XC energy contributions for Density Functional Perturbation Theory (DFPT).
-   Interfacing with external XC libraries (like LibXC via `m_libpaw_libxc`) or host-code specific XC routines (like those in ABINIT).
-   Managing different spin polarizations (non-polarized, collinear, non-collinear).
-   Providing utility functions to query characteristics of XC functionals (e.g., LDA/GGA/MGGA level, need for density gradients/Laplacian/kinetic energy density).

## Key Components

### Main Public XC Calculation Subroutines

-   **`pawxc(corexc, enxc, enxcdc, ..., rhor, ..., vxc, ...)`**:
    -   **Purpose**: Computes XC energy densities, potentials, and optionally kernels, using density input (`rhor`) provided on a full spherical grid (radial points `nrad`, angular points from `pawang%angl_size`).
    -   **Inputs**: Core density (`corexc`), input density moments (`rhor`), PAW setup (`pawang`, `pawrad`), XC choice (`ixc`, `xclevel`, `hyb_mixing`), options for calculation (`option`, `usecore`, `usexcnhat`), and various tolerances. For meta-GGAs, kinetic energy densities (`taur`, `coretau`) can be provided.
    -   **Outputs**: XC energy density (`exc` per point, then integrated to `enxc`), double-counting energy (`enxcdc`), XC potential (`vxc` on the angular grid), and optionally XC kernels (`kxc`, `k3xc`).
-   **`pawxcm(corexc, enxc, enxcdc, ..., rhor, ..., vxc, ...)`**:
    -   **Purpose**: Similar to `pawxc`, but uses density input (`rhor`) provided as spherical harmonic (LM) moments for each radial point. This is typically controlled by `pawxcdev > 0`.
    -   **Logic**: Internally, it reconstructs spherical densities from moments, calls `pawxcsph` to get spherical XC quantities, and then projects these back to LM moments if needed for output `vxc` or `kxc`.
-   **`pawxcpositron(...)` / `pawxcmpositron(...)`**:
    -   **Purpose**: Specialized versions for electron-positron correlation functionals, for grid-based and moment-based densities respectively.
-   **`pawxc_dfpt(...)` / `pawxcm_dfpt(...)`**:
    -   **Purpose**: Compute first-order changes in XC potential (\(V_{xc}^{(1)}\)) and contributions to the second-order XC energy (\(E_{xc}^{(2)}\)) for DFPT calculations, for grid-based and moment-based densities/potentials respectively.

### Utility Functions for XC Functional Properties

-   **`pawxc_get_nkxc(nkxc, nspden, xclevel)`**: Subroutine to determine the size (`nkxc`) of the XC kernel array based on spin polarization (`nspden`) and XC level (`xclevel`).
-   **`pawxc_get_xclevel(ixc)`**: Integer function, returns the level of an XC functional (1 for LDA, 2 for GGA/MGGA, 3 for TDDFT kernels) given its ABINIT `ixc` identifier or LibXC ID.
-   **`pawxc_get_usekden(ixc)`**: Integer function (0 or 1), returns whether a functional requires kinetic energy density (\(\tau\)).
-   **`pawxc_get_uselaplacian(ixc)`**: Integer function (0 or 1), returns whether a functional requires the Laplacian of the density (\(\nabla^2 \rho\)).
-   **`pawxc_is_tb09(ixc)`**: Logical function, checks if the functional is the Tran-Blaha 09 (modified Becke-Johnson) meta-GGA.

### Private Core XC Calculation Routines

-   **`pawxcsph(exc, exexch, ..., rho_updn, vxc, ...)`**:
    -   **Purpose**: Computes XC energy density, potential, and kernel for a spherically symmetric (or \(m=0\) component of) density provided as up and down spin components (`rho_updn`). This is the workhorse for `pawxcm`.
    -   **Logic**: Calls `pawxc_drivexc_wrapper` to interface with the actual XC library (ABINIT's or LibXC). Handles LDA and GGA.
-   **`pawxcsphpositron(...)`**: Similar to `pawxcsph` but for electron-positron correlation.
-   **`pawxcsph_dfpt(...)`**: Similar to `pawxcsph` but for DFPT quantities.

### Wrapper Subroutines for Host Code/LibXC

-   **`pawxc_drivexc_wrapper(...)`**: Acts as a direct wrapper to either ABINIT's `drivexc` routine or `libxc_functionals_getvxc` (from `m_libpaw_libxc_funcs` which calls C wrappers for LibXC). This is where the actual call to the XC functional library happens.
-   **`pawxc_mkdenpos_wrapper(...)`**: Wrapper for `mkdenpos` (ensures density is above a threshold).
-   **`pawxc_xcmult_wrapper(...)`**: Wrapper for `xcmult` (multiplies density gradients by \(dE_{xc}/d|\nabla\rho|\) for GGA).
-   **`pawxc_size_dvxc_wrapper(...)`**: Wrapper for `size_dvxc` (gets dimensions of XC derivative arrays).
-   **`pawxc_xcpositron_wrapper(...)`**: Wrapper for `xcpositron`.

### Non-Collinear Magnetism Helpers

-   **`pawxc_rotate_mag(rho_in, rho_out, mag, ...)`**: Projects a non-collinear density (input as total density + magnetization vector) onto a local magnetization direction `mag` to get effective up/down spin densities for collinear XC evaluation.
-   **`pawxc_rotate_back_mag(vxc_in, vxc_out, mag, ...)`**: Rotates a collinear XC potential (derived from projected densities) back into the non-collinear representation.
-   **`pawxc_rotate_back_mag_dfpt(...)`**: Similar for first-order potentials in DFPT.

## Important Variables/Constants

-   `rho_min`: `real(dp), parameter :: tol14`. A threshold below which density is considered zero, used by `pawxc_mkdenpos_wrapper`.
-   `ixc`: Integer input selecting the XC functional. Negative values usually indicate LibXC functionals.
-   `xclevel`: Integer indicating functional type (LDA, GGA, etc.).
-   `pawxcdev`: Integer option for how XC is treated on angular grids (0 for full grid calculation, 1 or 2 for LM moment expansion of density/potential).

## Usage Examples

The main routines `pawxc` or `pawxcm` are called by higher-level PAW modules (like `m_pawdij` when computing \(D_{ij}^{XC}\) or routines that calculate on-site XC energy).

```fortran
module m_example_pawxc_usage
    use m_pawxc
    use m_libpaw_defs, only: dp
    use m_pawang, only: pawang_type
    use m_pawrad, only: pawrad_type
    implicit none

    subroutine calculate_onsite_xc(radial_grid, angular_setup, &
                                   ixc_functional_id, xc_level_id, nspden_val, &
                                   density_moments, core_density_radial, &
                                   xc_energy_density, xc_potential_moments)
        type(pawrad_type), intent(in) :: radial_grid
        type(pawang_type), intent(in) :: angular_setup
        integer, intent(in) :: ixc_functional_id, xc_level_id, nspden_val
        real(dp), intent(in) :: density_moments(radial_grid%mesh_size, angular_setup%l_size_max**2, nspden_val)
        real(dp), intent(in) :: core_density_radial(radial_grid%mesh_size)

        real(dp), intent(out) :: xc_energy_density(radial_grid%mesh_size, angular_setup%angl_size) ! Example output structure
        real(dp), intent(out) :: xc_potential_moments(radial_grid%mesh_size, angular_setup%l_size_max**2, nspden_val)

        real(dp) :: enxc_total, enxcdc_total
        logical, dimension(angular_setup%l_size_max**2) :: lm_selection_map
        ! ... (Initialize lm_selection_map, other parameters like nhat, kxc, etc.) ...
        lm_selection_map = .true. ! Assume all moments are used for simplicity

        ! Using pawxcm (LM moment based)
        call pawxcm(corexc=core_density_radial, enxc=enxc_total, enxcdc=enxcdc_total, &
                    hyb_mixing=0.0_dp, ixc=ixc_functional_id, kxc=dummy_kxc, & ! dummy_kxc needs proper dimensioning if nkxc > 0
                    lm_size=size(lm_selection_map), lmselect=lm_selection_map, &
                    nhat=dummy_nhat, nkxc=0, nk3xc=0, non_magnetic_xc=.false., & ! dummy_nhat similar
                    nrad=radial_grid%mesh_size, nspden=nspden_val, option=0, & ! option 0: E_xc and V_xc
                    pawang=angular_setup, pawrad=radial_grid, pawxcdev=1, & ! pawxcdev=1 for moments
                    rhor=density_moments, usecore=1, usexcnhat=0, &
                    vxc=xc_potential_moments, xclevel=xc_level_id, xc_denpos=1e-14_dp)

        ! xc_potential_moments now contains Vxc(r, L, M)
        ! enxc_total contains the integrated on-site XC energy for this atom

contains
        real(dp) :: dummy_kxc(1,1,0), dummy_nhat(1,1,1) ! Minimal dummy arrays
    end subroutine calculate_onsite_xc

end module m_example_pawxc_usage
```

## Dependencies and Interactions

-   **`m_libpaw_defs`, `m_libpaw_tools`, `m_libpaw_memory`**: For basic definitions, utilities, and memory.
-   **`m_libpaw_libxc`**: Crucial for calculations using LibXC functionals (`ixc < 0`). It provides the Fortran interface to the C wrappers around LibXC.
-   **ABINIT specific modules (if `HAVE_LIBPAW_ABINIT` is defined)**:
    -   `m_xcpositron`: For electron-positron correlation.
    -   `m_drivexc`: ABINIT's main XC driver.
    -   `m_xc_noncoll`: For non-collinear magnetism density/potential rotations.
-   **`m_pawang`**: Defines `pawang_type`, providing angular grid information and \(S_{lm}\) values used in `pawxc`.
-   **`m_pawrad`**: Defines `pawrad_type`, providing radial grid information and routines like `nderiv_gen`, `simp_gen`.

This module serves as the primary engine for evaluating XC contributions within PAW spheres, offering flexibility in terms of XC functional choice (via LibXC or host code), density representation (grid or moments), and physical context (ground-state, DFPT, e-e, e-p). The wrappers ensure that `libPAW` can seamlessly integrate with different underlying XC calculation libraries.
