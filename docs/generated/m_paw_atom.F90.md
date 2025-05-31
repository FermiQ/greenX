# `m_paw_atom.F90`

## Overview

The Fortran module `m_paw_atom` provides subroutines for performing various calculations related to individual atoms within the Projector Augmented Wave (PAW) formalism. These routines are often referred to as "atompaw" operations and deal with quantities like PAW shape functions, Hartree potentials from compensation charges, and the construction of the \(D_{ij}\) matrix elements (both the frozen core part \(D_{ij}^0\) and the kinetic energy contribution \(K_{ij}\)).

This module relies on radial grid information (`pawrad_type`) and tabulated PAW setup data (`pawtab_type`) from other modules.

## Key Components

### Public Subroutines

-   **`atompaw_shpfun(ll, mesh, norm, pawtab, shapefunc)`**:
    -   **Purpose**: Computes the PAW shape function \(g_l(r)\) for a given angular momentum `ll`.
    -   **Inputs**:
        -   `ll`: Integer, the angular momentum quantum number.
        -   `mesh`: `type(pawrad_type)`, contains radial grid information.
        -   `pawtab`: `type(pawtab_type)`, contains PAW setup data, including shape function type and parameters (`rshp`, `shape_sigma`, `shape_lambda`, `shape_alpha`, `shape_q`).
    -   **Output**:
        -   `norm`: `real(dp)`, the normalization factor computed for the shape function (relevant for types -1, 1, 2).
        -   `shapefunc(:)`: `real(dp)`, array (intent inout), filled with the values of the shape function on the radial grid up to `pawtab%rshp`. If `pawtab%shape_type == -1`, this array is expected to contain the numerical shape function on input, and it's then normalized.
    -   **Logic**:
        -   Selects the shape function calculation based on `pawtab%shape_type`:
            -   Type -1: Uses a numerical shape function provided in `pawtab%shapefunc(:,ll+1)`.
            -   Type 1: \(g(r) = k(r) \cdot r^l\), where \(k(r) = \exp(-(r/\sigma)^\lambda)\).
            -   Type 2: \(g(r) = k(r) \cdot r^l\), where \(k(r) = [\sin(\pi r/r_{shp})/(\pi r/r_{shp})]^2\).
            -   Type 3: \(g(r) = \alpha_1 j_l(q_1 r) + \alpha_2 j_l(q_2 r)\) (Bessel function based).
        -   Normalizes the shape function (for types -1, 1, 2) such that \(\int_0^{r_{shp}} g_l(r) r^{l+2} dr = 1\). For type 3, `norm` is set to 1.

-   **`atompaw_shapebes(al, ql, ll, rc)`**:
    -   **Purpose**: Determines the parameters `al(2)` and `ql(2)` for a Bessel-type shape function: \(g_l(r) = al_1 j_l(ql_1 r) + al_2 j_l(ql_2 r)\).
    -   The parameters are found such that the shape function and its first two derivatives are zero at \(r=rc\), and its integral \(\int_0^{rc} g_l(r) r^{l+2} dr = 1\).
    -   **Inputs**:
        -   `ll`: Integer, angular momentum.
        -   `rc`: `real(dp)`, cut-off radius.
    -   **Outputs**:
        -   `al(2)`: `real(dp)`, the \(\alpha\) coefficients.
        -   `ql(2)`: `real(dp)`, the \(q\) factors (wave numbers).
    -   **Logic**: Uses `paw_solvbes` to find roots for \(q\), then solves a 2x2 linear system for \(\alpha\) coefficients based on boundary conditions and normalization.

-   **`atompaw_vhnzc(ncore, radmesh_core, vhnzc, znucl)`**:
    -   **Purpose**: Computes the Hartree potential \(V_H(n_{Zc})\) due to the PAW compensation charge density \(n_{Zc}(r)\). The compensation charge density is related to the core density `ncore` and the total nuclear charge `znucl`.
    -   **Inputs**:
        -   `ncore(:)`: `real(dp)`, atomic core density on `radmesh_core`.
        -   `radmesh_core`: `type(pawrad_type)`, radial grid for `ncore`.
        -   `znucl`: `real(dp)`, effective nuclear charge for the PAW sphere (valence + core charge considered inside PAW).
    -   **Output**:
        -   `vhnzc(:)`: `real(dp)`, the Hartree potential \(V_H(n_{Zc})\) on `radmesh_core`.
    -   **Logic**:
        1.  Forms \(r^2 \rho_{comp}(r)\) where \(\rho_{comp}\) is related to `ncore`.
        2.  Solves Poisson's equation for this charge using `poisson` subroutine.
        3.  Adjusts the potential: \((V_{Poisson} - z_{nucl}) / r\).
        4.  Ensures correct behavior at \(r=0\) using `pawrad_deducer0`.

-   **`atompaw_dij0(...)`**:
    -   **Purpose**: Computes the "frozen" core part of the PAW strength matrix elements, \(D_{ij}^0\). This typically includes kinetic energy terms and interactions with the frozen core density and the local part of the pseudopotential.
    -   \(D_{ij}^0 = K_{ij} + \langle \phi_i | V_H(n_{Zc}) | \phi_j \rangle - \langle \tilde{\phi}_i | V_H(\tilde{n}_{Zc}) | \tilde{\phi}_j \rangle - \int V_H(\tilde{n}_{Zc}(r)) \hat{Q}_{ij}(r) dr + \langle \phi_i | V_{extra} | \phi_j \rangle\). (The exact formula depends on the terms included).
    -   **Inputs**:
        -   `indlmn`, `lmnmax`: Describe mapping of (l,m,n) projector indices.
        -   `kij`: `real(dp)`, kinetic part of \(D_{ij}\).
        -   `ncore`: Core density.
        -   `opt_init`: Flag, if 0, shape function for Q compensation term is computed; if 1, uses `pawtab%qijl`.
        -   `pawtab`: PAW setup data (wavefunctions \(\phi, \tilde{\phi}\), `qijl`).
        -   `radmesh`, `radmesh_core`, `radmesh_vloc`: Radial grid information.
        -   `vhtnzc`: Local potential \(V_H(\tilde{n}_{Zc})\).
        -   `znucl`: Nuclear charge for \(V_H(n_{Zc})\).
    -   **Output**:
        -   `pawtab%dij0(:)`: `real(dp)`, filled with the computed \(D_{ij}^0\) values.
    -   **Logic**:
        1.  Initializes `pawtab%dij0` with `kij`.
        2.  Computes \(\langle \phi_i | V_H(n_{Zc}) | \phi_j \rangle\) term:
            -   Calculates \(V_H(n_{Zc})\) using `atompaw_vhnzc`.
            -   Integrates \( \int \phi_i(r) \phi_j(r) V_H(n_{Zc})(r) r^2 dr \).
        3.  Computes \(-\langle \tilde{\phi}_i | V_H(\tilde{n}_{Zc}) | \tilde{\phi}_j \rangle\) term:
            -   Interpolates `vhtnzc` onto `radmesh` if grids differ.
            -   Integrates \( -\int \tilde{\phi}_i(r) \tilde{\phi}_j(r) V_H(\tilde{n}_{Zc})(r) r^2 dr \).
        4.  Adds \(\langle \phi_i | V_{minus\_half} | \phi_j \rangle\) term if `pawtab%has_vminushalf == 1`.
        5.  Computes \(-\int V_H(\tilde{n}_{Zc}(r)) \hat{Q}_{ij}(r) dr\) term:
            -   \(\hat{Q}_{ij}(r) = (\phi_i \phi_j - \tilde{\phi}_i \tilde{\phi}_j)\) for \(l=0\) (isotropic part) or uses precomputed `pawtab%qijl` if `opt_init == 1`.
            -   The integral involves the shape function \(g_0(r)\) for the compensation charge: \(-\left( \int V_H(\tilde{n}_{Zc}) g_0(r) r^2 dr \right) \left( \int (\phi_i \phi_j - \tilde{\phi}_i \tilde{\phi}_j) r^2 dr \right)\).

-   **`atompaw_kij(..., kij)`**:
    -   **Purpose**: Deduces the kinetic energy contribution \(K_{ij}\) to the \(D_{ij}\) matrix, effectively by subtracting the potential energy terms from \(D_{ij}^0\).
    -   \(K_{ij} = D_{ij}^0 - \langle \phi_i | V_H(n_{Zc}) | \phi_j \rangle + \langle \tilde{\phi}_i | V_H(\tilde{n}_{Zc}) | \tilde{\phi}_j \rangle + \int V_H(\tilde{n}_{Zc}(r)) \hat{Q}_{ij}(r) dr - \langle \phi_i | V_{extra} | \phi_j \rangle\).
    -   **Inputs**: Similar to `atompaw_dij0`, including `pawtab%dij0` (which are the \(D_{ij}^0\) values). `opt_vhnzc` controls if the \(V_H(n_{Zc})\) term is included.
    -   **Output**:
        -   `kij(:)`: `real(dp)`, filled with the computed kinetic energy components.
    -   **Logic**: Essentially reverses the potential energy calculations performed in `atompaw_dij0` and subtracts them from `pawtab%dij0`.

## Important Variables/Constants

-   `pawrad_type`: Derived type for radial grid information.
-   `pawtab_type`: Derived type for PAW setup data (atomic wavefunctions, projector info, shape function parameters, etc.).
-   `shape_type`: Integer in `pawtab_type` indicating the type of shape function to use.
-   `phi`, `tphi`: Arrays in `pawtab_type` for all-electron and pseudo partial waves.
-   `dij0`: Array in `pawtab_type` storing the frozen \(D_{ij}^0\) matrix elements.

## Usage Examples

These subroutines are typically called during the setup phase of a PAW calculation, often when processing PAW dataset files (e.g., XML pseudopotential files) to prepare the necessary atomic quantities.

```fortran
module m_paw_setup_example
    use m_libpaw_defs, only: dp
    use m_paw_atom
    use m_pawrad, only: pawrad_type ! Assume these are initialized
    use m_pawtab, only: pawtab_type   ! Assume this is initialized

    implicit none

    subroutine setup_paw_data_for_atom(atomic_l, radial_grid, paw_setup_data, &
                                       core_density, vloc_potential, nuclear_charge, &
                                       kinetic_contrib, frozen_dij)
        integer, intent(in) :: atomic_l
        type(pawrad_type), intent(in) :: radial_grid
        type(pawtab_type), intent(inout) :: paw_setup_data ! dij0 will be modified
        real(dp), intent(in) :: core_density(:)
        real(dp), intent(in) :: vloc_potential(:) ! This is V_H(tilde_n_Zc)
        real(dp), intent(in) :: nuclear_charge    ! This is Z_nucl for V_H(n_Zc)
        real(dp), intent(in) :: kinetic_contrib(paw_setup_data%lmn2_size)
        real(dp), intent(out) :: frozen_dij(paw_setup_data%lmn2_size) ! To store computed Kij

        real(dp) :: shape_norm
        real(dp), dimension(radial_grid%mesh_size) :: g_shape_function

        ! For atompaw_dij0/kij, need indlmn, lmnmax, radmesh_core, radmesh_vloc
        ! These would be part of a fuller PAW setup. Simplified here.
        integer, dimension(6,paw_setup_data%lmn_size) :: example_indlmn
        integer :: example_lmnmax = paw_setup_data%lmn_size
        ! example_indlmn would map (l,m,n) indices for projectors

        ! 1. Compute a shape function
        call atompaw_shpfun(atomic_l, radial_grid, shape_norm, paw_setup_data, g_shape_function)
        ! paw_setup_data%shapefunc might be updated if it was type -1 and normalized

        ! 2. Compute D_ij^0 (modifies paw_setup_data%dij0)
        ! Assuming radmesh_core and radmesh_vloc are same as radial_grid for simplicity
        call atompaw_dij0(example_indlmn, kinetic_contrib, example_lmnmax, core_density, 0, &
                          paw_setup_data, radial_grid, radial_grid, radial_grid, &
                          vloc_potential, nuclear_charge)

        ! 3. Compute K_ij from D_ij^0
        call atompaw_kij(example_indlmn, frozen_dij, example_lmnmax, core_density, 0, 1, &
                         paw_setup_data, radial_grid, radial_grid, radial_grid, &
                         vloc_potential, nuclear_charge)

    end subroutine setup_paw_data_for_atom

end module m_paw_setup_example
```

## Dependencies and Interactions

-   **`m_libpaw_defs` (`USE_DEFS`)**: For `dp` and mathematical constants like `pi`, `four_pi`.
-   **`m_libpaw_tools` (`USE_MSG_HANDLING`)**: For error reporting (`LIBPAW_BUG`).
-   **`m_libpaw_memory` (`USE_MEMORY_PROFILING`)**: For memory allocation macros (`LIBPAW_ALLOCATE`, `LIBPAW_DEALLOCATE`).
-   **`m_paw_numeric`**: Provides numerical routines:
    -   `paw_jbessel`: For spherical Bessel functions (used in `atompaw_shpfun`, `atompaw_shapebes`).
    -   `paw_solvbes`: For finding roots related to Bessel functions (used in `atompaw_shapebes`).
    -   `paw_spline`, `paw_splint`: For spline interpolation (used in `atompaw_dij0`, `atompaw_kij` if grids differ).
-   **`m_pawrad`**: Defines `pawrad_type` and provides radial integration (`simp_gen`), Poisson solver (`poisson`), and utility functions like `pawrad_deducer0`, `bound_deriv`, `pawrad_ifromr`.
-   **`m_pawtab`**: Defines `pawtab_type` which contains essential input data like PAW partial waves (\(\phi, \tilde{\phi}\)), projector information, shape function parameters, and stores the output `dij0`.

This module encapsulates key atomic-centered calculations required for setting up the PAW method, particularly the construction of the \(D_{ij}\) matrix elements which are fundamental to the PAW transformation.
