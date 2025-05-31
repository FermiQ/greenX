# `m_pawtab.F90`

## Overview

The Fortran module `m_pawtab` defines and manages crucial derived data types for Projector Augmented Wave (PAW) calculations, primarily the `pawtab_type`. This extensive data structure serves as a container for all tabulated atomic data read or derived from a PAW pseudopotential (PSP) file for a single atom type. It holds a wide array of information, including radial wavefunctions, projectors, densities, potentials, pre-calculated matrix elements, and parameters defining the PAW setup.

The module also defines `wvlpaw_type` and its sub-type `wvlpaw_rholoc_type`, which are specifically for PAW implementations using wavelet basis sets, storing Gaussian fits to projectors and local density/potential information.

Public procedures are provided for lifecycle management (allocation, deallocation, nullification), setting status flags, printing, and MPI broadcasting of `pawtab_type` objects.

## Key Components

### Derived Type: `wvlpaw_rholoc_type`

-   **Purpose**: Stores local density \(\rho(r)\) and local potential \(V(r)\) and their second derivatives on a radial mesh, primarily for wavelet-based PAW.
-   **Members**:
    -   `msz`: Integer, mesh size.
    -   `d(msz, 4)`: `real(dp)`, allocatable. Stores:
        -   `d(:,1)`: Local density \(\rho(r)\).
        -   `d(:,2)`: Second derivative of local density \(\rho''(r)\).
        -   `d(:,3)`: Local potential \(V(r)\).
        -   `d(:,4)`: Second derivative of local potential \(V''(r)\).
    -   `rad(msz)`: `real(dp)`, allocatable, radial mesh points.
-   **Procedures**: `wvlpaw_rholoc_free`, `wvlpaw_rholoc_nullify`.

### Derived Type: `wvlpaw_type`

-   **Purpose**: Stores data specific to PAW calculations using wavelet basis sets, particularly Gaussian fits of projectors.
-   **Members**:
    -   `npspcode_init_guess`: Integer.
    -   `ptotgau`: Integer, total number of complex Gaussians for fitting projectors.
    -   `pngau(:)`: Integer, allocatable, number of complex Gaussians per basis element (projector).
    -   `parg(:,:)`: `real(dp)`, allocatable (`2, ptotgau`), exponents of Gaussian functions (real and imaginary parts).
    -   `pfac(:,:)`: `real(dp)`, allocatable (`2, ptotgau`), coefficients of Gaussian functions (real and imaginary parts).
    -   `rholoc`: `type(wvlpaw_rholoc_type)`, local density and potential data.
-   **Procedures**: `wvlpaw_allocate`, `wvlpaw_free`, `wvlpaw_nullify`.

### Derived Type: `pawtab_type`

-   **Purpose**: A comprehensive container for all tabulated PAW data for a single atom type. This is a central data structure in `libPAW`.
-   **Integer Scalars (Metadata, Sizes, Flags)**:
    -   `basis_size`: Number of \((l,n)\) PAW basis functions (partial waves/projectors).
    -   `has_*`: Numerous flags (e.g., `has_coretau`, `has_fock`, `has_kij`, `has_tproj`, `has_nabla`, `has_wvl`) indicating if optional arrays are allocated (1) or computed and stored (2).
    -   `ij_proj`: Number of \((i,j)\) pairs for orbitals on which DFT+U or local exact exchange is applied.
    -   `ij_size`: `basis_size*(basis_size+1)/2`.
    -   `lcut_size`: \(l_{max\_Gaunt}+1\), may be reduced by input parameter `pawlcutd`.
    -   `l_size`: \(2 \cdot l_{max\_psp} - 1\), related to max \(L\) in Gaunt coefficients.
    -   `lexexch`, `lpawu`: Angular momentum \(l\) for local exact exchange and DFT+U.
    -   `lmn_size`, `lmn2_size`: Number of \((l,m,n)\) channels and unique pairs.
    -   `lmnmix_sz`: Number of \(\rho_{ij}\) elements in SCF mixing.
    -   `mesh_size`, `partialwave_mesh_size`, `core_mesh_size`, etc.: Sizes of various radial meshes.
    -   `mqgrid`, `mqgrid_shp`: Number of q-points for reciprocal space data.
    -   `shape_lambda`, `shape_type`: Parameters for the PAW compensation charge shape function.
    -   `use*`: Flags like `usetcore` (use pseudo-core), `usexcnhat` (use \(\hat{n}\) in XC), `usepawu`, `useexexch`, `usespnorb`.
-   **Real(`dp`) Scalars**:
    -   `beta`: Integral \(\int (V_H[n_{Zc}] - V_H[\tilde{n}_{Zc}]) 4\pi r^2 dr\).
    -   `dncdq0`, `d2ncdq0`, `dnvdq0`, `dtaucdq0`: Derivatives of Fourier transformed core/valence densities/kinetic energy densities at \(q=0\).
    -   `ex_cc`: Core-core exact exchange energy.
    -   `exccore`: XC energy of the core density.
    -   `exchmix`: Mixing factor for exact exchange in hybrid functionals.
    -   `jpawu`, `upawu`: Hubbard J and U parameters.
    -   `lamb_shielding`: Lamb shielding correction for NMR.
    -   `rpaw`, `rshp`, `rcore`, `rcoretau`: Various PAW radii.
    -   `shape_sigma`: Parameter for Gaussian shape functions.
-   **Objects**:
    -   `wvl`: `type(wvlpaw_type), pointer`, data for wavelet basis.
-   **Integer Arrays (Allocatable)**:
    -   `indklmn(8, lmn2_size)`: Maps packed \((ij)\) index to various \((l,m,n)\) related indices.
    -   `indlmn(6, lmn_size)`: Maps full \((lmn)\) index to individual \(l,m,n\), angular momentum index `lm`, basis function index `ln`, and spin index `s`.
    -   `klmntomn(4, lmn2_size)`: Maps packed \((ij)\) index to \(m_i, m_j, n_i, n_j\).
    -   `kmix(lmnmix_sz)`: Indices of \(\rho_{ij}\) elements for SCF mixing.
    -   `lnproju(nproju)`: `ln` indices for orbitals involved in DFT+U/EXX.
    -   `orbitals(basis_size)`: \(l\) quantum number for each \((l,n)\) basis function.
-   **Real(`dp`) Arrays (Allocatable)**: A large number of arrays storing:
    -   Radial functions: `coredens`, `coretau` (AE core kin. energy density), `phi` (AE partial waves), `rhoij0` (initial \(\rho_{ij}\)), `shapefunc` (radial part of \(g_l S_{lm}\)), `tcoredens` (pseudo-core density), `tcoretau` (pseudo-core kin. energy density), `tphi` (pseudo partial waves), `tproj` (projectors \(\tilde{p}\)), `vhtnzc` (\(V_H[\tilde{n}_{Zc}]\)), `vhnzc` (\(V_H[n_{Zc}]\)), `vminushalf` (LDA-1/2 potential).
    -   Reciprocal space functions: `shapefncg` (\(g_l(q)\)), `tcorespl` (\(\tilde{n}_{core}(q)\)), `tcoretauspl` (\(\tilde{\tau}_{core}(q)\)), `tvalespl` (\(\tilde{n}_{val}(q)\)).
    -   Matrix elements: `dij0` (\(D_{ij}^0\)), `dltij` (normalization \(1/(1+\delta_{ij})\) for sums), `eijkl` & `eijkl_sr` (Hartree/Fock kernels \(\langle ij|kl \rangle\)), `euijkl` & `euij_fll` (DFT+U kernels), `fk` (Slater integrals \(F^k\)), `gammaij` (background \(D_{ij}\)), `gnorm` (shape function norms), `ex_cvij` (core-valence exact exchange), `kij` (\(K_{ij}\)), `nabla_ij` & `nabla_im_ij` (\(\langle \phi_i|\nabla|\phi_j\rangle - \dots\)), `sij` (\(S_{ij} = \langle \phi_i|\phi_j\rangle - \langle \tilde{\phi}_i|\tilde{\phi}_j\rangle\)).
    -   Auxiliary radial products: `phiphj` (\(\phi_i \phi_j\)), `phiphjint` (\(\int \phi_i \phi_j r^2 dr\)), `qijl` (moments \(\int Q_{ij}^{LM} r^L dr\)).
    -   `rad_for_spline`: Radial grid specifically for spline operations if shape functions are numerical.
    -   `vee`, `vex`: DFT+U and EXX interaction matrices on correlated orbitals.
    -   `zioneff`: Effective ionic charge seen by projectors.
    -   `nablaphi`, `tnablaphi`: \(d\phi/dr - \phi/r\) and \(d\tilde{\phi}/dr - \tilde{\phi}/r\).

### Public Procedures

-   **`pawtab_free(Pawtab)` (Interface for `_0D`, `_1D`)**: Deallocates all allocatable members of `pawtab_type` instance(s).
-   **`pawtab_nullify(Pawtab)` (Interface for `_0D`, `_1D`)**: Resets flags and scalar members to default/uninitialized states.
-   **`pawtab_get_lsize(Pawtab, l_size_atm, natom, typat, mpi_atmtab)`**: Fills `l_size_atm` with `Pawtab(typat(iatom))%lcut_size` for each atom.
-   **`pawtab_set_flags(Pawtab, has_coretau, ...)` (Interface for `_0D`, `_1D`)**: Sets the various `has_*` flags in `Pawtab` instance(s).
-   **`pawtab_print(Pawtab, header, unit, prtvol, mode_paral)`**: Prints a summary of the contents of `Pawtab` instance(s).
-   **`pawtab_bcast(pawtab, comm_mpi, only_from_file)`**: Broadcasts `pawtab_type` data using MPI. `only_from_file` controls if only data directly read from PSP is broadcast or all computed data. Involves extensive serialization/deserialization.

## Important Variables/Constants

-   The numerous `has_*` flags are critical for determining the state of the `pawtab_type` object and what data is available or needs to be computed.
-   The various `*_mesh_size` integers define the dimensions of the radial arrays.
-   `lmn_size`, `lmn2_size`, `ij_size` define the number of PAW channels and pairs.

## Usage Examples

A `pawtab_type` variable is typically populated by routines in `m_pawpsp` (e.g., `pawpsp_main` calls `pawpsp_7in` or `pawpsp_17in`, which in turn call `pawpsp_read` and `pawpsp_calc`).

```fortran
module m_example_pawtab_usage
    use m_pawtab
    use m_libpaw_defs, only: dp
    ! Assume paw_setup is a variable of type pawtab_type, already populated
    ! by a call to, e.g., pawpsp_main.

    subroutine check_paw_data(paw_setup)
        type(pawtab_type), intent(in) :: paw_setup

        call pawtab_print(paw_setup, header="PAW Setup for Current Atom Type", prtvol=1)

        if (paw_setup%has_kij == 2) then
            ! print *, "Kinetic energy matrix K_ij(1,1): ", paw_setup%kij(1)
        end if

        ! To free:
        ! type(pawtab_type) :: mutable_paw_setup
        ! mutable_paw_setup = paw_setup ! if needed to modify, or pass intent(inout)
        ! call pawtab_free(mutable_paw_setup)
    end subroutine check_paw_data

end module m_example_pawtab_usage
```

## Dependencies and Interactions

-   **`m_libpaw_defs`, `m_libpaw_tools`, `m_libpaw_mpi`, `m_libpaw_memory`**: For basic definitions, utilities, MPI, and memory management.
-   **`m_paw_numeric`**: `paw_derfc` is used by `screened_coul_kernel` (which seems to be part of `m_pawrad` but mentioned here, perhaps an old location or indirect use).
-   **`m_paw_sphharm`, `m_pawang`**: Provide data/routines for angular parts, which are implicitly needed to make sense of some radial functions stored here (e.g., the `l` quantum number for `shapefunc(:,l)`).
-   **`m_pawrad`**: Defines `pawrad_type` used for radial grids, and its routines like `pawrad_init`, `simp_gen`, `nderiv_gen` are used by `m_pawpsp` and `m_paw_atom` to compute many of the arrays stored in `pawtab_type`.
-   **`m_pawxmlps`**: For reading XML PAW files, whose parsed data is then used to populate `pawtab_type`.
-   **`m_pawxc`**: For on-site XC calculations that might use/generate data related to `pawtab_type`.
-   **`m_paw_atom`**: Routines like `atompaw_dij0`, `atompaw_kij` directly compute and store into `pawtab%dij0` and `pawtab%kij`.

The `pawtab_type` is arguably the most crucial data structure in `libPAW` for representing the fixed atomic data of a PAW dataset. Its comprehensive nature reflects the complexity of the PAW method. The MPI broadcast routine is vital for distributing this dataset in parallel calculations.
