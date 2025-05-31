# `m_paw_onsite.F90`

## Overview

The Fortran module `m_paw_onsite` is dedicated to computing on-site matrix elements of the nabla (gradient) operator, \(\nabla\), within the Projector Augmented Wave (PAW) formalism. These matrix elements are essential for various physical properties, particularly those involving momentum or current operators. The module provides routines to calculate both valence-valence and core-valence contributions to these matrix elements.

## Key Components

### Public Subroutines

-   **`pawnabla_init(mpsang, ntypat, pawrad, pawtab)`**:
    -   **Purpose**: Evaluates the valence-valence on-site contributions of the nabla operator in Cartesian coordinates. Specifically, it computes for each atom type:
        \[ \langle\phi_i|\nabla|\phi_j\rangle - \langle\tilde\phi_i|\nabla|\tilde\phi_j\rangle \]
        where \(\phi_i\) are the all-electron partial waves and \(\tilde\phi_i\) are the pseudo partial waves.
    -   **Inputs**:
        -   `mpsang`: Integer, \(1 + \text{maximum angular momentum } l\).
        -   `ntypat`: Integer, number of atom types.
        -   `pawrad(ntypat)`: `type(pawrad_type)`, PAW radial grid data for each atom type.
        -   `pawtab(ntypat)`: `type(pawtab_type)` (intent inout), PAW setup data. The results are stored within this structure.
    -   **Side Effects**:
        -   Updates `pawtab(itypat)%nabla_ij(3, lmn_size, lmn_size)` with the computed matrix elements for each atom type `itypat`. The first dimension corresponds to the Cartesian components (x,y,z).
        -   Sets `pawtab(itypat)%has_nabla = 2` upon completion.
    -   **Logic**:
        1.  Calls `setnabla_ylm` (from `m_paw_sphharm`) to get pre-calculated angular integrals involving \(\nabla\) and spherical harmonics. These are stored in `ang_phipphj`.
        2.  For each atom type:
            a.  Allocates `pawtab(itypat)%nabla_ij`.
            b.  Computes two types of radial integrals:
                -   `int1(iln,jln)`: \(\int [ \phi_{iln} \frac{d\phi_{jln}}{dr} - \frac{\phi_{iln}\phi_{jln}}{r} ] - [ \tilde\phi_{iln} \frac{d\tilde\phi_{jln}}{dr} - \frac{\tilde\phi_{iln}\tilde\phi_{jln}}{r} ] dr \)
                -   `int2(iln,jln)`: \(\int [ \frac{\phi_{iln}\phi_{jln}}{r} - \frac{\tilde\phi_{iln}\tilde\phi_{jln}}{r} ] dr \)
                (Note: The integrands are effectively \(r^2 \times (\text{terms})\) due to the \(r^2 dr\) volume element, but the formulas are written as if operating on \(P(r)=r\phi(r)\), so \(P \frac{d P'}{dr} - \frac{PP'}{r}\) and \(PP'/r\).)
            c.  Combines these radial integrals with the angular matrix elements `ang_phipphj` to form the Cartesian components of \(\langle\phi_i|\nabla|\phi_j\rangle - \langle\tilde\phi_i|\nabla|\tilde\phi_j\rangle\).
            d.  Symmetrizes/anti-symmetrizes the result: `nabla_ij(ii,ilmn,jlmn) = -nabla_ij(ii,jlmn,ilmn)`.

-   **`pawnabla_core_init(mpsang, ntypat, pawrad, pawtab, phi_cor, indlmn_cor, diracflag)`**:
    -   **Purpose**: Evaluates the core-valence on-site contributions of the nabla operator:
        \[ \langle\phi_i|\nabla|\phi_{core,j}\rangle \]
        (The tilde part \(\langle\tilde\phi_i|\nabla|\tilde\phi_{core,j}\rangle\) is implicitly zero as pseudo core states are typically not defined or don't contribute to this term). This routine can handle relativistic (Dirac) core wave functions.
    -   **Inputs**:
        -   `mpsang`, `ntypat`, `pawrad`, `pawtab`: Similar to `pawnabla_init`.
        -   `phi_cor(:,:)`: `real(dp)`, core wave functions for the current atom type.
        -   `indlmn_cor(:,:)`: Integer array, maps \((l,m,n,s,\kappa, ...)\) indices for core states.
        -   `diracflag` (optional): Integer. If 1, indicates Dirac relativistic core wave functions are used, requiring handling of complex spherical harmonics and Clebsch-Gordan coefficients.
    -   **Side Effects**:
        -   Updates `pawtab(itypat)%nabla_ij(3, lmn_size, lmncmax)` with core-valence matrix elements. `lmncmax` is the number of core states.
        -   If `diracflag=1`, also updates `pawtab(itypat)%nabla_im_ij` with the imaginary parts of the matrix elements.
        -   Sets `pawtab(itypat)%has_nabla = 3` (non-relativistic/scalar-relativistic) or `4` (Dirac relativistic).
    -   **Logic**:
        1.  Similar to `pawnabla_init`, calls `setnabla_ylm` for angular integrals.
        2.  For each atom type:
            a.  Allocates `pawtab(itypat)%nabla_ij` (and `nabla_im_ij` if Dirac).
            b.  Computes radial integrals `int1` and `int2`, but now between valence partial waves \(\phi_i\) and core wave functions \(\phi_{core,j}\):
                -   `int1(iln,jln_core)`: \(\int [ \phi_{iln} \frac{d\phi_{core,jln}}{dr} - \frac{\phi_{iln}\phi_{core,jln}}{r} ] dr \)
                -   `int2(iln,jln_core)`: \(\int [ \frac{\phi_{iln}\phi_{core,jln}}{r} ] dr \)
            c.  If `diracflag=1`:
                -   Uses Clebsch-Gordan coefficients (`cgc`) to couple spin and orbital angular momentum.
                -   Handles conversion from complex to real spherical harmonics for \(m \neq 0\) core states, calculating real and imaginary parts of the matrix elements separately.
            d.  Else (scalar-relativistic or non-relativistic):
                -   Combines radial and angular integrals directly, similar to `pawnabla_init`.

## Important Variables/Constants

-   `pawrad_type`: Defined in `m_pawrad`, provides radial grid information.
-   `pawtab_type`: Defined in `m_pawtab`, stores PAW setup data, including partial waves, projector info, and output arrays like `nabla_ij`.
-   `ang_phipphj`: Array storing pre-calculated angular integrals involving spherical harmonics and the components of the \(\nabla\) operator. These are computed by `setnabla_ylm`.
-   `indlmn`, `indlmn_cor`: Arrays mapping composite indices to \((l,m,n)\) quantum numbers for valence and core states, respectively.

## Usage Examples

These routines are typically called during the initialization phase of a PAW calculation, after PAW datasets have been read and processed.

```fortran
module m_example_onsite_nabla
    use m_libpaw_defs, only: dp
    use m_paw_onsite
    use m_pawrad, only: pawrad_type
    use m_pawtab, only: pawtab_type
    ! ... other necessary modules ...
    implicit none

    subroutine calculate_nabla_matrix_elements(num_atom_types, all_radial_grids, all_paw_setups, &
                                               core_wfs_type1, core_indices_type1)
        integer, intent(in) :: num_atom_types
        type(pawrad_type), intent(in) :: all_radial_grids(num_atom_types)
        type(pawtab_type), intent(inout) :: all_paw_setups(num_atom_types) ! Modified

        ! Example for one atom type's core data
        real(dp), intent(in) :: core_wfs_type1(:,:) ! (mesh_size, num_core_states)
        integer, intent(in) :: core_indices_type1(:,:) ! (8, num_core_states) for Dirac

        integer, parameter :: max_ang_mom_plus_1 = 4 ! Example: up to l=3

        ! Calculate valence-valence nabla contributions
        call pawnabla_init(mpsang=max_ang_mom_plus_1, ntypat=num_atom_types, &
                           pawrad=all_radial_grids, pawtab=all_paw_setups)

        ! Example: Calculate core-valence for the first atom type (assuming it has core states)
        ! This would typically be inside a loop over atom types if core_wfs/indices vary
        if (num_atom_types > 0 .and. size(core_wfs_type1,2) > 0) then
             call pawnabla_core_init(mpsang=max_ang_mom_plus_1, ntypat=1, &
                                    pawrad=all_radial_grids(1:1), pawtab=all_paw_setups(1:1), &
                                    phi_cor=core_wfs_type1, indlmn_cor=core_indices_type1, &
                                    diracflag=1) ! Example: Dirac core
        end if

        ! Now all_paw_setups(ityp)%nabla_ij contains the computed matrix elements.
        ! all_paw_setups(ityp)%has_nabla indicates the status.

    end subroutine calculate_nabla_matrix_elements

end module m_example_onsite_nabla
```

## Dependencies and Interactions

-   **`m_libpaw_defs` (`USE_DEFS`)**: For `dp`, `half`, `zero`, etc.
-   **`m_libpaw_tools` (`USE_MSG_HANDLING`)**: For error reporting (`LIBPAW_BUG`).
-   **`m_libpaw_memory` (`USE_MEMORY_PROFILING`)**: For memory allocation macros.
-   **`m_pawrad`**: Defines `pawrad_type` and provides radial integration (`simp_gen`) and derivative (`nderiv_gen`) routines.
-   **`m_pawtab`**: Defines `pawtab_type` which stores the input PAW functions and the output nabla matrix elements.
-   **`m_paw_sphharm`**: Provides `setnabla_ylm` which computes and tabulates the necessary angular matrix elements of the \(\nabla\) operator components.

This module calculates fundamental PAW matrix elements involving the gradient operator, which are building blocks for properties related to momentum, currents, or response to vector potentials. The handling of relativistic core states in `pawnabla_core_init` adds a layer of sophistication for heavy elements.
