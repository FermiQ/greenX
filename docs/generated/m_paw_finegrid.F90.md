# `m_paw_finegrid.F90`

## Overview

The Fortran module `m_paw_finegrid` is dedicated to calculations involving quantities defined on a fine real-space grid within the PAW atomic spheres. These quantities are essential for accurately representing functions and operators in the vicinity of atoms where the PAW method makes corrections to pseudopotential calculations. The module provides routines to:
-   Generate points on this fine grid relative to an atom's position (\(\vec{r}-\vec{R}\)).
-   Compute products of radial shape functions \(g_l(|\vec{r}-\vec{R}|)\) and real spherical harmonics \(Y_{lm}(\widehat{\vec{r}-\vec{R}})\), along with their Cartesian derivatives.
-   Compute the Fourier transform of these \(g_l Y_{lm}\) products.
-   Compute phase factors \(e^{i\vec{q}\cdot\vec{r}}\) on the fine grid for phonon calculations.

The module handles different analytical forms for shape functions and includes logic for their behavior near the origin and for numerical shape functions provided via splines.

## Key Components

### Public Subroutines

-   **`pawgylm(gylm, gylmgr, gylmgr2, lm_size, nfgd, optgr0, optgr1, optgr2, pawtab, rfgd)`**:
    -   **Purpose**: Computes \(g_l(|\vec{r}-\vec{R}|) Y_{lm}(\widehat{\vec{r}-\vec{R}})\) and optionally its first and second Cartesian derivatives on the fine grid points `rfgd`.
    -   **Inputs**:
        -   `lm_size`: Number of \( (l,m) \) components.
        -   `nfgd`: Number of fine grid points.
        -   `optgr0`, `optgr1`, `optgr2`: Integer flags (0 or 1) to control computation of the function, first derivatives, and second derivatives, respectively.
        -   `pawtab`: `type(pawtab_type)`, PAW setup data for the current atom type (contains shape function type, parameters, radial functions for splining, etc.).
        -   `rfgd(3,nfgd)`: `real(dp)`, Cartesian coordinates of fine grid points relative to the atom center (\(\vec{r}-\vec{R}\)).
    -   **Outputs**:
        -   `gylm(nfgd, lm_size)` (if `optgr0=1`): The values \(g_l Y_{lm}\).
        -   `gylmgr(3, nfgd, lm_size)` (if `optgr1=1`): First Cartesian derivatives \(\nabla (g_l Y_{lm})\).
        -   `gylmgr2(6, nfgd, lm_size)` (if `optgr2=1`): Second Cartesian derivatives (components xx, yy, zz, yz, xz, xy).
    -   **Logic**:
        1.  Initializes output arrays.
        2.  Determines `l_size` from `pawtab%lcut_size`.
        3.  Retrieves shape function parameters (`shape_type`, `sigma`, `lambda`, `pi_over_rshp`, `rcut`, `alpha`, `qq`) from `pawtab`.
        4.  If `shape_type == -1` (numerical), it splines the shape function and its derivatives from `pawtab` onto the provided `rnrm` values (radial distances of `rfgd` points).
        5.  Computes \(Y_{lm}\) and their derivatives (if needed) using `initylmr`.
        6.  Computes radial factors \(g_l(r)\), \(g_l'(r)/r\), and \((g_l''(r) - g_l'(r)/r)/r^2\) based on `shape_type` using internal helper functions (see below) or splined values. Special care is taken for \(r \to 0\) limits.
        7.  Combines radial factors with \(Y_{lm}\) and their derivatives to form `gylm`, `gylmgr`, `gylmgr2`.

-   **`pawgylmg(gprimd, gylmg, kg, kpg, kpt, lmax, nkpg, npw, ntypat, pawtab, ylm)`**:
    -   **Purpose**: Computes the Fourier transform of \(g_l(r)Y_{lm}(r)\) for each atom type, for a set of plane waves (G-vectors).
    -   \( (g_l Y_{lm})(|\vec{k}+\vec{G}|) = 4\pi i^l \int dr r^2 j_l(|\vec{k}+\vec{G}|r) g_l(r) \cdot Y_{lm}(\widehat{\vec{k}+\vec{G}}) \). The routine uses pre-calculated \(Y_{lm}(\widehat{\vec{k}+\vec{G}})\) and \(g_l(|\vec{k}+\vec{G}|)\) (Fourier transform of radial part).
    -   **Inputs**:
        -   `gprimd(3,3)`: Reciprocal space primitive vectors.
        -   `kg(3,npw)`: Integer G-vector coordinates.
        -   `kpg(npw,nkpg)`: \(|\vec{k}+\vec{G}|\) components (if `useylm=1`, not fully clear from context).
        -   `kpt(3)`: Reduced coordinates of the k-point.
        -   `lmax`: Max \(l+1\).
        -   `nkpg`: Second dimension of `kpg`.
        -   `npw`: Number of plane waves.
        -   `ntypat`: Number of atom types.
        -   `pawtab(ntypat)`: PAW setup data (contains pre-tabulated Fourier transforms of shape functions `shapefncg`).
        -   `ylm(npw, lmax**2)`: Real spherical harmonics \(Y_{lm}(\widehat{\vec{k}+\vec{G}})\).
    -   **Output**:
        -   `gylmg(npw, lmax**2, ntypat)`: Fourier transformed values.
    -   **Logic**:
        1.  Computes \(|\vec{k}+\vec{G}|\) norms (`kpgnorm`).
        2.  For each atom type and each \(l\):
            a.  Interpolates `pawtab(itypat)%shapefncg` (which is \( \int dr r^2 j_l(qr) g_l(r) \)) onto the `kpgnorm` values using `paw_uniform_splfit` to get `glg(:) = \tilde{g}_l(|\vec{k}+\vec{G}|)\).
            b.  Multiplies by `ylm` to get the final result: `gylmg = ylm * glg`.

-   **`pawrfgd_fft(ifftsph, gmet, n1, n2, n3, nfgd, rcut, rfgd, rprimd, ucvol, xred, fft_distrib, fft_index, me_fft)`**:
    -   **Purpose**: Determines the fine grid points \(\vec{r}\) within a sphere of radius `rcut` around an atom at `xred` (reduced coordinates) for an FFT grid, and computes \(\vec{r}-\vec{R}\).
    -   **Inputs**:
        -   `gmet(3,3)`: Reciprocal space metric tensor.
        -   `n1,n2,n3`: Dimensions of the full FFT grid.
        -   `rcut`: Cutoff radius for the sphere.
        -   `rprimd(3,3)`: Real space primitive vectors.
        -   `ucvol`: Unit cell volume.
        -   `xred(3)`: Reduced coordinates of the atom center \(\vec{R}\).
        -   `fft_distrib`, `fft_index`, `me_fft` (optional): MPI distribution info for FFT grid.
    -   **Outputs**:
        -   `ifftsph(:)`: Allocatable, FFT indices of points within the sphere.
        -   `nfgd`: Number of points found.
        -   `rfgd(:,:)`: Allocatable, Cartesian coordinates \(\vec{r}-\vec{R}\) for these points.
    -   **Logic**:
        1.  Defines a bounding box around the atom based on `rcut` and `gmet`.
        2.  Loops over FFT grid indices (i1,i2,i3) within this box.
        3.  Handles MPI distribution of FFT planes.
        4.  For each point \(\vec{r}\) (from i1,i2,i3), calculates \(\vec{r}-\vec{R}\) and its norm squared.
        5.  If \(|\vec{r}-\vec{R}|^2 \le rcut^2\), the point is selected, and its FFT index and \(\vec{r}-\vec{R}\) vector are stored.

-   **`pawrfgd_wvl(...)`**:
    -   **Purpose**: Similar to `pawrfgd_fft` but for wavelet grids (BigDFT context). Determines fine grid points around an atom.
    -   **Inputs**: Includes geometry code `geocode`, grid spacing `hh`, various wavelet grid parameters (`i3s`, `n1i`, `n2i`, `n3pi`, `shift`), and atom Cartesian coordinates `xcart`.
    -   **Logic**: Loops over a box defined by `cutoff` (related to `rloc`) around the atom, checks periodicity, and selects points if \(|\vec{r}-\vec{R}|^2 \le rcut^2\). Uses internal helper subroutines `my_ind_positions` and `my_ext_buffers` for wavelet grid indexing and boundary conditions.

-   **`pawexpiqr(expiqr, gprimd, nfgd, qphon, rfgd, xred)`**:
    -   **Purpose**: Computes \(e^{i\vec{q}\cdot\vec{r}}\) for each point \(\vec{r}\) on the fine grid. \(\vec{r}\) is effectively \((\vec{r}-\vec{R}) + \vec{R}\).
    -   **Inputs**:
        -   `gprimd(3,3)`: Reciprocal space primitive vectors.
        -   `nfgd`: Number of fine grid points.
        -   `qphon(3)`: Reduced coordinates of the phonon wavevector \(\vec{q}\).
        -   `rfgd(3,nfgd)`: Coordinates \(\vec{r}-\vec{R}\).
        -   `xred(3)`: Reduced coordinates of the atom \(\vec{R}\).
    -   **Output**:
        -   `expiqr(2,nfgd)`: Stores `cos(phase)` and `sin(phase)` where `phase = \vec{q} \cdot \vec{r}\).
    -   **Logic**: Converts \(\vec{q}_{reduced}\) to Cartesian \(\vec{q}_{cart}\). Computes phase as \(\vec{q}_{cart} \cdot (\vec{r}-\vec{R}) + 2\pi \vec{q}_{reduced} \cdot \vec{R}_{reduced}\).

### Internal Helper Functions (within `pawgylm`)

-   A series of small functions (`shapefunc1`, `shapefunc1_0`, `dshpfunc1`, `d2shpfunc1_ovr2_0_3`, etc.) are defined within `pawgylm`'s `CONTAINS` block. These provide analytical expressions for \(g_l(r)\) and its relevant derivatives (or combinations like \(g'(r)/r\)) for different `shape_type` values (1: Gaussian-like, 2: Sinc-squared-like, 3: Bessel-based) and their Taylor expansions for small \(r\).
-   Module-level private variables `lambda`, `pi_over_rshp`, `sigma`, `alpha`, `qq` are used by these internal functions, capturing parameters from `pawtab`.

## Important Variables/Constants

-   `pawtab_type`: From `m_pawtab`, provides PAW setup data like shape function parameters and tabulated functions.
-   `ffact`: Parameter array in `pawgylm` for \( (2l+1)!! \) like factors for Bessel shape functions.
-   `toldev`: Tolerance used in `pawgylm` to switch to Taylor expansions near \(r=0\).

## Usage Examples

```fortran
module m_example_finegrid_usage
    use m_paw_finegrid
    use m_pawtab, only: pawtab_type
    use m_libpaw_defs, only: dp
    ! ... other necessary modules ...
    implicit none

    subroutine calculate_on_fine_grid(atom_paw_setup, atom_pos_reduced, &
                                      fft_grid_dims, recip_metric, real_vectors, cell_volume, &
                                      cutoff_radius)
        type(pawtab_type), intent(in) :: atom_paw_setup
        real(dp), intent(in) :: atom_pos_reduced(3)
        integer, intent(in) :: fft_grid_dims(3)
        real(dp), intent(in) :: recip_metric(3,3), real_vectors(3,3), cell_volume
        real(dp), intent(in) :: cutoff_radius

        integer :: num_fine_grid_pts
        integer, allocatable :: fine_grid_fft_indices(:)
        real(dp), allocatable :: fine_grid_coords_rel(:,:) ! (3, num_fine_grid_pts)

        real(dp), allocatable :: gylm_values(:,:)     ! (num_fine_grid_pts, lm_size)
        integer, parameter :: lm_size_example = (2*atom_paw_setup%lcut_max + 1)**2 ! Example, lcut_max from pawtab

        ! 1. Determine fine grid points around the atom
        call pawrfgd_fft(ifftsph=fine_grid_fft_indices, gmet=recip_metric, &
                         n1=fft_grid_dims(1), n2=fft_grid_dims(2), n3=fft_grid_dims(3), &
                         nfgd=num_fine_grid_pts, rcut=cutoff_radius, rfgd=fine_grid_coords_rel, &
                         rprimd=real_vectors, ucvol=cell_volume, xred=atom_pos_reduced)

        if (num_fine_grid_pts > 0) then
            allocate(gylm_values(num_fine_grid_pts, lm_size_example))

            ! 2. Compute g_l(r)Y_lm(r) on these points
            call pawgylm(gylm=gylm_values, gylmgr=legacy_dummy_gr, gylmgr2=legacy_dummy_gr2, &
                         lm_size=lm_size_example, nfgd=num_fine_grid_pts, &
                         optgr0=1, optgr1=0, optgr2=0, &
                         pawtab=atom_paw_setup, rfgd=fine_grid_coords_rel)
            ! Note: legacy_dummy_gr/gr2 would need to be declared if optgr1/2 were 1.

            ! ... use gylm_values ...

            deallocate(gylm_values)
            deallocate(fine_grid_fft_indices, fine_grid_coords_rel)
        end if
    contains
        ! Dummy arrays for unused optional outputs if not needed
        real(dp) :: legacy_dummy_gr(3,0,0), legacy_dummy_gr2(6,0,0)
    end subroutine calculate_on_fine_grid

end module m_example_finegrid_usage
```

## Dependencies and Interactions

-   **`m_libpaw_defs` (`USE_DEFS`)**: For `dp` and constants like `pi`, `tol*`.
-   **`m_libpaw_tools` (`USE_MSG_HANDLING`)**: For error reporting (`LIBPAW_BUG`, `LIBPAW_ERROR`).
-   **`m_libpaw_memory` (`USE_MEMORY_PROFILING`)**: For memory allocation macros.
-   **`m_pawtab`**: Defines `pawtab_type`, which is a key input providing PAW dataset parameters.
-   **`m_paw_sphharm`**: Provides `initylmr` for calculating real spherical harmonics and their derivatives.
-   **`m_paw_numeric`**: Provides numerical utilities:
    -   `paw_jbessel`: Spherical Bessel functions.
    -   `paw_splint`: Spline interpolation.
    -   `paw_uniform_splfit`: Spline fitting for uniform grid to non-uniform points (used in `pawgylmg`).
    -   `paw_sort_dp`: Sorting double precision numbers.

This module is critical for mapping PAW quantities onto a real-space fine grid, which is often necessary for evaluating matrix elements involving plane waves or other basis functions, and for calculating densities and potentials in real space. The handling of different shape functions and their derivatives, especially near the origin, demonstrates the detailed numerical work involved.
