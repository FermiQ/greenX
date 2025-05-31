# `m_pawang.F90`

## Overview

The Fortran module `m_pawang` defines and manages the `pawang_type` derived data type. This data structure is crucial for Projector Augmented Wave (PAW) calculations as it encapsulates data related to the angular discretization within the PAW augmentation spheres. This includes the definition of the angular mesh, evaluation of real spherical harmonics (\(S_{lm}\)) and their gradients on this mesh, computation of Gaunt coefficients, "nabla Gaunt" coefficients (integrals involving \(\nabla S_{lm}\)), rotation matrices for \(S_{lm}\) under crystal symmetries, and matrix elements of the spin-orbit coupling operator (\(\vec{L} \cdot \vec{S}\)).

The module provides two primary public procedures: `pawang_init` for initializing and populating a `pawang_type` object, and `pawang_free` for deallocating its components.

## Key Components

### Derived Type: `pawang_type`

-   **Purpose**: Stores all data related to the angular aspects of PAW calculations for a given setup.
-   **Integer Scalars**:
    -   `angl_size`: Total number of points on the angular mesh (`ntheta * nphi`).
    -   `l_max`: \(l_{max}+1\), where \(l_{max}\) is the maximum angular momentum.
    -   `l_size_max`: \(2 \cdot l_{max} - 1\) (likely related to the maximum \(L\) in Gaunt coefficients like \(\langle Y_{LM} | Y_{l_1m_1} Y_{l_2m_2} \rangle\)).
    -   `ngnt`: Number of non-zero real Gaunt coefficients stored.
    -   `nnablagnt`: Number of non-zero "nabla real Gaunt" coefficients stored.
    -   `ntheta`, `nphi`: Number of \(\theta\) and \(\phi\) points for the angular mesh.
    -   `nsym`: Number of symmetry operations.
    -   `nabgnt_option`: Flag (0 or 1) controlling computation of nabla Gaunt coefficients.
    -   `gnt_option`: Flag (0, 1, or 2) controlling computation and extent of Gaunt coefficients.
    -   `use_ls_ylm`: Flag (0 or 1) indicating if \(\vec{L} \cdot \vec{S}\) matrix elements are computed.
    -   `ylm_size`: Size of the `ylmr` and `ylmrgr` arrays, typically \((l_{max\_for\_Ylm}+1)^2\).
-   **Integer Arrays (Allocatable)**:
    -   `gntselect(:,:)`: Stores indices of non-zero real Gaunt coefficients.
    -   `nablagntselect(:,:,:)`: Stores indices of non-zero nabla real Gaunt coefficients.
-   **Real(`dp`) Arrays (Allocatable)**:
    -   `anginit(3,angl_size)`: Cartesian coordinates \((x,y,z)\) of points on the unit sphere for the angular mesh.
    -   `angwgth(angl_size)`: Weights for each angular mesh point, for spherical integration.
    -   `ls_ylm(:,:,:)`: Stores matrix elements of the \(\vec{L} \cdot \vec{S}\) operator in the real spherical harmonics basis.
    -   `nablarealgnt(:)`: Stores values of non-zero nabla real Gaunt coefficients \(\int S_{LM} (\nabla S_{l_1m_1}) \cdot (\nabla S_{l_2m_2}) d\Omega\).
    -   `realgnt(:)`: Stores values of non-zero real Gaunt coefficients \(\int S_{LM} S_{l_1m_1} S_{l_2m_2} d\Omega\).
    -   `ylmr(ylm_size,angl_size)`: Real spherical harmonics \(S_{lm}\) evaluated at each angular mesh point.
    -   `ylmrgr(:,:,:)`: Cartesian gradients \(\nabla S_{lm}\) (and optionally second derivatives) evaluated at each angular mesh point.
    -   `zarot(:,:,:,:)`: Rotation matrices for real spherical harmonics corresponding to crystal symmetry operations.

### Public Subroutines

-   **`pawang_init(Pawang, gnt_option, nabgnt_option, lmax, nphi, ntheta, nsym, ngrad2_ylm, use_angular_grid, use_ylm, use_ls_ylm)`**:
    -   **Purpose**: Initializes a `pawang_type` variable (`Pawang`). It allocates and computes various components based on the input flags and parameters.
    -   **Inputs**:
        -   `gnt_option`, `nabgnt_option`: Control calculation of Gaunt and nabla Gaunt coefficients.
        -   `lmax`: Maximum angular momentum \(l\).
        -   `nphi`, `ntheta`: Dimensions for the angular grid.
        -   `nsym`: Number of symmetry operations.
        -   `ngrad2_ylm`: Order of derivatives for \(S_{lm}\) to compute (0, 1, or 2).
        -   `use_angular_grid`: Flag to compute `anginit` and `angwgth`.
        -   `use_ylm`: Flag to compute `ylmr` and `ylmrgr`.
        -   `use_ls_ylm`: Flag to compute `ls_ylm`.
    -   **Output/Side Effect**: The `Pawang` variable is populated.
    -   **Logic**:
        1.  Sets scalar members of `Pawang` (e.g., `l_max`, `nsym`).
        2.  If `use_angular_grid == 1`, calls `ylm_angular_mesh` to generate the grid.
        3.  If `use_ylm > 0` and grid exists, calls `initylmr` to compute \(S_{lm}\) and their derivatives on the grid.
        4.  If `gnt_option > 0`, calls `realgaunt` to compute Gaunt coefficients.
        5.  If `nabgnt_option == 1`, calls `nablarealgaunt` to compute nabla Gaunt coefficients.
        6.  If `use_ls_ylm > 0`, calls `lsylm` to compute \(\vec{L} \cdot \vec{S}\) matrix elements.
        7.  If `nsym > 0`, allocates `Pawang%zarot` (rotation matrices are computed elsewhere, likely by `setsym_ylm` which would be called by a higher-level routine).

-   **`pawang_free(Pawang)`**:
    -   **Purpose**: Deallocates all allocatable arrays within the `Pawang` variable and resets scalar members to default/uninitialized values (e.g., `angl_size = 0`, `l_max = -1`).

## Important Variables/Constants

-   The module relies on constants from `m_libpaw_defs` (e.g., `dp`).
-   The various `*_option` and `use_*` flags passed to `pawang_init` are critical for determining which components of `pawang_type` are computed and stored, allowing for tailored setups based on calculation needs.

## Usage Examples

```fortran
module m_example_pawang_usage
    use m_pawang
    use m_libpaw_defs, only: dp
    ! For a full example, one would also need m_paw_sphharm for the routines called by pawang_init
    implicit none

    subroutine setup_angular_data
        type(pawang_type) :: angular_data
        integer, parameter :: lmax_calc = 2  ! Max l for this calculation
        integer, parameter :: ntheta_val = 10
        integer, parameter :: nphi_val = 20
        integer, parameter :: num_symmetries = 1 ! e.g., identity
        integer, parameter :: calc_gaunt = 1     ! Compute Gaunt coefficients
        integer, parameter :: calc_nabla_gaunt = 0 ! Skip nabla Gaunt
        integer, parameter :: ylm_deriv_order = 1  ! Compute Ylm and first derivatives
        integer, parameter :: need_angular_grid = 1
        integer, parameter :: need_ylm_on_grid = 1
        integer, parameter :: need_ls_matrix = 1

        call pawang_init(Pawang=angular_data, &
                         gnt_option=calc_gaunt, nabgnt_option=calc_nabla_gaunt, &
                         lmax=lmax_calc, nphi=nphi_val, ntheta=ntheta_val, &
                         nsym=num_symmetries, ngrad2_ylm=ylm_deriv_order, &
                         use_angular_grid=need_angular_grid, &
                         use_ylm=need_ylm_on_grid, use_ls_ylm=need_ls_matrix)

        ! Now angular_data contains the computed quantities, e.g.:
        ! if (angular_data%angl_size > 0) then
        !    print *, "First angular grid point (cartesian): ", angular_data%anginit(:,1)
        !    print *, "Weight of first angular point: ", angular_data%angwgth(1)
        ! end if
        ! if (angular_data%ngnt > 0) then
        !    print *, "First non-zero Gaunt coefficient: ", angular_data%realgnt(1)
        ! end if

        call pawang_free(angular_data)

    end subroutine setup_angular_data

end module m_example_pawang_usage
```

## Dependencies and Interactions

-   **`m_libpaw_defs` (`USE_DEFS`)**: For `dp`.
-   **`m_libpaw_tools` (`USE_MSG_HANDLING`)**: For error reporting.
-   **`m_libpaw_memory` (`USE_MEMORY_PROFILING`)**: For memory allocation macros.
-   **`m_paw_sphharm`**: This is a major dependency. `pawang_init` calls several routines from `m_paw_sphharm` to perform the actual calculations:
    -   `ylm_angular_mesh`
    -   `initylmr`
    -   `realgaunt`
    -   `nablarealgaunt`
    -   `lsylm`
    -   (Note: `setsym_ylm` for `zarot` is typically called externally and the results passed or set into `Pawang%zarot`).
-   Other `libPAW` modules that require pre-calculated angular data (like \(S_{lm}\) on a grid, Gaunt coefficients for integrating products of functions, etc.) will use an initialized `pawang_type` variable. For example, `m_paw_an` initialization takes `Pawang` as input.

The `m_pawang` module serves as a constructor and container for essential angular data used in PAW calculations, centralizing their computation and storage. The options for initialization allow for flexibility in what data is actually generated, optimizing memory and computation time.
