# `m_pawrad.F90`

## Overview

The Fortran module `m_pawrad` is dedicated to defining and managing radial meshes (grids) used in Projector Augmented Wave (PAW) calculations. It introduces the `pawrad_type` derived data type to encapsulate all information related to a specific radial grid, such as grid points, integration weights, and mesh type.

The module provides a suite of public procedures for:
-   Initializing various types of radial grids (linear and different logarithmic schemes).
-   Allocating and deallocating memory associated with these grid structures.
-   Printing information about a grid.
-   Comparing and copying grid structures.
-   Performing numerical operations on functions defined on these grids, including:
    -   Numerical integration using Simpson's rule (`simp_gen`).
    -   Numerical differentiation (`nderiv_gen`, `nderiv_lin`).
    -   Solving Poisson's equation for a given radial charge distribution (`poisson`).
    -   Calculating the radial part of Slater integrals (`calc_slatradl`).
-   Utility functions like finding the grid index corresponding to a given radius (`pawrad_ifromr`), extrapolating function values to \(r=0\) (`pawrad_deducer0`), and evaluating screened Coulomb interaction kernels.
-   Broadcasting `pawrad_type` objects in MPI parallel environments.

## Key Components

### Derived Type: `pawrad_type`

-   **Purpose**: Stores all information defining a specific radial grid and associated integration weights.
-   **Integer Scalars**:
    -   `int_meshsz`: The number of grid points to be used for numerical integration (can be less than `mesh_size` if `r_for_intg` is specified).
    -   `mesh_size`: Total number of points in the radial mesh.
    -   `mesh_type`: Integer code specifying the type of radial grid:
        -   `1` (RMESH_LINEAR): Linear grid, \(r(i) = (i-1) \cdot AA\).
        -   `2` (RMESH_LOG1): Logarithmic grid, \(r(i) = AA \cdot (e^{BB \cdot (i-1)} - 1)\).
        -   `3` (RMESH_LOG2): Logarithmic grid, \(r(i>1) = AA \cdot e^{BB \cdot (i-2)}\), \(r(1)=0\).
        -   `4` (RMESH_LOG3): Logarithmic grid, \(r(i) = -AA \cdot \ln(1 - BB \cdot (i-1))\), where \(BB=1/N\).
        -   `5` (RMESH_NL): Non-linear grid, \(r(i) = AA \cdot i / (N-i)\) (N is `lstep` here).
-   **Real(`dp`) Scalars**:
    -   `lstep`: The \(BB\) parameter for logarithmic grids, or \(N\) for `mesh_type=5`.
    -   `rmax`: The radius of the last grid point, \(r(\text{mesh\_size})\).
    -   `rstep`: The \(AA\) parameter (radial step or scaling factor).
    -   `stepint`: Effective step size used for generalized integration/differentiation formulas.
-   **Real(`dp`) Arrays (Allocatable)**:
    -   `rad(mesh_size)`: Stores the radial coordinate \(r_i\) of each grid point.
    -   `radfact(mesh_size)`: Stores factors \(dr/d\xi\) used in integration/differentiation when transforming to a uniform auxiliary grid \(\xi\). For linear grids, this is 1. For log grids, it involves \(r_i\) or \(r_i+AA\).
    -   `simfact(mesh_size)`: Simpson's rule integration weights, pre-multiplied by `radfact` and \(h/3\) factors.

### Public Procedures

-   **`pawrad_init(mesh, mesh_size, mesh_type, rstep, lstep, r_for_intg)`**:
    -   Main constructor. Initializes a `pawrad_type` variable `mesh`. Sets up `rad`, `radfact`, `simfact` arrays based on the specified type and parameters. `r_for_intg` defines up to which radius integrations are performed, setting `int_meshsz`.
-   **`pawrad_free(Rmesh)` (Interface for `pawrad_free_0D`, `pawrad_free_1D`)**:
    -   `pawrad_free_0D(Rmesh)`: Deallocates allocatable arrays in a single `pawrad_type` object.
    -   `pawrad_free_1D(Rmesh(:))`: Calls `pawrad_free_0D` for each element in an array of `pawrad_type`.
-   **`pawrad_print(Rmesh, header, unit, prtvol, mode_paral)`**: Prints information about the radial mesh `Rmesh`.
-   **`pawrad_isame(Rmesh1, Rmesh2, hasameq, whichdenser)`**: Compares two meshes. `hasameq` is true if they have the same defining equation. `whichdenser` indicates which mesh is denser or if they are incompatible.
-   **`pawrad_copy(mesh1, mesh2)`**: Performs a deep copy from `mesh1` to `mesh2`.
-   **`pawrad_ifromr(radmesh, rr)`**: Integer function, returns the grid index `i` such that `radmesh%rad(i)` is closest to or appropriate for the given radius `rr`.
-   **`pawrad_deducer0(func, funcsz, radmesh)`**: Extrapolates the value of `func(1)` (i.e., at \(r=0\)) using a 3-point formula based on `func(2), func(3), func(4)`.
-   **`pawrad_bcast(pawrad, comm_mpi)`**: Broadcasts a `pawrad_type` object using MPI. Serializes data into integer and real buffers.
-   **`simp_gen(intg, func, radmesh, r_for_intg)`**: Performs numerical integration of `func` over `radmesh` using Simpson's rule up to `r_for_intg` (or `radmesh%int_meshsz`). Uses precomputed `simfact` or computes them if `r_for_intg` is specified.
-   **`nderiv_gen(der, func, radmesh, der2)`**: Computes first (`der`) and optionally second (`der2`) derivatives of `func` on `radmesh`. It transforms to an equivalent linear grid for differentiation via `nderiv_lin`.
-   **`nderiv_lin(hh, yy, zz, ndim, norder)`**: Computes 1st or 2nd derivative (`norder`) of `yy` on a linear grid with spacing `hh`. Uses finite difference formulas.
-   **`bound_deriv(func, mesh, nn, yp1, ypn)`**: Computes derivatives at the boundaries of the mesh (`yp1` at `r(1)`, `ypn` at `r(nn)`).
-   **`poisson(den, ll, radmesh, rv, screened_sr_separation, qq)`**: Solves the radial Poisson equation for a given multipolar component `ll` of the charge density `den` (which is expected to be \(4\pi r^2 \rho_l(r)\) or \(r^2 \rho_l(r)\)). Outputs \(r \cdot V_l(r)\) in `rv`. Can handle screened Coulomb interaction if `screened_sr_separation` is provided.
-   **`screened_coul_kernel(order, r1, r2, formula)`**: Real function, evaluates the kernel \(\frac{r_<^l}{r_>^{l+1}}\) (or its screened version) for the screened Coulomb interaction, based on Angyan et al., J. Phys. A 39, 8613 (2006). `order` is \(l\).
-   **`calc_slatradl(ll, mesh_size, ff1, ff2, Pawrad, integral)`**: Calculates the radial part of Slater integrals: \(\frac{4\pi}{2L+1} \int ff1(r_1) \frac{r_<^L}{r_>^{L+1}} ff2(r_2) dr_1 dr_2\). Internally uses `poisson`.

## Important Variables/Constants

-   `RMESH_*` parameters: Integer constants defining different mesh types.
-   `rad()`, `radfact()`, `simfact()`: Core arrays within `pawrad_type` defining the grid and integration weights.
-   `mesh_type`, `rstep`, `lstep`: Key parameters defining the grid generation.

## Usage Examples

```fortran
module m_example_pawrad_usage
    use m_pawrad
    use m_libpaw_defs, only: dp
    implicit none

    subroutine test_radial_grid_operations
        type(pawrad_type) :: my_grid
        real(dp), dimension(100) :: my_function, my_integral_values
        real(dp) :: integral_result
        integer :: num_points = 100

        ! 1. Initialize a linear radial grid
        call pawrad_init(mesh=my_grid, mesh_size=num_points, mesh_type=RMESH_LINEAR, &
                         rstep=0.02_dp, r_for_intg=my_grid%rad(num_points))

        ! call pawrad_print(my_grid, header="My Linear Grid")

        ! 2. Define a function on this grid (e.g., r^2 * exp(-r))
        do i = 1, num_points
            my_function(i) = my_grid%rad(i)**2 * exp(-my_grid%rad(i))
        end do

        ! 3. Integrate the function
        ! For simp_gen, the input function should be f(r), not r^2*f(r) or 4*pi*r^2*f(r)
        ! if radfact and simfact correctly incorporate r^2 and 4pi.
        ! Assuming func is the bare function f(r) and simfact includes r^2 and 4pi/ (2l+1) etc.
        ! The example here is conceptual. The exact form of func depends on what simp_gen expects.
        ! Based on poisson and calc_slatradl, simp_gen likely expects the full integrand including r^2 factors.
        ! Let's assume my_function is already f(r)*r^2
        call simp_gen(integral_result, my_function, my_grid)
        ! print *, "Integral of r^2*exp(-r) dr from 0 to rmax: ", integral_result

        ! 4. Free the grid
        call pawrad_free(my_grid)

    end subroutine test_radial_grid_operations

end module m_example_pawrad_usage
```

## Dependencies and Interactions

-   **`m_libpaw_defs` (`USE_DEFS`)**: For `dp` and various mathematical constants (`pi`, `tol*`).
-   **`m_libpaw_tools` (`USE_MSG_HANDLING`)**: For error reporting (`LIBPAW_ERROR`, `LIBPAW_BUG`).
-   **`m_libpaw_mpi` (`USE_MPI_WRAPPERS`)**: For `xmpi_comm_rank` and `xmpi_bcast` in `pawrad_bcast`.
-   **`m_libpaw_memory` (`USE_MEMORY_PROFILING`)**: For memory allocation macros.
-   **`m_paw_numeric`**: Provides `paw_derfc` used by `screened_coul_kernel`.
-   This module is fundamental for all other PAW modules that operate on radial grids, such as `m_pawtab` (which stores functions on these grids), `m_paw_atom` (for atomic calculations), `m_pawdij` (for \(D_{ij}\) calculations involving radial integrals).

The `m_pawrad` module provides the foundational tools for defining and working with radial discretizations in `libPAW`. The support for different grid types allows flexibility, and the provision of integration and differentiation routines is essential for subsequent PAW calculations. The Simpson's rule implementation (`simp_gen`) and its interaction with `radfact` and `simfact` are key for accurate radial integrals.
