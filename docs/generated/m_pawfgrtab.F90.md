# `m_pawfgrtab.F90`

## Overview

The Fortran module `m_pawfgrtab` defines and manages the `pawfgrtab_type` derived data type. This data structure is used in Projector Augmented Wave (PAW) calculations to store various atom-centered quantities that are tabulated on a fine real-space grid (denoted "Fine Grid Radial TABulated" or similar). This fine grid is typically a dense rectangular (Cartesian) grid within the PAW augmentation sphere of an atom.

The `pawfgrtab_type` holds data such as:
-   The coordinates of the fine grid points relative to the atom's center (\(\vec{r}_k - \vec{R}_a\)).
-   Mapping indices from these local fine grid points to a global FFT grid (if applicable).
-   Pre-calculated values of radial shape functions multiplied by real spherical harmonics, \(g_l(r)S_{lm}(\hat{r})\), and their Cartesian derivatives on this fine grid.
-   Phase factors like \(e^{i\vec{q}\cdot\vec{r}}\) used in phonon calculations.
-   First-order "frozen" compensation charge densities (\(\hat{n}_{fr}\)) and their gradients, relevant for response function calculations.

The module provides public procedures for initializing (`pawfgrtab_init`), deallocating (`pawfgrtab_free`), nullifying flags in (`pawfgrtab_nullify`), copying (`pawfgrtab_copy`), printing (`pawfgrtab_print`), and handling MPI distribution (`pawfgrtab_gather`, `pawfgrtab_redistribute`) of `pawfgrtab_type` arrays.

## Key Components

### Derived Type: `pawfgrtab_type`

-   **Purpose**: Stores atom-centered data defined on a local fine real-space grid.
-   **Integer Scalars**:
    -   `cplex`: 1 if real data, 2 if complex (though most arrays here seem real, this might be for consistency or future use).
    -   `expiqr_allocated`, `gylm_allocated`, `gylmgr_allocated`, `gylmgr2_allocated`, `nhatfr_allocated`, `nhatfrgr_allocated`, `rfgd_allocated`: Integer flags (0: not allocated/computed, 1: allocated, 2: allocated and computed).
    -   `itypat`: Atom type index.
    -   `l_size`: \(l_{max\_gaunt}+1\), maximum angular momentum for Gaunt coefficients related to this atom type.
    -   `nfgd`: Number of fine grid points within the PAW sphere for this atom.
    -   `nspden`: Number of spin-density components.
-   **Integer Arrays (Allocatable)**:
    -   `ifftsph(nfgd)`: Maps local fine grid points to indices on a global FFT grid.
-   **Real(`dp`) Arrays (Allocatable)**:
    -   `expiqr(2, nfgd)`: Stores \(\cos(\vec{q}\cdot\vec{r})\) and \(\sin(\vec{q}\cdot\vec{r})\) on the fine grid.
    -   `gylm(nfgd, l_size*l_size)`: Stores \(g_l(r)S_{lm}(\hat{r})\) values on the fine grid.
    -   `gylmgr(3, nfgd, l_size*l_size)`: Cartesian gradients \(\nabla (g_l S_{lm})\) on the fine grid.
    -   `gylmgr2(6, nfgd, l_size*l_size)`: Second Cartesian derivatives of \(g_l S_{lm}\) (xx, yy, zz, yz, xz, xy components).
    -   `nhatfr(nfgd, nspden)`: "Frozen part" of the first-order compensation charge density \(\hat{n}_{fr}^{(1)}\) on the fine grid (for response functions).
    -   `nhatfrgr(3, nfgd, nspden)`: Cartesian gradients of \(\hat{n}_{fr}^{(1)}\).
    -   `rfgd(3, nfgd)`: Cartesian coordinates \(\vec{r}_k - \vec{R}_a\) of the fine grid points.

### Public Procedures

-   **`pawfgrtab_init(Pawfgrtab, cplex, l_size_atm, nspden, typat, mpi_atmtab, comm_atom)`**:
    -   Initializes an array of `pawfgrtab_type`. Sets scalar members and allocates arrays to zero size (actual data like `rfgd` points are computed elsewhere, e.g., by `pawrfgd_fft` or `pawrfgd_wvl` from `m_paw_finegrid`). Sets `*_allocated` flags to 0.
-   **`pawfgrtab_free(Pawfgrtab)`**:
    -   Deallocates all allocatable arrays within each element of `Pawfgrtab` and resets `*_allocated` flags to 0.
-   **`pawfgrtab_nullify(Pawfgrtab)`**:
    -   Resets `nfgd` and all `*_allocated` flags to 0. Does not deallocate.
-   **`pawfgrtab_copy(pawfgrtab_in, pawfgrtab_cp, mpi_atmtab, comm_atom)`**:
    -   Performs a deep copy from `pawfgrtab_in` to `pawfgrtab_cp`. Handles MPI distribution (scatter/gather if distributions differ). Allocates arrays in `pawfgrtab_cp` based on sizes in `pawfgrtab_in`.
-   **`pawfgrtab_print(...)`**:
    -   Prints information about the `pawfgrtab_type` instances, such as dimensions and allocation statuses. Does not print the large arrays themselves unless `prtvol` is very high (and even then, typically not).
-   **`pawfgrtab_gather(pawfgrtab, pawfgrtab_gathered, comm_atom, istat, mpi_atmtab)`**:
    -   Gathers distributed `pawfgrtab` arrays from processes in `comm_atom` into `pawfgrtab_gathered` on all processes. Serializes data into buffers for MPI_Allgatherv.
-   **`pawfgrtab_redistribute(...)`**:
    -   Redistributes `pawfgrtab` arrays between different MPI distributions (`mpi_comm_in` to `mpi_comm_out`). Supports a brute-force (gather+scatter) or an asynchronous communication method using precomputed send/receive lists.

### Private Procedures for MPI Communication

-   **`pawfgrtab_isendreceive_fillbuffer(...)`**: Serializes data from `pawfgrtab_type` instances into integer and real buffers for MPI sending. Packs metadata and actual array data.
-   **`pawfgrtab_isendreceive_getbuffer(...)`**: Deserializes data from received MPI buffers into `pawfgrtab_type` instances. Allocates arrays based on received metadata.

## Important Variables/Constants

-   The `*_allocated` flags within `pawfgrtab_type` are crucial for tracking the status of potentially large arrays, determining if they need to be computed or can be used.
-   `nfgd`: The number of fine grid points, determining the size of most arrays. This can vary per atom.

## Usage Examples

```fortran
module m_example_pawfgrtab_usage
    use m_pawfgrtab
    use m_libpaw_defs, only: dp
    use m_paral_atom, only: get_my_natom ! For parallel example
    ! Assume an MPI module provides xmpi_world
    implicit none

    subroutine setup_fine_grid_tables(num_total_atoms, atom_types_map, l_max_values_per_atom)
        integer, intent(in) :: num_total_atoms
        integer, intent(in) :: atom_types_map(num_total_atoms)
        integer, intent(in) :: l_max_values_per_atom(num_total_atoms) ! l_size = l_max_val + 1 typically

        type(pawfgrtab_type), allocatable :: local_fgr_tables(:)
        integer :: num_local_atoms, i
        integer, parameter :: cplex_val = 1, nspden_val = 1

        ! Determine local atoms for current MPI process
        call get_my_natom(xmpi_world, num_local_atoms, num_total_atoms)

        if (num_local_atoms > 0) then
            allocate(local_fgr_tables(num_local_atoms))
        else
            allocate(local_fgr_tables(0))
        end if

        ! Initialize pawfgrtab_type for local atoms
        ! (mpi_atmtab would be needed if typat and l_size_atm were for all atoms)
        ! For simplicity, assume atom_types_map and l_max_values_per_atom are already subset for local atoms if parallel.
        call pawfgrtab_init(Pawfgrtab=local_fgr_tables, cplex=cplex_val, &
                            l_size_atm=l_max_values_per_atom(1:num_local_atoms), & ! This should be l_size, not l_max
                            nspden=nspden_val, typat=atom_types_map(1:num_local_atoms), &
                            comm_atom=xmpi_world)
                            ! Pass mpi_atmtab if typat/l_size_atm are global

        ! At this point, structures are initialized but arrays like rfgd, gylm are empty (size 0).
        ! Other routines (e.g., from m_paw_finegrid) would be called to populate them:
        ! e.g., call pawrfgd_fft(local_fgr_tables(i)%ifftsph, ..., local_fgr_tables(i)%rfgd, ...)
        ! then, local_fgr_tables(i)%rfgd_allocated would be set to 1 or 2.

        if (num_local_atoms > 0) then
            call pawfgrtab_print(local_fgr_tables(1:1), natom=1, unit=6, prtvol=1, mode_paral="COLL")
        end if

        call pawfgrtab_free(local_fgr_tables)
        if (allocated(local_fgr_tables)) deallocate(local_fgr_tables)

    end subroutine setup_fine_grid_tables

end module m_example_pawfgrtab_usage
```

## Dependencies and Interactions

-   **`m_libpaw_defs` (`USE_DEFS`)**: For `dp`.
-   **`m_libpaw_tools` (`USE_MSG_HANDLING`)**: For error reporting (`LIBPAW_BUG`).
-   **`m_libpaw_mpi` (`USE_MPI_WRAPPERS`)**: For MPI communication routines used in `pawfgrtab_gather` and `pawfgrtab_redistribute`.
-   **`m_paral_atom`**: Provides `get_my_atmtab` and `get_my_natom` for managing atom distribution in parallel operations.
-   **`m_paw_sphharm`**: Routines like `ylm_angular_mesh` and `initylmr` (though called by `m_pawang` which might then be used to fill `pawfgrtab`) are essential for computing some of the data stored in `pawfgrtab_type` (e.g., `gylm`, `gylmgr`). The `pawfgrtab_type` itself doesn't compute these, but stores them.
-   **`m_paw_finegrid`**: Routines like `pawrfgd_fft`, `pawgylm`, `pawexpiqr` are responsible for computing and filling the arrays within `pawfgrtab_type` instances.

The `m_pawfgrtab` module provides the data structures for holding pre-computed, atom-centered quantities on fine real-space grids. This pre-computation and storage are vital for the efficiency of PAW calculations, as these quantities are frequently reused. The MPI routines ensure this data can be handled correctly in parallel environments.
