# `m_paw_an.F90`

## Overview

The Fortran module `m_paw_an` defines and manages the `paw_an_type` derived data type. This data type is central to Projector Augmented Wave (PAW) method calculations, designed to store various physical quantities for a single atom that are represented either on an angular grid or as spherical harmonic (LM) moments. These quantities typically include potentials (Hartree, XC), XC kernels, and their tilde (pseudo) counterparts within the PAW formalism.

The module provides a suite of public procedures to:
-   Initialize (`paw_an_init`) arrays of `paw_an_type`.
-   Deallocate (`paw_an_free`) the allocatable members.
-   Nullify (`paw_an_nullify`) pointers and reset flags (though in Fortran, this primarily resets flags).
-   Copy (`paw_an_copy`) instances of the type, handling different parallel distributions.
-   Print (`paw_an_print`) information about the contents.
-   Gather (`paw_an_gather`) distributed `paw_an_type` arrays to a root process or all processes.
-   Redistribute (`paw_an_redistribute`) `paw_an_type` arrays between different MPI communicator layouts.
-   Reset status flags (`paw_an_reset_flags`) to force recomputation.

The module also includes private helper routines for serializing and deserializing `paw_an_type` data for MPI communication.

## Key Components

### Derived Type: `paw_an_type`

-   **Purpose**: Holds angular grid data or LM-moment expansions for various PAW quantities for a single atom.
-   **Integer Scalars**:
    -   `angl_size`: Dimension of the PAW angular mesh (`ntheta * nphi`).
    -   `cplex`: `1` for real data, `2` for complex data (for potentials/densities).
    -   `has_kxc`, `has_k3xc`, `has_vhartree`, `has_vxc`, `has_vxctau`, `has_vxcval`, `has_vxc_ex`: Integer flags (0: not allocated/used, 1: allocated and to be computed, 2: already computed). These flags track the status of different potential/kernel components.
    -   `itypat`: Atom type index.
    -   `lm_size`: `(l_size)**2`, where `l_size = 2*l_max_gaunt + 1`. `l_max_gaunt` is the max L for non-zero Gaunt coefficients.
    -   `mesh_size`: Dimension of the radial mesh for arrays in this type.
    -   `nkxc1`, `nk3xc1`: Number of independent components for XC kernels \(K_{xc}\) and \(K^{(3)}_{xc}\).
    -   `nspden`: Number of spin-density components (1 or 2).
-   **Logical Arrays**:
    -   `lmselect(lm_size)`: Allocatable. Selects non-zero LM-moments for one-center quantities.
-   **Real(`dp`) Arrays (Allocatable)**: These store the actual physical data. The first dimension is typically `cplex*mesh_size` (radial points, potentially complex). The second dimension is `lm_size` (for LM-moments) or `angl_size` (for values on angular grid), depending on `pawxcdev`. The third dimension is usually `nspden` or a kernel component index.
    -   `kxc1`, `kxct1`: XC kernel \(K_{xc}\) and pseudo XC kernel \(\tilde{K}_{xc}\).
    -   `k3xc1`, `k3xct1`: Derivatives of XC kernel \(K^{(3)}_{xc}\) and \(\tilde{K}^{(3)}_{xc}\).
    -   `vh1`, `vht1`: Hartree potential \(V_H\) and pseudo Hartree potential \(\tilde{V}_H\) (LM-moments).
    -   `vxc1`, `vxct1`: XC potential \(V_{xc}\) and pseudo XC potential \(\tilde{V}_{xc}\).
    -   `vxctau1`, `vxcttau1`: XC potential components related to kinetic energy density \(\tau\), \(V_{xc,\tau}\) and \(\tilde{V}_{xc,\tau}\).
    -   `vxc1_val`, `vxct1_val`: XC potential and pseudo XC potential from valence electrons only.
    -   `vxc_ex`: XC potential for local exact exchange.

### Public Procedures

-   **`paw_an_init(...)`**: Initializes an array of `paw_an_type`. Allocates members based on input parameters (dimensions, which potentials/kernels are needed via `has_*` flags). Handles atom parallelization using `mpi_atmtab` and `comm_atom`.
-   **`paw_an_free(Paw_an)`**: Deallocates all allocatable arrays within each element of the `Paw_an` array and resets `has_*` flags to 0.
-   **`paw_an_nullify(Paw_an)`**: Resets `has_*` flags to 0. (Note: Fortran pointers are nullified upon deallocation or explicitly with `nullify()`. This routine primarily resets status flags).
-   **`paw_an_copy(paw_an_in, paw_an_cpy, mpi_atmtab, comm_atom)`**: Copies data from `paw_an_in` to `paw_an_cpy`. Can handle cases where input is global and output is distributed, or vice-versa, or a simple serial copy. It reallocates arrays in `paw_an_cpy` according to `paw_an_in`'s metadata and copies data if the corresponding `has_*` flag in `paw_an_in` is 2 (computed).
-   **`paw_an_print(...)`**: Prints information about the `Paw_an` data structure, including sizes and status flags for each atom.
-   **`paw_an_gather(Paw_an_in, paw_an_gathered, master, comm_atom, mpi_atmtab)`**: Gathers distributed `Paw_an_in` arrays from all processes in `comm_atom` into `paw_an_gathered` on the `master` process, or all processes if `master == -1` (Allgather). This involves serializing the data into integer and real buffers, communicating, and then deserializing.
-   **`paw_an_redistribute(...)`**: Redistributes `paw_an` arrays from an input MPI distribution (`mpi_comm_in`) to an output distribution (`mpi_comm_out`).
    -   Supports two algorithms:
        1.  Brute-force: Allgather to all processes in `mpi_comm_in`, then each process in `mpi_comm_out` scatters/selects its required data.
        2.  Asynchronous: Uses non-blocking send/receive operations based on precomputed send/receive lists (`SendAtomProc`, `RecvAtomList`, etc.). This is more efficient for complex redistributions.
-   **`paw_an_reset_flags(Paw_an)`**: Sets all `has_*` status flags in `Paw_an` elements to 1 if they were previously > 0, forcing recomputation of the corresponding arrays.

### Private Procedures for MPI Communication

-   **`paw_an_isendreceive_fillbuffer(...)`**: Serializes data from `paw_an_type` instances (selected atoms) into integer (`buf_int`) and real (`buf_dp`) buffers for sending via MPI. It packs metadata (sizes, flags) and actual array data.
-   **`paw_an_isendreceive_getbuffer(...)`**: Deserializes data from received integer and real buffers into the appropriate `paw_an_type` instance. Allocates arrays based on received metadata and fills them.

## Important Variables/Constants

-   The module relies on constants and types from `m_libpaw_defs` (e.g., `dp`), `m_paral_atom` (for atom distribution logic), `m_pawang` (for `Pawang_type`, defining angular grid properties), and `m_pawtab` (for `Pawtab_type`, providing radial grid info per atom type).
-   The `has_*` integer flags within `paw_an_type` are crucial for controlling computation and data copying/communication (0=ignore, 1=allocate/compute, 2=computed/has_data).

## Usage Examples

```fortran
module m_example_paw_an_usage
    use m_libpaw_defs, only: dp
    use m_paw_an
    use m_pawang, only: pawang_type ! Assuming these are initialized elsewhere
    use m_pawtab, only: pawtab_type   ! Assuming these are initialized elsewhere
    use m_libpaw_mpi, only: xmpi_world, xmpi_comm_self ! For MPI context

    implicit none

    subroutine test_paw_an_operations(num_atoms_total, atom_types, paw_angular_grid, paw_radial_tables)
        integer, intent(in) :: num_atoms_total
        integer, intent(in) :: atom_types(num_atoms_total) ! Global types for all atoms
        type(pawang_type), intent(in) :: paw_angular_grid
        type(pawtab_type), intent(in) :: paw_radial_tables(:) ! Array, one per atom type

        type(paw_an_type), allocatable :: local_paw_an_data(:)
        integer :: local_num_atoms, i
        integer, parameter :: default_cplex = 1, default_nspden = 1
        integer, parameter :: default_nkxc1 = 3, default_nk3xc1 = 4
        integer, parameter :: default_pawxcdev = 1 ! 0 for grid, 1 for moments
        integer, optional_has_vxc = 2 ! Request Vxc and mark as "computed" (for example)

        ! Determine local atoms for current MPI process
        call get_my_natom(xmpi_world, local_num_atoms, num_atoms_total)

        if (local_num_atoms > 0) then
            allocate(local_paw_an_data(local_num_atoms))
        else
            allocate(local_paw_an_data(0)) ! Still allocate if no atoms locally
        end if

        ! Initialize paw_an_type for local atoms
        ! mpi_atmtab would be needed if local_paw_an_data was sized for all atoms
        ! but we are sizing it for local_num_atoms. For init, it maps local i to global atom index.
        ! A proper mpi_atmtab needs to be constructed if num_atoms_total > local_num_atoms
        ! For simplicity, assuming serial or that typat is already subset for local atoms if parallel.
        ! A full example would use get_my_atmtab to get the global indices for local atoms.

        ! Simplified: if serial, typat is atom_types. If parallel, typat needs to be the types
        ! of the *local* atoms. The paw_an_init will use mpi_atmtab if provided to map local
        ! Paw_an(i) to global atom typat(mpi_atmtab(i)).

        call paw_an_init(local_paw_an_data, natom=num_atoms_total, ntypat=size(paw_radial_tables), &
                         nkxc1=default_nkxc1, nk3xc1=default_nk3xc1, nspden=default_nspden, &
                         cplex=default_cplex, pawxcdev=default_pawxcdev, typat=atom_types, &
                         Pawang=paw_angular_grid, Pawtab=paw_radial_tables, &
                         has_vxc=optional_has_vxc, comm_atom=xmpi_world)
                         ! mpi_atmtab would be passed if local_paw_an_data was for a subset

        ! Print info (master prints info about the first local atom if any)
        if (local_num_atoms > 0) then
            call paw_an_print(local_paw_an_data(1:1), unit=6, mode_paral="COLL", natom=1)
        end if

        ! ... computations that might fill or use paw_an_data ...
        ! e.g., if has_vxc was 1, compute it and set Paw_an(i)%has_vxc = 2

        ! Clean up
        call paw_an_free(local_paw_an_data)
        if (allocated(local_paw_an_data)) deallocate(local_paw_an_data)

    end subroutine test_paw_an_operations

end module m_example_paw_an_usage
```

## Dependencies and Interactions

-   **`m_libpaw_defs` (`USE_DEFS`)**: For basic kinds (`dp`) and constants.
-   **`m_libpaw_tools` (`USE_MSG_HANDLING`)**: For error reporting (`LIBPAW_BUG`) and output (`wrtout`).
-   **`m_libpaw_mpi` (`USE_MPI_WRAPPERS`)**: For MPI communication primitives (`xmpi_sum`, `xmpi_allgather`, `xmpi_allgatherv`, etc.) used in `paw_an_gather` and `paw_an_redistribute`.
-   **`m_paral_atom`**: Provides `get_my_atmtab` and `get_my_natom` to manage atom distribution in parallel initializations and data transfers.
-   **`m_pawang`**: Defines `pawang_type`, which provides angular grid information needed by `paw_an_init`.
-   **`m_pawtab`**: Defines `pawtab_type`, which provides radial mesh information (per atom type) needed by `paw_an_init`.
-   **`m_libpaw_memory` (`USE_MEMORY_PROFILING`)**: For memory allocation macros (`LIBPAW_ALLOCATE`, `LIBPAW_DEALLOCATE`).

This module is crucial for managing atom-centered angular data within the PAW formalism, providing structured storage and robust mechanisms for parallel data handling. The serialization/deserialization routines (`paw_an_isendreceive_fillbuffer` and `paw_an_isendreceive_getbuffer`) are complex due to the need to pack/unpack various components of `paw_an_type`, including dynamically sized arrays and status flags.
