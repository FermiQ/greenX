# `m_paral_atom.F90`

## Overview

The Fortran module `m_paral_atom` provides routines to manage the parallel distribution of atomic sites (atoms) across different MPI processes. This is a common parallelization strategy where calculations related to specific atoms are assigned to specific processes. The module helps in determining which atoms are local to the current process and how atoms are generally mapped to processes.

## Key Components

### Public Procedures

-   **`get_my_natom(comm_atom, my_natom, natom)`**:
    -   **Purpose**: Calculates the number of atoms (`my_natom`) that are treated by the current MPI process.
    -   **Inputs**:
        -   `comm_atom`: Integer, the MPI communicator used for distributing atoms.
        -   `natom`: Integer, the total number of atoms in the system.
    -   **Output**:
        -   `my_natom`: Integer, the number of atoms assigned to the current process.
    -   **Logic**: If not a serial run (`comm_atom` is not `xmpi_comm_self` or `xmpi_comm_null`), it divides `natom` by the number of processes (`nproc` in `comm_atom`). The remainder (`mod(natom, nproc)`) is distributed one by one to the processes with rank `0` to `mod(natom, nproc)-1`.

-   **`get_my_atmtab(comm_atom, my_atmtab, my_atmtab_allocated, paral_atom, natom, my_natom_ref)`**:
    -   **Purpose**: Determines and allocates a table (`my_atmtab`) containing the global indices (1-based) of atoms treated by the current MPI process.
    -   **Inputs**:
        -   `comm_atom`: MPI communicator for atom parallelization.
        -   `paral_atom`: Logical (intent inout), indicates if atom parallelization is active. It can be set to `.false.` by this routine if `comm_atom` implies a serial run or `nproc` is 1.
        -   `natom`: Total number of atoms.
        -   `my_natom_ref` (optional): Integer, a reference value for the local number of atoms (for checking consistency).
    -   **Outputs**:
        -   `my_atmtab`: Integer pointer array, allocated to the size of `my_natom` and filled with global atom indices.
        -   `my_atmtab_allocated`: Logical, set to `.true.` if `my_atmtab` is allocated.
        -   `paral_atom`: Updated based on actual parallelism.
    -   **Logic**:
        -   Only proceeds if `paral_atom` is initially true and `comm_atom` indicates a parallel environment.
        -   Calculates `my_natom` (local number of atoms) similar to `get_my_natom`.
        -   If `my_natom > 0` and `my_atmtab` is not already associated, it allocates `my_atmtab`.
        -   Distributes atoms contiguously: each process gets a block of atoms. The first `mod(natom, nproc)` processes get `natom/nproc + 1` atoms, and the rest get `natom/nproc` atoms.
        -   Calculates `natom_bef` (number of atoms on preceding processes) to determine the starting global index for the current process.
        -   Fills `my_atmtab` with `iatom + natom_bef`.
        -   If `my_natom_ref` is provided, it checks if it matches the computed `size(my_atmtab)`.

-   **`free_my_atmtab(my_atmtab, my_atmtab_allocated)`**:
    -   **Purpose**: Deallocates the `my_atmtab` pointer array and resets `my_atmtab_allocated` to `.false.`.
    -   **Inputs/Outputs**:
        -   `my_atmtab`: Integer pointer array to be deallocated.
        -   `my_atmtab_allocated`: Logical flag to be updated.

-   **`get_proc_atmtab(iproc, atmtab, natom_out, natom, comm_atom_size)`**:
    -   **Purpose**: For a given processor rank `iproc` (0-based), determines the list of global atom indices (`atmtab`) it treats and the count (`natom_out`).
    -   **Inputs**:
        -   `iproc`: Integer, the rank of the processor for which the atom table is requested.
        -   `natom`: Integer, total number of atoms.
        -   `comm_atom_size`: Integer, the total number of processes in the atom communicator.
    -   **Outputs**:
        -   `atmtab`: Integer allocatable array, filled with global atom indices for processor `iproc`.
        -   `natom_out`: Integer, number of atoms on processor `iproc`.
    -   **Logic**: Similar atom distribution logic as in `get_my_atmtab`, but generalized for any `iproc`.

-   **`get_atm_proc(atom_list, natom, nproc, proc_list)`**:
    -   **Purpose**: For a given list of global atom indices (`atom_list`), determines the MPI process rank (`proc_list`) that treats each atom.
    -   **Inputs**:
        -   `atom_list`: Integer array, list of global atom indices.
        -   `natom`: Integer, total number of atoms.
        -   `nproc`: Integer, total number of processes in the atom communicator.
    -   **Output**:
        -   `proc_list`: Integer array, same size as `atom_list`, filled with the rank of the process responsible for each corresponding atom.
    -   **Logic**: Inverse of the logic in `get_my_atmtab` / `get_proc_atmtab`. It calculates which process rank would be assigned a given global atom index based on the contiguous block distribution.
        -   `nmod = mod(natom, nproc)`
        -   `dn = natom / nproc` (base number of atoms per proc)
        -   `dn1 = natom / nproc + 1` (number of atoms for procs getting an extra one)
        -   `jproclim = nmod - 1` (highest rank of procs that get `dn1` atoms)
        -   `natomlim = dn1 * (jproclim + 1)` (total atoms covered by procs with `dn1` atoms)
        -   For each `iatom` in `atom_list`:
            -   If `iatom <= natomlim`, it belongs to a proc in the `0..jproclim` range: `proc_list(iatm) = (iatom - 1) / dn1`.
            -   Else, it belongs to a later proc: `proc_list(iatm) = jproclim + 1 + (iatom - 1 - natomlim) / dn`.

## Important Variables/Constants

-   This module does not define its own public constants but uses those from `m_libpaw_defs` (implicitly via `USE_DEFS`) and `m_libpaw_mpi` (e.g., `xmpi_comm_self`, `xmpi_comm_null`, `xmpi_comm_size`, `xmpi_comm_rank`).

## Usage Examples

```fortran
module m_example_atom_paral
    use m_paral_atom
    use m_libpaw_defs, only: dp
    use m_libpaw_mpi, only: xmpi_world, xmpi_comm_rank, xmpi_comm_size
    implicit none

    subroutine distribute_atomic_tasks
        integer :: total_atoms = 100
        integer :: my_rank, num_procs
        integer :: local_natoms
        integer, pointer :: local_atom_indices(:)
        logical :: is_parallel_over_atoms = .true.
        logical :: atoms_allocated = .false.
        integer, allocatable :: atoms_on_proc0(:)
        integer :: natoms_on_proc0
        integer, dimension(3) :: query_atoms = [5, 50, 95]
        integer, dimension(3) :: proc_for_atoms

        my_rank = xmpi_comm_rank(xmpi_world)
        num_procs = xmpi_comm_size(xmpi_world)

        ! Get number of atoms for current process
        call get_my_natom(xmpi_world, local_natoms, total_atoms)
        print *, "Process ", my_rank, " has ", local_natoms, " atoms."

        ! Get table of atoms for current process
        call get_my_atmtab(xmpi_world, local_atom_indices, atoms_allocated, &
                           is_parallel_over_atoms, total_atoms)

        if (atoms_allocated) then
            print *, "Process ", my_rank, " atoms: ", local_atom_indices
            call free_my_atmtab(local_atom_indices, atoms_allocated)
        else if (local_natoms > 0 .and. .not. is_parallel_over_atoms) then
             ! Serial case, all atoms are local implicitly
             print *, "Process ", my_rank, " (serial) has all ", total_atoms, " atoms."
        end if

        ! Get atoms for a specific process (e.g., process 0)
        if (my_rank == 0) then ! Only one proc needs to do this
            call get_proc_atmtab(0, atoms_on_proc0, natoms_on_proc0, total_atoms, num_procs)
            print *, "Process 0 is expected to have atoms: ", atoms_on_proc0
            if (allocated(atoms_on_proc0)) deallocate(atoms_on_proc0)
        end if

        ! Find which process handles specific atoms
        call get_atm_proc(query_atoms, total_atoms, num_procs, proc_for_atoms)
        if (my_rank == 0) then
            do i = 1, size(query_atoms)
                print *, "Atom ", query_atoms(i), " is on process ", proc_for_atoms(i)
            end do
        end if

    end subroutine distribute_atomic_tasks

end module m_example_atom_paral
```

## Dependencies and Interactions

-   **`libpaw.h`**: Included for preprocessor directives.
-   **`m_libpaw_defs` (`USE_DEFS`)**: Provides basic type kinds (e.g., `dp`).
-   **`m_libpaw_tools` (`USE_MSG_HANDLING`)**: For error reporting (`LIBPAW_BUG`).
-   **`m_libpaw_mpi` (`USE_MPI_WRAPPERS`)**: Provides MPI communication primitives like `xmpi_comm_size`, `xmpi_comm_rank`, `xmpi_comm_self`, `xmpi_comm_null`.
-   **`m_libpaw_memory` (`USE_MEMORY_PROFILING`)**: Provides memory allocation macros like `LIBPAW_POINTER_ALLOCATE`, `LIBPAW_POINTER_DEALLOCATE`.

This module is fundamental for any `libPAW` calculation that aims to distribute work over atoms in an MPI parallel environment. The consistency of the atom distribution logic across these routines is key to correctly assigning and retrieving atom-specific data.
