# `m_pawrhoij.F90`

## Overview

The Fortran module `m_pawrhoij` defines and manages the `pawrhoij_type` derived data type. This data structure is fundamental to the Projector Augmented Wave (PAW) method, as it is designed to store the on-site PAW density matrices (also referred to as PAW augmentation occupancies or \(\rho_{ij}\) matrices). For a given atom `a`, these are defined as \(\rho_{ij}^a = \sum_{nk}^{occ} \langle C_{nk} | \tilde{p}_i^a \rangle \langle \tilde{p}_j^a | C_{nk} \rangle\), where \(|\tilde{p}_i^a\rangle\) are the PAW projector functions and \(|C_{nk}\rangle\) are the Kohn-Sham wavefunctions.

The module provides functionalities for:
-   Storing \(\rho_{ij}\) in both full (`rhoij_`) and packed (`rhoijp` with `rhoijselect`) formats.
-   Storing gradients of \(\rho_{ij}\) (`grhoij`), e.g., with respect to atomic positions or strain.
-   Storing residuals of \(\rho_{ij}\) (`rhoijres`) for Self-Consistent Field (SCF) convergence acceleration.
-   Initializing, allocating, deallocating, copying, and managing these data structures.
-   Performing various MPI operations (gather, broadcast, sum, redistribute) for parallel calculations.
-   Reading and writing \(\rho_{ij}\) data, potentially handling different file formats or versions.

## Key Components

### Derived Type: `pawrhoij_type`

-   **Purpose**: Stores the on-site PAW density matrix \(\rho_{ij}\) and related quantities for a single atom.
-   **Integer Scalars**:
    -   `cplex_rhoij`: 1 if \(\rho_{ij}\) are real, 2 if complex (e.g., due to spin-orbit coupling or if `pawcpxocc=2`).
    -   `itypat`: Atom type index.
    -   `lmn_size`: Number of \((l,m,n)\) projector channels.
    -   `lmn2_size`: `lmn_size * (lmn_size + 1) / 2` (for symmetric packed storage).
    -   `lmnmix_sz`: Number of \(\rho_{ij}\) elements included in SCF mixing.
    -   `ngrhoij`: Number of gradient components stored for \(\rho_{ij}\) (e.g., 3 for atomic displacements).
    -   `nrhoijsel`: Number of non-zero (or significant) elements stored in the packed format (`rhoijp`).
    -   `nspden`: Number of spin-density components for \(\rho_{ij}\) (1 for non-magnetic, 2 for collinear spin-polarized, 4 for non-collinear).
    -   `nspinor`: Number of spinor components of the wavefunctions.
    -   `nsppol`: Number of independent spin-polarization components.
    -   `qphase`: 1 if no \(e^{-i\vec{q}\cdot\vec{R}}\) phase, 2 if such a phase is included (for response functions at \(\vec{q} \neq 0\)).
    -   `use_rhoij_`, `use_rhoijp`, `use_rhoijres`: Flags (0 or 1) indicating if the corresponding array (`rhoij_`, `rhoijp`/`rhoijselect`, `rhoijres`) is allocated and in use.
-   **Integer Arrays (Allocatable)**:
    -   `kpawmix(lmnmix_sz)`: Indices of \(\rho_{ij}\) elements that are part of an SCF mixing scheme.
    -   `rhoijselect(lmn2_size)`: Stores the packed indices \(k\) corresponding to the full \( (i,j) \) pair index for non-zero \(\rho_{ij}\) elements.
-   **Real(`dp`) Arrays (Allocatable)**:
    -   `grhoij(ngrhoij, cplex_rhoij*qphase*lmn2_size, nspden)`: Gradients of \(\rho_{ij}\). Stored in non-packed format.
    -   `rhoij_(cplex_rhoij*qphase*lmn2_size, nspden)`: \(\rho_{ij}\) stored in a non-packed (full matrix, though still using `lmn2_size` for symmetry) format.
    -   `rhoijp(cplex_rhoij*qphase*nrhoijsel, nspden)`: \(\rho_{ij}\) stored in packed format (only non-zero elements).
    -   `rhoijres(cplex_rhoij*qphase*lmn2_size, nspden)`: Residuals of \(\rho_{ij}\) for SCF mixing. Stored in non-packed format.
    -   *Storage convention for complex and q-phase*: The first dimension combines these. If complex, real and imaginary parts are interleaved. If q-phased, components related to \(\cos(\vec{q}\cdot\vec{R})\) and \(\sin(\vec{q}\cdot\vec{R})\) are typically stored contiguously.
    -   *Storage convention for spin*: The last dimension `nspden` handles spin. For `nspden=4` (non-collinear), components are typically total density, x-magnetization, y-magnetization, z-magnetization.

### Public Procedures

-   **`pawrhoij_alloc(...)`**: Initializes and allocates members of `pawrhoij_type` array elements.
-   **`pawrhoij_free(pawrhoij)`**: Deallocates all allocatable arrays within `pawrhoij` elements and resets flags.
-   **`pawrhoij_nullify(pawrhoij)`**: Resets scalar flags to their default (uninitialized) state.
-   **`pawrhoij_copy(...)`**: Performs a deep copy, with options to preserve certain dimensionalities of the target if they differ from the source (e.g., `keep_cplex`, `keep_qphase`). Handles MPI distribution.
-   **`pawrhoij_gather(...)`**: Gathers distributed `pawrhoij_type` arrays to all processes in an MPI communicator. Allows selective gathering of components (e.g., with/without `grhoij`).
-   **`pawrhoij_bcast(...)`**: Broadcasts `pawrhoij_type` data from a master process.
-   **`pawrhoij_redistribute(...)`**: Redistributes `pawrhoij_type` arrays between different MPI distributions.
-   **`pawrhoij_io(...)`**: Handles reading/writing of `pawrhoij_type` data from/to a file unit. Supports different file formats (formatted, unformatted binary, NetCDF) and can handle different header versions.
-   **`pawrhoij_unpack(rhoij)`**: Converts \(\rho_{ij}\) from packed storage (`rhoijp`, `rhoijselect`) to full storage (`rhoij_`). Allocates `rhoij_` if not already.
-   **`pawrhoij_init_unpacked(rhoij)`**: Allocates and zeros the `rhoij_` array for full storage.
-   **`pawrhoij_free_unpacked(rhoij)`**: Deallocates the `rhoij_` array.
-   **`pawrhoij_filter(...)`**: Filters elements from a full `rhoij_input` (or `rhoij` itself) to populate the packed `rhoijp` and `rhoijselect` arrays, keeping only elements above `tol_rhoij`.
-   **`pawrhoij_inquire_dim(...)`**: Helper to determine the expected `cplex_rhoij`, `qphase_rhoij`, `nspden_rhoij` based on calculation parameters like `cpxocc`, `qpt`, `spnorb`.
-   **`pawrhoij_print_rhoij(...)`**: Prints the \(\rho_{ij}\) matrix (from packed or unpacked storage) using `pawio_print_ij`.
-   **`pawrhoij_symrhoij(...)`**: Symmetrizes \(\rho_{ij}\) (and its gradients if `choice > 1`) using crystal symmetries. It can also compute residuals. Takes unsymmetrized data from `pawrhoij_unsym`.
-   **`pawrhoij_mpisum_unpacked(pawrhoij, comm1, comm2)` (Interface for 1D/2D arrays)**: Performs an MPI sum reduction on the `rhoij_` (unpacked) component of `pawrhoij_type` arrays.

### Private Procedures for MPI Communication

-   `pawrhoij_isendreceive_fillbuffer(...)`: Serializes `pawrhoij_type` data into integer and real buffers for MPI sending.
-   `pawrhoij_isendreceive_getbuffer(...)`: Deserializes data from received MPI buffers into `pawrhoij_type`.

## Important Variables/Constants

-   The various `use_*` flags and `*_allocated` members in `pawrhoij_type` are crucial for managing which arrays are active and valid.
-   `tol_rhoij`: A tolerance (parameter `tol10` from `m_libpaw_defs`) used in `pawrhoij_filter`.

## Usage Examples

```fortran
module m_example_pawrhoij_usage
    use m_pawrhoij
    use m_libpaw_defs, only: dp
    use m_pawtab, only: pawtab_type ! For lmn_size
    ! Assume an MPI module provides xmpi_world
    implicit none

    subroutine manage_rhoij_data(num_local_atoms, num_total_atoms, &
                                 local_atom_types, local_lmn_sizes, &
                                 cplex_val, qphase_val, nspden_val, nspinor_val, nsppol_val)
        integer, intent(in) :: num_local_atoms, num_total_atoms
        integer, intent(in) :: local_atom_types(num_local_atoms)
        integer, intent(in) :: local_lmn_sizes(num_local_atoms) ! LMN size for each local atom
        integer, intent(in) :: cplex_val, qphase_val, nspden_val, nspinor_val, nsppol_val

        type(pawrhoij_type), allocatable :: local_rhoij_data(:)
        integer :: i

        allocate(local_rhoij_data(num_local_atoms))

        ! Allocate based on lmn_sizes which would come from pawtab for each itypat
        ! This is a simplified allocation call; typically pawtab would be passed.
        call pawrhoij_alloc(pawrhoij=local_rhoij_data, cplex_rhoij=cplex_val, &
                            nspden=nspden_val, nspinor=nspinor_val, nsppol=nsppol_val, &
                            typat=local_atom_types, lmnsize=local_lmn_sizes, &
                            qphase=qphase_val, use_rhoijp=1, use_rhoijres=1)
                            ! comm_atom and mpi_atmtab if natom > num_local_atoms

        ! ... populate pawrhoij_data(i)%rhoijp and %rhoijselect ...
        ! For example, after computing from wavefunction projections.

        ! Symmetrize (assuming pawrhoij_unsym is same as pawrhoij for in-place)
        ! call pawrhoij_symrhoij(pawrhoij=local_rhoij_data, pawrhoij_unsym=local_rhoij_data, &
        !                        choice=1, ... other sym args ... )

        ! Print for one atom
        if (num_local_atoms > 0) then
            call pawrhoij_print_rhoij(local_rhoij_data(1)%rhoijp, &
                                      cplex=local_rhoij_data(1)%cplex_rhoij, &
                                      qphase=local_rhoij_data(1)%qphase, &
                                      iatom=1, natom=num_local_atoms, &
                                      rhoijselect=local_rhoij_data(1)%rhoijselect)
        end if

        call pawrhoij_free(local_rhoij_data)
        deallocate(local_rhoij_data)

    end subroutine manage_rhoij_data

end module m_example_pawrhoij_usage
```

## Dependencies and Interactions

-   **`m_libpaw_defs`, `m_libpaw_tools`, `m_libpaw_mpi`, `m_libpaw_memory`**: For basic definitions, utilities, MPI, and memory.
-   **`m_paw_io`**: `pawio_print_ij` is used by `pawrhoij_print_rhoij`.
-   **`m_pawang`**: `pawang_type` is used by `pawrhoij_symrhoij` to access symmetry rotation matrices (`zarot`).
-   **`m_pawtab`**: `pawtab_type` is used by `pawrhoij_alloc` (if `lmnsize` is not given) and `pawrhoij_symrhoij` (for `indlmn`).
-   **`m_paral_atom`**: For atom distribution logic in MPI operations.
-   **NetCDF library**: Optionally used by `pawrhoij_io` if `LIBPAW_HAVE_NETCDF` is defined and `form="netcdf"`.

This module is essential for managing the PAW on-site density matrices, which are key intermediates in calculating various physical quantities and updating potentials within the PAW method. The support for packed storage, gradients, residuals, and extensive MPI operations makes it versatile for different types of PAW-based simulations.
