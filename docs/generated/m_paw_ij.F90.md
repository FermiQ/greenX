# `m_paw_ij.F90`

## Overview

The Fortran module `m_paw_ij` is central to `libPAW` as it defines and manages the `paw_ij_type` derived data type. This type is designed to store various PAW \(D_{ij}\) matrix elements for a single atom. The \(D_{ij}\) matrix elements are crucial components of the PAW method, representing \(\langle \tilde{p}_i | \hat{H} - \hat{H}_{pseudo} | \tilde{p}_j \rangle\) effectively, broken down into different contributions (e.g., kinetic, local potential, XC, Hartree, etc.) or the total value. The indices \(i\) and \(j\) refer to PAW projector/partial-wave channels.

The module provides a suite of public procedures to:
-   Initialize (`paw_ij_init`) arrays of `paw_ij_type`, allocating specific \(D_{ij}\) components based on input flags.
-   Deallocate (`paw_ij_free`) the allocatable members.
-   Nullify (`paw_ij_nullify`) status flags.
-   Copy (`paw_ij_copy`) instances, handling different parallel distributions.
-   Print (`paw_ij_print`) detailed information about the stored \(D_{ij}\) components.
-   Gather (`paw_ij_gather`) distributed `paw_ij_type` arrays.
-   Redistribute (`paw_ij_redistribute`) `paw_ij_type` arrays between different MPI communicator layouts.
-   Reset status flags (`paw_ij_reset_flags`) to indicate that components need recomputation.

The module also includes private helper routines for serializing and deserializing `paw_ij_type` data for efficient MPI communication.

## Key Components

### Derived Type: `paw_ij_type`

-   **Purpose**: Holds various \(D_{ij}\) matrix elements and related data for a single atom.
-   **Integer Scalars (Metadata and Flags)**:
    -   `cplex_dij`: 1 if \(D_{ij}\) are real, 2 if complex (e.g., spin-orbit, non-collinear magnetism).
    -   `has_dij`, `has_dij0`, `has_dijexxc`, `has_dijfock`, `has_dijfr`, `has_dijhartree`, `has_dijhat`, `has_dijnd`, `has_dijso`, `has_dijU`, `has_dijxc`, `has_dijxc_hat`, `has_dijxc_val`, `has_exexch_pot`, `has_pawu_occ`: Integer flags (0: not used/allocated, 1: allocated/to be computed, 2: computed/has data). These control which \(D_{ij}\) components are active.
    -   `itypat`: Atom type index.
    -   `lmn_size`: Number of \((l,m,n)\) projector/partial-wave channels.
    -   `lmn2_size`: `lmn_size * (lmn_size + 1) / 2` (for symmetric storage of \(D_{ij}\)).
    -   `ndij`: Number of spin components for \(D_{ij}\) (e.g., 1 for non-spin-polarized, 2 for spin-polarized collinear, 4 for non-collinear/spinors).
    -   `nspden`: Number of spin-density components (1 or 2).
    -   `nsppol`: Number of independent spin-polarizations.
    -   `qphase`: 1 if no phase, 2 if \(D_{ij}\) includes an \(e^{-i\vec{q}\cdot\vec{r}}\) phase (for response function calculations at \(\vec{q} \neq 0\)).

-   **Real(`dp`) Arrays (Allocatable)**: These store the \(D_{ij}\) matrix elements. The first dimension typically runs up to `cplex_dij * qphase * lmn2_size`, and the second dimension up to `ndij`.
    -   `dij(:,:)`: Total \(D_{ij}\) term.
    -   `dij0(:)`: Atomic (frozen core) part of \(D_{ij}\), read from PAW dataset. Real and spin-independent.
    -   `dijexxc(:,:)`: On-site matrix elements for local exact exchange (Fock operator).
    -   `dijfock(:,:)`: Total Fock exchange contribution to \(D_{ij}\).
    -   `dijfr(:,:)`: Frozen part of \(D_{ij}\) for response function calculations (q-dependent).
    -   `dijhartree(:)`: Hartree contribution to \(D_{ij}\). Spin-independent.
    -   `dijhat(:,:)`: \(\sum_{LM} \int Q_{ij}^{LM} V_{trial}\) term.
    -   `dijnd(:,:)`: Contribution from nuclear dipole moment.
    -   `dijso(:,:)`: Spin-orbit coupling contribution \(\langle \phi_i | L \cdot S | \phi_j \rangle\).
    -   `dijU(:,:)`: DFT+U contribution to \(D_{ij}\).
    -   `dijxc(:,:)`: XC contribution: \(\langle \phi_i | V_{xc}[n_{AE}] | \phi_j \rangle - \langle \tilde{\phi}_i | V_{xc}[\tilde{n}_{PS}] | \tilde{\phi}_j \rangle\).
    -   `dijxc_hat(:,:)`: \(\sum_{LM} \int Q_{ij}^{LM} V_{xc}\) term.
    -   `dijxc_val(:,:)`: Valence-only XC contribution.
    -   `noccmmp(:,:,:,:)`: PAW+U occupation matrix \(n^{spin}_{m,m'}\).
    -   `nocctot(:)`: Trace of the PAW+U occupation matrix.
    -   `vpawx(:,:,:)`: Exact exchange potential component for PAW+(local exact exchange).

### Public Procedures

-   **`paw_ij_init(...)`**: Initializes an array of `paw_ij_type`. Allocates specific `dij*` arrays based on input `has_*` flags and system parameters (spin, SO coupling, DFT+U, etc.). Handles atom parallelization.
-   **`paw_ij_free(Paw_ij)`**: Deallocates all allocatable arrays within `Paw_ij` and resets `has_*` flags to 0.
-   **`paw_ij_nullify(Paw_ij)`**: Resets all `has_*` flags to 0.
-   **`paw_ij_copy(paw_ij_in, paw_ij_cpy, mpi_atmtab, comm_atom)`**: Deep copies data from `paw_ij_in` to `paw_ij_cpy`. Handles allocation in `paw_ij_cpy`. Manages MPI scatter/gather if `comm_atom` implies different distributions.
-   **`paw_ij_print(...)`**: Prints detailed information about the various \(D_{ij}\) components stored in `Paw_ij`, including their values. Output formatting depends on `pawprtvol`, `enunit` (Ha or eV), etc.
-   **`paw_ij_gather(paw_ij_in, paw_ij_gathered, master, comm_atom)`**: Gathers distributed `paw_ij_in` arrays to `paw_ij_gathered` on a `master` process or all processes. Involves serialization/deserialization.
-   **`paw_ij_redistribute(...)`**: Redistributes `paw_ij` arrays between different MPI communicator layouts (`mpi_comm_in` to `mpi_comm_out`). Offers a brute-force (gather+scatter) or an asynchronous communication algorithm.
-   **`paw_ij_reset_flags(Paw_ij, all, dijhartree, self_consistent)`**: Sets specified `has_*` flags to 1, indicating that these components need recomputation. Allows selective or full reset.

### Private Procedures for MPI Communication

-   **`paw_ij_isendreceive_fillbuffer(...)`**: Serializes `paw_ij_type` data (metadata and arrays) into integer and real buffers for MPI sending.
-   **`paw_ij_isendreceive_getbuffer(...)`**: Deserializes data from received MPI buffers into a `paw_ij_type` structure.

## Important Variables/Constants

-   The `has_*` integer flags within `paw_ij_type` are crucial for managing which \(D_{ij}\) components are active, computed, or need allocation.
-   `lmn_size` and `lmn2_size` define the number of PAW channels and unique pairs, respectively.
-   `cplex_dij` and `qphase` determine the complexity and structure of the stored \(D_{ij}\) arrays.

## Usage Examples

```fortran
module m_example_paw_ij_usage
    use m_libpaw_defs, only: dp
    use m_paw_ij
    use m_pawtab, only: pawtab_type ! Assuming these are initialized elsewhere
    use m_libpaw_mpi, only: xmpi_world, xmpi_comm_self ! For MPI context

    implicit none

    subroutine setup_paw_dij_data(num_atoms_total, atom_types, paw_setups)
        integer, intent(in) :: num_atoms_total
        integer, intent(in) :: atom_types(num_atoms_total)
        type(pawtab_type), intent(in) :: paw_setups(:) ! Array, one per atom type

        type(paw_ij_type), allocatable :: local_paw_dij(:)
        integer :: local_num_atoms

        ! Parameters for initialization
        integer, parameter :: cplex_val = 1 ! For q=0
        integer, parameter :: nspinor_val = 1, nsppol_val = 1, nspden_val = 1
        integer, parameter :: pawspnorb_val = 0
        integer, parameter :: request_dij = 1, request_dij0 = 1, request_dijxc = 1

        ! Determine local atoms for current MPI process
        call get_my_natom(xmpi_world, local_num_atoms, num_atoms_total)

        if (local_num_atoms > 0) then
            allocate(local_paw_dij(local_num_atoms))
        else
            allocate(local_paw_dij(0))
        end if

        call paw_ij_init(Paw_ij=local_paw_dij, cplex=cplex_val, &
                         nspinor=nspinor_val, nsppol=nsppol_val, nspden=nspden_val, &
                         pawspnorb=pawspnorb_val, natom=num_atoms_total, ntypat=size(paw_setups), &
                         typat=atom_types, Pawtab=paw_setups, &
                         has_dij=request_dij, has_dij0=request_dij0, has_dijxc=request_dijxc, &
                         comm_atom=xmpi_world)
                         ! mpi_atmtab would be passed if local_paw_dij was for a subset

        ! ... operations using local_paw_dij ...
        ! e.g., compute Dij components and set has_dij=2 etc.

        if (local_num_atoms > 0) then
            call paw_ij_print(local_paw_dij(1:1), unit=6, pawprtvol=1, mode_paral="COLL", natom=1)
        end if

        call paw_ij_free(local_paw_dij)
        if (allocated(local_paw_dij)) deallocate(local_paw_dij)

    end subroutine setup_paw_dij_data

end module m_example_paw_ij_usage
```

## Dependencies and Interactions

-   **`m_libpaw_defs` (`USE_DEFS`)**: For `dp`, `tol8`, etc.
-   **`m_libpaw_tools` (`USE_MSG_HANDLING`)**: For error reporting (`LIBPAW_BUG`) and output (`wrtout`).
-   **`m_libpaw_mpi` (`USE_MPI_WRAPPERS`)**: For MPI communication used in `paw_ij_gather` and `paw_ij_redistribute`.
-   **`m_paral_atom`**: Provides `get_my_atmtab` and `get_my_natom` for managing atom distribution.
-   **`m_pawtab`**: Defines `pawtab_type`, which provides PAW dataset information (like `lmn_size`, `lmn2_size`) necessary for sizing arrays in `paw_ij_type`.
-   **`m_paw_io`**: `pawio_print_ij` is used by `paw_ij_print` to format the output of \(D_{ij}\) matrices.
-   **`m_libpaw_memory` (`USE_MEMORY_PROFILING`)**: For memory allocation macros.

The `paw_ij_type` and its associated routines in `m_paw_ij` are fundamental for storing and manipulating the PAW \(D_{ij}\)補償項 (compensation term) matrix elements, which are essential for calculating various physical properties. The module's MPI capabilities allow these potentially large data structures to be handled efficiently in parallel computations.
