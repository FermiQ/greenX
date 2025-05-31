# `m_pawcprj.F90`

## Overview

The Fortran module `m_pawcprj` defines and manages the `pawcprj_type` derived data type. This data structure is essential in Projector Augmented Wave (PAW) calculations as it stores the scalar projections of wavefunctions \(|C_{nk}\rangle\) onto the PAW projector functions \(|\tilde{p}_{lmn}^a\rangle\) for a specific atom `a`. These projections, denoted as \(c_{nk}^{a,lmn} = \langle \tilde{p}_{lmn}^a | C_{nk} \rangle\), are fundamental for reconstructing all-electron quantities from pseudo-wavefunctions within the PAW augmentation spheres. The module also supports storing derivatives of these projections, e.g., with respect to atomic positions or cell parameters, which are needed for force calculations or response functions.

The module provides a comprehensive set of public procedures for:
-   Allocating (`pawcprj_alloc`) and deallocating (`pawcprj_free`) memory for `pawcprj_type` arrays.
-   Basic operations like setting to zero (`pawcprj_set_zero`), copying (`pawcprj_copy`), conjugation (`pawcprj_conjg`), and linear combinations (`pawcprj_axpby`, `pawcprj_zaxpby`, `pawcprj_lincom`).
-   Applying symmetry operations (`pawcprj_symkn`).
-   Input/Output, primarily for debugging (`pawcprj_output`, `pawcprj_get`, `pawcprj_put`).
-   MPI communication: gathering, broadcasting, exchanging, summing, and redistributing these projection data across parallel processes (`pawcprj_mpi_*` routines).
-   Packing and unpacking data into simple buffers for efficient communication or storage (`pawcprj_pack`, `pawcprj_unpack`).
-   Calculating on-site overlap contributions (`paw_overlap`).

## Key Components

### Derived Type: `pawcprj_type`

-   **Purpose**: Stores the complex projection coefficients \(\langle \tilde{p}_{lmn} | C_{nk} \rangle\) for one atom and one wavefunction (or band/k-point/spin combination), and optionally their gradients.
-   **Integer Scalars**:
    -   `ncpgr`: Integer, number of gradients stored (e.g., 3 for Cartesian derivatives \(\nabla_{\vec{R}}\)).
    -   `nlmn`: Integer, number of \((l,m,n)\) projector channels for the atom.
-   **Real(`dp`) Arrays (Allocatable)**:
    -   `cp(2, nlmn)`: Stores the projection coefficients. `cp(1,:)` for real parts, `cp(2,:)` for imaginary parts.
    -   `dcp(2, ncpgr, nlmn)`: Stores derivatives of the projection coefficients. `dcp(1,:,:)` for real parts, `dcp(2,:,:)` for imaginary parts.

### Public Procedures

A large number of public procedures are provided, covering various aspects of managing `pawcprj_type` data:

**Lifecycle & Basic Operations:**
-   `pawcprj_alloc(cprj, ncpgr, nlmn)`: Allocates `cp` and `dcp` arrays for each element in `cprj(:,:)`.
-   `pawcprj_free(cprj)`: Deallocates `cp` and `dcp`.
-   `pawcprj_set_zero(cprj)`: Sets `cp` and `dcp` to zero.
-   `pawcprj_copy(cprj_in, cprj_out, icpgr)`: Copies data. `icpgr` allows copying a specific gradient component or only the projections `cp`.
-   `pawcprj_conjg(cprj)`: Conjugates the complex values in `cp` and `dcp` (negates the imaginary part).

**Mathematical Operations:**
-   `pawcprj_axpby(alpha, beta, cprjx, cprjy)`: Performs \(Y \leftarrow \alpha X + \beta Y\) where \(\alpha, \beta\) are real scalars.
-   `pawcprj_zaxpby(alpha, beta, cprjx, cprjy)`: Performs \(Y \leftarrow \alpha X + \beta Y\) where \(\alpha, \beta\) are complex scalars (passed as `real(dp) alpha(2)`).
-   `pawcprj_projbd(alpha, cprjx, cprjy)`: Performs \(Y_j \leftarrow Y_j + \sum_i \alpha_{ij} X_i\). Seems to be a block/matrix version of ZAXPBY.
-   `pawcprj_lincom(alpha, cprj_in, cprj_out, nn)`: Computes a linear combination \(C_{out} = \sum_i \alpha_i C_{in,i}\).
-   `paw_overlap(cprj1, cprj2, typat, pawtab, spinor_comm)`: Real function, computes the on-site contribution to the overlap \(\langle C_1 | S_{PAW} | C_2 \rangle\) between two states represented by `cprj1` and `cprj2`. \(S_{PAW} = 1 + \sum_{ij} |\tilde{p}_i\rangle (Q_{ij} - \delta_{ij}) \langle\tilde{p}_j|\). This function calculates \(\sum_{a,ij} c_{1}^{a,i*} (Q_{ij}^a - \delta_{ij}) c_{2}^{a,j}\).

**Symmetry & Reordering:**
-   `pawcprj_symkn(...)`: Constructs `cprj_fkn` (projections for a k-point) from `cprj_ikn` (projections at a symmetry-related k-point in the IBZ) using symmetry operations.
-   `pawcprj_reorder(cprj, atm_indx)`: Reorders the first dimension (atom index) of the `cprj` array according to `atm_indx`.

**I/O & Debugging:**
-   `pawcprj_output(cprj, prtgrads)`: Prints the content of `cprj` (projections and optionally gradients) to standard output. Useful for debugging.
-   `pawcprj_get(...)`, `pawcprj_put(...)`: Read/write `cprj` data for a given k-point/band from/to memory or a temporary file. Handles MPI distribution and data reordering.

**MPI Communication:**
-   `pawcprj_mpi_allgather(...)`: Performs an `MPI_ALLGATHER` operation.
-   `pawcprj_bcast(...)`: Broadcasts `cprj` data from a master process.
-   `pawcprj_transpose(...)`: Transposes a `cprj` data structure to change parallel distribution (e.g., from atom-distributed to band-distributed).
-   `pawcprj_gather_spin(...)`: Collects spin-distributed `cprj` data.
-   `pawcprj_mpi_exch(...)`: Exchanges `cprj` data between two specific MPI processes.
-   `pawcprj_mpi_send(...)`, `pawcprj_mpi_recv(...)`: Point-to-point send/receive operations.
-   `pawcprj_mpi_sum(...)`: Performs an element-wise sum (MPI_ALLREDUCE with MPI_SUM) of `cprj` data across an MPI communicator.

**Utilities:**
-   `pawcprj_getdim(...)`: Helper function to determine the `nlmn` (number of projectors) for each atom based on atom type and sorting mode.
-   `pawcprj_pack(nlmn, cprj, buffer, buffer_gr)`, `pawcprj_unpack(...)`: Pack `cprj` data into flat real buffers for communication/storage, and unpack back.

## Important Variables/Constants

-   The `pawcprj_type` itself is the main "variable" structure managed by this module.
-   Its components `cp` (coefficients) and `dcp` (derivatives of coefficients) store the core physical data.
-   `nlmn` (number of \(l,m,n\) projectors) and `ncpgr` (number of gradient components) define the dimensions of these arrays for each atom.

## Usage Examples

```fortran
module m_example_pawcprj_usage
    use m_libpaw_defs, only: dp
    use m_pawcprj
    use m_pawtab, only: pawtab_type ! For Pawtab to get nlmn
    ! Assume an MPI module provides xmpi_world, xmpi_comm_rank, etc.
    implicit none

    subroutine process_paw_projections(num_atoms, num_bands, num_kpts, current_k_idx, &
                                       atom_types_global, paw_setups, &
                                       all_wavefunction_projections)
        integer, intent(in) :: num_atoms, num_bands, num_kpts, current_k_idx
        integer, intent(in) :: atom_types_global(num_atoms)
        type(pawtab_type), intent(in) :: paw_setups(:) ! Array for each atom type

        ! Assume all_wavefunction_projections is distributed over k-points,
        ! and for current_k_idx, it holds data for all bands and atoms.
        type(pawcprj_type), intent(in) :: all_wavefunction_projections(num_atoms, num_bands)

        type(pawcprj_type), allocatable :: cprj_one_band(:,:), cprj_another_band(:,:)
        integer, dimension(num_atoms) :: nlmn_for_atoms
        integer :: i, num_gradients = 0 ! Example: no gradients stored/needed
        real(dp) :: overlap_val(2)

        ! Get nlmn for each atom
        call pawcprj_getdim(nlmn_for_atoms, num_atoms, nattyp_dummy, size(paw_setups), &
                            atom_types_global, paw_setups, "Random") ! Use "Ordered" if sorted by type

        allocate(cprj_one_band(num_atoms,1), cprj_another_band(num_atoms,1))
        call pawcprj_alloc(cprj_one_band, num_gradients, nlmn_for_atoms)
        call pawcprj_alloc(cprj_another_band, num_gradients, nlmn_for_atoms)

        ! Copy data for band 1 and band 2 (assuming nspinor=1 for simplicity)
        call pawcprj_copy(all_wavefunction_projections(:,1), cprj_one_band)
        call pawcprj_copy(all_wavefunction_projections(:,2), cprj_another_band)

        ! Compute onsite overlap between these two bands
        overlap_val = paw_overlap(cprj_one_band, cprj_another_band, atom_types_global, paw_setups)
        ! print *, "Onsite overlap for k=", current_k_idx, " bands 1&2 :", overlap_val(1), "(real part)"

        call pawcprj_free(cprj_one_band)
        call pawcprj_free(cprj_another_band)
        deallocate(cprj_one_band, cprj_another_band)
        ! (nattyp_dummy would be needed for pawcprj_getdim if sort_mode="Ordered")
contains
        integer, dimension(size(paw_setups)) :: nattyp_dummy ! Dummy, not used for random sort
    end subroutine process_paw_projections

end module m_example_pawcprj_usage
```

## Dependencies and Interactions

-   **`m_libpaw_defs` (`USE_DEFS`)**: For `dp`, `zero`, `tol16`, etc.
-   **`m_libpaw_tools` (`USE_MSG_HANDLING`)**: For error reporting (`LIBPAW_ERROR`, `LIBPAW_BUG`).
-   **`m_libpaw_mpi` (`USE_MPI_WRAPPERS`)**: For all MPI communication routines (`xmpi_*`).
-   **`m_libpaw_memory` (`USE_MEMORY_PROFILING`)**: For memory allocation macros.
-   **`m_pawtab`**: Defines `pawtab_type`, which is used by `paw_overlap` to get \(S_{ij}\) matrix elements (\(Q_{ij} - \delta_{ij}\)) and by `pawcprj_getdim` to get `lmn_size` per atom type.

This module is fundamental for PAW calculations, as the \(\langle \tilde{p} | C \rangle\) projections are the primary quantities used to calculate expectation values and perform other PAW-specific operations. The extensive MPI support highlights its role in enabling large-scale parallel simulations.
