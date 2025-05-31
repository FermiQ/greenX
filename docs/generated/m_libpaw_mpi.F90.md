# `m_libpaw_mpi.F90`

## Overview

The Fortran module `m_libpaw_mpi` serves as a wrapper for MPI (Message Passing Interface) library calls within the `libPAW` ecosystem. Its primary purpose is to provide a consistent interface to MPI functionalities, handling cases where MPI might not be available (falling back to serial behavior) and abstracting differences between MPI versions or specific MPI calls. This module is crucial for enabling parallel execution of `libPAW` in distributed memory environments.

The module defines:
-   MPI related constants (e.g., `xpaw_mpi_world`).
-   Wrappers for basic MPI routines (e.g., `xpaw_mpi_comm_rank`, `xpaw_mpi_barrier`).
-   Generic interfaces for various MPI collective and point-to-point communication patterns, specialized for different data types (integer, double precision) and array dimensions (1D, 2D, 3D).

The inclusion of MPI functionalities is controlled by preprocessor directives like `HAVE_MPI`, `HAVE_MPI1`, `HAVE_MPI2`, which are typically set during the build process based on the system's MPI availability.

## Key Components

### MPI Constants

-   `xpaw_mpi_world`: Integer parameter, defaults to `MPI_COMM_WORLD` if MPI is available, else `0`.
-   `xpaw_mpi_comm_self`: Integer parameter, defaults to `MPI_COMM_SELF` if MPI is available, else `0`.
-   `xpaw_mpi_comm_null`: Integer parameter, defaults to `MPI_COMM_NULL` if MPI is available, else `0`.

### Public Procedures (Wrappers for MPI routines)

-   **`xpaw_mpi_abort(comm, mpierr, msg, exit_status)`**: Wraps `MPI_ABORT`. Terminates MPI execution across the specified communicator.
-   **`xpaw_mpi_comm_rank(comm)`**: Integer function, wraps `MPI_COMM_RANK`. Returns the rank of the calling process in the communicator. Returns `0` if no MPI.
-   **`xpaw_mpi_comm_size(comm)`**: Integer function, wraps `MPI_COMM_SIZE`. Returns the total number of processes in the communicator. Returns `1` if no MPI.
-   **`xpaw_mpi_barrier(comm)`**: Subroutine, wraps `MPI_BARRIER`. Blocks until all processes in the communicator have reached this routine. No-op if no MPI or communicator size is 1.
-   **`xpaw_mpi_wait(request, mpierr)`**: Subroutine, wraps `MPI_WAIT`. Waits for a non-blocking operation (identified by `request`) to complete.
-   **`xpaw_mpi_waitall(array_of_requests, mpierr)`**: Subroutine, wraps `MPI_WAITALL`. Waits for all specified non-blocking operations to complete.
-   **`xpaw_mpi_iprobe(source, tag, mpicomm, flag, mpierr)`**: Subroutine, wraps `MPI_IPROBE`. Checks for an incoming message without actually receiving it.

### Generic Interfaces for MPI Communication Patterns

The module defines several generic interfaces (e.g., `xpaw_mpi_allgather`, `xpaw_mpi_bcast`) that map to specific subroutines based on the data type and rank of the arrays being communicated. These specific subroutines then call the corresponding MPI primitives.

Examples of wrapped MPI routines through these interfaces:
-   `MPI_ALLGATHER` / `MPI_ALLGATHERV`: (e.g., `xpaw_mpi_allgather_dp1d`, `xpaw_mpi_allgatherv_int1d`)
-   `MPI_SCATTERV`: (e.g., `xpaw_mpi_scatterv_dp2d`)
-   `MPI_ALLTOALL` / `MPI_ALLTOALLV`: (e.g., `xpaw_mpi_alltoall_dp1d`, `xpaw_mpi_alltoallv_dp2d`)
-   `MPI_BCAST`: (e.g., `xpaw_mpi_bcast_int`, `xpaw_mpi_bcast_dp3d`)
-   `MPI_GATHER` / `MPI_GATHERV`: (e.g., `xpaw_mpi_gather_dp1d`, `xpaw_mpi_gatherv_int1d`)
-   `MPI_RECV` / `MPI_IRECV` (Receive / Non-blocking Receive): (e.g., `xpaw_mpi_recv_dp2d`, `xpaw_mpi_irecv_dp1d`)
-   `MPI_SEND` / `MPI_ISEND` (Send / Non-blocking Send): (e.g., `xpaw_mpi_send_dp3d`, `xpaw_mpi_isend_dp1d`)
-   `xpaw_mpi_exch`: A custom exchange routine likely implemented using `MPI_SEND` and `MPI_RECV`.
-   `MPI_ALLREDUCE` (with `MPI_SUM` operation): (e.g., `xpaw_mpi_sum_int`, `xpaw_mpi_sum_dp2d`). Adapts for `MPI_IN_PLACE` if `HAVE_MPI2_INPLACE` is defined.

Each specific subroutine (e.g., `xpaw_mpi_bcast_dp1d`) handles the actual MPI call for its corresponding data type and array structure. If MPI is not available (`#else` branch of `#if defined HAVE_MPI`), these routines typically perform a direct copy for operations that would involve data transfer (like collectives) or do nothing for synchronization primitives.

### Private Procedures

-   **`xpaw_mpi_get_tag_ub(comm)`**: Integer function. Retrieves the `MPI_TAG_UB` attribute from a communicator, which is the upper bound for message tags. Provides a default if MPI is not present. This is used by send/recv routines to ensure tag values are valid.

## Important Variables/Constants

-   **Preprocessor Macros**: `HAVE_MPI`, `HAVE_MPI1`, `HAVE_MPI2`, `HAVE_MPI2_INPLACE`. These are critical and determine whether actual MPI calls are made or if serial fallbacks are used. They are defined in `libpaw.h` or a `config.h` file.
-   **MPI Datatype Constants**: `MPI_INTEGER`, `MPI_DOUBLE_PRECISION` are used in MPI calls.
-   **MPI Operation Constants**: `MPI_SUM` for reduction operations.
-   **MPI Status**: `MPI_STATUS_IGNORE` or `status(MPI_STATUS_SIZE)` for message status.

## Usage Examples

Modules within `libPAW` that need to perform parallel communication would `use m_libpaw_mpi`.

```fortran
module m_some_parallel_paw_module
    use m_libpaw_defs, only: dp
    use m_libpaw_mpi
    implicit none

    subroutine parallel_sum_example(local_value, global_sum)
        real(dp), intent(in) :: local_value
        real(dp), intent(out) :: global_sum
        integer :: ier, comm_world

        comm_world = xpaw_mpi_world ! Use the global communicator constant
        global_sum = local_value

        ! Perform a sum reduction across all processes in comm_world
        ! This is a simplified call; the actual xpaw_mpi_sum interface
        ! would take an array. For a scalar, one might wrap it in a 1-element array.
        ! Let's assume a 1D array version for this example:

        real(dp), dimension(1) :: send_buf, recv_buf
        send_buf(1) = local_value
        recv_buf(1) = local_value ! Important for MPI_IN_PLACE or if master only receives

        call xpaw_mpi_sum(recv_buf, comm_world, ier) ! Generic interface
        ! Or more specifically:
        ! call xpaw_mpi_sum_dp1d(recv_buf, comm_world, ier)

        global_sum = recv_buf(1)

        if (ier /= 0) then
            ! Handle MPI error
            call xpaw_mpi_abort(comm_world, ier, "Error in xpaw_mpi_sum")
        end if
    end subroutine parallel_sum_example

    subroutine broadcast_data_example(data_array, root_node, current_rank)
        real(dp), dimension(:), intent(inout) :: data_array
        integer, intent(in) :: root_node, current_rank
        integer :: ier, comm_world

        comm_world = xpaw_mpi_world

        ! Broadcast data_array from root_node to all other nodes
        call xpaw_mpi_bcast(data_array, root_node, comm_world, ier) ! Generic interface
        ! Or more specifically:
        ! call xpaw_mpi_bcast_dp1d(data_array, root_node, comm_world, ier)

        if (ier /= 0) then
             call xpaw_mpi_abort(comm_world, ier, "Error in xpaw_mpi_bcast")
        end if
    end subroutine broadcast_data_example

end module m_some_parallel_paw_module
```

## Dependencies and Interactions

-   **`libpaw.h`**: This header is included and provides the `HAVE_MPI*` preprocessor definitions that control the behavior of this module.
-   **`m_libpaw_defs`**: Provides `dp` (double precision kind) and standard I/O unit numbers (`std_out`, `ab_out`).
-   **MPI Library**: If `HAVE_MPI` is defined, this module directly calls MPI library routines. It may include `mpi.f90` (for `use mpi`) or `mpif.h` depending on the MPI version and configuration.
-   This module is fundamental for any parallel execution strategy within `libPAW`. Other `libPAW` modules that implement parallel algorithms will depend on `m_libpaw_mpi` for communication.
-   The design with generic interfaces allows for type-safe MPI calls for common data types and array structures while keeping the call sites in user code clean.

The module provides a critical abstraction layer, promoting portability and simplifying the use of MPI within `libPAW` by centralizing MPI calls and handling the serial case transparently.
