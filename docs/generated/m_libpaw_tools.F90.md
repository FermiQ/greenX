# `m_libpaw_tools.F90`

## Overview

The Fortran module `m_libpaw_tools` provides a suite of utility functions essential for the `libPAW` library. These tools cover a range of functionalities including message printing with MPI awareness, error handling, string manipulation, and basic I/O operations. Many of these routines are adapted from ABINIT's utility modules, ensuring consistent behavior when `libPAW` is used outside of the ABINIT environment. The module's behavior can be modified by preprocessor directives defined in `libpaw.h` (e.g., `HAVE_YAML`, `LIBPAW_HAVE_NETCDF`).

## Key Components

### Message Handling and Output

-   **`libpaw_wrtout(unit, msg, mode_paral)`**:
    -   Organizes sequential or parallel output of messages.
    -   `unit`: Fortran unit number for output.
    -   `msg`: Character string, the message to write.
    -   `mode_paral` (optional):
        -   `'COLL'` (default): All MPI processes call with the same message; only the master process prints.
        -   `'PERS'`: Each process calls with potentially different messages, or a single process calls.
        -   `'INIT'`: Used to set which rank is the master for `'COLL'` mode (master rank passed in `unit`).
-   **`libpaw_msg_hndl(msg, level, mode_paral, file, line)`**:
    -   Basic error and message handler.
    -   `msg`: String, information about the event.
    -   `level`: String, type of event (`'COMMENT'`, `'WARNING'`, `'ERROR'`, `'BUG'`).
    -   `mode_paral`: `'COLL'` or `'PERS'`.
    -   `file` (optional): Source file name where the event occurred.
    -   `line` (optional): Line number in the source file.
    -   Formats the message and uses `libpaw_wrtout`. For `'ERROR'` or `'BUG'`, it calls `libpaw_leave` to terminate.
-   **`libpaw_flush(unit)`**:
    -   Wrapper for the Fortran `flush` intrinsic, if available (controlled by `HAVE_FC_FLUSH` or `HAVE_FC_FLUSH_`).
-   **`libpaw_spmsg_getcount(ncomment, nwarning, nexit)`**:
    -   Retrieves the counts of `'COMMENT'` and `'WARNING'` messages printed, and an exit flag.
-   **`libpaw_spmsg_mpisum(mpicomm)`**:
    -   Reduces the special message counts (`COMMENT`, `WARNING`, `EXIT`) over an MPI communicator using `xpaw_mpi_sum`.
-   **`libpaw_write_comm_set(new_write_comm)`**:
    -   Sets the MPI communicator (`LIBPAW_WRITE_COMM`) to be used by `libpaw_wrtout`.
-   **`libpaw_log_flag_set(log_flag)`**:
    -   Sets a flag (`LIBPAW_HAS_LOG_FILE`) to enable/disable writing to `std_out`.
-   **`libpaw_netcdf_check(ncerr, msg, file, line)`**:
    -   Error handler for NetCDF calls. If `ncerr` indicates an error, it formats a message (including `nf90_strerror`) and calls `libpaw_msg_hndl` with level `'ERROR'`.
-   **`libpaw_leave(mode_paral, exit_status)`**:
    -   Routine for a "clean" exit, especially in parallel. Calls `xpaw_mpi_abort`.
-   **`libpaw_die(message, file, line)`**:
    -   Stops execution in case of unexpected severe errors, reporting file, line, and MPI rank. Calls `libpaw_leave`.

### String Handling

-   **`libpaw_basename(istr)`**: Pure function, returns the final component of a pathname string.
-   **`libpaw_to_upper(istr)`**: Pure function, converts an input string to uppercase.
-   **`libpaw_lstrip(istr)`**: Pure function, removes leading spaces from a string.
-   **`libpaw_indent(istr)`**: Pure function, indents a string (adds spaces at the beginning of each line, where lines are delimited by `char(10)`).

### I/O Tools

-   **`libpaw_get_free_unit()`**: Integer function, finds and returns an available Fortran I/O unit number. Returns -1 if none are found in the range (10 to 1024, or 10 to 64 for NAG compiler).
-   **`libpaw_lock_and_write(filename, string)`**:
    -   Writes a string to a specified file using a simple lock file mechanism (`filename.lock`) to prevent race conditions (though this is a basic form of locking).

### Private Variables

-   `LIBPAW_WRITE_COMM`: Integer, stores the MPI communicator handle for `libpaw_wrtout`. Default is `xpaw_mpi_world`.
-   `LIBPAW_COMMENT_COUNT`, `LIBPAW_WARNING_COUNT`: Integers, count printed comments and warnings.
-   `LIBPAW_EXIT_FLAG`: Integer, flag set if an exit is requested.
-   `LIBPAW_HAS_LOG_FILE`: Logical, controls output to `std_out`. Default is `.TRUE.`.
-   `LIBPAW_MPIABORTFILE`: Character parameter, name of a file used in `libpaw_msg_hndl` for MPI abort synchronization.

## Important Variables/Constants

-   `std_out`, `std_err`: Integer constants from `m_libpaw_defs` for standard output/error units.
-   `ch10`: Character constant for newline from `m_libpaw_defs`.
-   The module uses preprocessor directives from `libpaw.h` such as `HAVE_YAML` (to enable/disable YAML-formatted comments via `yaml_comment`) and `LIBPAW_HAVE_NETCDF` (to enable NetCDF features like `nf90_strerror`).

## Usage Examples

```fortran
module m_example_using_tools
    use m_libpaw_tools
    use m_libpaw_defs, only: std_out, dp
    implicit none

    subroutine my_routine(condition_is_bad, data_path)
        logical, intent(in) :: condition_is_bad
        character(len=*), intent(in) :: data_path
        character(len=100) :: base_name_val
        integer :: unit_num

        call libpaw_wrtout(std_out, "Starting my_routine.", "COLL") ! All MPI procs, master prints

        base_name_val = libpaw_basename(data_path)
        call libpaw_wrtout(std_out, "Basename: " // trim(base_name_val), "PERS") ! Each proc prints (if data_path differs)

        if (condition_is_bad) then
            call libpaw_msg_hndl("A bad condition was met.", "ERROR", "COLL", &
                                 __FILE__, __LINE__) ! __FILE__, __LINE__ are compiler specific for source info
            ! Execution will likely stop here
        end if

        unit_num = libpaw_get_free_unit()
        if (unit_num /= -1) then
            ! Use unit_num for I/O
            ! close(unit_num)
        else
            call libpaw_msg_hndl("No free I/O unit available!", "WARNING", "COLL")
        end if

        call libpaw_wrtout(std_out, "Finished my_routine.", "COLL")
        call libpaw_flush(std_out)

    end subroutine my_routine

end module m_example_using_tools
```

## Dependencies and Interactions

-   **`libpaw.h`**: Provides preprocessor definitions.
-   **`m_libpaw_defs`**: Provides basic constants like `std_out`, `ch10`, `dp`.
-   **`m_libpaw_mpi` (via `USE_MPI_WRAPPERS`)**: Provides MPI communication wrappers like `xpaw_mpi_comm_size`, `xpaw_mpi_comm_rank`, `xpaw_mpi_sum`, `xpaw_mpi_abort`. This is crucial for parallel-aware output and error handling.
-   **YAML library (`yaml_output`)**: If `HAVE_YAML` is defined, `libpaw_write_lines` will attempt to format comments for YAML output using `yaml_comment`.
-   **NetCDF library (`netcdf`)**: If `LIBPAW_HAVE_NETCDF` is defined, `libpaw_netcdf_check` uses `nf90_strerror`.
-   **Compiler specifics**:
    -   The `flush` capability depends on `HAVE_FC_FLUSH` or `HAVE_FC_FLUSH_`.
    -   The behavior of `exit` in `xpaw_mpi_abort` (called by `libpaw_leave`) can depend on compiler flags like `FC_NAG` or `HAVE_FC_EXIT`.
    -   Use of `__FILE__` and `__LINE__` for source information in `libpaw_msg_hndl` relies on compiler support for these common preprocessor macros.

This module provides essential low-level utilities that promote robust logging, error reporting, and basic operations within `libPAW`, adapting to different compilation environments and parallel execution modes.
