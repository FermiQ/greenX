# `libpaw.h`

## Overview

`libpaw.h` is a C preprocessor header file crucial for the `libPAW` library. Its primary role is to adapt `libPAW`'s behavior and dependencies based on the host code it's being compiled with. It achieves this by defining a series of macros that abstract functionalities like error handling, memory allocation, MPI communication, and access to external libraries (libXC, NetCDF).

The file is structured with conditional compilation blocks (`#if defined HAVE_LIBPAW_ABINIT`, `#elif defined HAVE_LIBPAW_BIGDFT`, `#else`) to provide specific configurations for:
- **ABINIT**: A plane-wave DFT code.
- **BigDFT**: A wavelet-based DFT code.
- **Default**: A generic setup if neither ABINIT nor BigDFT is specified.

Additionally, it includes a common definitions section, for instance, for NetCDF error checking.

## Key Components

The file doesn't define functions or classes in the traditional sense but rather a set of **macros** that Fortran source files (presumably using a preprocessor like `fpp` or `cpp`) will use. These macros then expand to host-specific Fortran code.

### Configuration Blocks:

1.  **ABINIT Specific (`HAVE_LIBPAW_ABINIT`)**:
    *   Includes `abi_common.h` (ABINIT's common header).
    *   Defines `USE_DEFS` to use `defs_basis`.
    *   Defines `USE_MPI_WRAPPERS` to use `m_xmpi`.
    *   Defines message handling macros (`LIBPAW_COMMENT`, `LIBPAW_WARNING`, etc.) to map to ABINIT's versions (`ABI_COMMENT`, etc.).
    *   Defines memory allocation macros (`LIBPAW_ALLOCATE`, etc.) to map to ABINIT's profiling-aware allocators (`ABI_MALLOC`, `ABI_FREE`).
    *   Conditionally defines `LIBPAW_HAVE_LIBXC` and `LIBPAW_HAVE_NETCDF` based on ABINIT's configuration.
    *   Defines `LIBPAW_CONTIGUOUS` to `ABI_CONTIGUOUS`.

2.  **BigDFT Specific (`HAVE_LIBPAW_BIGDFT`)**:
    *   Defines `USE_DEFS` to use `m_libpaw_defs`.
    *   Defines `USE_MPI_WRAPPERS` to use `m_libpaw_mpi`.
    *   Defines message handling macros to use `libpaw_msg_hndl` from `m_libpaw_tools`.
    *   Defines memory allocation macros to use BigDFT's `f_malloc`, `f_free`, etc. It distinguishes between basic types and user-defined types for allocation.
    *   Conditionally defines `LIBPAW_HAVE_LIBXC`.
    *   Disables NetCDF and FoX support for BigDFT within `libPAW`.
    *   Defines `LIBPAW_CONTIGUOUS` (as empty) and enables `LIBPAW_ISO_C_BINDING`.

3.  **Default Configuration**:
    *   Defines `USE_DEFS` to use `m_libpaw_defs`.
    *   Defines `USE_MPI_WRAPPERS` to use `m_libpaw_mpi`.
    *   Defines message handling macros to use `libpaw_msg_hndl` from `m_libpaw_tools`.
    *   Defines memory allocation macros to use standard Fortran `allocate` and `deallocate`.
    *   Disables libXC, NetCDF, and FoX support by default.
    *   Defines `LIBPAW_CONTIGUOUS` (as empty) and enables `LIBPAW_ISO_C_BINDING`.

### Common Macros (Examples, behavior depends on configuration block):

*   **`USE_DEFS`**: Specifies the Fortran module providing basic definitions.
*   **`USE_MPI_WRAPPERS`**: Specifies the Fortran module for MPI communication.
*   **`USE_MSG_HANDLING`**: Specifies the Fortran module for message and error handling.
*   **`LIBPAW_COMMENT(msg)`**: Macro for issuing a comment message.
*   **`LIBPAW_WARNING(msg)`**: Macro for issuing a warning message.
*   **`LIBPAW_ERROR(msg)`**: Macro for issuing an error message (likely terminating).
*   **`LIBPAW_BUG(msg)`**: Macro for issuing a bug message (likely terminating).
*   **`LIBPAW_ALLOCATE(ARR,SIZE)`**: Macro for allocating memory for an array.
*   **`LIBPAW_DEALLOCATE(ARR)`**: Macro for deallocating memory.
    *   Similar variants exist for pointers (`LIBPAW_POINTER_ALLOCATE`), user-defined types (`LIBPAW_DATATYPE_ALLOCATE`), and arrays with explicit bounds (`LIBPAW_BOUND1_ALLOCATE`).
*   **`BOUNDS(LBND,UBND)`**: Macro to specify array bounds, adapting to host code syntax (e.g., `:` for ABINIT/Default, `.to.` for BigDFT).
*   **`LIBPAW_HAVE_LIBXC`**: Defined if libXC support is enabled.
*   **`LIBPAW_HAVE_NETCDF`**: Defined if NetCDF support is enabled.
*   **`LIBPAW_HAVE_FOX`**: Defined if FoX (Fortran XML library) support is enabled (currently always undefined in this file).
*   **`LIBPAW_CONTIGUOUS`**: Macro for specifying Fortran `contiguous` attribute.
*   **`LIBPAW_ISO_C_BINDING`**: Flag indicating if ISO C Binding features are used.

### Common Definitions:

*   **`NCF_CHECK(ncerr)`**: A macro to check NetCDF error codes, calling `libpaw_netcdf_check` (presumably defined in `m_libpaw_tools` or a host-specific module) if an error occurs.

## Important Variables/Constants

This file does not declare variables or constants in the typical sense. The preprocessor symbols like `HAVE_LIBPAW_ABINIT`, `HAVE_LIBPAW_BIGDFT`, `HAVE_LIBXC`, `HAVE_NETCDF` are critical. These are typically defined at compile time (e.g., via compiler flags `-D`) based on the build configuration of the host code or `libPAW` itself.

## Usage Examples

This header file is not directly "used" like a typical C header that's `#include`d in C source files for function declarations. Instead, it's meant to be processed by a C preprocessor during the compilation of Fortran files.

A Fortran file within `libPAW` might contain:

```fortran
! In some m_paw_*.F90 file
#include "libpaw.h" ! This line is processed by fpp/cpp

module m_some_paw_module
    USE_DEFS
    USE_MSG_HANDLING
    implicit none

contains

    subroutine_some_paw_subroutine(data, size_data)
        real(dp), dimension(:), allocatable :: data
        integer :: size_data
        ! ...
        LIBPAW_ALLOCATE(data, (size_data)) ! Macro expands to host-specific allocation
        ! ...
        if (some_error_condition) then
            LIBPAW_ERROR("An error occurred in some_paw_subroutine") ! Macro for error
        end if
        ! ...
        LIBPAW_DEALLOCATE(data) ! Macro for deallocation
    end subroutine some_paw_subroutine

end module m_some_paw_module
```

The C preprocessor would replace `USE_DEFS`, `LIBPAW_ALLOCATE`, `LIBPAW_ERROR`, `LIBPAW_DEALLOCATE` etc., with the specific Fortran code defined in `libpaw.h` corresponding to the active configuration (ABINIT, BigDFT, or Default).

## Dependencies and Interactions

- **Host Code Configuration**: The content and behavior of `libPAW` Fortran files are heavily influenced by the macros defined in this file, which in turn depend on compile-time definitions like `HAVE_LIBPAW_ABINIT` or `HAVE_LIBPAW_BIGDFT`.
- **`config.h`**: If `HAVE_CONFIG_H` is defined, `config.h` is included. This file usually contains auto-generated configuration details from a build system (like Autotools), including definitions for `HAVE_LIBXC`, `HAVE_NETCDF`, etc.
- **ABINIT specific**:
    - `abi_common.h`: Provides ABINIT's common macros and definitions.
    - `defs_basis`, `m_xmpi`, `m_errors`, `m_abicore`, `m_profiling_abi`: ABINIT Fortran modules.
- **BigDFT specific**:
    - `m_libpaw_defs`, `m_libpaw_mpi`, `m_libpaw_tools`: `libPAW`'s own Fortran modules for when used with BigDFT (or default).
    - `dynamic_memory`: BigDFT's Fortran module for memory allocation.
- **Default configuration**:
    - Relies on `libPAW`'s internal modules (`m_libpaw_defs`, `m_libpaw_mpi`, `m_libpaw_tools`) and standard Fortran features.
- **External Libraries**:
    - **libXC**: The `LIBPAW_HAVE_LIBXC` macro controls conditional compilation of code that uses libXC.
    - **NetCDF**: The `LIBPAW_HAVE_NETCDF` macro controls conditional compilation for NetCDF usage, and `NCF_CHECK` provides a common error checking mechanism.
- **Fortran Preprocessor**: This file is designed to be processed by a C preprocessor before the Fortran compiler sees the code. This is a common way to achieve conditional compilation and macro-like functionality in Fortran.

In summary, `libpaw.h` is a central configuration file that uses C preprocessor directives to tailor `libPAW` to its execution environment, ensuring compatibility and leveraging host-specific features or providing default implementations.
