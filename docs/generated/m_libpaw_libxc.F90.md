# `m_libpaw_libxc.F90`

## Overview

This Fortran source file defines two modules related to the usage of the LibXC library within `libPAW`:
1.  **`m_libpaw_libxc_funcs`**: This is the core module that provides a comprehensive Fortran Application Programming Interface (API) to the LibXC library. It achieves this by defining a derived type, `libpaw_libxc_type`, to manage XC functional data, and a suite of procedures to initialize, query properties of, and evaluate these functionals. These procedures often make calls to C wrapper functions (defined in `libpaw_libxc.c`) using ISO C Binding for direct interaction with LibXC.
2.  **`m_libpaw_libxc`**: This module serves as an abstraction layer or a switch. Depending on the compilation context (specifically, if `libPAW` is compiled as part of ABINIT via `HAVE_LIBPAW_ABINIT`), it will either:
    *   `use libxc_functionals` (ABINIT's native LibXC interface module).
    *   `use m_libpaw_libxc_funcs` and rename its public entities to match the expected names (e.g., `libxc_functionals_check` becomes `libpaw_libxc_check`). This provides a consistent interface for `libPAW` regardless of whether it's using ABINIT's LibXC handling or its own.

The file heavily uses C preprocessor directives (`#include "libpaw.h"`, `#if defined LIBPAW_HAVE_LIBXC && defined LIBPAW_ISO_C_BINDING`) for conditional compilation, particularly for enabling LibXC support and ISO C Binding features.

## Key Components of `m_libpaw_libxc_funcs`

### Derived Type: `libpaw_libxc_type`

-   **Purpose**: Encapsulates all necessary information for a LibXC functional (or a pair, often exchange and correlation).
-   **Members**:
    -   `id`: Integer, LibXC numerical identifier for the functional.
    -   `family`: Integer, family of the functional (LDA, GGA, MGGA, etc.), using `LIBPAW_XC_FAMILY_*` constants.
    -   `kind`: Integer, kind of functional (exchange, correlation, etc.), using `LIBPAW_XC_*` constants.
    -   `nspin`: Integer, number of spin components (1 for unpolarized, 2 for polarized).
    -   `abi_ixc`: Integer, corresponding ABINIT `ixc` identifier.
    -   `has_exc`, `has_vxc`, `has_fxc`, `has_kxc`: Logical flags indicating availability of energy, potential, and 2nd/3rd derivatives.
    -   `needs_laplacian`: Logical, true if the functional requires the Laplacian of the density.
    -   `is_hybrid`: Logical, true if it's a hybrid functional.
    -   `hyb_mixing`, `hyb_mixing_sr`, `hyb_range`: `real(dp)`, parameters for hybrid functionals (mixing coefficient, short-range mixing, range separation omega).
    -   `temperature`: `real(dp)`, electronic temperature (if functional is temperature-dependent).
    -   `xc_tb09_c`: `real(dp)`, 'c' parameter for Tran-Blaha 09 MGGA functional.
    -   `sigma_threshold`: `real(dp)`, threshold for density gradient (sigma) for certain functionals.
    -   `conf`: `type(C_PTR), pointer`, C pointer to the underlying LibXC `xc_func_type` object.

### Public Constants (Loaded by `libpaw_libxc_constants_load`)

-   `LIBPAW_XC_FAMILY_*`: Integer constants for functional families (e.g., `LIBPAW_XC_FAMILY_LDA`, `LIBPAW_XC_FAMILY_GGA`).
-   `LIBPAW_XC_FLAGS_*`: Integer constants for functional properties/capabilities (e.g., `LIBPAW_XC_FLAGS_HAVE_EXC`, `LIBPAW_XC_FLAGS_NEEDS_LAPLACIAN`).
-   `LIBPAW_XC_*`: Integer constants for functional kinds (e.g., `LIBPAW_XC_EXCHANGE`, `LIBPAW_XC_CORRELATION`).
-   `LIBPAW_XC_SINGLE_PRECISION`: Integer, 1 if LibXC is single precision, 0 otherwise.

### Main Public Procedures

-   **`libpaw_libxc_check(stop_if_error)`**: Logical function, checks if `libPAW` was compiled with LibXC support and if LibXC is double precision. Optionally stops execution if checks fail.
-   **`libpaw_libxc_init(ixc, nspden, xc_functionals, el_temp, xc_tb09_c)`**: Subroutine to initialize one or two XC functionals based on an ABINIT-like `ixc` code. Populates the `xc_functionals` array (or the global `paw_xc_global` if `xc_functionals` is absent).
    -   `ixc`: Combined LibXC ID (e.g., `- (ID_X * 1000 + ID_C)`).
    -   `nspden`: Number of spin density components (1 or 2).
    -   `xc_functionals` (optional): Array of two `libpaw_libxc_type` to store functional info.
    -   `el_temp` (optional): Electronic temperature.
    -   `xc_tb09_c` (optional): 'c' parameter for TB09 functional.
-   **`libpaw_libxc_end(xc_functionals)`**: Subroutine to finalize and deallocate LibXC functional data stored in `xc_functionals` (or `paw_xc_global`).
-   **`libpaw_libxc_fullname(xc_functionals)`**: Character function, returns the full name(s) of the initialized functional(s) (e.g., "PBE+LDA").
-   **`libpaw_libxc_getid(xcname)`**: Integer function, returns the LibXC ID for a given functional name string.
-   **`libpaw_libxc_family_from_id(xcid)`**: Integer function, returns the family constant for a given LibXC ID.
-   **`libpaw_libxc_ixc(xc_functionals)`**: Integer function, returns the original `ixc` value used for initialization.
-   **`libpaw_libxc_islda(xc_functionals)`, `libpaw_libxc_isgga(...)`, `libpaw_libxc_ismgga(...)`**: Logical functions, check if the current functional set is LDA, GGA, or MGGA, respectively.
-   **`libpaw_libxc_is_tb09(xc_functionals)`**: Logical function, checks if any active functional is Tran-Blaha 09.
-   **`libpaw_libxc_set_c_tb09(xc_tb09_c, xc_functionals)`**: Subroutine to set the 'c' parameter for TB09.
-   **`libpaw_libxc_needs_laplacian(xc_functionals)`**: Logical function, checks if any active functional needs density Laplacian.
-   **`libpaw_libxc_needs_temperature(xc_functionals)`**: Logical function, checks if any active functional is temperature-dependent.
-   **`libpaw_libxc_set_temperature(temperature, xc_functionals)`**: Subroutine to set the electronic temperature for T-dependent functionals.
-   **`libpaw_libxc_has_kxc(xc_functionals)`, `libpaw_libxc_has_k3xc(...)`**: Logical functions, check availability of 2nd (fxc) and 3rd (kxc) derivatives from LibXC.
-   **`libpaw_libxc_nspin(xc_functionals)`**: Integer function, returns the number of spin components (1 or 2).
-   **`libpaw_libxc_is_hybrid(xc_functionals)`**: Logical function, checks if any active functional is hybrid.
-   **`libpaw_libxc_is_hybrid_from_id(xcid)`**: Logical function, checks if a functional ID corresponds to a hybrid.
-   **`libpaw_libxc_get_hybridparams(hyb_mixing, hyb_mixing_sr, hyb_range, xc_functionals)`**: Subroutine to retrieve parameters of a hybrid functional.
-   **`libpaw_libxc_set_hybridparams(hyb_mixing, hyb_mixing_sr, hyb_range, xc_functionals)`**: Subroutine to set parameters for specific hybrid functionals (PBE0, HSE).
-   **`libpaw_libxc_gga_from_hybrid(gga_id, hybrid_id, xc_functionals)`**: Logical function, attempts to find the underlying GGA functional IDs for a given hybrid functional ID.
-   **`libpaw_libxc_getvxc(...)`**: Subroutine to get XC energy density (`exc`), potential (`vxc`), and higher-order derivatives based on input density (`rho`), gradient (`grho2`), Laplacian (`lrho`), and kinetic energy density (`tau`). This is the main computational routine that calls the C wrappers for `xc_lda`, `xc_gga`, `xc_mgga`.

### Private Procedures

-   `libpaw_libxc_constants_load()`: Loads LibXC constants by calling C interface functions.
-   `char_f_to_c(f_string)` / `char_c_to_f(c_string, f_string)`: Utility functions for converting Fortran strings to C strings and vice-versa, handling null termination. Essential for C interoperability.
-   Other private helpers for specific tasks like TB09 parameter computation, temperature setting, etc.

### Global Variable

-   `paw_xc_global(2)`: A `save`d array of `libpaw_libxc_type`, acting as a default global storage for XC functional information if the user doesn't provide a local one to `libpaw_libxc_init` and other routines.

### C Interface Declarations

-   The module contains numerous `interface` blocks that declare the C functions (from `libpaw_libxc.c`) that these Fortran routines will call. These use `bind(C, name="...")` and `iso_c_binding` types.

## Key Components of `m_libpaw_libxc`

-   **Purpose**: To provide a stable interface name (`libxc_functionals_*`) for `libPAW`'s internal code, regardless of whether it's using ABINIT's LibXC layer or `libPAW`'s own (`m_libpaw_libxc_funcs`).
-   **Logic**:
    -   If `HAVE_LIBPAW_ABINIT` is defined: It `use libxc_functionals` (from ABINIT).
    -   Else: It `use m_libpaw_libxc_funcs` and renames all its public entities. For example, `libpaw_libxc_init` is made available as `libxc_functionals_init`.

## Important Variables/Constants

-   The constants defined in `m_libpaw_libxc_funcs` (e.g., `LIBPAW_XC_FAMILY_LDA`) are crucial for interpreting functional properties.
-   `libpaw_sigma_filtered(:)`: An array of character strings listing hybrid functionals that require special sigma thresholding due to (at the time of writing) limitations in LibXC v6.
-   `libpaw_sigma_threshold_def`: Default value for this sigma threshold.

## Usage Examples

```fortran
! Example of using the m_libpaw_libxc (which in turn uses m_libpaw_libxc_funcs)
module my_paw_calculation
    use m_libpaw_defs, only: dp
    use m_libpaw_libxc ! This will provide libxc_functionals_* routines
    implicit none

    subroutine setup_and_calculate_xc(ixc_abinit_style, npts_val, nspin_val, rho_data, exc_out, vxc_out)
        integer, intent(in) :: ixc_abinit_style, npts_val, nspin_val
        real(dp), intent(in) :: rho_data(npts_val, nspin_val)
        real(dp), intent(out) :: exc_out(npts_val)
        real(dp), intent(out) :: vxc_out(npts_val, nspin_val)

        ! Optional: define a local variable for XC functionals
        ! type(libpaw_libxc_type) :: local_xc_funcs(2)

        if (.not. libxc_functionals_check(.true.)) then
            ! Error message already handled by libxc_functionals_check
            return
        end if

        ! Initialize using the global paw_xc_global
        call libxc_functionals_init(ixc_abinit_style, nspin_val)
        ! Alternatively, to use a local variable:
        ! call libxc_functionals_init(ixc_abinit_style, nspin_val, xc_functionals=local_xc_funcs)

        print *, "Initialized XC functional: ", trim(libxc_functionals_fullname())
        ! print *, "Is it LDA? ", libxc_functionals_islda()

        ! For GGA/MGGA, grho2 etc. would be needed
        ! call libxc_functionals_getvxc(ndvxc, nd2vxc, npts, nspden, order, rho, exc, vxc, &
        !                               grho2, vxcgr, lrho, vxclrho, tau, vxctau, dvxc, d2vxc)
        call libxc_functionals_getvxc(0, 0, npts_val, nspin_val, 0, rho_data, exc_out, vxc_out)
                                      ! Simplified call for LDA, assuming order 0 (energy+potential)
                                      ! and no optional gradient/laplacian/tau inputs.

        ! Clean up (using global here)
        call libxc_functionals_end()
        ! Or for local:
        ! call libxc_functionals_end(xc_functionals=local_xc_funcs)

    end subroutine setup_and_calculate_xc

end module my_paw_calculation
```

## Dependencies and Interactions

-   **`libpaw.h`**: Provides preprocessor definitions like `LIBPAW_HAVE_LIBXC` and `LIBPAW_ISO_C_BINDING`.
-   **`m_libpaw_defs`**: Provides `dp` and other basic definitions.
-   **`m_libpaw_tools` (via `USE_MSG_HANDLING`)**: For error messages and output (`wrtout`).
-   **`libpaw_libxc.c`**: The C functions declared in the `interface` blocks are implemented in this C file. This forms the bridge to the actual LibXC library.
-   **LibXC Library**: The ultimate dependency for all XC calculations.
-   **ABINIT's `libxc_functionals`**: If `HAVE_LIBPAW_ABINIT` is true, this module defers to ABINIT's own LibXC interface.

This module is central to integrating advanced DFT capabilities (various XC functionals) into `libPAW`, providing a robust Fortran layer over the C-based LibXC library, and adapting to different host environments like ABINIT.
