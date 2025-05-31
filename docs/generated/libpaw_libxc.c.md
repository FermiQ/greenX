# `libpaw_libxc.c`

## Overview

This C source file, `libpaw_libxc.c`, acts as an interface or wrapper layer between the `libPAW` library and the LibXC library. LibXC is a widely used library that provides a collection of exchange-correlation (XC) functionals for density functional theory (DFT) calculations.

The main purposes of this file are:
1.  To provide a consistent C API for `libPAW` (specifically its Fortran parts) to access LibXC functionalities, abstracting away direct calls to LibXC.
2.  To handle compatibility issues across different versions of LibXC. The code contains numerous preprocessor directives (`#if (XC_MAJOR_VERSION < ... )`) to adapt function calls and constant definitions based on the LibXC version being used during compilation.
3.  To expose LibXC constants (like family types, flags, kinds of functionals) and functional evaluation routines (LDA, GGA, MGGA) to `libPAW`.
4.  To manage the memory for LibXC's `xc_func_type` objects.

This file is only compiled and used if `LIBPAW_HAVE_LIBXC` is defined (which is set in `libpaw.h`).

## Key Components

The file consists of several C functions, all prefixed with `libpaw_xc_`. These functions are designed to be callable from Fortran, likely using ISO C Binding.

### Constant Retrieval Functions:

These functions retrieve predefined constants from LibXC and pass them back to the caller (Fortran) via pointers. This is essential because these constants might have different values or might not exist across LibXC versions.

-   **`libpaw_xc_get_singleprecision_constant(int *xc_cst_singleprecision)`**:
    Determines if LibXC was compiled in single precision (`FLOAT` being `float`) or double precision (`FLOAT` being `double`). Sets `*xc_cst_singleprecision` to 1 for single, 0 for double.
-   **`libpaw_xc_get_family_constants(...)`**:
    Retrieves integer constants for different XC functional families (e.g., `XC_FAMILY_LDA`, `XC_FAMILY_GGA`, `XC_FAMILY_HYB_GGA`). Handles missing `XC_FAMILY_HYB_LDA` for LibXC versions before 6.
-   **`libpaw_xc_get_flags_constants(...)`**:
    Retrieves integer constants for XC functional flags (e.g., `XC_FLAGS_HAVE_EXC`, `XC_FLAGS_HAVE_VXC`, `XC_FLAGS_NEEDS_LAPLACIAN`). Handles `XC_FLAGS_NEEDS_LAPLACIAN` for LibXC versions before 4.
-   **`libpaw_xc_get_kind_constants(...)`**:
    Retrieves integer constants for the kind of functional (e.g., `XC_EXCHANGE`, `XC_CORRELATION`).

### Memory Management for `xc_func_type`:

-   **`XC(func_type) * libpaw_xc_func_type_malloc()`**:
    Allocates memory for an `xc_func_type` object (which holds information and parameters for a specific XC functional) and returns a pointer to it.
-   **`void libpaw_xc_func_type_free(XC(func_type) **xc_func)`**:
    Frees the memory previously allocated for an `xc_func_type` object and sets the pointer to NULL.

### LibXC Functional Wrappers:

These functions wrap the core LibXC routines for evaluating XC energies and potentials. They adapt the call signatures for different LibXC versions (e.g., LibXC v5.0 introduced more output arguments for higher-order derivatives, which are passed as NULL here if not used by `libPAW`).

-   **`libpaw_xc_get_lda(...)`**: Wrapper for `xc_lda` (LDA functional evaluation).
-   **`libpaw_xc_get_gga(...)`**: Wrapper for `xc_gga` (GGA functional evaluation).
-   **`libpaw_xc_get_mgga(...)`**: Wrapper for `xc_mgga` (meta-GGA functional evaluation).
    *Input arguments typically include `xc_func` (the functional pointer), `np` (number of points), `rho` (density), `sigma` (gradient of density for GGA/MGGA), `lapl` (Laplacian of density for MGGA), `tau` (kinetic energy density for MGGA).*
    *Output arguments include `zk` (XC energy density), `vrho` (potential derivative w.r.t. rho), `vsigma`, `vlapl`, `vtau`, and higher-order derivatives.*

### Information Retrieval from `xc_func_type`:

These functions provide access to metadata stored within an initialized `xc_func_type` object, with compatibility layers for older LibXC versions where direct struct member access was used.

-   **`libpaw_xc_get_info_name(XC(func_type) *xc_func)`**: Returns the name of the functional.
-   **`libpaw_xc_get_info_flags(XC(func_type) *xc_func)`**: Returns the flags of the functional.
-   **`libpaw_xc_get_info_kind(XC(func_type) *xc_func)`**: Returns the kind of the functional.
-   **`libpaw_xc_get_info_refs(XC(func_type) *xc_func, const int *number)`**: Returns a pointer to the character string of a specific reference for the functional.

### Parameter Setting and Retrieval for Functionals:

-   **`libpaw_xc_func_set_params(XC(func_type) *xc_func, double *ext_params, int n_ext_params)`**:
    Wrapper to set external parameters for some functionals (e.g., mixing parameter in PBE0). This function has extensive version-specific code due to changes in how LibXC handles external parameters and specific functionals like `XC_HYB_GGA_XC_PBEH` or `XC_MGGA_X_TB09`.
-   **`libpaw_xc_func_get_params_name(XC(func_type) *xc_func, const int *number)`**:
    Gets the name of an external parameter by its index. Available from LibXC v5.0.
-   **`libpaw_xc_func_get_params_description(XC(func_type) *xc_func, int *number)`**:
    Gets the description of an external parameter by its index. Available from LibXC v5.0.
-   **`libpaw_xc_func_set_params_name(XC(func_type) *xc_func, const char *name, double *par)`**:
    Sets an external parameter by its name. Available from LibXC v5.0.

### Threshold Setting for Functionals:

-   **`libpaw_xc_func_set_density_threshold(XC(func_type) *xc_func, double *dens_threshold)`**:
    Sets the density threshold below which the density is considered zero. Available from LibXC v4.0.
-   **`libpaw_xc_func_set_sig_threshold(XC(func_type) *xc_func, double *sigma_threshold)`**:
    Sets the sigma (density gradient) threshold. Available from LibXC v5.0.

### Utility Function:

-   **`libpaw_xc_func_is_hybrid_from_id(int func_id)`**:
    Checks if a functional (given its ID) is a hybrid functional by querying its family type. Adapts for `XC_FAMILY_HYB_LDA` in LibXC v6.0+.

## Important Variables/Constants

-   **`XC_MAJOR_VERSION`, `XC_MINOR_VERSION`**: Preprocessor constants from LibXC's `xc_version.h` used extensively for conditional compilation.
-   **`FLOAT`**: A type definition that is `double` by default but changes to `float` if LibXC is compiled in single precision (relevant for LibXC < v4). The `XC(func)` macro (e.g., `XC(func_type)`) resolves to `xc_func_type` or `xc_s_func_type` depending on this.
-   Various LibXC constants like `XC_FAMILY_LDA`, `XC_FLAGS_HAVE_EXC`, `XC_EXCHANGE`, etc., are used internally.

## Usage Examples

These C functions are not typically called directly by an end-user of `libPAW`. Instead, they are called by the Fortran module `m_libpaw_libxc.F90`, which provides a Fortran API to this C layer.

A conceptual call from Fortran (via `m_libpaw_libxc.F90`) might look like this:

```fortran
! In m_libpaw_libxc.F90 or a module using it
type(c_ptr) :: xc_func_ptr
integer :: functional_id, np, info
real(dp), dimension(np) :: rho, sigma, zk, vrho, vsigma
! ...
! Initialize functional
call paw_xc_func_init(xc_func_ptr, functional_id, PAW_SPIN_UNPOLARIZED)
! ...
! Call GGA evaluation
call paw_xc_gga(xc_func_ptr, np, rho, sigma, zk, vrho, vsigma)
! which internally calls libpaw_xc_get_gga via C binding
! ...
! Free functional
call paw_xc_func_end(xc_func_ptr)
```

The C functions in `libpaw_libxc.c` handle the actual interaction with the LibXC library.

## Dependencies and Interactions

-   **`libpaw.h`**: This header is included first and provides the `LIBPAW_HAVE_LIBXC` macro which determines if this entire C file's content is active.
-   **LibXC library**: This entire file is a wrapper around LibXC. It requires LibXC headers (`xc.h`, `xc_version.h`, `xc_funcs.h`, and potentially `xc_config.h` for older versions) and linking against the LibXC library.
-   **C Standard Library**: Uses `stdlib.h` (for `malloc`, `free`, `abort`) and `stdio.h` (for `fprintf`, `stderr`).
-   **Fortran part of `libPAW`**: The C functions defined here are exposed to Fortran, primarily through the `m_libpaw_libxc.F90` module, which uses ISO C Binding to call them.

This file is a critical component for enabling `libPAW` to use the wide variety of exchange-correlation functionals provided by LibXC, while maintaining robustness across different LibXC versions. The version-specific preprocessor directives make the code somewhat complex but are necessary for this backward and forward compatibility.
