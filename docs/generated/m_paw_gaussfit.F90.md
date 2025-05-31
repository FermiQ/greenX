# `m_paw_gaussfit.F90`

## Overview

The Fortran module `m_paw_gaussfit` provides routines to fit radial functions, specifically PAW non-local projectors, to a sum of complex Gaussian functions. This is often done to simplify calculations involving these projectors, for example, by enabling analytical or more efficient evaluation of matrix elements. The module uses the Levenberg-Marquardt algorithm for the non-linear least-squares fitting procedure. It supports several predefined functional forms for the sum of Gaussians and includes MPI parallelization to distribute the workload when trying different numbers of Gaussian terms in the fit.

## Key Components

### Public Subroutines

-   **`gaussfit_projector(basis_size, mparam, nparam_array, nterm_bounds, orbitals, param, pawrad, rpaw, tproj, comm_mpi)`**:
    -   **Purpose**: Main public routine to fit each PAW non-local projector `tproj` to a sum of Gaussians.
    -   **Inputs**:
        -   `basis_size`: Integer, number of projectors to fit.
        -   `mparam`: Integer, maximum number of parameters for the Gaussian sum.
        -   `nterm_bounds(2)`: Integer array, min and max number of Gaussian terms to try.
        -   `orbitals(basis_size)`: Integer array, `l` quantum number for each projector.
        -   `pawrad`: `type(pawrad_type)`, input radial grid for `tproj`.
        -   `rpaw`: `real(dp)`, PAW radius.
        -   `tproj(:,:)`: `real(dp)`, array (`pawrad%mesh_size`, `basis_size`) of projectors to be fitted.
        -   `comm_mpi` (optional): Integer, MPI communicator for parallel execution.
    -   **Outputs**:
        -   `nparam_array(basis_size)`: Integer array, number of parameters found for the best fit for each projector.
        -   `param(mparam, basis_size)`: `real(dp)`, fitted Gaussian parameters for each projector.
    -   **Logic**:
        1.  For each projector `ibasis` from 1 to `basis_size`:
            a.  Prepares the projector data: extracts \( \tilde{p}_l(r) / r^l / r \) (divides by \(r^l\) due to spherical harmonic part and by an additional \(r\) due to a convention in some PAW formulations).
            b.  Splines this modified projector onto a new, finer, and potentially extended linear radial grid (`mesh_tmp`) to ensure Gaussians decay to zero.
            c.  Calls `gaussfit_main` to perform the fitting for the current projector on `mesh_tmp`.
        2.  Manages MPI communication if `comm_mpi` is provided.

-   **`gaussfit_main(mparam, nparam_out, nterm_bounds, nr, param_out, pawrad, option, outfile, rpaw, y, comm_mpi)`**:
    -   **Purpose**: Core routine to fit a single function `y(nr)` to a sum of Gaussians, trying different numbers of terms.
    -   **Inputs**:
        -   `mparam`: Max number of parameters.
        -   `nterm_bounds(2)`: Min/max number of Gaussian terms.
        -   `nr`: Number of points in the radial grid `pawrad%rad`.
        -   `pawrad`: `type(pawrad_type)` for the input function `y`.
        -   `option`: Integer, selects the functional form of the Gaussian sum (see below).
        -   `outfile`: Character string, base name for output files (fit details).
        -   `rpaw`: `real(dp)`, PAW radius (used for initial guess generation).
        -   `y(nr)`: `real(dp)`, the function to be fitted.
        -   `comm_mpi` (optional): MPI communicator.
    -   **Outputs**:
        -   `nparam_out`: Integer, number of parameters in the best fit.
        -   `param_out(mparam)`: `real(dp)`, parameters of the best fit.
    -   **Logic**:
        1.  Handles MPI distribution: If parallel, `gaussfit_mpi_main` is called by the master process to distribute the fitting tasks (different `nterm` values) among processes.
        2.  Each process (or the master in serial) iterates `nterm` from `nterm_bounds(1)` to `nterm_bounds(2)` for its assigned tasks:
            a.  Sets initial guess for parameters using `gaussfit_set_param*` routines based on `option`.
            b.  Initializes constraints for parameters using `gaussfit_constrains_init`.
            c.  Calls `gaussfit_fit` to perform the Levenberg-Marquardt optimization.
            d.  Stores the resulting \(\chi^2\) (`chisq`).
        3.  If parallel, results (\(\chi^2\) values) are gathered to the master.
        4.  The master process identifies the `nterm` (`minterm`) that gave the lowest \(\chi^2\).
        5.  The master process re-runs the fit for this `minterm` with more iterations (`maxiter=1000`) to refine the solution.
        6.  The best parameters and `nparam_out` are broadcast to all processes if parallel.

### Core Fitting and MPI Helper Subroutines (Private)

-   **`gaussfit_fit(...)`**: Implements the Levenberg-Marquardt iterations. Calls `gaussfit_rlsf`. If `verbosity >= 1`, writes the original and fitted function to `outfile`.
-   **`gaussfit_rlsf(...)`**: The Levenberg-Marquardt algorithm itself. Iteratively adjusts parameters to minimize \(\chi^2\). Uses `gaussfit_chisq_alpha_beta` to calculate \(\chi^2\), Jacobian (\(\beta\)), and Hessian approximation (\(\alpha\)). Calls `dgetrf` and `dgetri` (LAPACK) for matrix inversion.
-   **`gaussfit_chisq_alpha_beta(...)`**: Calculates \(\chi^2\), \(\alpha_{kl} = \sum_i \frac{\partial y_i}{\partial a_k} \frac{\partial y_i}{\partial a_l}\), and \(\beta_k = \sum_i (y_i - y_{fit}(a)_i) \frac{\partial y_i}{\partial a_k}\).
-   **`gaussfit_calc_deriv_r`, `_c`, `_c2`, `_c3`, `_c4`**: Subroutines to calculate \(y_{fit}(x)\) and its derivatives \(\partial y_{fit} / \partial a_k\) for the different `option` values (functional forms):
    -   `option=1`: \( \sum_j a_{1j} \cos(a_{2j} x^2) + a_{3j} \sin(a_{4j} x^2) \)
    -   `option=2`: \( \sum_j a_{1j} e^{-a_{2j} x^2} (a_{3j} \cos(a_{4j} x^2) + a_{5j} \sin(a_{6j} x^2)) \)
    -   `option=3`: \( \sum_j a_{1j} \cos(k_j x^2) + a_{2j} \sin(k_j x^2) \) (where \(k_j\) are fixed, e.g., `sep**j`)
    -   `option=4`: \( \sum_j a_{1j} e^{-a_{2j} x^2} (a_{3j} \cos(k_j x^2) + a_{4j} \sin(k_j x^2)) \) (where \(k_j\) are fixed)
-   **`gaussfit_set_param1` to `gaussfit_set_param5`**: Provide initial guesses for the parameters `a_k` for the different options.
-   **`gaussfit_constrains_init(...)`, `gaussfit_apply_constrains(...)`**: Initialize and apply constraints to parameters during fitting (e.g., positivity for widths).
-   **MPI Helper Routines (`gaussfit_mpi_*`)**:
    -   `gaussfit_mpi_main`: Orchestrates the distribution of `nterm` values across MPI processes, aiming for load balancing based on estimated computational weight (`gaussfit_mpi_set_weight`).
    -   `gaussfit_mpi_set_weight`, `_remove_item`, `_add_item`, `_calc_deviation`, `_swap`, `_assign`: Utilities for the load balancing algorithm.

### Constants for Constraints

-   `positive=2`, `restricted=3`, `restricted_and_positive=4`: Integer parameters for `gaussfit_constrains_init`.

## Important Variables/Constants

-   `option`: Integer determining the functional form of the Gaussian sum used in the fit.
-   `nterm_bounds`: Determines the range of Gaussian terms explored. The optimal number of terms is chosen based on \(\chi^2\).
-   `rpaw`: PAW radius, used for scaling initial guesses for Gaussian parameters.
-   `param_out`: Array storing the fitted parameters (amplitudes, exponents, frequencies of Gaussians).

## Usage Examples

The primary public interface is `gaussfit_projector`.

```fortran
module m_example_gaussfit_usage
    use m_paw_gaussfit
    use m_libpaw_defs, only: dp
    use m_pawrad, only: pawrad_type
    ! Assuming mpi_comm_world is available from an MPI module if used
    implicit none

    subroutine fit_my_projectors(num_proj, l_values, radial_grid_data, &
                                 paw_cutoff_radius, projector_functions, &
                                 fitted_params, num_fitted_params)
        integer, intent(in) :: num_proj
        integer, intent(in) :: l_values(num_proj)
        type(pawrad_type), intent(in) :: radial_grid_data
        real(dp), intent(in) :: paw_cutoff_radius
        real(dp), intent(in) :: projector_functions(size(radial_grid_data%rad,1), num_proj)

        integer, intent(out) :: num_fitted_params(num_proj)
        real(dp), intent(out) :: fitted_params(24, num_proj) ! Max params, e.g., 6 params/term * 4 terms_max

        integer, parameter :: mparam_val = 24 ! Max parameters, e.g. option 2, max_terms=4 => 6*4=24
        integer, dimension(2) :: nterm_bnds = [2, 4] ! Try fitting with 2, 3, or 4 Gaussian terms (option 2)

        ! Optional: MPI communicator from an MPI module
        ! integer :: mpi_comm = mpi_comm_world

        call gaussfit_projector(basis_size=num_proj, mparam=mparam_val, &
                                nparam_array=num_fitted_params, nterm_bounds=nterm_bnds, &
                                orbitals=l_values, param=fitted_params, &
                                pawrad=radial_grid_data, rpaw=paw_cutoff_radius, &
                                tproj=projector_functions)
                                ! , comm_mpi=mpi_comm) ! Uncomment for MPI

    end subroutine fit_my_projectors

end module m_example_gaussfit_usage
```

## Dependencies and Interactions

-   **`m_libpaw_defs` (`USE_DEFS`)**: For `dp` and numerical constants like `tol8`, `tol10`.
-   **`m_libpaw_tools` (`USE_MSG_HANDLING`)**: For error reporting (`LIBPAW_BUG`) and output (`wrtout`).
-   **`m_libpaw_mpi` (`USE_MPI_WRAPPERS`)**: For MPI communication (`xmpi_comm_rank`, `xmpi_comm_size`, `xmpi_bcast`, `xmpi_gatherv`).
-   **`m_libpaw_memory` (`USE_MEMORY_PROFILING`)**: For memory allocation macros.
-   **`m_paw_numeric`**: Provides `paw_splint` and `paw_spline` for spline interpolation.
-   **`m_pawrad`**: Defines `pawrad_type` and provides `pawrad_init`, `pawrad_deducer0`, `pawrad_free`, `pawrad_ifromr` for radial grid operations.
-   **LAPACK**: `dgetrf` and `dgetri` (via system LAPACK libraries) are used in `gaussfit_rlsf` for solving linear systems arising in the Levenberg-Marquardt algorithm.

This module provides a sophisticated tool for approximating radial functions with sums of Gaussians, which can be very beneficial for performance in subsequent calculations. The MPI parallelization strategy allows for efficient exploration of the optimal number of Gaussian terms.
