# `lapack_interfaces.f90`

## Overview

The `lapack_interfaces` module in the `GX-common` component provides explicit Fortran interface declarations for a set of commonly used BLAS (Basic Linear Algebra Subprograms) and LAPACK (Linear Algebra PACKage) subroutines. Defining these interfaces is crucial for ensuring type safety and enabling the compiler to verify the correctness of calls to these external, highly optimized numerical libraries.

The routines interfaced are primarily double-precision versions.

## Key Components

The module consists of several `interface` blocks, each declaring the signature of one or more external procedures.

### Interface: `blas3`

- **`subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)`**:
    - **Purpose**: Performs double-precision general matrix-matrix multiplication.
    - Operation: `C := alpha*op(A)*op(B) + beta*C`, where `op(X)` is `X` or `X^T`.
    - Key arguments:
        - `transa`, `transb`: Specify if matrices A and B are transposed.
        - `m`, `n`, `k`: Dimensions of the matrices.
        - `alpha`, `beta`: Scaling factors.
        - `a`, `b`, `c`: The matrices involved.
        - `lda`, `ldb`, `ldc`: Leading dimensions of the matrices.

### Interface: `svd`

- **`subroutine dgesdd(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info)`**:
    - **Purpose**: Computes the singular value decomposition (SVD) of a double-precision M-by-N matrix A: `A = U * SIGMA * V^T`. It uses a divide-and-conquer algorithm.
    - Key arguments:
        - `jobz`: Specifies whether to compute singular vectors (U, VT) and how many.
        - `m`, `n`: Dimensions of matrix A.
        - `a`: The input matrix A, overwritten on exit.
        - `s`: Output array of singular values.
        - `u`, `vt`: Output arrays for left and right singular vectors.
        - `work`, `lwork`, `iwork`: Workspace arrays.
        - `info`: Output error status.

### Interface: `presicion` (Note: "precision" is misspelled in the Fortran code)

- **`double precision function dlamch(cmach)`**:
    - **Purpose**: Determines double-precision floating-point machine parameters.
    - Key arguments:
        - `cmach`: Character input specifying which parameter to return (e.g., 'E' for epsilon, 'S' for safe minimum, 'B' for base, 'P' for precision, 'O' for overflow, 'U' for underflow).

### Interface: `diag`

- **`subroutine dsyevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, iwork, ifail, info)`**:
    - **Purpose**: Computes selected eigenvalues and, optionally, eigenvectors of a real symmetric matrix A.
    - Key arguments:
        - `jobz`: 'N' for eigenvalues only, 'V' for eigenvalues and eigenvectors.
        - `range`: 'A' for all, 'V' for values in `[vl, vu]`, 'I' for i-th through j-th.
        - `uplo`: 'U' if upper triangle of A is stored, 'L' if lower.
        - `n`: Order of matrix A.
        - `a`: Input symmetric matrix.
        - `vl`, `vu`, `il`, `iu`: Define range of eigenvalues to find.
        - `abstol`: Absolute error tolerance for eigenvalues.
        - `m`: Output, total number of eigenvalues found.
        - `w`: Output, selected eigenvalues.
        - `z`: Output, selected eigenvectors.
        - `info`: Output error status.

### Interface: `unitary`

- **`subroutine dlaset(uplo, m, n, alpha, beta, a, lda)`**:
    - **Purpose**: Initializes a double-precision M-by-N matrix A to `beta` on the off-diagonals and `alpha` on the diagonal.
    - Key arguments:
        - `uplo`: 'U' for upper triangle, 'L' for lower; specifies which part to set. If 'F' (Full), all elements are set.
        - `m`, `n`: Dimensions of matrix A.
        - `alpha`: Value for diagonal elements.
        - `beta`: Value for off-diagonal elements.
        - `a`: The matrix to be initialized.

## Usage Examples

These interfaces are not called directly by end-users of the GX-TimeFrequency library but are used internally by other modules that perform numerical computations. For example, `minimax_grids.F90` uses `dgemm` and `dgesdd`.

```fortran
! Example of how dgemm might be used in another module
module some_computation
  use kinds, only: dp
  use lapack_interfaces, only: blas3 ! Or directly 'use lapack_interfaces, only: dgemm'
  implicit none

  subroutine matrix_multiply_example(A, B, C, m, n, k)
    real(kind=dp), dimension(:,:), intent(in) :: A, B
    real(kind=dp), dimension(:,:), intent(inout) :: C
    integer, intent(in) :: m, n, k

    ! Assuming A(m,k), B(k,n), C(m,n)
    ! Perform C = 1.0 * A * B + 0.0 * C
    call dgemm('N', 'N', m, n, k, 1.0_dp, A, m, B, k, 0.0_dp, C, m)
  end subroutine matrix_multiply_example

end module some_computation
```

## Dependencies and Interactions

- **`kinds`**: All interfaces within this module `use kinds, only: dp` to ensure that real variables are declared with double precision, matching the 'd' prefix in the LAPACK/BLAS routine names (e.g., `dgemm`).
- **External Libraries (LAPACK/BLAS)**: This module is critically dependent on the availability of a LAPACK and BLAS library at link time. The object code from these libraries provides the actual implementations of the subroutines.
- Other modules in GreenX (e.g., `minimax_grids`, other numerical algorithm modules) will `use lapack_interfaces` to make these routines callable in a type-safe manner.

By providing these explicit interfaces, the module promotes robust and portable use of standard numerical linear algebra libraries, which are fundamental for high-performance scientific computing.
