# `m_paw_sphharm.F90`

## Overview

The Fortran module `m_paw_sphharm` provides a comprehensive suite of routines for computations involving complex spherical harmonics (\(Y_{lm}\)), real spherical harmonics (\(S_{lm}\)), Associated Legendre Polynomials, and related angular momentum algebra. These tools are fundamental for Projector Augmented Wave (PAW) calculations and other quantum mechanical methods that utilize spherical harmonic basis sets or expansions.

The module includes functionalities to:
-   Evaluate spherical harmonics and their derivatives.
-   Compute matrix elements of angular momentum operators (\(\vec{L}\), \(\vec{L} \cdot \vec{S}\)).
-   Transform representations between different spherical harmonic bases (complex vs. real, \(m_l m_s\) vs. \(j m_j\)).
-   Generate angular meshes and associated weights for numerical integration.
-   Calculate rotation matrices for spherical harmonics under symmetry operations.
-   Compute Gaunt coefficients (integrals of products of three spherical harmonics) and similar integrals involving gradients.

## Key Components

### Public Procedures

-   **`ylmc(il, im, kcart)`**: Complex function, returns the complex spherical harmonic \(Y_{lm}(\theta, \phi)\) for \(l \le 3\), given Cartesian coordinates `kcart`.
-   **`ylmcd(il, im, kcart, dth, dphi)`**: Subroutine, computes first derivatives of complex \(Y_{lm}\) (\(l \le 3\)) with respect to \(\theta\) (`dth`) and \(\phi\) (`dphi`).
-   **`ylm_cmplx(lx, ylm, xx, yy, zz)`**: Subroutine, calculates all complex spherical harmonics \(Y_{lm}\) up to `lx` (\(\le 4\)), storing them in `ylm`.
-   **`initylmr(mpsang, normchoice, npts, nrm, option, rr, ylmr, ylmr_gr)`**: Subroutine, calculates real spherical harmonics \(S_{lm}\) and optionally their Cartesian gradients and second derivatives on a set of vectors `rr`.
    -   `mpsang`: \(1 + l_{max}\).
    -   `option`: Controls calculation of values (1), values+gradients (2), or values+gradients+second derivatives (3).
-   **`ys(l2, m2, l1, m1, ys_val)`**: Subroutine, computes the matrix element \(\langle Y_{l_2 m_2} | S_{l_1 m_1} \rangle\).
-   **`lxyz(lp, mp, idir, ll, mm, lidir)`**: Subroutine, computes matrix element \(\langle Y_{l'm'} | L_{idir} | Y_{lm} \rangle\) where \(L_{idir}\) is a component of the angular momentum operator.
-   **`slxyzs(lp, mp, idir, ll, mm, sls_val)`**: Subroutine, computes matrix element \(\langle S_{l'm'} | L_{idir} | S_{lm} \rangle\).
-   **`lsylm(ls_ylm, lmax)`**: Subroutine, computes matrix elements of the \(\vec{L} \cdot \vec{S}\) operator in the real spherical harmonics basis coupled with spin: \(\langle \sigma, S_{lm_1} | \vec{L} \cdot \vec{S} | S_{lm_2}, \sigma' \rangle\).
-   **`plm_coeff(blm, mpsang, xx)`**: Subroutine, computes coefficients related to Associated Legendre Polynomials \(P_{lm}(x)\) and their derivatives, used for second derivatives of \(Y_{lm}\).
-   **`ass_leg_pol(l, m, xarg)`**: Real function, computes the Associated Legendre Polynomial \(P_l^m(x)\).
-   **`plm_dphi(ll, mm, xx)`**: Real function, computes \(m P_{lm}(x) / \sqrt{1-x^2}\).
-   **`plm_dtheta(ll, mm, xx)`**: Real function, computes \(-\sqrt{1-x^2} \frac{d}{dx} P_{lm}(x)\).
-   **`plm_d2theta(mpsang, plm_d2t, xx)`**: Subroutine, computes \(\frac{d^2}{d\theta^2} P_{lm}(\cos\theta)\).
-   **`pl_deriv(mpsang, pl_d2, xx)`**: Subroutine, computes \(\frac{d^2}{dx^2} P_l(x)\) (Legendre polynomial second derivative).
-   **`ylm_angular_mesh(ntheta, nphi, angl_size, cart_coord, ang_wgth)`**: Subroutine, generates points (`cart_coord`) and weights (`ang_wgth`) for a Gauss-Legendre angular mesh on a unit sphere.
-   **`mat_mlms2jmj(lcor, mat_mlms, mat_jmj, ndij, option, optspin, prtvol, unitfi, wrt_mode)`**: Subroutine, transforms a matrix from the \(|l,s,m_l,m_s\rangle\) basis (`mat_mlms`) to the \(|l,s,j,m_j\rangle\) basis (`mat_jmj`) or vice-versa.
    -   `option=1`: \(m_l m_s \to j m_j\).
    -   `option=2`: \(j m_j \to m_l m_s\).
-   **`mat_slm2ylm(lcor, mat_inp_c, mat_out_c, ndij, option, optspin, prtvol, unitfi, wrt_mode)`**: Subroutine, transforms a matrix between real spherical harmonic (\(S_{lm}\)) and complex spherical harmonic (\(Y_{lm}\)) bases.
    -   `option=1`: \(S_{lm} \to Y_{lm}\).
    -   `option=2`: \(Y_{lm} \to S_{lm}\).
-   **`setsym_ylm(gprimd, lmax, nsym, pawprtvol, rprimd, sym, zarot)`**: Subroutine, computes rotation matrices (`zarot`) for real spherical harmonics under crystal symmetry operations `sym`.
-   **`setnabla_ylm(ang_phipphj, mpsang)`**: Subroutine, computes and tabulates various angular integrals involving spherical harmonics and components of the gradient operator. These are used in `m_paw_onsite`. The integrals are of the form \(\int S_i S_j \sin\theta \cos\phi d\Omega\), \(\int S_i \frac{dS_j}{d\theta} \cos\theta \cos\phi d\Omega\), etc.
-   **`gaunt(ll, mm, l1, m1, l2, m2)`**: Real function, computes Gaunt coefficients \(\int Y_{l_1 m_1}^* Y_{l_2 m_2} Y_{ll mm} d\Omega\). (Note: The docstring says \(\sqrt{4\pi} \int Y_{l_i m_i}^* Y_{ll mm}^* Y_{l_j m_j} d\Omega \), which is more standard).
-   **`realgaunt(l_max, ngnt, gntselect, realgnt)`**: Subroutine, computes "real Gaunt coefficients" \(\int S_{lm} S_{l'm'} S_{l''m''} d\Omega\). Stores non-zero coefficients and their indices.
-   **`nablarealgaunt(l_max, l_max_ij, nnablagnt, nabgauntselect, nablagaunt)`**: Subroutine, computes integrals of the form \(\int S_{lm} (\nabla S_{l'm'}) \cdot (\nabla S_{l''m''}) d\Omega\). Stores non-zero values and their indices. Some values for low \(l, l_{ij}\) are hardcoded; higher values are computed numerically via `initylmr` and angular quadrature.

### Private Helper Functions

-   `create_slm2ylm`, `create_mlms2jmj`: Generate transformation matrices.
-   `mkeuler`: Determines Euler angles from a rotation matrix.
-   `dbeta`: Calculates Wigner d-matrix elements \(d^l_{m'm}(\beta)\).
-   `phim`: Computes basis functions for azimuthal angle dependence in real spherical harmonics.
-   `gauleg`: Computes Gauss-Legendre quadrature points and weights.
-   `perms`: Computes \(N!/(N-k)!\).
-   `rfactorial`: Computes \(N!\) as a real number.

## Important Variables/Constants

-   The module relies on constants from `m_libpaw_defs` (e.g., `dp`, `pi`, `sqrt2`, `tol*`).
-   Mathematical precision and algorithm choices (e.g., series expansions vs. recurrence relations for Bessel functions in `ylmc`) are embedded in the routines.

## Usage Examples

These routines are fundamental tools used by various other modules within `libPAW` that deal with angular dependencies, spherical expansions, or symmetry operations.

```fortran
module m_example_sphharm_usage
    use m_paw_sphharm
    use m_libpaw_defs, only: dp, pi, czero, cone
    implicit none

    subroutine test_sphharm_calcs
        complex(dpc) :: y21_val, y21_dth, y21_dphi
        real(dp) :: kvec(3) = [1.0_dp, 1.0_dp, 1.0_dp]
        real(dp) :: real_ylm_grid( (2+1)**2, 100 ) ! For lmax=2, 100 grid points
        real(dp) :: grid_coords(3,100), grid_norms(100)
        integer :: npts_val = 100, lmax_val = 2

        ! 1. Calculate a complex spherical harmonic Y_2,1
        y21_val = ylmc(2, 1, kvec)
        ! print *, "Y_2,1 at (1,1,1) = ", y21_val

        ! 2. Calculate its derivatives
        call ylmcd(2, 1, kvec, y21_dth, y21_dphi)
        ! print *, "dY_2,1/dtheta = ", y21_dth, ", dY_2,1/dphi = ", y21_dphi

        ! 3. Initialize real spherical harmonics on a dummy grid
        ! (grid_coords and grid_norms would need to be properly set up)
        ! call initylmr(lmax_val+1, 1, npts_val, grid_norms, 1, grid_coords, real_ylm_grid)

        ! 4. Compute Gaunt coefficient
        ! real(dp) :: gaunt_val = gaunt(l1=1,m1=0, l2=1,m2=0, ll=2,mm=0) ! <Y10 Y10 | Y20>
        ! print *, "Gaunt <1010|20> = ", gaunt_val / sqrt(4.0_dp*pi) ! To match common definition

    end subroutine test_sphharm_calcs

end module m_example_sphharm_usage
```

## Dependencies and Interactions

-   **`m_libpaw_defs` (`USE_DEFS`)**: For `dp` and mathematical constants.
-   **`m_libpaw_tools` (`USE_MSG_HANDLING`)**: For error reporting (`LIBPAW_ERROR`, `LIBPAW_BUG`).
-   **`m_libpaw_memory` (`USE_MEMORY_PROFILING`)**: For memory allocation macros.
-   Many other `libPAW` modules that require manipulation of angular quantities or PAW functions in spherical coordinates will depend on `m_paw_sphharm`. For example, `m_paw_finegrid` uses `initylmr`, and `m_paw_onsite` uses `setnabla_ylm`.

This module encapsulates a significant amount of the mathematical machinery required for spherical harmonic manipulations in an atomic or PAW context. Its correctness and efficiency are vital for the overall performance and accuracy of `libPAW`. The hardcoded values in `setnabla_ylm` and `nablarealgaunt` for low \(l\) values suggest these are pre-calculated for efficiency, with numerical integration used for higher \(l\) if needed.
