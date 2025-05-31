# `m_pawpsp.F90`

## Overview

The Fortran module `m_pawpsp` is a critical component of `libPAW` for handling Projector Augmented Wave (PAW) pseudopotential (PSP) datasets. Its primary function is to read PAW atomic data from files, process this data, and compute various derived quantities necessary for PAW calculations. The module supports both traditional formatted PAW PSP files (e.g., "pspcod=7") and XML-based PAW PSP files (e.g., "pspcod=17", often from ATOMPAW).

Key operations include:
-   Reading header information and raw data (radial functions for partial waves, projectors, core densities, local potentials) from PSP files.
-   Initializing PAW data structures (`pawtab_type`, `pawrad_type`).
-   Computing derived quantities like projector form factors (\(f_l(q)\)) and the q-space local potential (\(V_{loc}(q)\)).
-   Handling specific data transformations if wavelets are used (e.g., fitting projectors to Gaussians).
-   Performing consistency checks on the read data.

## Key Components

### Derived Type: `pawpsp_header_type`

-   **Purpose**: Stores common header information read from a PAW PSP file.
-   **Members**:
    -   `basis_size`: Integer, number of \((l,n)\) PAW basis functions (partial waves/projectors).
    -   `l_size`: Integer, \(l_{max\_gaunt}+1\).
    -   `lmn_size`: Integer, total number of \((l,m,n)\) PAW basis channels.
    -   `mesh_size`: Integer, dimension of the main radial mesh used in the PSP.
    -   `pawver`: Integer, version number of the PAW PSP format.
    -   `shape_type`: Integer, type of shape function used for compensation charge.
    -   `rpaw`: `real(dp)`, PAW radius.
    -   `rshp`: `real(dp)`, cut-off radius of the shape function.

### Main Public Subroutines

-   **`pawpsp_main(...)`**:
    -   **Purpose**: Orchestrates the reading and processing of a PAW PSP file. This is the main entry point for this module.
    -   **Inputs**: Filename (`filpsp`), flags for wavelet usage (`usewvl`), Coulomb treatment (`icoulomb`), XC functional info (`ixc`, `xclevel`, `pawxcdev`), compensation charge handling (`usexcnhat_in`), q-space grids (`qgrid_ff`, `qgrid_vl`), and optionally pre-read XML data (`psxml`), MPI communicator (`comm_mpi`), and XC tolerances.
    -   **Outputs**: Populated `pawrad` and `pawtab` structures, \(f_l(q)\) (`ffspl`), \(V_{loc}(q)\) (`vlspl`), integrated local potential term (`epsatm`), and XC core correction radius (`xcccrc`).
    -   **Logic**:
        1.  Checks if the PSP file is XML or traditional format.
        2.  If not XML, opens the file and calls `pawpsp_read_header` to get basic info.
        3.  If XML and `psxml` is provided, uses `pawpsp_read_header_xml` and `pawpsp_read_pawheader`.
        4.  Performs consistency checks on header data vs. input parameters (`pawpsp_consistency`).
        5.  Calls format-specific routines:
            -   `pawpsp_7in` for traditional "pspcod=7" files.
            -   `pawpsp_17in` for XML "pspcod=17" files (which uses `paw_setuploc` from XML parsing).
        6.  If MPI is used, broadcasts the processed data (`pawpsp_bcast`).
        7.  If wavelets are enabled (`usewvl=1` or `icoulomb > 0`), calls `pawpsp_wvl` to fit projectors to Gaussians.

-   **`pawpsp_7in(...)`**:
    -   **Purpose**: Handles reading of traditional "pspcod=7" formatted PAW files and subsequent calculations.
    -   **Logic**: Calls `pawpsp_read` to parse the file, then `pawpsp_calc` to compute derived quantities. If wavelets are used, calls `pawpsp_calc_d5` and `pawpsp_wvl_calc`.
-   **`pawpsp_17in(...)`**:
    -   **Purpose**: Handles reading of XML-based "pspcod=17" PAW files (using data pre-parsed into `paw_setup_t` by `m_pawxmlps`) and subsequent calculations.
    -   **Logic**: Similar to `pawpsp_7in`, it initializes radial meshes and PAW data from `paw_setuploc` (derived from `psxml`), then calls `pawpsp_calc`. If wavelets are used, calls `pawpsp_calc_d5` and `pawpsp_wvl_calc`.

### Core Data Reading and Processing Subroutines

-   **`pawpsp_read_header(funit, lloc, lmax, mmax, pspcod, pspxc, r2well, zion, znucl)`**: Reads the first 3 lines of a traditional PAW PSP file.
-   **`pawpsp_read_header_2(funit, pspversion, basis_size, lmn_size)`**: Reads specific header parts (version, basis sizes) for traditional files.
-   **`pawpsp_read_header_xml(...)`**: Extracts header-like information from a pre-parsed XML structure (`psxml` of `paw_setup_t`).
-   **`pawpsp_read_pawheader(...)`**: Similar to `pawpsp_read_header_xml`, populates `pawpsp_header_type` from `psxml`.
-   **`pawpsp_read(core_mesh, ...)`**: The main routine for parsing the content of a traditional (pspcod=7) PAW PSP file after the initial header lines. Reads radial meshes, partial waves, projectors, core densities, local potential, \(D_{ij}^0\), \(\rho_{ij}^0\), and shape functions.
-   **`pawpsp_read_corewf(...)`**: Reads core wavefunctions from a separate file (e.g., `corewf.abinit`, `corewf.xml`, `corewf.dat`). Handles different formats and potential relativistic effects.
-   **`pawpsp_calc(...)`**: Central calculation routine. After raw data is read (by `pawpsp_read` or from `psxml` via `pawpsp_17in`), this routine:
    -   Performs consistency checks.
    -   Computes \(D_{ij}^0\) if not directly provided (or if `vlocopt=0` requires recalculation with \(V_H(\tilde{n}_{Zc})\)).
    -   Computes kinetic energy components \(K_{ij}\) if needed.
    -   Calculates \(f_l(q)\) (`ffspl`) using `pawpsp_nl`.
    -   Calculates \(V_{loc}(q)\) (`vlspl`) and \(\epsilon_{atm}\) using `pawpsp_lo`.
    -   Calculates \(n_{core}(q)\) (`pawtab%tcorespl`) and \(\tau_{core}(q)\) (`pawtab%tcoretauspl`) using `pawpsp_cg`.
    -   Calculates XC energy of the core density (`pawtab%exccore`).
    -   May spline radial functions onto finer/different grids if needed for consistency or numerical accuracy.
-   **`pawpsp_nl(ffspl, ...)`**: Computes projector form factors \(f_l(q) = \int j_l(2\pi qr) u_l(r) r dr\).
-   **`pawpsp_lo(epsatm, ..., q2vq, ...)`**: Computes \(q^2 V_{loc}(q)\) and \(\epsilon_{atm} = 4\pi \int [r^2 V_{loc}(r) + rZ_v] dr\).
-   **`pawpsp_cg(dnqdq0, d2nqdq0, ..., nq, ...)`**: Computes \(n(q) = 4\pi \int \frac{\sin(2\pi qr)}{2\pi qr} r^2 n(r) dr\).
-   **`pawpsp_calc_d5(mesh, mesh_size, tcoredens)`**: Calculates up to 5th derivatives of `tcoredens(:,1)` using splines and stores them in `tcoredens(:,2:6)`.
-   **`pawpsp_vhar2rho(radmesh, rho, vv)`**: Solves Poisson's equation \(\nabla^2 V = -4\pi\rho\) to get \(\rho\) from \(V\). (Note: the formula in comments is for \(\nabla^2 V = 4\pi\rho\), usually it's \(-\rho/\epsilon_0\) or \(-4\pi\rho\) in atomic units).
-   **`pawpsp_wvl_calc(...)`**: Performs wavelet-specific calculations, primarily setting up `pawtab%wvl%rholoc` (local density and potential for wavelets) by calling `pawpsp_vhar2rho`.
-   **`pawpsp_wvl_sin2gauss(...)`**: Converts a function represented as \(\sum (a \sin(kx^2) + b \cos(kx^2))\) or similar into a sum of complex Gaussians. Used for wavelet basis.
-   **`pawpsp_rw_atompaw(...)`**: Reads an ATOMPAW PSP file and writes out a new one with Gaussian projector information if `pawtab%wvl` is populated.
-   **`pawpsp_bcast(...)`**: Broadcasts processed PAW data (`pawrad`, `pawtab`, `ffspl`, `vlspl`, etc.) to all MPI processes.

## Important Variables/Constants

-   `pawtab_type`: Defined in `m_pawtab`. This is the primary data structure populated by this module, holding all PAW setup information for an atom type.
-   `pawrad_type`: Defined in `m_pawrad`. Describes radial grids.
-   `pspversion`: Integer indicating the format version of the PAW PSP file.
-   `pspcod`: Integer, code for pseudopotential type (7 for PAW traditional, 17 for PAW XML).
-   `usexcnhat`: Flag controlling whether the compensation charge density \(\hat{n}\) is included in the XC functional evaluation within PAW spheres.
-   `vlocopt`: Option for how the local potential \(V_{loc}\) is defined or read from the PSP file (0: \(V_{bare}\), 1: \(V_H(\tilde{n}_{Zc})\) with \(\hat{n}\) in XC, 2: \(V_H(\tilde{n}_{Zc})\) without \(\hat{n}\) in XC).

## Usage Examples

The main entry point `pawpsp_main` is typically called by a code like ABINIT or Quantum Espresso when it needs to load PAW dataset information for a particular atomic species.

```fortran
! Simplified example of how pawpsp_main might be invoked
module m_dft_simulation
    use m_pawpsp
    use m_libpaw_defs, only: dp
    ! ... other necessary modules ...
    implicit none

    subroutine load_paw_dataset(filename, atom_z, valence_z, xc_choice, &
                                lmax_psp, paw_setup_data, radial_grid_data, &
                                form_factors, vloc_qspace, eps_atm_val, xc_core_rad)
        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: atom_z, valence_z
        integer, intent(in) :: xc_choice, lmax_psp

        type(pawrad_type), intent(inout) :: radial_grid_data
        type(pawtab_type), intent(inout) :: paw_setup_data
        real(dp), intent(out) :: form_factors(:,:,:) ! (mqgrid_ff, 2, lnmax)
        real(dp), intent(out) :: vloc_qspace(:,:)   ! (mqgrid_vl, 2)
        real(dp), intent(out) :: eps_atm_val, xc_core_rad

        ! Dummy q-grids for example
        integer, parameter :: mqgrid_ff_val = 100, mqgrid_vl_val = 100
        real(dp), dimension(mqgrid_ff_val) :: q_grid_ff_vals = [(0.01_dp * real(i), i=1,mqgrid_ff_val)]
        real(dp), dimension(mqgrid_vl_val) :: q_grid_vl_vals = [(0.01_dp * real(i), i=1,mqgrid_vl_val)]

        ! Other parameters would be set based on input/defaults
        integer :: use_wavelets = 0, coulomb_treat = 0, paw_xc_dev = 1, use_comp_charge_xc = 1
        integer :: xc_level_val = 1 ! LDA
        real(dp) :: hybrid_mix = 0.0_dp

        call pawpsp_main(pawrad=radial_grid_data, pawtab=paw_setup_data, &
                         filpsp=filename, usewvl=use_wavelets, icoulomb=coulomb_treat, &
                         hyb_mixing=hybrid_mix, ixc=xc_choice, xclevel=xc_level_val, &
                         pawxcdev=paw_xc_dev, usexcnhat_in=use_comp_charge_xc, &
                         qgrid_ff=q_grid_ff_vals, qgrid_vl=q_grid_vl_vals, &
                         ffspl=form_factors, vlspl=vloc_qspace, &
                         epsatm=eps_atm_val, xcccrc=xc_core_rad, &
                         zionpsp=valence_z, znuclpsp=atom_z)
                         ! Optional args like psxml, comm_mpi, xc_denpos, etc. omitted for brevity

    end subroutine load_paw_dataset

end module m_dft_simulation
```

## Dependencies and Interactions

-   **`m_libpaw_defs`, `m_libpaw_tools`, `m_libpaw_mpi`, `m_libpaw_memory`**: For basic definitions, utilities, MPI, and memory management.
-   **`m_libpaw_libxc`**: Used by `pawpsp_calc` to get `libxc_functionals_nspin` if `ixc < 0`.
-   **`m_pawxmlps`**: If an XML PAW dataset is being read, `rdpawpsxml_core` (for core wavefunctions) and other routines from `m_pawxmlps` (like `pawpsxml2ab`, though `pawpsp_17in` seems to directly use the `paw_setup_t` type) are used to parse the XML file into a `paw_setup_t` structure before `pawpsp_17in` processes it.
-   **`m_pawang`, `m_pawtab`, `m_pawrad`**: Define the core data structures (`pawang_type`, `pawtab_type`, `pawrad_type`) that are populated or used by this module.
-   **`m_paw_numeric`**: Provides numerical tools like spline interpolation (`paw_spline`, `paw_splint`), Bessel functions (`paw_jbessel_4spline`), etc., used in `pawpsp_nl`, `pawpsp_calc`, `pawpsp_calc_d5`.
-   **`m_paw_atom`**: Provides atomic PAW calculations like `atompaw_shapebes`, `atompaw_vhnzc`, `atompaw_dij0`, `atompaw_kij`, used in `pawpsp_calc`.
-   **`m_pawxc`**: Provides `pawxc` and `pawxcm` for calculating on-site XC energy contributions, used in `pawpsp_calc`.
-   **`m_paw_gaussfit`**: Provides `gaussfit_projector` used by `pawpsp_wvl` if wavelets are enabled.
-   **`fox_sax`**: (Fortran XML library) Used if `LIBPAW_HAVE_FOX` is defined, likely for parsing XML PSP files if `m_pawxmlps` relies on it.

This module acts as the primary interface for ingesting PAW atomic datasets, making it a foundational step for any subsequent PAW calculation. It bridges the file format data with the internal `pawtab_type` representation and computes several key derived quantities.
