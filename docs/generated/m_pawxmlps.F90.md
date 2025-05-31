# `m_pawxmlps.F90`

## Overview

The Fortran module `m_pawxmlps` is designed to read Projector Augmented Wave (PAW) pseudopotential (PSP) data from files formatted in XML. This is in contrast to older, traditional formatted PAW files. The module can interface with the FoX (Fortran XML) library for parsing if `libPAW` is compiled with FoX support. Alternatively, it contains pure Fortran routines to parse specific information from the XML files, particularly for header data and core wavefunctions, by manually processing lines from the XML file.

The module defines several derived data types (like `radial_grid_t`, `radialfunc_t`, `state_t`, and the main container `paw_setup_t`) that mirror the hierarchical structure of the PAW XML data. These types are used to temporarily store the information as it is read from the file before being processed and transferred into the main `libPAW` data structures (like `pawtab_type`).

## Key Components

### Derived Data Types (Mirroring XML Structure)

These types are primarily used internally by the FoX SAX parser or the manual Fortran parser to hold data as it's extracted from the XML file.

-   **`radial_grid_t`**: Stores attributes of a `<radial_grid>` tag (equation type `eq`, parameters `aa`, `bb`, `dd`, `nn`, start/end indices `istart`, `iend`, and grid ID `id`).
-   **`radialfunc_t`**: Stores a radial function, including the ID of the grid it's defined on (`grid`), an optional state identifier (`state`), and the actual `data(:)`.
-   **`gaussian_expansion_t`**: Stores parameters for a Gaussian expansion of a function, including `ngauss`, `state` ID, and arrays for `factors` and `expos`(nents) of the Gaussians.
-   **`shape_function_t`**: Describes a PAW shape function, including its type (`gtype`), cutoff radius (`rc`), grid ID, `lamb` parameter (for Gaussian types), and numerical `data(:,:)` if applicable.
-   **`state_t`**: Describes an atomic state (valence or core) with quantum numbers `nn` (principal), `ll` (angular), `kk` (kappa, for relativistic), occupation `ff`, cutoff radius `rc`, energy `ee`, and state ID `id`.
-   **`valence_states_t`**: Container for an array of `state_t` describing valence states.
-   **`generator_t`**: Information about the code/method used to generate the PAW dataset.
-   **`xc_functional_t`**: Describes the XC functional used (type and name).
-   **`atom_t`**: Basic atomic information (symbol, Z, core charge, valence charge).
-   **`paw_setup_t`**: The top-level container that aggregates all the above information read from a PAW XML file. It includes members for version, atomic data, XC functional, generator info, valence states, multiple radial grids, shape function, various radial functions (core densities, partial waves, projectors, local potentials, kinetic energy densities, etc.), and Gaussian fits for projectors.
    -   Public procedures associated: `paw_setup_free`, `paw_setup_copy`.

### Public Procedures

-   **`rdpawpsxml(filename, paw_setup)`**:
    -   **Purpose**: Main routine to read a complete PAW dataset from an XML file into the `paw_setup` structure. This version is primarily intended for use with the FoX XML library.
    -   **Logic (Conceptual, relies on FoX SAX events)**: Initializes `paw_setuploc` (a module-level variable of `paw_setup_t`), then calls `processFile` (a FoX routine) which uses `paw_begin_element1`, `paw_end_element1`, and `pawdata_chunk` as callbacks to parse the XML and populate `paw_setuploc`. Finally, copies `paw_setuploc` to the output `paw_setup`.
-   **`rdpawpsxml_header(ecut_tmp, filename, paw_setup)`**:
    -   **Purpose**: Reads only the header portion of a PAW XML file using pure Fortran string parsing (line-by-line reading and keyword searching via `paw_rdfromline`). It populates parts of the `paw_setup` structure, primarily header-like information and recommended energy cutoffs `ecut_tmp`.
    -   **Note**: This is a fallback or alternative to a full FoX parse, focusing on extracting essential initial information.
-   **`rdpawpsxml_core(energy_cor, filename, lcor, ncor, nphicor, pawrad, phi_cor, kappacor)`**:
    -   **Purpose**: Reads core wavefunctions (\(\phi_{core}\)) and their properties (energies `energy_cor`, quantum numbers `lcor`, `ncor`, `kappacor`) from a PAW XML file. Uses pure Fortran string parsing.
    -   **Inputs**: `filename`, `pawrad` (for target radial grid).
    -   **Outputs**: `energy_cor`, `lcor`, `ncor`, `nphicor` (number of core states), `phi_cor` (core wavefunctions splined onto `pawrad`), and optionally `kappacor` (kappa quantum numbers for relativistic states).
    -   **Logic**: Parses the XML file for `<core_states>` and `<ae_core_wavefunction>` sections, reads the data, and splines it onto the provided `pawrad` if the grids differ.

### FoX SAX Event Handlers (Private, but public due to FoX requirements)

These are only active if `LIBPAW_HAVE_FOX` is defined.
-   **`paw_begin_element1(namespaceURI, localName, name, attributes)`**: Callback for FoX when an XML start tag (`<name ...>`) is encountered. It inspects `name` and `attributes` to identify the data type being defined and prepares module-level variables (like `rp`, `valstate`, `grids`) to receive data.
-   **`paw_end_element1(namespaceURI, localName, name)`**: Callback for FoX when an XML end tag (`</name>`) is encountered. It finalizes the data stored in module-level variables (e.g., copies from `valstate` to `paw_setuploc%valence_states%state`).
-   **`pawdata_chunk(chunk)`**: Callback for FoX when character data (content between tags) is encountered. It parses the numerical data from `chunk` and stores it into the currently active array pointed to by `rp%data`.

### Private Helper Procedures

-   **`paw_rdfromline(keyword, line, output, ierr)`**:
    -   A utility function for the pure Fortran XML parsing. It searches for `keyword="value"` within a given `line` and extracts `value` into `output`.

### Module-level Variables

-   `ipsp2xml(:)`: Integer array, likely a map.
-   `npsp_pawxml`: Integer, count of XML PSPs.
-   `paw_setup(:)`: Array of `paw_setup_t`, potentially for multiple atom types if read by FoX in one go.
-   `paw_setuploc`: Single `paw_setup_t` instance, used by the FoX parser to accumulate data for the current PSP being parsed, and also by the Fortran readers.
-   Several private variables (`in_valenceStates`, `ndata`, `rp`, etc.) are used by the FoX SAX handlers to manage parsing state.

## Important Variables/Constants

-   `dpxml`: Real kind parameter `selected_real_kind(14)`, used for variables storing data read from XML to maintain precision.
-   `XML_RECL`: Integer parameter, maximum record length for sequential file access.
-   `NGAUSSIAN_MAX`: Integer parameter, maximum number of Gaussians in a fit for `gaussian_expansion_t`.

## Usage Examples

The primary use is internal to `libPAW`, specifically by `m_pawpsp` when an XML file is provided.

```fortran
! Conceptual usage within another module (e.g., m_pawpsp)
module m_process_xml_psp
    use m_pawxmlps
    use m_libpaw_defs, only: dp
    implicit none

    subroutine load_paw_setup_from_xml(psp_filename)
        character(len=*), intent(in) :: psp_filename
        type(paw_setup_t) :: loaded_paw_data
        real(dp) :: energy_cutoffs(3,2) ! Example structure for ecut

#if defined LIBPAW_HAVE_FOX
        ! Using FoX SAX parser (simplified representation)
        ! paw_setuploc is populated by FoX callbacks (begin_element, end_element, data_chunk)
        ! when processFile is called on psp_filename.
        ! After parsing, paw_setuploc would be copied to loaded_paw_data.
        ! This is typically orchestrated by a higher-level routine that calls FoX.
        ! For example, rdpawpsxml calls FoX's processFile.
        call rdpawpsxml(psp_filename, loaded_paw_data)
#else
        ! Using pure Fortran reader for header info
        call rdpawpsxml_header(energy_cutoffs, psp_filename, loaded_paw_data)
        ! And then potentially for core wavefunctions if needed separately
        ! call rdpawpsxml_core(...)
#endif
        ! Now loaded_paw_data contains the information from the XML file
        ! (fully if FoX was used, partially if only rdpawpsxml_header was called)

        ! ... process loaded_paw_data, transfer to pawtab_type ...

        call paw_setup_free(loaded_paw_data) ! Clean up

    end subroutine load_paw_setup_from_xml

end module m_process_xml_psp
```

## Dependencies and Interactions

-   **`libpaw.h`**: For preprocessor directives, especially `LIBPAW_HAVE_FOX`.
-   **`m_libpaw_defs` (`USE_DEFS`)**: For basic kinds (`dp`).
-   **`m_libpaw_tools` (`USE_MSG_HANDLING`)**: For error reporting (`LIBPAW_ERROR`, `LIBPAW_BUG`).
-   **`m_libpaw_memory` (`USE_MEMORY_PROFILING`)**: For memory allocation macros.
-   **`m_pawrad`**: Defines `pawrad_type` and provides `pawrad_init`, `pawrad_free`, `pawrad_ifromr`, `bound_deriv`. These are used when `rdpawpsxml_core` needs to spline data onto a different radial grid.
-   **`m_paw_numeric`**: Provides `paw_spline`, `paw_splint` used by `rdpawpsxml_core`.
-   **FoX SAX library**: If `LIBPAW_HAVE_FOX` is defined, this module heavily relies on FoX for XML parsing through the SAX interface (`processFile`, `getValue`, and the callback routines).
-   **`m_pawpsp`**: This module is primarily used by `m_pawpsp` (specifically `pawpsp_17in`) to handle XML formatted PAW datasets. The `paw_setup_t` structure populated here is then used by `m_pawpsp` to fill the main `pawtab_type` structure.

This module provides the bridge between PAW dataset files in XML format and the internal data structures of `libPAW`. The dual parsing capability (FoX or native Fortran) offers flexibility depending on available libraries.
