# BigDFT2Wannier.f90

## Overview

*This file provides a utility to interface BigDFT with the Wannier90 code. It handles the conversion of BigDFT wavefunction data into a format suitable for Wannier90 (e.g., .amn, .mmn, .eig files) and prepares inputs for a Wannier90 calculation. It can also perform a pre-check to select optimal virtual orbitals.*

## Key Components

*List of the primary functions, classes, or modules defined in the file, with a short description of each.*

*   **PROGRAM BigDFT2Wannier:** The main program that orchestrates the BigDFT to Wannier90 interface.
*   **SUBROUTINE read_inter_header:** Reads the header of the `input.inter` file which controls various parameters for the interface.
*   **SUBROUTINE read_inter_list:** Reads the list of virtual orbitals from `input.inter` if pre-check is not performed.
*   **SUBROUTINE read_nnkp_int_alloc:** Reads integer values from the Wannier90 `.nnkp` file for initial allocations.
*   **SUBROUTINE read_nnkp:** Reads detailed information from the Wannier90 `.nnkp` file (lattice vectors, k-points, projection centers, etc.).
*   **SUBROUTINE angularpart:** Calculates the angular part of spherical harmonics.
*   **SUBROUTINE radialpart:** Calculates the radial part of spherical harmonics.
*   **SUBROUTINE write_functions:** Writes cube files for spherical harmonics, their radial parts, or angular parts if requested.
*   **SUBROUTINE write_inter:** Writes/updates the `input.inter` file, potentially adding the list of chosen virtual orbitals after a pre-check.
*   **SUBROUTINE write_amn:** Writes the `.amn` file (overlap matrix between Bloch states and trial orbitals) for Wannier90.
*   **SUBROUTINE write_mmn:** Writes the `.mmn` file (overlap matrix between cell-periodic parts of Bloch states at neighboring k-points) for Wannier90.
*   **SUBROUTINE write_unk_bin:** Writes `UNK` files (Bloch wavefunctions on a real-space grid) for Wannier90, if requested.
*   **SUBROUTINE allocate_initial, mmnk_calculation_allocation, deallocate_projectors, deallocate_amnk_calculation, final_deallocations:** Various routines for memory management (allocation and deallocation of arrays).
*   **SUBROUTINE split_vectors_for_parallel:** Distributes orbital/vector data across parallel processes.
*   **MODULE BigDFT_API, bigdft_run, Poisson_Solver, etc.:** Utilizes various modules from the BigDFT codebase for core functionalities like MPI communication, input/output, orbital manipulation, and physical calculations.

## Important Variables/Constants

*Descriptions of any critical variables or constants defined in the file that affect its behavior.*

*   **seedname (character):** The base name for Wannier90 input/output files (e.g., `wannier.nnkp`, `wannier.amn`). Read from `input.inter`.
*   **filetype (character):** Format of BigDFT wavefunctions (e.g., "ETSF", "BIN"). Read from `input.inter`.
*   **n_occ (integer):** Number of occupied orbitals. Read from `input.inter`.
*   **n_virt (integer):** Number of virtual orbitals to be used for Wannierization. Read from `input.inter`.
*   **n_virt_tot (integer):** Total number of virtual orbitals available from BigDFT calculation (used in pre-check).
*   **pre_check (logical):** If true, performs a pre-check to determine the optimal virtual orbitals.
*   **residentity (logical):** If true, uses a resolution of the identity approach to construct virtual states/projectors.
*   **input (type(input_variables)):** Stores input parameters for the BigDFT run.
*   **atoms (type(atoms_data)):** Stores atomic structure information.
*   **orbs, orbsp, orbsv, orbsb (type(orbitals_data)):** Data structures to hold information about occupied, projector, virtual, and combined orbitals, respectively.
*   **amnk (real, allocatable):** Stores the overlap matrix <psi_nk | g_mn> (Bloch states vs. localized trial orbitals).
*   **mmnk_re, mmnk_im (real, allocatable):** Stores the real and imaginary parts of the overlap matrix <u_nk | u_mk+b>.
*   **kpts (real, allocatable):** Array storing k-point coordinates.
*   **ctr_proj, x_proj, y_proj, z_proj (real(kind=8), allocatable):** Arrays storing the centers and orientation of projection functions for Wannier90.
*   **l, mr, rvalue, zona (integer/real, allocatable):** Parameters defining the type (s, p, d, f, hybrids) and radial extent of the projection functions.

## Usage Examples

*If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.*

*The program is typically run after a BigDFT calculation. It requires an `input.inter` file and a `seedname.nnkp` file (Wannier90 format) to be present in the working directory or specified path.*

A typical command-line invocation might look like:
```bash
BigDFT2Wannier
```
(The program reads `input.inter` and `log.yaml` by default to get necessary parameters and BigDFT run information).

Alternatively, command-line options can be used (though the primary control is often via `input.inter` and `seedname.nnkp`):
```bash
BigDFT2Wannier -i specific_input.yaml
```

**input.inter (example structure):**
```
mywanniersystem  # seedname
ETSF             # filetype for wavefunctions
F F 4            # ResIdentity(F/T), WriteResId(F/T), n_occ
T 10 5           # Pre_check(F/T), n_virt_tot, n_virt_to_use
F F F F          # Write UNK, Write Sph Harm, Write Ang Part, Write Rad Part
<nx> <ny> <nz>   # Grid dimensions (can be read from BigDFT output)
data-mywannier   # Directory of BigDFT run output (containing wavefunctions)
# Optional: list of virtual bands if pre_check is False
```

## Dependencies and Interactions

*Notes on any dependencies the file has on other parts of the system or external libraries, and how it interacts with other components.*

*   **BigDFT Output Files:** Depends on output files from a previous BigDFT calculation, primarily wavefunction files (e.g., `wavefunction.etsf`, `virtuals.etsf`) and potentially log files (`log.yaml`) for run parameters. The location of these is often specified in `input.inter` or inferred.
*   **Wannier90 Input Files:**
    *   Requires a `seedname.nnkp` file, which is a standard Wannier90 input defining the trial orbitals, k-point mesh, and other parameters.
    *   Reads an `input.inter` file (custom to BigDFT2Wannier) that controls its operational mode, specifies input wavefunction formats, number of bands, etc.
*   **Wannier90 Output Files:** Generates files required by Wannier90:
    *   `seedname.amn`: Overlaps between Bloch states and trial orbitals.
    *   `seedname.mmn`: Overlaps between cell-periodic parts of Bloch states at neighboring k-points.
    *   `seedname.eig`: Eigenvalues of the Bloch states.
    *   Optionally `UNK*.s` files: Bloch states on a real-space grid.
*   **Wannier90 Executable/Library:** While this program *prepares* files for Wannier90, it does not directly call the Wannier90 executable itself. The generated files are then used as input for a separate Wannier90 run.
*   **Internal BigDFT Modules:** Heavily relies on various BigDFT modules for tasks such as:
    *   Reading BigDFT input/output.
    *   Handling atomic structures and orbital data.
    *   Performing wavefunction manipulations (e.g., reading, transforming from Daubechies basis to real space).
    *   MPI communications for parallel execution.
*   **Fortran Libraries:** Uses standard Fortran libraries and MPI for parallel processing.
```
