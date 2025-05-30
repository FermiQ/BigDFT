# BigDFT.f90

## Overview

*This file contains the main program for the BigDFT electronic structure code. It handles the overall workflow, including input parsing, initialization, calculation execution, and output generation.*

## Key Components

*List of the primary functions, classes, or modules defined in the file, with a short description of each.*

*   **PROGRAM BigDFT:** The main executable unit of the BigDFT code. It orchestrates the entire calculation process.
*   **MODULE module_base:** Provides fundamental types and utilities used throughout the BigDFT code. (Referenced via `use module_base`)
*   **MODULE bigdft_run:** Contains the core routines for performing the DFT calculations. (Referenced via `use bigdft_run`)
*   **MODULE public_keys:** Defines public keys or constants used in the code, like `SKIP_RUN`. (Referenced via `use public_keys, only: SKIP_RUN`)
*   **SUBROUTINE f_lib_initialize():** Initializes the underlying Fortran library environment. (External, called)
*   **SUBROUTINE bigdft_command_line_options(options):** Parses command-line arguments and options provided to the BigDFT executable. (Called)
*   **SUBROUTINE bigdft_init(options):** Initializes the BigDFT environment, including MPI setup and run naming conventions, based on the parsed options. (Called)
*   **SUBROUTINE run_objects_init(runObj, run):** Initializes the `run_objects` type, likely setting up data structures for a specific calculation run. (Called within a loop)
*   **SUBROUTINE init_state_properties(outs, natoms):** Initializes data structures for storing the properties of the system state. (Called)
*   **SUBROUTINE bigdft_get_run_properties(run, posinp_id):** Retrieves specific properties or settings for the current run. (Called)
*   **SUBROUTINE process_run(posinp_id, runObj, outs):** The main routine that processes a given calculation run, likely performing the SCF cycle and other computational steps. (Called)
*   **SUBROUTINE deallocate_state_properties(outs):** Deallocates memory used for storing state properties. (Called)
*   **SUBROUTINE free_run_objects(runObj):** Deallocates memory associated with the `run_objects` type. (Called)
*   **SUBROUTINE bigdft_finalize(ierr):** Finalizes the BigDFT execution, handling cleanup and MPI finalization. (Called)
*   **SUBROUTINE f_lib_finalize():** Finalizes the underlying Fortran library environment. (External, called)

## Important Variables/Constants

*Descriptions of any critical variables or constants defined in the file that affect its behavior.*

*   **runObj (type(run_objects)):** A derived-type variable that likely encapsulates all data and parameters for a single BigDFT calculation run.
*   **options (type(dictionary), pointer):** A pointer to a dictionary-like structure holding the command-line options and input parameters.
*   **runs (type(dictionary), pointer):** A pointer to a dictionary structure that seems to manage multiple runs or configurations defined in the input.
*   **run (type(dictionary), pointer):** A pointer used to iterate through individual run configurations within the `runs` dictionary.
*   **ierr (integer):** Standard integer variable used for error status returns from subroutines.
*   **posinp_id (character(len=60)):** Stores an identifier related to the input, possibly for naming or referencing specific parts of the input.
*   **outs (type(state_properties)):** A derived-type variable holding output properties of the calculation.
*   **SKIP_RUN (public_keys):** A constant imported from `public_keys` module, used to indicate if a particular run should be skipped.

## Usage Examples

*If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.*

```
bigdft < input_file.yaml
```

## Dependencies and Interactions

*Notes on any dependencies the file has on other parts of the system or external libraries, and how it interacts with other components.*

*   **Internal Module Dependencies:**
    *   `module_base`: For basic definitions and utilities.
    *   `bigdft_run`: For core calculation routines.
    *   `public_keys`: For shared constants/keys.
*   **External Library Dependencies:** Implicitly depends on an underlying Fortran library (e.g., for `f_lib_initialize`, `f_lib_finalize`) and an MPI library for parallel execution.
*   **Input Files:** Interacts with input files, typically in YAML format, which define the system to be calculated and various parameters. The `bigdft_command_line_options` and subsequent processing of the `options` dictionary handle this.
*   **Output Files:** Generates various output files, including log files detailing the calculation progress, files containing wavefunctions, energies, and other computed properties. The specifics are handled within `process_run` and its sub-calls.
*   **Interactions:** The main program `BigDFT` coordinates the calls to various subroutines and modules to perform the sequence of operations: initialization, option parsing, loop over calculation runs (if any), processing each run, and finalization.
```
