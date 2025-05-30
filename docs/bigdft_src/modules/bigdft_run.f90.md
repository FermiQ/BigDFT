# modules/bigdft_run.f90

## Overview

*This file defines the Fortran module `bigdft_run`, which is central to orchestrating calculations within the BigDFT program. It provides data structures and routines to manage a single computational run or a series of calculations. This includes handling input parameters, atomic structures, restart information, calculation states (like energy and forces), and coordinating different calculation phases such as ground state SCF, geometry optimization, or molecular dynamics.*

## Key Components

*List of the primary functions, classes, or modules defined in the file, with a short description of each.*

*   **MODULE bigdft_run:** The main module containing all definitions and procedures.
*   **TYPE QM_restart_objects:** A derived data type for storing information needed to restart a Quantum Mechanics (QM) calculation or for post-processing. This includes Kohn-Sham wavefunctions (`KSwfn`), support functions for linear scaling (`tmb`), old atomic positions, and GPU pointers.
*   **TYPE MM_restart_objects:** A derived data type for supplementary information related to Molecular Mechanics (MM) or QM/MM calculations.
*   **TYPE run_objects:** A central derived data type that encapsulates all information for a BigDFT run. It holds pointers to input parameters (`inputs`), atomic data (`atoms`), QM restart objects (`rst`), MM restart objects (`mm_rst`), user input dictionaries (`user_inputs`), and potentially sections for multi-stage runs.
*   **TYPE state_properties:** A derived data type to store the results of a calculation, such as total energy (`energy`), energy components (`energs`), atomic forces (`fxyz`), and stress tensor (`strten`).
*   **SUBROUTINE bigdft_init(options, with_taskgroups):** Initializes the BigDFT environment, including MPI setup and task group configurations, based on command-line options or an input dictionary. It identifies the runs to be performed by the current BigDFT instance.
*   **SUBROUTINE bigdft_command_line_options(options):** Parses command-line arguments and populates an options dictionary.
*   **SUBROUTINE run_objects_init(runObj, run_dict, source):** Initializes a `run_objects` variable, setting up input parameters, atomic structures, and restart information based on a run dictionary and optionally a source `run_objects` for restart.
*   **SUBROUTINE set_run_objects(runObj):** Parses the input dictionary within `runObj%user_inputs` to populate `runObj%inputs` and `runObj%atoms`.
*   **SUBROUTINE bigdft_state(runObj, outs, infocode):** The core routine that performs a single state calculation (e.g., SCF for QM, energy/force for MM). It takes a `run_objects` as input and populates a `state_properties` with results.
*   **SUBROUTINE process_run(id, runObj, outs):** A higher-level routine that manages a calculation run. It calls `bigdft_state` and, if requested, proceeds to geometry optimization (`geopt`) or molecular dynamics (`md`).
*   **SUBROUTINE quantum_mechanical_state(runObj, outs, infocode):** Specifically handles the quantum mechanical part of a calculation by calling the `cluster` routine (the main SCF driver).
*   **SUBROUTINE geopt(runObj, outs, nproc, iproc, ncount_bigdft):** (Interface, implementation in `geopt/geometry.f90`) Drives the geometry optimization process.
*   **SUBROUTINE md(runObj, outs, nproc, iproc):** (Interface, implementation likely in an MD module) Drives molecular dynamics simulations.
*   **SUBROUTINE init_state_properties / deallocate_state_properties:** Allocate and deallocate `state_properties` structures.
*   **SUBROUTINE free_run_objects / release_run_objects:** Deallocate and release `run_objects` structures, managing reference counting for shared components.
*   **Various accessor functions:** `bigdft_nat`, `bigdft_norb`, `bigdft_get_eval`, `bigdft_get_cell`, `bigdft_get_rxyz`, `bigdft_set_rxyz`, etc., provide interfaces to query information from `run_objects`.
*   **Signal handling routines:** `run_objects_type_init`, `process_run_type_init` for setting up signals used for plugins or advanced workflow control.

## Important Variables/Constants

*Descriptions of any critical variables or constants defined in the file that affect its behavior.*

*   **runObj%inputs (type(input_variables), pointer):** Stores all input parameters parsed from the input file(s).
*   **runObj%atoms (type(atoms_data), pointer):** Stores atomic coordinates, species information, pseudopotentials, etc.
*   **runObj%rst (type(QM_restart_objects), pointer):** Contains data for restarting QM calculations, including wavefunctions.
*   **runObj%run_mode (type(f_enumerator), pointer):** Specifies the type of calculation (e.g., `QM_RUN_MODE`, `LENNARD_JONES_RUN_MODE`, `MULTI_RUN_MODE`).
*   **outs%energy (real(gp)):** Total energy of the system.
*   **outs%fxyz (real(gp), dimension(:,:), pointer):** Atomic forces.
*   **outs%strten (real(gp), dimension(6)):** Stress tensor components.
*   **INPUT_POLICY_SCRATCH, INPUT_POLICY_MEMORY, INPUT_POLICY_DISK (integer parameters):** Define policies for how a calculation should start (e.g., from scratch, from memory, from disk files).

## Usage Examples

*Provide a conceptual Fortran code snippet from the main BigDFT program showing how this module might be used:*

```fortran
PROGRAM main_bigdft_driver
  USE bigdft_run
  USE module_types ! For state_properties
  USE dictionaries   ! For options dictionary
  IMPLICIT NONE

  TYPE(dictionary), POINTER :: opts
  TYPE(run_objects) :: current_run
  TYPE(state_properties) :: current_state_props
  INTEGER :: ierr, i_run

  CALL f_lib_initialize() ! Initialize Fortran environment

  ! 1. Get command line options
  CALL bigdft_command_line_options(opts)

  ! 2. Initialize BigDFT environment (MPI, task groups, identify runs)
  CALL bigdft_init(opts)

  ! 3. Loop over identified runs
  DO i_run = 0, bigdft_nruns(opts) - 1
    ! 3a. Initialize run_objects for the current run
    ! (run_dict is obtained from opts // 'BigDFT' // i_run)
    CALL run_objects_init(current_run, opts // 'BigDFT' // i_run)

    ! Skip run if indicated (e.g. logfile already exists and 'skip' is true)
    IF (current_run%user_inputs .get. SKIP_RUN%name, .FALSE.) CYCLE

    ! 3b. Initialize state_properties
    CALL init_state_properties(current_state_props, bigdft_nat(current_run))

    ! 3c. Process the run (e.g., SCF, geopt, MD)
    CALL process_run(bigdft_run_id_toa(current_run), current_run, current_state_props)

    ! 3d. Clean up for the current run
    CALL deallocate_state_properties(current_state_props)
    CALL free_run_objects(current_run)
  END DO

  CALL dict_free(opts)
  CALL bigdft_finalize(ierr) ! Finalize MPI etc.
  CALL f_lib_finalize()

END PROGRAM main_bigdft_driver
```

## Dependencies and Interactions

*Notes on any dependencies the file has on other parts of the system or external libraries, and how it interacts with other components.*

*   **Core Data Types (`module_defs`, `module_types`, `module_atoms`):** Depends on these for fundamental data structures like `input_variables`, `atoms_data`, `DFT_wavefunction`, `energy_terms`, and various enumerators.
*   **Input System (`module_input_dicts`, `module_input_keys`, `dictionaries`, `yaml_parse`):** Heavily interacts with the input system to parse command-line options, read input files (YAML format), and manage run parameters through dictionaries.
*   **Calculation Engines:**
    *   **QM Calculations (`cluster` routine - external):** The `quantum_mechanical_state` subroutine calls the main `cluster` routine (likely in `cluster.f90` or similar) to perform the core self-consistent DFT calculations.
    *   **MM Calculations (various modules like `module_lj`, `module_tersoff`, etc.):** `bigdft_state` can dispatch to various classical potential routines if `run_mode` is set accordingly.
    *   **Geometry Optimization (`geopt/geometry.f90`):** `process_run` calls `geopt` (defined in `geopt/geometry.f90`) if geometry optimization is requested.
    *   **Molecular Dynamics (MD modules):** `process_run` can call `md` if MD is requested.
*   **Force Calculation (`forces.f90` implicitly):** Relies on force calculations (performed within `cluster` or MM routines) which populate `outs%fxyz`, used by geometry optimization or MD.
*   **MPI and Communication (`communications_base`):** Uses MPI for parallel execution and data distribution, managed via `bigdft_mpi` global object.
*   **Logging and Output (`yaml_output`, `module_io` implicitly):** Uses YAML for structured logging and standard I/O routines for writing files (e.g., atomic coordinates, restart files).
*   **Plugin System (`module_f_objects`, `f_utils`):** Supports a plugin architecture through signals and dynamic library loading, allowing extensibility.
*   **Role as Coordinator:** This module acts as a high-level coordinator. The `run_objects` type serves as a container passing all necessary information between different parts of the code (input parsing, SCF engine, geometry optimization, MD, I/O). `process_run` and `bigdft_state` are central routines that manage the flow of a typical BigDFT calculation.
```
