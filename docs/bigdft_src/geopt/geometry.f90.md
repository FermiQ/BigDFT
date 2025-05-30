# geopt/geometry.f90

## Overview

*This file provides a suite of Fortran subroutines for performing geometry optimization of molecular and periodic systems within the BigDFT code. It includes a main driver routine (`geopt`) that manages the optimization process and calls various specific optimization algorithms. It also contains routines for initializing optimization parameters, checking convergence, and handling atomic coordinates and constraints.*

## Key Components

*List of the primary functions, classes, or modules defined in the file, with a short description of each.*

*   **SUBROUTINE geopt_init:** Initializes default parameters for geometry optimization algorithms (e.g., tolerances, maximum iterations) stored in the `minpar` module (implicitly used, not defined in this file).
*   **SUBROUTINE geopt_set_verbosity(verb):** Sets the verbosity level for the geometry optimization output.
*   **SUBROUTINE geopt(runObj, outs, nproc, iproc, ncount_bigdft):** The main driver routine for geometry optimization. It initializes the process, selects the optimization method based on input, and iteratively calls the chosen algorithm until convergence or a maximum number of steps is reached.
*   **SUBROUTINE geopt_conditional(outs, it, fluct, frac_fluct, forcemax, check):** A generic routine to check for convergence based on maximum force component (`fmax`), force fluctuations (`fluct`), and specified thresholds.
*   **SUBROUTINE loop(runObj, outs, nproc, iproc, ncount_bigdft, fail):** A simple loop mode that repeatedly calls the BigDFT state calculation, useful for generating trajectories or when external tools handle the geometry updates.
*   **SUBROUTINE ab6md(runObj, outs, nproc, iproc, ncount_bigdft, fail):** Implements molecular dynamics using an interface to ABINIT's molecular dynamics routines (`abi_moldyn`).
*   **SUBROUTINE timeleft(tt):** Calculates the remaining CPU time against a specified limit.
*   **SUBROUTINE convcheck(fmax, fluctfrac_fluct, forcemax, check):** Checks if the geometry optimization has converged based on maximum force and fluctuation criteria.
*   **SUBROUTINE fnrmandforcemax(ff, fnrm, fmax, nat):** Calculates the L2-norm (`fnrm`) and the maximum component (`fmax`) of the forces `ff`.
*   **SUBROUTINE updatefluctsum(fnoise, fluct):** Updates a running average of force fluctuations.
*   **SUBROUTINE transforce(at, fxyz, sumx, sumy, sumz):** Calculates the total translational force on the system.
*   **SUBROUTINE rundiis(runObj, outs, nproc, iproc, ncount_bigdft, fail):** Implements a DIIS (Direct Inversion in the Iterative Subspace) based geometry optimization algorithm.
*   **SUBROUTINE fire(runObj, outs, nproc, iproc, ncount_bigdft, fail):** Implements the FIRE (Fast Inertial Relaxation Engine) algorithm, a damped molecular dynamics approach for geometry optimization.
*   **SUBROUTINE geometry_output(method, ncount_bigdft, it, fmax, fnrm, fluct, extra_tree):** Handles YAML output for geometry optimization steps, logging key metrics.
*   **SUBROUTINE f2fslave(runObj, outs, nproc, iproc, ncount_bigdft, fail):** Implements a client mode for geometry optimization via socket communication with an external master program (e.g., i-PI). Includes sub-routines for socket communication like `send_status`, `get_init`, `get_data`, `send_data`.
*   **SUBROUTINE rxyz_cart2int / rxyz_int2cart:** Convert atomic coordinates between Cartesian and internal (fractional with respect to lattice vectors) representations.
*   **SUBROUTINE invertmat:** Inverts a 3x3 matrix (or general NxN using LAPACK for N>3).

## Important Variables/Constants

*Descriptions of any critical variables or constants defined in the file that affect its behavior.*

*   **runObj (type(run_objects)):** A key derived type that encapsulates most of the BigDFT run information, including input parameters (`runObj%inputs`), atomic structure (`runObj%atoms`), etc.
*   **outs (type(state_properties)):** A derived type holding output properties from a BigDFT calculation, such as total energy (`outs%energy`) and atomic forces (`outs%fxyz`).
*   **parmin (MODULE minpar):** This module (defined elsewhere) stores various parameters controlling the minimization algorithms, e.g.:
    *   `parmin%approach`: String specifying the optimization algorithm (e.g., 'LBFGS', 'FIRE').
    *   `parmin%iter`: Current iteration number.
    *   `parmin%GTOL`, `parmin%XTOL`, `parmin%FTOL`: Convergence tolerances.
    *   `parmin%STPMIN`, `parmin%STPMAX`: Min/max step sizes for line searches.
*   **ncount_bigdft (integer):** Counter for the number of BigDFT self-consistent calculations performed.
*   **forcemax (real(gp)):** The input convergence criterion for the maximum force component.
*   **frac_fluct (real(gp)):** Factor used with force fluctuations for convergence check.
*   **dt, dtmax (real(gp)):** Timestep and maximum timestep for MD-like optimizers (FIRE, ab6md).
*   **alpha (real(gp)):** Damping parameter for the FIRE algorithm.
*   **velcur (real(gp) array):** Stores current atomic velocities in MD-based optimizers.
*   `ugeopt` (integer parameter): Fortran unit number for geometry optimization log file (`geopt.mon`).

## Usage Examples

*If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.*

*The `geopt` subroutine is the main entry point. It's typically called after an initial SCF calculation. The specific optimization algorithm is chosen based on input parameters.*

```fortran
USE bigdft_run ! For run_objects, state_properties
USE geopt_module ! (Conceptual - geopt routines are directly available or used via geopt)
USE minpar       ! For optimization parameters

IMPLICIT NONE

TYPE(run_objects) :: current_run
TYPE(state_properties) :: current_state
INTEGER :: num_procs, my_proc_id, total_scf_calls

! ... Initialize current_run with input parameters, atomic structure ...
! ... Perform an initial SCF calculation to populate current_state (energy, forces) ...
! CALL bigdft_state(current_run, current_state, initial_infocode)
! total_scf_calls = 1

! Set geometry optimization parameters (can be read from input)
CALL geopt_init() ! Initialize defaults
! runObj%inputs%geopt_approach = "FIRE" ! Or "LBFGS", "DIIS", etc.
! runObj%inputs%forcemax = 1.0E-4
! runObj%inputs%ncount_cluster_x = 100 ! Max SCF calls for geometry part

IF (my_proc_id == 0) THEN
  CALL yaml_open_stream('geopt_run.yaml')
END IF

! Perform geometry optimization
CALL geopt(current_run, current_state, num_procs, my_proc_id, total_scf_calls)

IF (my_proc_id == 0) THEN
  CALL yaml_close_stream()
  PRINT *, "Geometry optimization finished."
  PRINT *, "Total SCF calls:", total_scf_calls
  PRINT *, "Final energy:", current_state%energy
  ! ... Output final geometry, forces ...
END IF

! ...
```

## Dependencies and Interactions

*Notes on any dependencies the file has on other parts of the system or external libraries, and how it interacts with other components.*

*   **BigDFT Core (`bigdft_run`, `scfloop_API`):** The `geopt` routine and its sub-algorithms repeatedly call `bigdft_state` (or similar via `scfloop_API` for `ab6md`) to recompute energies and forces for new atomic configurations.
*   **Input Parameters (`module_input_keys`, `minpar`):** Relies on input parameters to select the optimization algorithm, set convergence criteria, define maximum iterations, and control algorithm-specific settings (e.g., timestep for FIRE). The `minpar` module holds many of these parameters.
*   **Atomic Structure (`module_atoms`):** Directly modifies atomic positions (`runObj%atoms%astruct%rxyz`) based on calculated forces and the chosen optimization algorithm.
*   **Force Calculations (e.g., `forces.f90` content):** The forces calculated by routines in (conceptually) `forces.f90` (like `outs%fxyz`) are the primary input driving the geometry updates.
*   **Output and Logging (`yaml_output`, `module_io` implicitly):** Uses YAML for structured logging of the optimization progress (`geopt.mon`). Atomic coordinate files are written at various stages.
*   **External Libraries (Potentially):**
    *   LAPACK (`gesv` called in `rundiis`): For solving linear systems in the DIIS algorithm.
    *   Socket Communication (`f90sockets` used in `f2fslave`): For client-server based optimization with external drivers like i-PI.
    *   ABINIT MD (`m_ab6_moldyn` used in `ab6md`): For performing molecular dynamics runs using ABINIT's MD engine.
*   **Optimization Algorithms:** This file implements or interfaces with various local optimization algorithms:
    *   LBFGS (via `lbfgsdriver`, not fully shown but implied by `minpar` usage)
    *   Conjugate Gradient (via `conjgrad`)
    *   Steepest Descent variants (VSSD via `vstepsd`)
    *   FIRE (`fire`)
    *   DIIS (`rundiis`)
    *   SQNM/SBFGS (`sqnm`)
*   **Constraints:** While not extensively detailed in the provided snippet, geometry optimization often involves handling constraints (e.g., fixed atoms, fixed bond lengths). The `keep_internal_coordinates_constraints` routine in `forces.f90` suggests such capabilities exist and would be managed alongside geometry updates.
```
