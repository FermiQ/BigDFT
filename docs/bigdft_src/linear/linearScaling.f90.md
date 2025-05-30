# linear/linearScaling.f90

## Overview

*This file provides the core Fortran subroutines for performing linear-scaling Density Functional Theory (DFT) calculations within the BigDFT framework. The central routine, `linearScaling`, manages the self-consistent cycle, which involves optimizing localized support functions (often referred to as "minimal basis" or TMB in the code) and the density kernel. This approach is designed to tackle very large systems where traditional cubic-scaling DFT becomes computationally prohibitive.*

## Key Components

*List of the primary functions, classes, or modules defined in the file, with a short description of each.*

*   **SUBROUTINE linearScaling(...):** The main driver for linear-scaling DFT calculations. It orchestrates the overall self-consistent field (SCF) procedure, including:
    *   Initialization of parameters and data structures.
    *   Outer loop for self-consistency, potentially with multiple accuracy levels (low and high).
    *   Inner loops for:
        *   Optimization of support functions (calls `getLocalizedBasis`).
        *   Optimization of the density kernel (coefficients of the support functions) (calls `get_coeff` and manages the SCF cycle for the kernel via `scf_kernel`).
    *   Calculation of total energy and other properties.
    *   Potential post-processing like charge analysis or wavefunction output.
*   **SUBROUTINE scf_kernel(nit_scc, remove_coupling_terms, update_phi):** (Internal to `linearScaling`) Manages the self-consistent iterations for optimizing the density kernel. It repeatedly calls `get_coeff` to update the kernel and `sumrho_for_TMBs` and `updatePotential` to update density and potential until convergence.
*   **SUBROUTINE check_inputguess():** (Internal to `linearScaling`) Checks the quality of an initial guess for wavefunctions/kernel, potentially deciding on the strategy for starting the SCF (e.g., number of low/high accuracy iterations).
*   **SUBROUTINE check_for_exit():** (Internal to `linearScaling`) Determines if the outer SCF loop should terminate based on convergence criteria or maximum iteration counts.
*   **SUBROUTINE printSummary():** (Internal to `linearScaling`) Prints a summary of values for the current SCF iteration (e.g., energy, delta, kernel method).
*   **SUBROUTINE print_info(final):** (Internal to `linearScaling`) Prints a more detailed summary of the SCF cycle, including information about support function and kernel optimization convergence.
*   **SUBROUTINE intermediate_forces():** (Internal to `linearScaling`, experimental/debug) Calculates forces at an intermediate step of the SCF cycle.
*   **SUBROUTINE get_lagrange_mult(cdft_it, vgrad):** (Internal to `linearScaling`) Handles Lagrange multiplier updates for constrained DFT (CDFT) calculations.
*   **SUBROUTINE output_fragment_rotations(...):** Calculates and outputs rotation information between fragments, likely used in fragment-based calculations.
*   **SUBROUTINE get_boundary_weight(...):** Analyzes the weight of support functions at the boundaries of their localization regions.
*   **SUBROUTINE allocate_local_arrays / deallocate_local_arrays:** (Internal to `linearScaling`) Utility routines for memory management.

## Important Variables/Constants

*Descriptions of any critical variables or constants defined in the file that affect its behavior.*

*   **tmb (type(DFT_wavefunction)):** Derived type representing the "minimal basis" or localized support functions. It holds coefficients, localization information, and associated matrices (overlap, Hamiltonian, density kernel).
*   **KSwfn (type(DFT_wavefunction)):** Derived type representing Kohn-Sham wavefunctions, often used here to define the space or provide initial guesses.
*   **tmb%linmat%kernel_ (type(matrices)):** Stores the density kernel matrix in the basis of the support functions.
*   **tmb%linmat%ovrlp_ (type(matrices)):** Stores the overlap matrix of the support functions.
*   **tmb%linmat%ham_ (type(matrices)):** Stores the Hamiltonian matrix in the basis of the support functions.
*   **tmb%psi (real(wp) array):** Coefficients of the localized support functions in the underlying wavelet basis.
*   **tmb%coeff (real(wp) array):** Coefficients used to construct Kohn-Sham orbitals from support functions (if diagonalization is performed).
*   **input%lin%...:** Various input parameters prefixed with `lin` control aspects of the linear-scaling calculation, such as:
    *   `input%lin%scf_mode`: Defines the SCF strategy (e.g., `LINEAR_FOE`, `LINEAR_DIRECT_MINIMIZATION`, `LINEAR_MIXDENS_SIMPLE`).
    *   `input%lin%nit_lowaccuracy`, `input%lin%nit_highaccuracy`: Number of iterations for different accuracy levels.
    *   `input%lin%lowaccuracy_conv_crit`, `input%lin%highaccuracy_conv_crit`: Convergence criteria for density/potential.
    *   `input%lin%support_functions_converged`: Convergence criterion for support function optimization.
    *   `input%lin%order_taylor`: Order of Taylor expansion for matrix functions (e.g., in FOE).
    *   `input%lin%constrained_dft`: Logical flag to enable constrained DFT.
*   **denspot (type(DFT_local_fields)):** Stores the electron density (`rhov`) and local potentials (`V_ext`, `V_hartree`, `V_xc`).
*   **energs (type(energy_terms)):** Stores different components of the total energy.
*   **target_function (integer):** Specifies the objective function for support function optimization (e.g., trace minimization, energy minimization).

## Usage Examples

*Provide a conceptual Fortran code snippet illustrating the main steps in a linear-scaling calculation:*

```fortran
USE linear_scaling_module_or_routines ! (linearScaling.f90 provides subroutines directly)
USE module_types          ! For DFT_wavefunction, atoms_data, input_variables etc.
USE module_hamiltonian    ! For DFT_local_fields, DFT_PSP_projectors
USE module_energy         ! For energy_terms

IMPLICIT NONE

INTEGER :: ip, np, info_code
TYPE(DFT_wavefunction) :: tmb_wfn, ks_wfn_guess
TYPE(atoms_data) :: atom_info
TYPE(input_variables) :: inp
REAL(GP), DIMENSION(:,:), ALLOCATABLE :: positions_xyz, forces_pulay_out
TYPE(DFT_local_fields) :: density_potential
REAL(GP), ALLOCATABLE :: old_density_potential_array(:)
TYPE(DFT_PSP_projectors) :: psp_data
TYPE(GPU_pointers) :: gpu_ptrs
TYPE(energy_terms) :: energy_components
REAL(GP) :: total_energy_out
TYPE(system_fragment), POINTER :: reference_fragments(:) => NULL()
TYPE(cdft_data) :: constrained_dft_params

! ... Initialize MPI (ip, np), atom_info, inp (input parameters), ...
! ... positions_xyz, psp_data, gpu_ptrs. ...
! ... Setup initial guess for ks_wfn_guess and tmb_wfn (support functions). ...
! ... Allocate and initialize density_potential, old_density_potential_array, ...
! ... energy_components, forces_pulay_out. ...
! ... Setup constrained_dft_params if doing CDFT.

! Call the main linear scaling routine
CALL linearScaling(ip, np, ks_wfn_guess, tmb_wfn, atom_info, inp, positions_xyz, &
                   density_potential, old_density_potential_array, psp_data, gpu_ptrs, &
                   energy_components, total_energy_out, forces_pulay_out, info_code, &
                   reference_fragments, constrained_dft_params, &
                   fdisp_dummy, fion_dummy) ! Pass dummy or actual dispersion/ionic forces

IF (info_code == 0) THEN
  PRINT *, "Linear scaling calculation converged."
  PRINT *, "Total Energy:", total_energy_out
  ! ... Process forces_pulay_out and other forces if calculated within linearScaling ...
ELSE
  PRINT *, "Linear scaling calculation did NOT converge properly. Infocode:", info_code
END IF

! ... Clean up allocated memory ...
```

## Dependencies and Interactions

*Notes on any dependencies the file has on other parts of the system or external libraries, and how it interacts with other components.*

*   **Core DFT Modules (`module_base`, `module_types`, `module_atoms`, `rhopotential`, `module_hamiltonian` implicitly):** Relies on these for fundamental data structures, atomic information, density and potential updates, Hamiltonian construction, and energy calculations.
*   **Support Function Optimization (`get_basis`, `diis_sd_optimization`):** Uses routines from these modules to optimize the localized support functions (TMBs). This involves iterative methods like DIIS or steepest descent.
*   **Density Kernel Optimization (`get_kernel`, `foe_base`, `sparsematrix` family):** Employs various techniques to optimize the density kernel, including:
    *   Direct minimization methods.
    *   Fermi Operator Expansion (FOE).
    *   Diagonalization (for smaller systems or specific purposes).
    *   Sparse matrix algebra for efficient handling of large matrices.
*   **Communication (`communications_base`, `communications`):** Heavily dependent on MPI-based communication routines for distributing data (wavefunctions, matrices, density, potential) and parallelizing computations across multiple processors.
*   **Grid Operations (`locreg_operations`, `locregs_init`):** Uses wavelet-based grid operations for representing functions and applying operators.
*   **Input/Output (`io`, `yaml_output`):** For reading input parameters, writing output data (wavefunctions, matrices, energies, logs), and YAML-based logging.
*   **Post-processing (`postprocessing_linear`):** May call routines for analyzing results, such as Loewdin charge analysis or constructing Kohn-Sham orbitals from the optimized support functions and kernel.
*   **Fragment and CDFT (`module_fragments`, `constrained_dft`):** Interacts with modules for fragment-based calculations and constrained DFT.
*   **Forces (`forces.f90` - implicitly):** While this file focuses on the electronic structure, the resulting density kernel and support functions are essential for subsequent force calculations in a linear-scaling context (often handled by routines in `forces_linear.f90`).
*   **Alternative to Cubic Scaling:** Provides the framework for O(N) DFT calculations, offering a significant advantage for very large systems where the O(N^3) scaling of traditional DFT becomes a bottleneck.
```
