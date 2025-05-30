# forces.f90

## Overview

*This file contains a collection of Fortran subroutines dedicated to calculating atomic forces and related quantities (like stress tensor) within the BigDFT code. These are crucial for performing geometry optimizations, molecular dynamics simulations, and understanding the interactions within the simulated system. The calculations include local contributions (from local potential and core charge density) and non-local contributions (from pseudopotential projectors).*

## Key Components

*List of the primary functions, classes, or modules defined in the file, with a short description of each.*

*   **SUBROUTINE calculate_forces:** The main driver routine that orchestrates the calculation of total atomic forces. It calls other subroutines to compute different contributions (local, non-local, core, ionic, dispersion, Pulay) and can also calculate the stress tensor.
*   **SUBROUTINE local_forces:** Calculates the local part of the atomic forces arising from the interaction of the electron density with the local potential (e.g., Hartree and local pseudopotential). It also computes the local contribution to the stress tensor.
*   **SUBROUTINE nonlocal_forces:** Calculates the non-local part of the atomic forces originating from the pseudopotential projectors acting on the wavefunctions. It handles both cubic and linear scaling modes and contributes to the non-local stress tensor.
*   **SUBROUTINE rhocore_forces:** Calculates the contribution to the forces due to the interaction of the core charge density (if Non-Linear Core Correction - NLCC is used) with the exchange-correlation potential.
*   **SUBROUTINE soft_PCM_forces:** Calculates force contributions specific to the Polarizable Continuum Model (PCM), accounting for the interaction between the quantum system and a dielectric continuum.
*   **SUBROUTINE atomic_density_matrix_delta:** Calculates the difference between the target and calculated atomic density matrix, likely used in schemes with occupancy control.
*   **SUBROUTINE calc_coeff_derproj:** Calculates coefficients for the derivatives of projectors, used in force calculations (details specific to implementation).
*   **SUBROUTINE elim_torque_reza:** Eliminates rotational forces with respect to the center of mass, ensuring that only translational forces remain for isolated systems.
*   **SUBROUTINE cross:** Computes the cross product of two 3D vectors.
*   **SUBROUTINE moment_of_inertia:** Calculates the moment of inertia tensor for a set of atoms.
*   **SUBROUTINE normalizevector:** Normalizes a given vector.
*   **SUBROUTINE symm_stress:** Symmetrizes the stress tensor based on the system's symmetries.
*   **SUBROUTINE symmetrise_forces:** Symmetrizes the atomic forces based on the system's symmetries.
*   **SUBROUTINE local_hamiltonian_stress:** Calculates the stress tensor contribution from the local part of the Hamiltonian (kinetic energy).
*   **SUBROUTINE internal_forces:** Transforms forces from Cartesian to internal coordinates and applies constraints.
*   **SUBROUTINE constraints_internal / keep_internal_coordinates_constraints:** Apply constraints to atomic positions based on internal coordinates.

## Important Variables/Constants

*Descriptions of any critical variables or constants defined in the file that affect its behavior.*

*   **fxyz (real(gp), dimension(3,nat), intent(inout)):** Array storing the Cartesian components of the forces on each atom. This is the primary output of many force routines.
*   **strten (real(gp), dimension(6), intent(out)):** Array storing the components of the stress tensor (Voigt notation: xx, yy, zz, yz, xz, xy).
*   **psi (real(wp), dimension(...), intent(in)):** Array containing the wavefunctions (coefficients in a basis).
*   **rho (real(wp), dimension(...), intent(in)):** Array containing the electron density on the grid.
*   **pot (real(wp), dimension(...), intent(in)):** Array containing the total local potential (Hartree + XC + local pseudopotential) on the grid.
*   **potxc (real(wp), dimension(...), intent(in)):** Array containing the exchange-correlation potential on the grid.
*   **nlpsp (type(DFT_PSP_projectors)):** Derived type data structure holding information about the non-local pseudopotential projectors.
*   **atoms (type(atoms_data)):** Derived type data structure holding atomic information (positions, types, pseudopotential parameters, etc.).
*   **Glr (type(locreg_descriptors)):** Derived type describing the local region grids and wavelet descriptors.
*   **ob (type(orbital_basis)):** Derived type describing the orbital basis.
*   **tmb (type(DFT_wavefunction)):** Derived type for DFT wavefunctions, likely used in linear scaling context.
*   **fion (real(gp), dimension(3,nat), intent(in)):** Forces due to ion-ion interactions (Ewald sum).
*   **fdisp (real(gp), dimension(3,nat), intent(in)):** Dispersion forces (e.g., from DFT-D methods).
*   **fpulay (real(gp), dimension(3,nat), intent(in)):** Pulay forces, arising from the basis set's dependence on atomic positions.

## Usage Examples

*If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.*

*Within the BigDFT workflow, after a self-consistent field (SCF) calculation has converged, forces can be computed as follows (conceptual):*

```fortran
USE module_forces ! (Note: forces.f90 doesn't declare a single module, but its routines would be used)
USE module_types  ! For atoms_data, DFT_PSP_projectors, locreg_descriptors, etc.
USE module_io     ! For reading/writing data
IMPLICIT NONE

! ... Declare variables for atoms, wavefunctions (psi_data), density (rho_data), ...
! ... potentials (potential_data, potential_xc_data), projectors (psp_projectors), ...
! ... grid descriptors (grid_lr), force arrays (forces_total, forces_ionic, forces_dispersion, forces_pulay), etc. ...
! ... and other necessary inputs for 'calculate_forces' like rxyz, hx, hy, hz, etc.

REAL(GP), DIMENSION(3, atoms%astruct%nat) :: total_forces, ionic_forces, dispersion_forces, pulay_forces
REAL(GP), DIMENSION(6) :: stress_tensor, ewald_stress, hartree_stress, xc_stress
LOGICAL :: compute_stress
INTEGER :: iproc, nproc, psolver_grp_size, i3start, n3plane_loc, nspin_val, imode_val
! ... (Initialize these variables based on the BigDFT run) ...

! Ensure all input data (density, potentials, wavefunctions, projectors) are up-to-date.
! Projectors might need to be refilled if atomic positions changed:
CALL fill_projectors(grid_lr, [hx,hy,hz], atoms%astruct, orbital_basis_data, rxyz, psp_projectors, 0)

! Initialize total_forces to zero or with specific components
total_forces = 0.0_GP

! Call the main force calculation routine
CALL calculate_forces(iproc, nproc, psolver_grp_size, grid_lr, atoms, orbital_basis_data, &
                      psp_projectors, rxyz, hx, hy, hz, &
                      density_potential_box_distribution, &
                      i3start, n3plane_loc, nspin_val, &
                      .TRUE., ! refill_proj flag
                      gather_array_info, rho_data, potential_data, potential_xc_data, &
                      size_of_psi_array, psi_data, ionic_forces, dispersion_forces, total_forces, &
                      compute_stress, ewald_stress, hartree_stress, xc_stress, stress_tensor, &
                      pressure_output, psp_offset_value, imode_val, tmb_data, pulay_forces)

! 'total_forces' now contains the calculated atomic forces.
! 'stress_tensor' (if compute_stress was .TRUE.) contains the stress.

! ... Further processing like geometry optimization step or writing output ...
```

## Dependencies and Interactions

*Notes on any dependencies the file has on other parts of the system or external libraries, and how it interacts with other components.*

*   **Core Data Modules (`module_base`, `module_types`, `module_atoms`):** Depends on these for fundamental data types (e.g., `atoms_data`, `DFT_PSP_projectors`, `locreg_descriptors`, `orbital_basis`, `DFT_wavefunction`), MPI communication structures, and basic parameters.
*   **Grid and Wavelet Modules (`locregs`, `module_dpbox`):** Uses these for grid information, wavelet descriptors, and iterating over distributed density/potential data.
*   **Pseudopotential Projectors (`psp_projectors_base`, `orbitalbasis`):** Relies on modules that manage and apply pseudopotential projectors, which are essential for the non-local force contributions.
*   **Hamiltonian and Potential (`module_xc`, `module_hamiltonian` - implicitly):** The forces are derivatives of the energy, which involves terms from the Hamiltonian, including exchange-correlation potentials.
*   **Input Data:** Requires converged electron density (`rho`), wavefunctions (`psi`), total local potential (`pot`), XC potential (`potxc`), atomic positions (`rxyz`), and pseudopotential information (`nlpsp`, `atoms`).
*   **Output to Other Components:**
    *   The calculated forces (`fxyz`) are the primary input for geometry optimization algorithms (e.g., BFGS, conjugate gradient) that aim to find energy minima.
    *   In molecular dynamics simulations, these forces are used to propagate atomic trajectories according to Newton's equations of motion.
    *   The stress tensor (`strten`) can be used for cell optimization under constant pressure or for analyzing mechanical properties.
*   **Linear Scaling (`forces_linear`):** Includes specific routines (e.g., `nonlocal_forces_linear`) for calculating forces in the context of linear scaling DFT methods, which use different data structures and algorithms.
*   **Symmetry (`m_ab6_symmetry`):** Uses symmetry routines to symmetrize forces and the stress tensor, which can be important for maintaining system symmetry during simulations and for calculations with special k-points.
```
