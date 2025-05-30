# BZ.py

## Overview

*This Python file, `BZ.py`, is part of the PyBigDFT suite and provides tools for handling Brillouin Zone (BZ) representations, k-point paths, and band structure data. It leverages libraries such as NumPy, ASE (Atomistic Simulation Environment), and spglib for tasks like symmetry analysis, k-point generation, and data interpolation. This module is essential for setting up and analyzing calculations involving periodic systems, particularly for band structure visualization.*

## Key Components

*List of the primary functions, classes, or modules defined in the file, with a short description of each.*

*   **`get_ev(ev, keys=None, ikpt=1)` (Function):**
    *   A utility function to extract eigenvalue energies from a dictionary-like structure (`ev`) which typically represents eigenvalue data from a BigDFT logfile. It can filter by specific keys, k-point index, and spin channel.

*   **`CLASS BandArray(numpy.ndarray)`:**
    *   A custom NumPy array subclass designed to store band energy data for a single k-point. It handles data for different spin channels and associates k-point information (index, coordinates, weight) with the energy array.
    *   **`__new__(cls, *args, **kwargs)`:** Constructor that can take raw data or parse eigenvalue data from a BigDFT log output (via `logdata` kwarg) to populate the array.
    *   **`__init__(self, *args, **kwargs)`:** Initializes k-point metadata (`ikpt`, `kpt`, `kwgt`).
    *   **`set_kpt(self, ikpt, kpt, kwgt=1.0)`:** Sets or updates the k-point metadata.
    *   **`__add__(self, b)`:** Defines addition for `BandArray` objects, ensuring compatibility of k-point information.

*   **`CLASS BZPath()`:**
    *   Represents a path in the Brillouin Zone, typically used for band structure calculations.
    *   **`__init__(self, lattice, path, special_points, npts=50)`:**
        *   Initializes a k-point path using ASE's `get_bandpath` functionality.
        *   `lattice`: The reciprocal lattice vectors.
        *   `path`: A list of k-point labels (e.g., `['G', 'X', 'M', 'G']`) or explicit k-point coordinates defining the segments of the path.
        *   `special_points`: A dictionary mapping labels to k-point coordinates (e.g., `{'G': [0,0,0], 'X': [0.5,0,0]}`).
        *   `npts`: Number of points to generate along each segment of the path.
    *   Stores the generated k-points (`self.path`), corresponding x-axis distances for plotting (`self.xaxis`), and labels for high-symmetry points (`self.xlabel`, `self.symbols`).

*   **`CLASS BrillouinZone()`:**
    *   The main class for representing and operating on the Brillouin Zone of a crystal structure.
    *   **`__init__(self, astruct, mesh, evals, fermi_energy)`:**
        *   `astruct`: A dictionary containing the crystal structure information (cell vectors, atomic positions).
        *   `mesh`: The k-point mesh (e.g., `[nx, ny, nz]`) used in the DFT calculation.
        *   `evals`: A list of `BandArray` objects (or similar) containing eigenvalues for all k-points in the mesh.
        *   `fermi_energy`: The Fermi energy of the system.
        *   Uses `spglib` to get symmetry information and the irreducible reciprocal mesh.
        *   Uses `ase.dft.kpoints` to identify special points and paths for the given lattice type.
        *   Creates an interpolator (`scipy.interpolate.interpnd.LinearNDInterpolator`) for the band structure data across the BZ, handling data duplication for periodic boundary conditions.
    *   **`plot(self, path=None, npts=50)`:**
        *   Plots the interpolated band structure along a specified `BZPath` or a default path.
        *   Uses `matplotlib.pyplot` for plotting.
        *   Draws lines for the Fermi energy and high-symmetry k-points.

## Important Variables/Constants

*Descriptions of any critical variables or constants defined in the file that affect its behavior.*

*   **`BandArray.ikpt` (int):** Index of the k-point.
*   **`BandArray.kpt` (tuple/list of 3 floats):** Fractional coordinates of the k-point.
*   **`BandArray.kwgt` (float):** Weight of the k-point.
*   **`BZPath.path` (numpy.ndarray):** Array of k-point coordinates along the path.
*   **`BZPath.xaxis` (numpy.ndarray):** Accumulated distances along the k-point path, for x-axis of band structure plots.
*   **`BZPath.xlabel` (list):** Positions of high-symmetry points on the x-axis.
*   **`BZPath.symbols` (list):** Labels for the high-symmetry points (e.g., "$\Gamma$", "X").
*   **`BrillouinZone.lattice` (numpy.ndarray):** The reciprocal lattice vectors.
*   **`BrillouinZone.special_points` (dict):** Dictionary of high-symmetry k-points and their coordinates.
*   **`BrillouinZone.special_paths` (list of lists):** Defines standard paths through the BZ for common lattice types.
*   **`BrillouinZone.fermi_energy` (float):** The Fermi energy.
*   **`BrillouinZone.interpolator` (scipy.interpolate.interpnd.LinearNDInterpolator):** The object used to interpolate band energies throughout the BZ.

## Usage Examples

*If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.*

```python
from PyBigDFT.BigDFT import BZ as PyBZ # Assuming BZ.py is accessible this way
from PyBigDFT.IO import Logfiles # To read data from BigDFT outputs
import numpy

# 1. Load data from a BigDFT output (e.g., logfile)
# log = Logfiles(<path_to_logfile>)
# astruct = log.get_posinp() # Atomic structure
# evals_raw = log.evals # Eigenvalue data from log
# fermi_energy = log.fermi_energy

# --- Example data (replace with actual data loading) ---
astruct_example = {
    'Cell': [10.0, 10.0, 10.0], # Orthorhombic for simplicity
    'Positions': [{'Si': [0.0, 0.0, 0.0]}, {'Si': [0.5, 0.5, 0.5]}]
}
# Create dummy BandArray instances for evals (replace with actual parsed data)
dummy_kpt1_data = PyBZ.BandArray(data=[[-2.0, -1.0, 0.5, 1.0],[ -1.5, 0.0, 0.8, 1.2]],
                               ikpt=1, kpt=[0.0, 0.0, 0.0], kwgt=0.125)
dummy_kpt2_data = PyBZ.BandArray(data=[[-1.8, -0.9, 0.6, 1.1],[-1.4, 0.1, 0.9, 1.3]],
                               ikpt=2, kpt=[0.5, 0.0, 0.0], kwgt=0.25)
evals_example = [dummy_kpt1_data, dummy_kpt2_data] # Needs a full mesh for real BZ class
fermi_energy_example = 0.0
kpoint_mesh_example = [2, 1, 1] # Simplified mesh for example
# --- End of Example data ---

# 2. Create a BrillouinZone object (assuming you have the full k-point mesh data)
# For a real case, 'evals_example' would contain BandArray objects for all k-points in 'kpoint_mesh_example'
# bz = PyBZ.BrillouinZone(astruct_example, kpoint_mesh_example, evals_example, fermi_energy_example)

# 3. Define a k-point path using labels (if special points are known for the lattice)
#    (Using dummy special_points and path for this standalone example)
special_points_dummy = {'G': [0.,0.,0.], 'X': [0.5,0.,0.], 'M': [0.5,0.5,0.]}
path_segments_dummy = [['G','X'], ['X','M'], ['M','G']] # Corresponds to bz.special_paths[0] for cubic

# Create a BZPath object
# Note: BrillouinZone.lattice would be the reciprocal lattice. Here, using direct for simplicity of example.
direct_lattice_example = numpy.diag(astruct_example['Cell'])
k_path_obj = PyBZ.BZPath(direct_lattice_example, path_segments_dummy[0], special_points_dummy, npts=30)

# 4. Get the k-points along the path (these would be fed into a band structure calculation)
# kpoints_for_calc = k_path_obj.path
# print("K-points for band structure:", kpoints_for_calc)

# 5. Plot the band structure (after BZ object is properly initialized with full mesh data)
# bz.plot(path=k_path_obj)
# Or plot along a default path:
# bz.plot()

# Example of creating and using BandArray directly
ba = PyBZ.BandArray(data=[[-1.0, 0.0, 1.0], [-0.5, 0.5, 1.5]],
                  ikpt=1, kpt=[0.1, 0.2, 0.3], kwgt=0.5)
print("BandArray k-point:", ba.kpt)
print("BandArray energies (spin 1):", ba[0])
if ba.shape[0] > 1:
    print("BandArray energies (spin 2):", ba[1])
```

## Dependencies and Interactions

*   **NumPy:** Extensively used for numerical arrays and operations (e.g., `BandArray` inherits from `numpy.ndarray`).
*   **ASE (Atomistic Simulation Environment):** Used for DFT k-point utilities, specifically:
    *   `ase.dft.kpoints.get_special_points`: To get coordinates of high-symmetry k-points.
    *   `ase.dft.kpoints.parse_path_string`: To parse path definitions.
    *   `ase.dft.kpoints.get_bandpath`: To generate k-points along a specified path in the BZ.
*   **spglib:** Used for symmetry analysis of the crystal structure, specifically `spglib.get_spacegroup` and `spglib.get_ir_reciprocal_mesh` to find irreducible k-points.
*   **SciPy:** Used for interpolation, specifically `scipy.interpolate.interpnd.LinearNDInterpolator` to interpolate band energies across the Brillouin Zone.
*   **Matplotlib:** Used by the `BrillouinZone.plot()` method for visualizing band structures.
*   **PyBigDFT Modules:**
    *   Typically interacts with modules that parse BigDFT output files (like `Logfiles` from `PyBigDFT.IO`) to obtain atomic structure, eigenvalue data, and Fermi energy.
    *   The k-points generated by this module can be used to set up input for BigDFT calculations, particularly for band structure computations.
*   **BigDFT Calculations:**
    *   This module processes outputs from BigDFT (eigenvalues, structure).
    *   It helps in preparing inputs for specific types of BigDFT calculations (e.g., generating a list of k-points for a band structure run).
```
