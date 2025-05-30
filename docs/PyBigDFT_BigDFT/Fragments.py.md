# Fragments.py

## Overview

*This Python file, `Fragments.py`, is part of the PyBigDFT suite and provides a collection of classes and functions for defining, manipulating, and analyzing molecular fragments and systems composed of such fragments. It includes tools for handling atomic coordinates, geometric transformations (rotations, translations), XYZ file I/O, and operations like finding centroids, dipoles, and decomposing a system into reference fragments. This module is particularly useful for setting up and analyzing calculations based on fragmentation methods or for studying subsystems within larger molecular arrangements.*

## Key Components

*List of the primary functions, classes, or modules defined in the file, with a short description of each.*

*   **`CLASS XYZfile()`:**
    *   **Purpose:** A utility class to facilitate the creation and writing of atomic coordinates in the XYZ file format.
    *   **Key Methods:**
        *   `__init__(self, filename=None, units='atomic')`: Initializes with an optional output filename and units for coordinates.
        *   `append(self, array, basename='', names=None, attributes=None)`: Adds a set of atomic positions to the internal list to be written.
        *   `dump(self, position='w')`: Writes the stored atomic coordinates to the specified file or standard output.

*   **Helper functions for XYZ files:**
    *   `open_xyz(filename, nat, unit, comment, position='a')`: Opens an XYZ file for writing.
    *   `close_xyz(f, filename)`: Closes an XYZ file.
    *   `dump_xyz_positions(f, array, basename='', names=None)`: Writes atomic positions to an already open file.
    *   `dump_xyz(array, basename='', units='atomic', names=None, filename=None, position='a')`: A convenience function to dump an array of positions to an XYZ file.

*   **`CLASS Lattice()`:**
    *   **Purpose:** Represents crystal lattice vectors and provides methods for generating points on the lattice.
    *   **Key Methods:**
        *   `__init__(self, vectors)`: Initializes with a list or array of lattice vectors.
        *   `grid(self, origin=[0.0,0.0,0.0], extremes=None, radius=None)`: Generates a set of translation vectors based on integer multiples of the lattice vectors, starting from an `origin`, within specified `extremes` or a `radius`.

*   **`CLASS RotoTranslation()`:**
    *   **Purpose:** Defines a combined rotation and translation transformation. It can determine the optimal transformation to align one set of points (`pos1`) onto another (`pos2`) using an external `wahba` module (presumably solving Wahba's problem).
    *   **Key Methods:**
        *   `__init__(self, pos1, pos2)`: Calculates the rotation matrix `self.R`, translation vector `self.t`, and a cost `self.J` for the transformation.
        *   `dot(self, pos)`: Applies the stored rototranslation to a given set of positions.
        *   `invert(self)`: Inverts the transformation.

*   **`CLASS Translation(RotoTranslation)`:**
    *   **Purpose:** A specialized `RotoTranslation` representing only a translation.
    *   `__init__(self, t)`: Initializes with a translation vector `t`.

*   **`CLASS Rotation(RotoTranslation)`:**
    *   **Purpose:** A specialized `RotoTranslation` representing only a rotation.
    *   `__init__(self, R)`: Initializes with a rotation matrix `R`.

*   **`CLASS Fragment()`:**
    *   **Purpose:** Represents a molecular fragment, storing a list of its constituent atoms, their coordinates, and an identifier.
    *   **Key Methods:**
        *   `__init__(self, atomlist=None, id='Unknown', units='AU')`: Initializes a fragment with a list of atoms (dictionaries), an ID, and units for input coordinates.
        *   `append(self, atom=None, sym=None, positions=None)`: Adds an atom to the fragment.
        *   `element(self, atom)`: Returns the chemical symbol of an atom in the fragment.
        *   `rxyz(self, atom)`: Returns the coordinates of an atom as a NumPy array.
        *   `__positions(self)`: Internal method to update the matrix of atomic positions.
        *   `centroid(self)`: Calculates the geometric center (centroid) of the fragment.
        *   `center_of_charge(self, zion)`: Calculates the center of charge, given atomic numbers `zion`.
        *   `transform(self, Rt)`: Applies a `RotoTranslation` object `Rt` to the fragment's atomic positions.
        *   `line_up(self)`: Aligns the fragment's principal axis of inertia with the coordinate axes and moves its centroid to the origin.
        *   `ellipsoid(self, center=0.0)`: Calculates the inertia-like tensor (ellipsoid of inertia).
        *   `q0(self, atom)`, `q1(self, atom)`: Retrieves atomic charges (q0) and dipoles (q1) if present in atom's dictionary.
        *   `Q(self, atom=None)`: Calculates the total charge of the fragment (or a specific element type within it) from atomic charges.
        *   `d0(self, center=None)`, `d1(self, center=None)`: Calculates the fragment's dipole moment based on atomic charges (`d0`) or atomic charges and atomic dipoles (`d1`).
        *   `xyz(self, filename=None, units='atomic')`: Writes the fragment's coordinates to an XYZ file.
        *   `dict(self)`: Converts fragment information to a list of dictionaries suitable for BigDFT's external potential format.

*   **`CLASS System()`:**
    *   **Purpose:** Represents a collection of `Fragment` objects, managing them as a whole system.
    *   **Key Methods:**
        *   `__init__(self, mp_dict=None, xyz=None, nat_reference=None, units='AU', transformations=None, reference_fragments=None)`: Initializes the system. Can populate from an XYZ file, a multipole dictionary (`mp_dict`), or by applying transformations to reference fragments. `nat_reference` helps in splitting concatenated XYZ files into multiple fragments.
        *   `fill_from_xyz(self, file, nat_reference)`: Populates fragments from an XYZ file.
        *   `fill_from_mp_dict(self, mpd, nat_reference)`: Populates fragments from a list of multipole dictionaries.
        *   `append(self, frag)`: Adds a `Fragment` object to the system.
        *   `pop(self, ifrag)`: Removes and returns a fragment from the system.
        *   `centroid(self)`: Calculates the centroid of the entire system (average of fragment centroids).
        *   `central_fragment(self)`: Finds the fragment closest to the system's centroid.
        *   `fragment_transformation(self, frag1, frag2)`: Determines the `RotoTranslation` to map `frag1` onto `frag2`.
        *   `decompose(self, reference_fragments)`: Decomposes the system's fragments by finding the best matching reference fragment and the corresponding transformation for each. Stores results in `self.decomposition`.
        *   `recompose(self, transformations=None, reference_fragments=None)`: Rebuilds the system's fragments based on a list of transformations and reference fragments.
        *   `Q(self)`: Calculates the total charge of the system.
        *   `xyz(self, filename=None, units='atomic')`: Writes all fragments in the system to a single XYZ file.
        *   `dict(self)`: Converts the entire system into a BigDFT dictionary format.

*   **`frag_average(ref, flist, clean_monopole=True)` (Function):**
    *   Averages atomic attributes (like multipole coefficients `q0`, `q1`, `q2`) from a list of fragments (`flist`) onto a reference fragment structure (`ref`). Can neutralize the resulting averaged fragment.

*   **`distance(i, j)` (Function):**
    *   Calculates the distance between the centroids of two `Fragment` objects.

## Important Variables/Constants

*Descriptions of any critical variables or constants defined in the file that affect its behavior.*

*   **`AU_to_A` (float module-level constant):** Conversion factor from Bohr (Atomic Units of length) to Angstroms.
*   **`Debye_to_AU` (float module-level constant):** Conversion factor from Debye to Atomic Units of dipole moment.
*   **`Fragment.atoms` (list of dicts):** List where each dictionary represents an atom (e.g., `{'C': [x,y,z], 'q0': [charge]}`).
*   **`Fragment.id` (str):** An identifier for the fragment.
*   **`Fragment.positions` (numpy.mat):** A NumPy matrix holding the Cartesian coordinates of atoms in the fragment (in AU).
*   **`Fragment.protected_keys` (list of str):** Keys in atom dictionaries that are not atom symbols (e.g., 'q0', 'q1', 'sigma').
*   **`System.fragments` (list of Fragment objects):** The list of fragments comprising the system.
*   **`System.CMs` (list of numpy.ndarray):** List of centroids of the fragments in the system.
*   **`System.decomposition` (list of dicts):** Stores the result of decomposing the system, where each dict contains the transformation (`'RT'`) and reference fragment (`'ref'`) for a system fragment.
*   **`RotoTranslation.R` (numpy.mat):** Rotation matrix.
*   **`RotoTranslation.t` (numpy.mat):** Translation vector.
*   **`RotoTranslation.J` (float):** Cost function value from Wahba's problem solution (measures goodness of fit for the transformation).

## Usage Examples

*Provide Python code snippets for common tasks:*

```python
from PyBigDFT.BigDFT import Fragments
import numpy as np

# 1. Define a single fragment
atom_list_water = [
    {'O': [0.000000, 0.000000, 0.117300]},
    {'H': [0.000000, 0.757200, -0.469200]},
    {'H': [0.000000, -0.757200, -0.469200]}
]
water_frag = Fragments.Fragment(atom_list_water, id="Water1", units="A")
print("Water fragment centroid (AU):", water_frag.centroid())
water_frag.xyz("water.xyz", units="angstroem")

# 2. Create a system of fragments
# Assume another fragment, water_frag2, is defined similarly or loaded
# water_frag2 = ...
system = Fragments.System(units="A")
system.append(water_frag)
# system.append(water_frag2) # If water_frag2 exists

# Get the system's centroid
# print("System centroid (AU):", system.centroid())

# Write the whole system to an XYZ file
system.xyz("system_of_fragments.xyz", units="angstroem")

# 3. Load a system from an XYZ file (e.g., two water molecules)
# Assuming "two_waters.xyz" contains 6 atoms, 3 per water molecule
# system_from_file = Fragments.System(xyz="two_waters.xyz", nat_reference=3, units="A")
# print(f"Loaded {len(system_from_file.fragments)} fragments.")
# if len(system_from_file.fragments) == 2:
# print("Distance between fragments:", Fragments.distance(system_from_file.fragments[0], system_from_file.fragments[1]))

# 4. Fragment transformation (conceptual)
# frag_a = system_from_file.fragments[0]
# frag_b = system_from_file.fragments[1]
# rt_ab = Fragments.RotoTranslation(frag_a.positions, frag_b.positions) # Find transform from A to B
# print("Transformation cost J:", rt_ab.J)
# transformed_a_positions = rt_ab.dot(frag_a.positions)
# print("RMSD after transformation:", np.sqrt(np.mean(np.square(transformed_a_positions - frag_b.positions))))

# 5. Decompose and Recompose (conceptual)
# reference_water = Fragments.Fragment(atom_list_water, units="A") # A single reference
# system_to_decompose = Fragments.System(xyz="multiple_waters.xyz", nat_reference=3, units="A")
# system_to_decompose.decompose(reference_fragments=[reference_water])
# for i, item in enumerate(system_to_decompose.decomposition):
#     print(f"Fragment {i}: best ref id {item['id']}, transform cost J {item['RT'].J}")
# system_to_decompose.recompose() # Rebuilds fragments from decomposition data
# system_to_decompose.xyz("recomposed_system.xyz")

```

## Dependencies and Interactions

*   **NumPy:** Heavily used for all numerical operations, especially for handling atomic coordinates as matrices and vectors.
*   **PyYAML (implicitly for `__str__` of Fragment, though not explicitly imported at top level):** Used if string representation of a Fragment were to dump YAML (the current `__str__` is commented out).
*   **`wahba` module (external, not part of PyBigDFT):** Required by the `RotoTranslation` class to solve Wahba's problem for finding optimal rigid transformations. If not available, `RotoTranslation` initialization will fail.
*   **PyBigDFT Modules:**
    *   The dictionaries produced by `Fragment.dict()` and `System.dict()` are formatted for use as `external_potential` or `posinp` in `PyBigDFT.Inputfiles.Inputfile`.
    *   Could be used to prepare systems for QM/MM calculations or fragment-based DFT methods (e.g., SIESTA-style fragmentation, ONETEP-style NGSCF optimization if such interfaces exist in PyBigDFT).
*   **BigDFT Calculations:**
    *   This module is primarily for pre-processing (setting up fragmented systems) or post-processing (analyzing fragment properties from a larger calculation's outputs, like atomic charges/multipoles).
    *   The fragment definitions can be used to construct specialized BigDFT inputs, for example, by defining different regions for QM or MM treatment, or by defining charge/spin constraints on specific fragments.
```
