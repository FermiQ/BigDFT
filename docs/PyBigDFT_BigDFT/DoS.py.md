# DoS.py

## Overview

*This Python file, `DoS.py`, from the PyBigDFT suite, is designed for calculating, processing, and visualizing the electronic Density of States (DoS). It provides classes to handle eigenvalue data, apply broadening (smearing), and plot the resulting DoS curves. It can also handle spatially decomposed DoS (sdos) for more detailed analysis.*

## Key Components

*List of the primary functions, classes, or modules defined in the file, with a short description of each.*

*   **`CLASS DiracSuperposition()`:**
    *   **Purpose:** Represents a raw density of states as a collection of Dirac delta functions located at the energy eigenvalues. It provides methods to generate a broadened curve (e.g., Gaussian or Lorentzian) from these discrete eigenvalues.
    *   **Key Methods:**
        *   **`__init__(self, dos, wgts=[1.0])`:** Initializes with an array of eigenvalues (`dos`) and corresponding weights (`wgts`). It determines the energy range (`self.xlim`).
        *   **`curve(self, xs, sigma, wgts=None)`:** Generates a broadened DoS curve over a given energy grid `xs` using a specified smearing width `sigma`. It can apply additional weights.
        *   **`peak(self, omega, e, sigma)`:** Calculates the value of a single broadened peak (currently Gaussian) centered at energy `e` evaluated at `omega` with width `sigma`.
        *   **`peaks(self, xs, dos, norms, sigma)`:** Sums up multiple broadened peaks to form the overall DoS curve.

*   **`CLASS DoS()`:**
    *   **Purpose:** The main class for managing and plotting one or more Density of States datasets. It can process eigenvalues from various sources, apply smearing, and generate interactive plots.
    *   **Key Methods:**
        *   **`__init__(self, bandarrays=None, energies=None, evals=None, units='eV', label='1', sigma=0.2/AU_eV, npts=2500, fermi_level=None, norm=1.0, sdos=None)`:**
            *   Initializes a DoS object. Can take input from:
                *   `bandarrays`: A list of `BandArray` objects (from `PyBigDFT.BigDFT.BZ`).
                *   `energies`: A raw NumPy array of energies.
                *   `evals`: A list of dictionaries (e.g., parsed from BigDFT logfiles).
            *   `units`: Units of the input energies ('eV' or 'AU').
            *   `sigma`: Default smearing width (in eV, converted from AU if necessary).
            *   `npts`: Number of points for the energy grid in plots.
            *   `fermi_level`: Value of the Fermi energy to be indicated on plots.
            *   `norm`: Normalization factor or weights for the DoS.
            *   `sdos`: Data for spatially decomposed DoS.
        *   **`append_from_bandarray(self, bandarrays, label)`:** Adds a DoS dataset from a list of `BandArray` objects, handling k-point weights and spin channels.
        *   **`append_from_dict(self, evals, label)`:** Adds a DoS dataset from a list of eigenvalue dictionaries.
        *   **`append(self, energies, label=None, units='eV', norm=1.0)`:** Adds a DoS dataset from a raw array of energies.
        *   **`conversion_factor(self, units)`:** Returns the conversion factor to eV from the given units.
        *   **`fermi_level(self, fermi_level, units='eV')`:** Sets the Fermi level, converting units if necessary.
        *   **`_set_range(self, npts=None)`:** Determines the energy range for plotting based on the min/max eigenvalues across all datasets.
        *   **`curve(self, dos, norm, sigma=None)`:** (Deprecated or internal style) Generates a broadened curve for a single dataset (functionality now mainly in `DiracSuperposition`).
        *   **`dump(self, sigma=None)`:** Prints the DoS data (energy vs. intensity) to standard output, suitable for piping to a file for gnuplot.
        *   **`plot(self, sigma=None, legend=False)`:** Generates an interactive plot of the DoS using `matplotlib`. Includes a slider to adjust the smearing width (`sigma`) dynamically. Also handles plotting of spatially decomposed DoS (sdos) if data is available.
        *   **`_embed_sdos(self, sdos)`:** Prepares sdos data for plotting.
        *   **`_set_sdos_selector(self)` / `_set_sdos_sliders(self, cmin, cmax)` / `_update_sdos(self, val)`:** Internal methods for managing interactive sdos plotting (coordinate selection, range sliders).
        *   **`update(self, val)`:** Callback function for the smearing slider to replot the DoS with the new sigma.
        *   **`get_ev(self, ev, keys=None)`:** Utility to extract eigenvalues from dictionary (similar to the module-level `get_ev`).

*   **`_bandarray_to_data(jspin, bandarrays)` (Internal helper function):**
    *   Extracts and organizes energy data and k-point weights from a list of `BandArray` objects for a specific spin channel (`jspin`).

## Important Variables/Constants

*Descriptions of any critical variables or constants defined in the file that affect its behavior.*

*   **`AU_eV` (float module-level constant):** Conversion factor from Hartrees (Atomic Units of energy) to electron-Volts (27.31138386 eV/Hartree).
*   **`DoS.ens` (list of dicts):** Stores the different DoS datasets. Each dictionary contains a `DiracSuperposition` object (`'dos'`).
*   **`DoS.labels` (list of str):** Labels for each DoS dataset, used in plot legends.
*   **`DoS.sigma` (float):** Current smearing width in eV.
*   **`DoS.npts` (int):** Number of points used for generating the energy grid for plots.
*   **`DoS.ef` (float):** Fermi level in eV.
*   **`DoS.range` (numpy.ndarray):** The common energy grid (x-axis) used for plotting all DoS curves.
*   **`DoS.sdos` (list):** Stores data for spatially decomposed DoS.
*   **`DiracSuperposition.dos` (numpy.ndarray):** Array of eigenvalues (energies of the Dirac deltas).
*   **`DiracSuperposition.norm` (list/array):** Weights associated with the eigenvalues or k-points.
*   **`DiracSuperposition.xlim` (tuple):** Minimum and maximum energy range for this specific set of deltas.

## Usage Examples

*Provide Python code snippets for common tasks:*

```python
from PyBigDFT.BigDFT import DoS, Logfiles, BZ # Assuming BZ.py for BandArray
import numpy as np

# --- Option 1: From a BigDFT Logfile ---
# log = Logfiles.Logfile("path/to/bigdft_log.yaml")
# evals_data = log.evals # This would be a list of BandArray objects if k-points are present
# fermi_e = log.fermi_energy
# dos_calc = DoS.DoS(bandarrays=evals_data, fermi_level=fermi_e, units='AU')

# --- Option 2: From raw energy eigenvalues (e.g., for a molecule) ---
# Assuming energies are already in eV
molecular_energies = np.array([-5.0, -4.5, -1.0, 0.5, 1.2, 1.5])
fermi_e_mol = 0.0 # Example Fermi level

dos_calc = DoS.DoS(energies=molecular_energies, fermi_level=fermi_e_mol, units='eV', label="MoleculeX")
dos_calc.sigma = 0.05 # Set smearing width in eV

# Add another set of energies for comparison
another_set_energies = molecular_energies * 0.9 - 0.5
dos_calc.append(another_set_energies, label="MoleculeY_Shifted", units='eV', norm=0.8)

# Get DoS data for the first dataset (index 0)
# This would involve calling the internal DiracSuperposition.curve method
# energies_grid = dos_calc.range
# dos_values_ds1 = dos_calc.ens[0]['dos'].curve(energies_grid, sigma=dos_calc.sigma)[1]

# Dump data to console (can be redirected to a file)
print("# Energy(eV) DoS_MoleculeX DoS_MoleculeY_Shifted")
dos_calc.dump(sigma=0.1) # Use a specific sigma for dumping

# Plot DoS interactively
dos_calc.plot(legend=True)

# --- Example with BandArray (conceptual if BZ.py is available) ---
# Assuming 'evals_from_bz' is a list of BZ.BandArray objects
# dos_from_bands = DoS.DoS(bandarrays=evals_from_bz, fermi_level=0.0, units='AU', label="SolidState")
# dos_from_bands.plot()
```

## Dependencies and Interactions

*   **NumPy:** Essential for all numerical operations, especially for handling arrays of energies and DoS values.
*   **Matplotlib:** Used for plotting the DoS curves and for the interactive smearing slider (`matplotlib.pyplot`, `matplotlib.widgets.Slider`, `matplotlib.widgets.RadioButtons`).
*   **`futile.Utils.write`:** Used for utility print statements.
*   **PyBigDFT Modules:**
    *   Can take `BandArray` objects from `PyBigDFT.BigDFT.BZ` as input, which is useful for k-point dependent DoS from solid-state calculations.
    *   Often used with `PyBigDFT.Logfiles` to extract eigenvalue data and Fermi levels directly from BigDFT output files.
*   **BigDFT Calculations:** This module primarily processes the output of BigDFT calculations, specifically the eigenvalues. It does not directly launch or control BigDFT runs.
*   **Workflow:** Typically used as a post-processing tool. After a BigDFT calculation finishes and eigenvalues are available (either from a logfile or other sources), `DoS.py` is used to generate, visualize, and analyze the density of states.
```
