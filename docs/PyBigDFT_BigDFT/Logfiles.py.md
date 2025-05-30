# Logfiles.py

## Overview

*This Python file, `Logfiles.py`, is a crucial part of the PyBigDFT suite, providing the `Logfile` class and associated helper functions to parse, analyze, and extract data from BigDFT's output log files. Since BigDFT typically produces detailed logs in YAML format, this module acts as the primary interface for accessing calculation results, structural information, convergence details, and other diagnostic data programmatically.*

## Key Components

*List of the primary functions, classes, or modules defined in the file, with a short description of each.*

*   **`CLASS Logfile()`:**
    *   **Purpose:** Represents one or more BigDFT log file documents. It parses the YAML content and provides convenient attributes and methods to access the stored data. It can handle single-point calculations, geometry optimizations (multiple documents), or lists of log data.
    *   **Key Methods:**
        *   **`__init__(self, *args, **kwargs)`:**
            *   Constructor that can load data from:
                *   Positional arguments (`*args`): One or more file paths to YAML log files.
                *   `archive` (kwarg): Path to a tar archive (.tgz) containing log files.
                *   `member` (kwarg): Specific member(s) to extract if `archive` is used or if YAML files contain multiple documents.
                *   `dictionary` (kwarg): A raw Python dictionary or a list of dictionaries representing parsed log content.
                *   `label` (kwarg): A label for the logfile instance.
            *   Initializes attributes by parsing the YAML data using predefined paths in the `BUILTIN` dictionary.
            *   If multiple documents/files are loaded (e.g., from a geometry optimization), it stores individual `Logfile` instances for each step in `self._instances` and often sets the main attributes to reflect the document with the lowest energy.
        *   **`__getitem__(self, index)`:** Allows accessing individual `Logfile` instances from a multi-document series using list-like indexing (e.g., `log[0]`).
        *   **`__str__(self)`:** Returns a string summarizing key information from the log file.
        *   **`__len__(self)`:** Returns the number of documents/steps if it's a series, otherwise 0.
        *   **`_initialize_class(self, d)`:** Internal method to parse a single YAML document (dictionary `d`) and set the object's attributes.
        *   **`_fermi_level_from_evals(self, evals)`:** Internal method to estimate the Fermi level if not explicitly found but eigenvalues are present.
        *   **`_sdos_line_to_orbitals(self, sorbs)`:** Internal helper for processing spatially decomposed DoS data.
        *   **`_get_bz(self, ev, kpts)`:** Internal method to process eigenvalues and k-points, likely creating `BandArray` objects (from `BZ.py`).
        *   **`get_dos(self, label=None, npts=2500)`:** Returns a `DoS` object (from `PyBigDFT.BigDFT.DoS`) initialized with eigenvalue data from the log file, ready for DoS calculation and plotting.
        *   **`get_brillouin_zone(self)`:** Returns a `BrillouinZone` object (from `PyBigDFT.BigDFT.BZ`) for band structure analysis, using k-point and eigenvalue data from the log.
        *   **`wfn_plot(self)`:** Plots the wavefunction convergence history (gnrm vs. iteration) using matplotlib.
        *   **`geopt_plot(self)`:** For a series of log files (geometry optimization), plots energy and maximum force vs. iteration number.
        *   **`_print_information(self)`:** Internal method to generate the string summary for `__str__`.

*   **`get_logs(files, select_document=None)` (Function):**
    *   Loads one or more YAML files specified in the `files` list.
    *   Uses `futile.Yaml.load` for robust loading of potentially multi-document YAML streams.
    *   `select_document`: Can be used to pick specific documents if a file contains multiple YAML entries.
    *   Returns a list of dictionaries, where each dictionary is a parsed YAML document.

*   **`floatify(scalar)` (Function):**
    *   Converts a string representation of a float (which might use 'd' or 'D' for exponents, common in Fortran output) to a Python float.

*   **`document_quantities(doc, to_extract)` (Function):**
    *   Extracts a predefined set of quantities (specified in `to_extract` dictionary, which maps desired names to paths within the YAML structure) from a single parsed YAML document (`doc`).
    *   Uses the `BUILTIN` dictionary for common predefined paths.

*   **`perform_operations(variables, ops, debug=False)` (Function):**
    *   Executes arbitrary Python code strings (`ops`) within a local namespace populated by `variables`. Potentially unsafe if `ops` comes from untrusted sources.

*   **`process_logfiles(files, instructions, debug=False)` (Function):**
    *   A more advanced processing utility that iterates through log files, extracts quantities using `document_quantities`, and then executes Python code snippets defined in the `instructions` dictionary (e.g., for evaluation or setup).

*   **`find_iterations(log)` (Function):**
    *   Navigates through a parsed log dictionary to find and extract the history of wavefunction convergence norms (gnrm) for each SCF iteration.

*   **`plot_wfn_convergence(wfn_it, gnrm_cv)` (Function):**
    *   Uses `matplotlib` to plot the wavefunction convergence data obtained from `find_iterations`.

*   **`BUILTIN` (Module-level dictionary):**
    *   A dictionary defining standard paths to access common quantities within the BigDFT YAML log structure. Keys are human-readable names (e.g., "energy", "fermi_level", "forces"), and values specify one or more potential list of keys (paths) to locate that data in the YAML.

## Important Variables/Constants

*Descriptions of any critical variables or constants defined in the file that affect its behavior.*

*   **`Logfile.log` (dict):** The primary dictionary holding the parsed YAML content of the current document.
*   **`Logfile.label` (str):** A label for the Logfile instance, often derived from the filename.
*   **`Logfile.srcdir` (str):** The source directory from which the logfile was loaded.
*   **`Logfile._instances` (list of Logfile):** If the input contained multiple YAML documents (e.g., a geometry optimization), this list holds `Logfile` objects for each individual document/step.
*   **`Logfile.reference_log` (int):** Index into `_instances` pointing to the log document used as the primary reference (often the one with the lowest energy).
*   **Attributes from `BUILTIN`:** The `_initialize_class` method dynamically creates attributes on the `Logfile` instance based on the keys in `BUILTIN` and the data found in the log. Common examples include:
    *   `Logfile.nat` (int): Number of atoms.
    *   `Logfile.energy` (float): Total energy (e.g., EKS or FKS).
    *   `Logfile.fermi_level` (float): Fermi energy.
    *   `Logfile.astruct` (dict): Atomic structure information.
    *   `Logfile.evals` (list of `BandArray` or similar): Eigenvalues.
    *   `Logfile.kpts` (list): K-point information.
    *   `Logfile.forces` (list/array): Atomic forces.
    *   `Logfile.forcemax` (float): Maximum force component.
    *   `Logfile.pressure` (float): Calculated pressure.
    *   `Logfile.sdos` (list): Data for spatially decomposed DoS.
    *   `Logfile.symmetry` (str): Space group symbol.

## Usage Examples

*Provide Python code snippets for common tasks:*

```python
from PyBigDFT.BigDFT import Logfiles
# Assuming other necessary PyBigDFT modules like DoS, BZ might be imported for specific methods

# 1. Load a single logfile
try:
    log_single = Logfiles.Logfile("path/to/single_calc_log.yaml")
    print(str(log_single)) # Print summary
    if hasattr(log_single, 'energy'):
        print(f"Energy: {log_single.energy} Ha")
    if hasattr(log_single, 'fermi_level'):
        print(f"Fermi Level: {log_single.fermi_level} Ha")

    # Get Density of States object
    # dos_obj = log_single.get_dos()
    # dos_obj.plot()

except ValueError as e:
    print(f"Error loading logfile: {e}")
except IOError:
    print("Logfile not found.")


# 2. Load multiple logfiles (e.g., from a geometry optimization)
#    or a single YAML file containing multiple documents.
try:
    log_series = Logfiles.Logfile("path/to/geopt_step1.yaml",
                                 "path/to/geopt_step2.yaml",
                                 "path/to/geopt_final.yaml",
                                 label="MyGeopt")
    # Or from a single multi-document file:
    # log_series = Logfiles.Logfile("path/to/all_geopt_steps.yaml", label="MyGeopt")

    print(f"Number of steps in series: {len(log_series)}")

    # Access data for the lowest energy structure (default for main attributes)
    print(f"Lowest energy in series: {log_series.energy} Ha")
    if hasattr(log_series, 'forcemax'):
        print(f"Max force at lowest energy point: {log_series.forcemax} Ha/Bohr")

    # Access data for a specific step (e.g., the first step, index 0)
    # if len(log_series) > 0:
    #     first_step_log = log_series[0]
    #     print(f"Energy of first step: {first_step_log.energy}")

    # Plot geometry optimization convergence
    # log_series.geopt_plot()

except ValueError as e:
    print(f"Error loading logfile series: {e}")
except IOError:
    print("One or more logfiles not found.")

# 3. Extract specific quantities using document_quantities (lower-level)
# raw_logs = Logfiles.get_logs(["path/to/single_calc_log.yaml"])
# if raw_logs:
#     my_quantities = {"total_energy": ["Last Iteration", "EKS"],
#                      "num_atoms": ["Atomic System Properties", "Number of atoms"]}
#     extracted_data = Logfiles.document_quantities(raw_logs[0], my_quantities)
#     print("Extracted Total Energy:", extracted_data.get("total_energy"))
#     print("Extracted Num Atoms:", extracted_data.get("num_atoms"))

```

## Dependencies and Interactions

*   **PyYAML:** Essential for parsing the YAML formatted log files. The code uses `yaml.load` and potentially `yaml.CLoader`.
*   **futile.Utils, futile.Yaml:** Used for utility functions like `write` and robust YAML loading (`Yaml.load`).
*   **NumPy:** Used for handling numerical data extracted from log files, especially for arrays of eigenvalues, forces, etc.
*   **Matplotlib:** Required for plotting methods like `wfn_plot()`, `geopt_plot()`, and implicitly by `DoS.plot()` and `BZ.plot()` called through this module.
*   **PyBigDFT Modules:**
    *   **`PyBigDFT.BigDFT.BZ`:** The `Logfile.evals` attribute is often processed into `BandArray` objects from `BZ.py`. The `get_brillouin_zone()` method returns a `BZ.BrillouinZone` object.
    *   **`PyBigDFT.BigDFT.DoS`:** The `get_dos()` method returns a `DoS.DoS` object.
*   **BigDFT Output:** This module is entirely dependent on the structure and content of BigDFT's YAML log files. Changes in BigDFT's logging format could potentially affect the parsing paths defined in `BUILTIN`.
*   **Role in Workflow:** `Logfiles.py` is a cornerstone for any post-processing or analysis workflow involving BigDFT. It provides the primary means to programmatically access detailed results from calculations, which can then be used for:
    *   Generating plots (DoS, band structures, convergence).
    *   Extracting specific numerical data for further analysis or reporting.
    *   Providing input data to other analysis tools or PyBigDFT modules.
    *   Checking the status and convergence of calculations.
```
