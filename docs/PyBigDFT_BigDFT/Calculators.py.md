# Calculators.py

## Overview

*This Python file, `Calculators.py`, provides classes for performing BigDFT calculations. It offers two distinct approaches: `GIBinding` which uses GObject Introspection to directly interface with the BigDFT library, and `SystemCalculator` which executes the `bigdft` program as an external system command. These calculators serve as a bridge between Python scripting and the core BigDFT functionalities.*

## Key Components

*List of the primary functions, classes, or modules defined in the file, with a short description of each.*

*   **`CLASS GIBinding()`:**
    *   **Purpose:** Provides a calculator interface to BigDFT by using its GObject Introspection bindings. This allows for more direct interaction with the BigDFT library from within Python, potentially offering finer control and avoiding the overhead of system calls.
    *   **Key Methods:**
        *   **`__init__(self)`:** Initializes the BigDFT library via its GObject Introspection interface and sets up the MPI environment if applicable.
        *   **`set(self, inputfile=None)`:** Initializes or resets the BigDFT run object (`self.runObj`) with a new set of input parameters provided as a dictionary (`inputfile`).
        *   **`update(self, inputfile)`:** Updates the current BigDFT run object with new or modified input parameters. This is useful for sequential calculations where only parts of the input change (e.g., atomic positions during a geometry scan).
        *   **`run(self)`:** Executes the BigDFT calculation using the current `runObj` and returns the output results (e.g., energies, forces stored in `self.out`).
        *   **`__del__(self)`:** Finalizes the BigDFT library when the `GIBinding` object is deleted.

*   **`CLASS SystemCalculator()`:**
    *   **Purpose:** Provides a calculator interface that runs the `bigdft` executable as an external command. This is a more traditional way of interfacing with compiled simulation codes.
    *   **Key Methods:**
        *   **`__init__(self, omp, mpi)`:** Initializes the calculator, storing the number of OMP threads (`omp`) and MPI processes (`mpi`) to be used. It also constructs the base command for running `bigdft`.
        *   **`run(self, name='', outdir='', run_name='', skip=False)`:** Executes the `bigdft` command as a system call.
            *   `name`: Corresponds to the `-n` option (job name).
            *   `outdir`: Corresponds to the `-d` option (output directory).
            *   `run_name`: Corresponds to the `-r` option (run name, often for specific input files).
            *   `skip`: Corresponds to the `-s Yes` option (to skip calculation if logfile indicates completion).
            It sets the `OMP_NUM_THREADS` environment variable before running.

## Important Variables/Constants

*Descriptions of any critical variables or constants defined in the file that affect its behavior.*

*   **`GIBinding.runObj`:** An instance of the BigDFT run object from the GObject Introspection bindings. This object holds all the input parameters and internal state for a calculation.
*   **`GIBinding.out`:** Stores the results returned by `runObj.calculate()`.
*   **`GIBinding.iproc`, `GIBinding.nproc`:** MPI process ID and total number of MPI processes.
*   **`SystemCalculator.omp` (str):** Number of OMP threads to use for the calculation.
*   **`SystemCalculator.mpi` (str):** Number of MPI processes to use for the calculation.
*   **`SystemCalculator.command` (str):** The base command string used to execute the `bigdft` executable (e.g., `mpirun -np <mpi_procs> $BIGDFT_ROOT/bigdft`).

## Usage Examples

*If applicable, provide examples or code snippets demonstrating how to use the functions, classes, or modules in the file.*

**Using `GIBinding`:**
```python
from PyBigDFT.BigDFT import Calculators # Assuming BZ.py is in the same directory for this relative import
import yaml

# Example input dictionary (can also be loaded from a YAML file)
input_dict = {
    "dft": {"ixc": "LDA", "itermax": 50},
    "posinp": {
        "positions": [{"H": [0.0, 0.0, 0.0]}, {"H": [0.0, 0.0, 0.74]}],
        "units": "angstroem"
    }
}

# Initialize the GIBinding calculator
calculator = Calculators.GIBinding()

# Set the input for the first calculation
calculator.set(input_dict)

# Run the calculation
results = calculator.run()
if calculator.iproc == 0: # Typically, only print/process results on the master process
    print("Energy:", results.eKS) # Access results (attributes depend on BigDFT bindings)

# Example of updating geometry and rerunning
new_posinp = {
    "posinp": {
        "positions": [{"H": [0.0, 0.0, 0.0]}, {"H": [0.0, 0.0, 0.80]}],
        "units": "angstroem"
    }
}
calculator.update(new_posinp)
results_updated = calculator.run()
if calculator.iproc == 0:
    print("Updated Energy:", results_updated.eKS)

# Calculator will automatically call BigDFT.lib_finalize() upon deletion
del calculator
```

**Using `SystemCalculator`:**
```python
from PyBigDFT.BigDFT import Calculators
from PyBigDFT.Inputfiles import Inputfile # For creating input files

# Create an Inputfile object and set parameters
inp = Inputfile()
inp.set_xc("PBE")
inp.add_posinp_xyz([[0.0, 0.0, 0.0, "O"], [0.757, 0.586, 0.0, "H"], [-0.757, 0.586, 0.0, "H"]])
inp.write_yaml("input_h2o.yaml") # Write the input to a file

# Initialize the SystemCalculator
# Assuming BIGDFT_ROOT is set in the environment
sys_calc = Calculators.SystemCalculator(omp=2, mpi=4)

# Run the calculation (BigDFT will read input_h2o.yaml if 'name' matches)
sys_calc.run(name="h2o", run_name="input_h2o")

# After the run, results would be in output files (e.g., log-h2o.yaml)
# and need to be parsed separately using PyBigDFT.Logfiles or similar.
```

## Dependencies and Interactions

*   **`GIBinding` Dependencies:**
    *   **GObject Introspection Bindings for BigDFT (`gi.repository.BigDFT`):** This is the core dependency, allowing direct calls to BigDFT library functions.
    *   **MPI:** Implicitly uses MPI through the BigDFT library bindings.
*   **`SystemCalculator` Dependencies:**
    *   **BigDFT Executable:** Requires a compiled `bigdft` executable located via `$BIGDFT_ROOT` environment variable.
    *   **MPI Executable (e.g., `mpirun`):** Used to launch the `bigdft` executable in parallel.
    *   **OS Environment:** Relies on `os.environ` for `BIGDFT_ROOT` and setting `OMP_NUM_THREADS`, and `os.system` for running commands.
    *   **`futile.Utils.write`:** Used for printing the command being executed.
*   **PyBigDFT Modules:**
    *   Can be used with `PyBigDFT.Inputfiles` to prepare input dictionaries or files.
    *   Results from calculations (especially from `SystemCalculator` which writes files, or `GIBinding` if results are saved) would typically be parsed using `PyBigDFT.Logfiles`.
*   **ASE (Atomic Simulation Environment):** While these specific calculator classes are **not** direct ASE-compatible calculators (they don't inherit from `ase.calculators.calculator.Calculator` or implement its full API like `get_potential_energy`, `get_forces`), they provide mechanisms to run BigDFT. The results could then be used to update ASE `Atoms` objects, or these classes could serve as a basis for creating a fully ASE-compliant calculator.
*   **Workflow:**
    *   `GIBinding` is suited for workflows where Python maintains more direct control over the BigDFT library, potentially for complex, iterative procedures or when fine-grained access to library functions is needed.
    *   `SystemCalculator` is suited for workflows where BigDFT is treated more like a black box, executed for a given input file, with results parsed from output files. This is a common pattern for many simulation codes.
```
