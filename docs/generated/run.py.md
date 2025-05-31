# `run.py`

## Overview

The `run.py` module, part of the `pygreenx` Python package, provides a utility class `BinaryRunner` for executing compiled binaries, particularly those associated with the GreenX library. It offers a structured way to run these executables with various configurations, including serial, OpenMP, MPI, and hybrid OpenMP/MPI parallelism. The module handles the construction of command-line arguments and environment setup for these different execution modes.

## Key Components

### `BuildType(enum.Enum)`
An enumeration that defines the type of build or execution mode for the binary.
- **Members**:
    - `serial`: For serial execution.
    - `pure_mpi`: For pure MPI parallel execution.
    - `pure_omp`: For pure OpenMP parallel execution.
    - `omp_mpi`: For hybrid OpenMP and MPI execution.

### `default_run_cmd` (Module variable)
A dictionary that maps `BuildType` enum members to default command prefixes for running the binary.
- Example: `BuildType.serial` maps to `'./'`, `BuildType.pure_mpi` maps to `'mpirun'`.

### `ProcessResults` Class
A simple data class to store and provide access to the results of an executed process.
- **Attributes**:
    - `stdout`: The standard output from the executed process.
    - `stderr`: The standard error output from the executed process.
    - `return_code`: The integer return code of the process.
    - `success`: Boolean, true if `return_code` is 0, false otherwise.

### `BinaryRunner` Class
The main class responsible for configuring and running an executable.
- **`__init__(self, binary: str, build_type: BuildType, omp_num_threads: Optional[int] = 1, mpi_processes: Optional[int] = 1, run_cmd: Optional[Union[str, dict]] = default_run_cmd, args: Optional[List[str]] = None)`**:
    - Constructor for the `BinaryRunner`.
    - **Parameters**:
        - `binary`: Path to the executable file.
        - `build_type`: A `BuildType` enum member indicating how the binary should be run.
        - `omp_num_threads`: (Optional) Number of OpenMP threads to use (default: 1).
        - `mpi_processes`: (Optional) Number of MPI processes to launch (default: 1).
        - `run_cmd`: (Optional) Either a string (e.g., 'mpirun') or a dictionary mapping `BuildType` to a command string. Defaults to `default_run_cmd`.
        - `args`: (Optional) A list of command-line arguments to pass to the binary.
    - Raises `FileNotFoundError` if the binary does not exist.
    - Asserts that `omp_num_threads` and `mpi_processes` are greater than 0.

- **`get_run_list(self) -> list`**:
    - Constructs the list of arguments that will be passed to `subprocess.run` to execute the binary.
    - For serial or pure OpenMP runs, it's typically `[binary_path] + args`.
    - For MPI or hybrid runs, it prepends the `run_cmd` (e.g., `mpirun`) and MPI process count arguments (e.g., `['-np', mpi_processes]`).
    - Raises `ValueError` if the `build_type` is not recognized.

- **`run(self, run_args: Optional[list] = None) -> ProcessResults`**:
    - Executes the configured binary.
    - **Parameters**:
        - `run_args`: (Optional) If provided, this list is used directly for `subprocess.run`. Otherwise, `get_run_list()` is called to generate the arguments.
    - Sets the `OMP_NUM_THREADS` environment variable for the subprocess according to `self.omp_num_threads`.
    - Uses `subprocess.run` to execute the command, capturing stdout, stderr, and the return code.
    - Returns a `ProcessResults` object containing the execution results.

## Important Variables/Constants
- `default_run_cmd`: Defines the standard way to launch executables based on their `BuildType`.

## Usage Examples

```python
import os
from pygreenx.run import BinaryRunner, BuildType, ProcessResults

# Assume 'my_program' is a compiled executable in the current directory
# and pygreenx is in the PYTHONPATH

# Example 1: Running a serial program
try:
    runner_serial = BinaryRunner(binary='./my_program', build_type=BuildType.serial, args=['--input', 'input.txt'])
    results_serial: ProcessResults = runner_serial.run()
    if results_serial.success:
        print("Serial program stdout:")
        print(results_serial.stdout.decode())
    else:
        print(f"Serial program failed with code {results_serial.return_code}:")
        print(results_serial.stderr.decode())
except FileNotFoundError as e:
    print(e)


# Example 2: Running an OpenMP program with 4 threads
try:
    # For OpenMP, OMP_NUM_THREADS is set by BinaryRunner
    runner_omp = BinaryRunner(binary='./my_program_omp',
                              build_type=BuildType.pure_omp,
                              omp_num_threads=4)
    results_omp: ProcessResults = runner_omp.run()
    if results_omp.success:
        print("OMP program stdout:")
        print(results_omp.stdout.decode())
except FileNotFoundError as e:
    print(e)


# Example 3: Running an MPI program with 2 processes
# (Requires mpirun to be in PATH and my_program_mpi to be compiled with MPI)
try:
    runner_mpi = BinaryRunner(binary='./my_program_mpi',
                              build_type=BuildType.pure_mpi,
                              mpi_processes=2)
    # The run_cmd for pure_mpi defaults to 'mpirun'
    results_mpi: ProcessResults = runner_mpi.run()
    if results_mpi.success:
        print("MPI program stdout:")
        print(results_mpi.stdout.decode())
    else:
        print(f"MPI program failed with code {results_mpi.return_code}:")
        print(results_mpi.stderr.decode())
except FileNotFoundError as e:
    print(e)
except ValueError as e: # If mpirun is not found by subprocess.run via default_run_cmd
    print(e)


# Example 4: Hybrid MPI/OpenMP
try:
    runner_hybrid = BinaryRunner(binary='./my_program_hybrid',
                                 build_type=BuildType.omp_mpi,
                                 mpi_processes=2,
                                 omp_num_threads=4,
                                 args=['--some-arg'])
    results_hybrid: ProcessResults = runner_hybrid.run()
    if results_hybrid.success:
        print("Hybrid program stdout:")
        print(results_hybrid.stdout.decode())
except FileNotFoundError as e:
    print(e)

```

## Dependencies and Interactions

- **Standard Python Libraries**:
    - `os`: Used for accessing environment variables (`os.environ`).
    - `pathlib`: Used for path manipulation and checking if the binary file exists (`pathlib.Path`).
    - `enum`: Used for creating the `BuildType` enumeration.
    - `subprocess`: Used for running the external executable (`subprocess.run`).
    - `typing`: Used for type hints (`Optional`, `Union`, `List`).
- This module is intended to be used as part of the `pygreenx` package to script and manage runs of GreenX executables.
- The actual executables (e.g., `./my_program`) are external to this Python script and are assumed to be compiled GreenX codes.
- For MPI execution, an MPI implementation (like OpenMPI, MPICH) and its `mpirun` (or equivalent) command must be available in the system's PATH.
```
