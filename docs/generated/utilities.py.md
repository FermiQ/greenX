# `utilities.py`

## Overview

The `utilities.py` module, as part of the `pygreenx` Python package, provides helper functions that can be used across various Python scripts and tests within the GreenX library ecosystem. Currently, it contains a function for locating test binaries within the project structure.

## Key Components

### Functions:

- **`find_test_binary(root, binary_name: str) -> Path`**:
    - **Purpose**: To recursively search for a specific binary file starting from a given `root` directory and traversing its subdirectories. This is particularly useful for test scripts that need to locate compiled test executables.
    - **Parameters**:
        - `root`: A string or `pathlib.Path` object representing the root directory from which the search should begin.
        - `binary_name`: A string representing the name of the binary file to find (e.g., `"test_program.exe"` or `"my_test_binary"`).
    - **Returns**:
        - A `pathlib.Path` object representing the full path to the found binary file.
    - **Raises**:
        - `AssertionError`: If zero or more than one file matching `binary_name` is found within the search scope. The error message includes the list of found paths and the starting root directory.
    - **Logic**:
        1. Uses `Path(root).rglob(binary_name)` to perform a recursive glob search for any file or directory matching `binary_name`.
        2. Filters the results to include only actual files using `path.is_file()`.
        3. Asserts that the number of found files is exactly one.
        4. Returns the `Path` object of the uniquely identified binary.
    - **Note**: The function's docstring includes a performance consideration, suggesting that for a large number of tests or a very deep directory structure, this recursive search might become inefficient.

## Important Variables/Constants

(No module-level constants are defined in this file.)

## Usage Examples

```python
from pathlib import Path
from pygreenx.utilities import find_test_binary

# Assume the project structure:
# /project_root/
#   tests/
#     bin/
#       my_test_executable
#   scripts/
#     run_my_test.py

# In run_my_test.py:
project_root_dir = Path(__file__).parent.parent # Assuming script is in /project_root/scripts/
binary_to_find = "my_test_executable"
tests_bin_directory = project_root_dir / "tests" # More specific root for searching

try:
    executable_path = find_test_binary(tests_bin_directory, binary_to_find)
    print(f"Found binary at: {executable_path}")
    # Now executable_path can be used with BinaryRunner or subprocess
except AssertionError as e:
    print(f"Error locating binary: {e}")
except Exception as e:
    print(f"An unexpected error occurred: {e}")

# Example of what might cause an error:
# try:
#     # Searching from a broader root where duplicates might exist or file might be missing
#     find_test_binary(project_root_dir, "non_existent_binary")
# except AssertionError as e:
#     print(e)
```

## Dependencies and Interactions

- **Standard Python Libraries**:
    - `pathlib`: Used for all path manipulations and the recursive file search (`Path.rglob`, `Path.is_file`).
- This module is likely used by test automation scripts or other utility scripts within the `pygreenx` package that need to dynamically locate executables.
- It interacts with the file system to search for files.

The `find_test_binary` function provides a convenient way to locate necessary executables for testing or other scripted operations, reducing the need for hardcoded paths.
