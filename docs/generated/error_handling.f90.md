# `error_handling.f90`

## Overview

The `error_handling` module in the `GX-common` component provides a basic framework for registering and accessing error messages within the GreenX library. It allows different parts of the library to report errors in a standardized format, storing the most recent error message for later retrieval.

## Key Components

### Module Variable:

- **`error_message__`**:
    - Type: `character(len=err_len), public, protected`
    - Default: `"No Error reported so far"`
    - Description: This is the primary variable that stores the formatted error message. `err_len` is imported from the `constants` module (typically 1024).
    - The `public` attribute allows other modules to read this variable (e.g., via `gx_get_error_message` in `api_utilites`).
    - The `protected` attribute ensures that it can only be modified by procedures within the `error_handling` module itself (i.e., by `register_exc`).

### Public Subroutines:

- **`register_exc(msg, filename, line_number)`**:
    - **Purpose**: To register an error. It formats the provided information and updates the global `error_message__` string.
    - **Inputs**:
        - `msg`: Character string (`character(len=*)`), intent(in). The main error description.
        - `filename` (Optional): Character string (`character(len=*)`), intent(in). The name of the source file where the error occurred. Defaults to "Unknown File" if not provided.
        - `line_number` (Optional): Integer, intent(in). The line number in the source file where the error occurred. Defaults to 0 if not provided.
    - **Logic**:
        1. Initializes local `my_filename` to "Unknown File" and `line_num` to 0.
        2. If `filename` is present, `my_filename` is updated with its trimmed value.
        3. If `line_number` is present, `line_num` is updated.
        4. A formatted error string is constructed using `write` into `error_message__`. The format includes:
            ```
            --- ! ERROR
            src_file: <filename>
            src_line: <line_number>
            message: |
            <msg>
            ...
            ```
            Newline characters (`ch10` from `constants`) are used for structuring the message. The actual message `msg` and `filename` are trimmed of trailing spaces.

## Usage Examples

```fortran
module example_error_reporting
  use error_handling
  use constants, only: err_len
  implicit none

  character(len=err_len) :: current_error
  integer :: status

  subroutine do_something(input_val)
    integer, intent(in) :: input_val
    if (input_val < 0) then
      call register_exc("Input value cannot be negative.", &
                        filename=__FILE__, line_number=__LINE__)
      status = -1 ! Indicate an error
      return
    end if
    status = 0 ! Indicate success
    ! ... normal operations ...
  end subroutine do_something

  ! Elsewhere, to check and retrieve the error:
  call do_something(-5)
  if (status /= 0) then
     current_error = error_message__ ! Accessing the public variable
     print *, trim(current_error)
     ! Or, if using an API wrapper:
     ! call gx_get_error_message(current_error)
     ! print *, trim(current_error)
  end if

end module example_error_reporting
```
*(Note: `__FILE__` and `__LINE__` are preprocessor macros typically available if the Fortran source is preprocessed, often indicated by `.F90` extension.)*

## Dependencies and Interactions

- **`constants`**: This module uses `ch10` (newline character) and `err_len` (length of error string) from the `constants` module.
- **`api_utilites`**: The `gx_get_error_message` subroutine in `api_utilites` directly reads the `error_message__` variable from this module to provide error messages to the API user.
- **Macros (e.g., `_REGISTER_EXC`)**: Files like `gx_common.h` often define macros (e.g., `_REGISTER_EXC(msg)`) that wrap calls to `register_exc`, potentially automatically passing `__FILE__` and `__LINE__` if the compiler/preprocessor supports it. This simplifies error reporting in other parts of the codebase.

The `error_handling` module provides a simple, global error messaging system. Only the last error registered is stored. It's a common pattern for libraries where a detailed error string can be fetched after a function indicates failure via an error code.
