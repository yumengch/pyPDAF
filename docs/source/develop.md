# Developer Guide for pyPDAF

The following guide explains the structure, implementation details, and mechanisms used in `pyPDAF`, a Python interface to the PDAF (Parallel Data Assimilation Framework) Fortran library. This guide is aimed at developers who wish to understand the existing framework, make modifications, or contribute to its development. Please refer to (installation from source code)[install] if you'd like to further develop pyPDAF.

---

## Overview

`pyPDAF` bridges Python and Fortran by leveraging the `Cython` library for seamless interaction. As Python is implemented in C, any interaction between Python and Fortran is effectively handled as C-to-Fortran communication. The Fortran 2003 standard facilitates this through the `iso_c_binding` module, which ensures compatibility between Fortran and C datatypes.

`pyPDAF` relies on wrapper functions and automated code generation to handle interactions between Python and Fortran efficiently. Contributions to the library can include raising issues, suggesting features, or submitting pull requests with code enhancements.

---

## Implementation Details

### Source Code Structure

The core Python code of `pyPDAF` resides in the `src/pyPDAF` directory. To simplify the development process and avoid manually writing numerous subroutine wrappers, the library includes automated code generation scripts located in the `tools` directory:

- `write_PDAF_funcs.py`
- `write_user_funcs.py`
- `write_PDAF_pyi.py`

#### Purpose of the Scripts
- **Analyse Fortran Subroutines:** Parse subroutine names, arguments, and properties from Fortran files in `src/fortran`.
- **Extract Documentation:** Retrieve comments for each argument to generate API documentation.
- **Generate Wrapper Code:** Write Python-compatible code for `pyPDAF`.

> **Note:** These scripts are tailored specifically for `pyPDAF` and are not general-purpose Fortran-to-Cython parsers. Outputs must be carefully reviewed to ensure correctness.

---

### Fortran Subroutines and Wrappers

#### Interoperability with `bind(c)`
Fortran subroutines use the `bind(c)` keyword for compatibility with C. Since PDAF does not use this keyword, `pyPDAF` provides its own wrapper subroutines that:
1. Begin with the prefix `c__` to denote compatibility.
2. Include argument-level comments above the argument definition formatted using Sphinx syntax for API documentation.
3. Follow specific naming conventions and formatting rules.

#### Example Wrapper Subroutine
```fortran
subroutine c__PDAFomi_set_doassim(i_obs, doassim) bind(c)
   ! Index of observation types
   integer(c_int), intent(in) :: i_obs
   ! Flag to determine assimilation (0: no, 1: yes)
   integer(c_int), intent(in) :: doassim
   thisobs(i_obs)%doassim = doassim
end subroutine c__PDAFomi_set_doassim
```

#### Additional Features
- **Setter Subroutines:** For PDAF-OMI modules, setters for derived types (`obs_f`, `obs_l`) are implemented to handle complex structures.
- **Explicit Interfaces:** User-defined callback subroutines are declared in `src/Fortran/U_PDAF_interface_c_binding.F90` to ensure type safety and compatibility.

---

### Cython Integration

#### Cython Declarations
To call Fortran subroutines in Python, `pyPDAF` uses Cython declarations defined in `.pxd` files (e.g., `src/pyPDAF/PDAF.pxd`). Example:
```cython
cdef extern void c__pdaf_eofcovar (
    int* dim_state, int* nstates, int* nfields, int* dim_fields,
    int* offsets, int* remove_mstate, int* do_mv, double* states,
    double* stddev, double* svals, double* svec, double* meanstate,
    int* verbose, int* status) noexcept;
```

#### Python Wrappers
These Cython declarations are wrapped into Python functions to make them accessible to users. Wrappers simplify argument handling by:
1. Removing arguments that can be inferred (e.g., array dimensions).
2. Removing arguments that have `intent(out)` property.
3. Returning all `intent(out)` or `intent(inout)` arguments in Python-friendly structures.

Example:
```cython
def eofcovar(int[::1] dim_fields, int[::1] offsets, int remove_mstate, int do_mv, 
             double[:,:] states, double[::1] meanstate, int verbose):
    # Convert inputs to Fortran-compatible format
    cdef double[::1] states_f = np.asfortranarray(states).ravel(order="F")
    cdef int dim_state = states.shape[0]
    cdef int nstates = states.shape[1]
    cdef int nfields = dim_fields.shape[0]

    # Allocate output arrays
    cdef double[::1] stddev = np.zeros((nfields), dtype=np.float64)
    cdef double[::1] svals = np.zeros((nstates), dtype=np.float64)
    cdef double[::1] svec = np.zeros((dim_state, nstates), dtype=np.float64)
    cdef int status

    # Call the Fortran subroutine
    c__pdaf_eofcovar(&dim_state, &nstates, &nfields, &dim_fields[0], &offsets[0], 
                     &remove_mstate, &do_mv, &states_f[0], &stddev[0], &svals[0], 
                     &svec[0], &meanstate[0], &verbose, &status)

    # Return results in Python-friendly format
    return (
        np.asarray(states).reshape((dim_state, nstates), order='F'),
        np.asarray(stddev),
        np.asarray(svals),
        np.asarray(svec),
        np.asarray(meanstate),
        status
    )
```

---

### Handling Callback Functions

Callback functions allow users to define custom routines for data assimilation. However, Fortran expects these routines to follow specific interfaces. To resolve this, `pyPDAF` provides:
- **Fortran Interface Declarations:** Defined in `U_PDAF_interface_c_binding.F90`.
- **Python-to-C Bridge:** Implemented in `src/pyPDAF/UserFunc.pyx`.

Example:
```cython
state_p_np, uinv_np, ens_p_np, flag[0] = (<object>init_ens_pdaf)(
    filtertype[0], dim_p[0], dim_ens[0], state_p_np.base, uinv_np.base, 
    ens_p_np.base, flag[0]
)
```
where `init_ens_pdaf` is defined in `src/pyPDAF/UserFunc.pxd` as a pointer: `cdef void*  init_ens_pdaf = NULL;`. The pointer is associated in `src/pyPDAF/PDAF.pyx`:
```cython
from . cimport UserFunc as c__PDAFcython
c__PDAFcython.init_ens_pdaf = <void*>py__init_ens_pdaf
```
where `py__init_ens_pdaf` is the Python call-back function.

### Caveats
1. **Pass-by-Reference in Fortran:** Fortran passes arguments by reference, while Python’s behavior depends on the object type.
2. **Maintaining Consistency:** When Python functions modify arguments, ensure the original reference is preserved.

Example Safety Check:
```python
    cdef double[::1] state_p_new
    if state_p != &state_p_np[0]:
        state_p_new = np.asarray(<double[:dim_p[0]]> state_p)
        state_p_new[...] = state_p_np
        warnings.warn("The memory address of state_p is changed in c__init_ens_pdaf." 
         "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)
```

---

### Special Considerations for PDAF-OMI

- **Allocatable Arrays:** Explicit sizes are provided in wrappers to avoid interoperability issues, as Fortran 2018 support for assumed-size arrays is still evolving.
- **Derived Types:** Setter routines are implemented for fields in derived types like `obs_f` and `obs_l` to simplify usage.

---

This guide provides an in-depth understanding of `pyPDAF`’s structure, implementation details, and best practices for contributing. Developers are encouraged to verify all generated code and follow the documented conventions for consistency and reliability.