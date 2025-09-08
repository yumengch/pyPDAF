# Developer Guide

The following guide explains the structure, implementation details, and mechanisms
used in `pyPDAF`. This guide is aimed at developers who wish to understand the
existing framework, make modifications, or contribute to its development.

---

## Overview

`pyPDAF` bridges Python and Fortran by leveraging the `Cython` library.
As Python is implemented in C, any interaction between Python and Fortran is
effectively handled as C-to-Fortran communication.

`Cython` automatically converts its module to `C` source code. The interoperability
with `Fortran` is achieved by the Fortran 2003 feature `iso_c_binding` module.

Contributions to the library can include raising issues, suggesting features,
or submitting pull requests with code enhancements.

---

## Adding a new Fortran function in pyPDAF

### Fortran Subroutines and Wrappers
The Fortran subroutine wrappers that are interoperable with C functions are
given in [src/fortran](src/fortran).

#### Interoperability with `bind(c)`
Fortran subroutines use the `bind(c)` keyword for compatibility with C.

Since PDAF does not use this keyword, `pyPDAF` provides its own wrapper subroutines that:
1. Subroutine names begin with the prefix `c__` to denote compatibility.
2. Arguments use corresponding C types.
3. User-supplied functions must be declared with [specified interface](src/fortran/pdaf_c_cb_interface.f90)
4. Bind(c) user-supplied functions must be converted to Fortran subroutines by
   pointers and wrapper subroutines in [src/fortran/pdaf_c_f_interface.f90](src/fortran/pdaf_c_f_interface.f90).
   This is a requirement in `flang`.
5. We do not use features that interoperable with derived types. This is because
   current standard does not support allocatable arrays in derived types.
   Therefore, `obs_f` and `obs_l` in PDAFomi are allocated in Fortran referenced
   by indices, and getting and setter functions.

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

---

### Cython Integration

#### Cython Declarations
To call Fortran subroutines in Python, `pyPDAF` uses Cython declarations
defined in `.pxd` files (e.g., `src/pyPDAF/PDAF.pxd`). Example:
```cython
cdef extern void c__pdaf_eofcovar (
    int* dim_state, int* nstates, int* nfields, int* dim_fields,
    int* offsets, int* remove_mstate, int* do_mv, double* states,
    double* stddev, double* svals, double* svec, double* meanstate,
    int* verbose, int* status) noexcept;
```

#### Python Wrappers
These Cython declarations are wrapped into Python functions to make them
accessible to users. Wrappers should return all `intent(out)` or `intent(inout)`
arguments in Python-friendly structures.

Example:
```cython
def eofcovar(int  dim, int  nstates, int  nfields, int [::1] dim_fields,
    int [::1] offsets, int  remove_mstate, int  do_mv,
    double [::1,:] states, double [::1] meanstate, int  verbose):
    """
    EOF analysis of an ensemble of state vectors by singular value decomposition.

    Typically, this function is used with :func:`pyPDAF.PDAF.SampleEns`
    to generate an ensemble of a chosen size (up to the number of EOFs plus one).

    Here, the function performs a singular value decomposition
    of the ensemble anomaly of the input matrix,
    which is usually an ensemble formed by state vectors
    at multiple time steps.
    The singular values and corresponding singular vectors
    can be used to construct a covariance matrix.
    This can be used as the initial error covariance for the initial ensemble.

    A multivariate scaling can be performed to ensure that all fields
    in the state vectors have unit variance.

    It can be useful to store more EOFs than one finally
    might want to use to have the flexibility
    to carry the ensemble size.


    See Also
    --------
    `PDAF webpage <https://pdaf.awi.de/trac/wiki/EnsembleGeneration>`_

    Parameters
    ----------
    dim : int
        Dimension of state vector
    nstates : int
        Number of state vectors
    nfields : int
        Number of fields in state vector
    dim_fields : ndarray[np.intc, ndim=1]
        Size of each field
        Array shape: (nfields)
    offsets : ndarray[np.intc, ndim=1]
        Start position of each field
        Array shape: (nfields)
    remove_mstate : int
        1: subtract mean state from states
    do_mv : int
        1: Do multivariate scaling; 0: no scaling
    states : ndarray[np.float64, ndim=2]
        State perturbations
        Array shape: (dim, nstates)
    meanstate : ndarray[np.float64, ndim=1]
        Mean state (only changed if remove_mstate=1)
        Array shape: (dim)
    verbose : int
        Verbosity flag

    Returns
    -------
    states : ndarray[np.float64, ndim=2]
        State perturbations
        Array shape: (dim, nstates)
    stddev : ndarray[np.float64, ndim=1]
        Standard deviation of field variability
        Array shape: (nfields)
    svals : ndarray[np.float64, ndim=1]
        Singular values divided by sqrt(nstates-1)
        Array shape: (nstates)
    svec : ndarray[np.float64, ndim=2]
        Singular vectors
        Array shape: (dim, nstates)
    meanstate : ndarray[np.float64, ndim=1]
        Mean state (only changed if remove_mstate=1)
        Array shape: (dim)
    status : int
        Status flag
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] states_np = np.asarray(states, dtype=np.float64, order="F")
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] stddev_np = np.zeros((nfields), dtype=np.float64, order="F")
    cdef double [::1] stddev = stddev_np
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] svals_np = np.zeros((nstates), dtype=np.float64, order="F")
    cdef double [::1] svals = svals_np
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] svec_np = np.zeros((dim, nstates), dtype=np.float64, order="F")
    cdef double [::1,:] svec = svec_np
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] meanstate_np = np.asarray(meanstate, dtype=np.float64, order="F")
    cdef int  status
    with nogil:
        c__pdaf_eofcovar(&dim, &nstates, &nfields, &dim_fields[0],
                         &offsets[0], &remove_mstate, &do_mv, &states[0,0],
                         &stddev[0], &svals[0], &svec[0,0], &meanstate[0],
                         &verbose, &status)

    return states_np, stddev_np, svals_np, svec_np, meanstate_np, status
```

---

#### Handling Callback Functions

Callback functions allow users to provide information for data assimilation.
However, Fortran expects these routines to follow specific interfaces.
These are handled in `src/pyPDAF/pdaf_c_cb_interface.pxd` and corresponding
``src/pyPDAF/pdaf_c_cb_interface.pyx`.

Example:
```cython
cdef void c__init_ens_pdaf(int* filtertype, int* dim_p, int* dim_ens,
    double* state_p, double* uinv, double* ens_p, int* flag) noexcept with gil:
    """Fill the ensemble array that is provided by PDAF with an initial ensemble of model states.

    This function is called by :func:`pyPDAF.PDAF.init`. The initialised
    ensemble array will be distributed to model by :func:`pyPDAF.PDAF.init_forecast`.

    Parameters
    ----------
    filtertype : int
            filter type given in PDAF_init
    dim_p : int
            PE-local state dimension given by PDAF_init
    dim_ens : int
            number of ensemble members
    state_p : ndarray[np.float64, ndim=1]
            PE-local model state
            This array must be filled with the initial
            state of the model for SEEK, but it is not
            used for ensemble-based filters.
            One can still make use of this array within
            this function.
            Array shape: (dim_p)
    uinv : ndarray[np.float64, ndim=2]
            This array is the inverse of matrix
            formed by right singular vectors of error
            covariance matrix of ensemble perturbations.
            This array has to be filled in SEEK, but it is
            not used for ensemble-based filters.
            Nevertheless, one can still make use of this
            array within this function e.g.,
            for generating an initial ensemble perturbation
            from a given covariance matrix.
            Dimension of this array is determined by the
            filter type.
            * (dim_ens, dim_ens) for (L)ETKF, (L)NETF, (L)KNETF, and SEEK
            * (dim_ens - 1, dim_ens - 1) for (L)SEIK, (L)ESTKF, and 3DVar using ensemble
            * (1, 1) for (L)EnKF, particle filters and gen_obs
            Array shape: (dim_ens - 1, dim_ens-1)
    ens_p : ndarray[np.float64, ndim=2]
            PE-local ensemble
            Array shape: (dim_p, dim_ens)
    flag : int
            pdaf status flag

    Returns
    -------
    state_p : ndarray[np.float64, ndim=1]
            PE-local model state
            This array must be filled with the initial
            state of the model for SEEK, but it is not
            used for ensemble-based filters.
            One can still make use of this array within
            this function.
            Array shape: (dim_p)
    uinv : ndarray[np.float64, ndim=2]
            This array is the inverse of matrix
            formed by right singular vectors of error
            covariance matrix of ensemble perturbations.
            This array has to be filled in SEEK, but it is
            not used for ensemble-based filters.
            Nevertheless, one can still make use of this
            array within this function e.g.,
            for generating an initial ensemble perturbation
            from a given covariance matrix.
            Dimension of this array is determined by the
            filter type.
            * (dim_ens, dim_ens) for (L)ETKF, (L)NETF, (L)KNETF, and SEEK
            * (dim_ens - 1, dim_ens - 1) for (L)SEIK, (L)ESTKF, and 3DVar using ensemble
            * (1, 1) for (L)EnKF, particle filters and gen_obs
            Array shape: (dim_ens - 1, dim_ens-1)
    ens_p : ndarray[np.float64, ndim=2]
            PE-local ensemble
            Array shape: (dim_p, dim_ens)
    flag : int
            pdaf status flag
    """
    cdef size_t uinv_len = max(dim_ens[0]-1, 1)
    cdef double[::1] state_p_np = np.asarray(<double[:dim_p[0]:1]> state_p, order="F")
    cdef double[::1,:] uinv_np = np.asarray(<double[:uinv_len:1,:uinv_len]> uinv, order="F")
    cdef double[::1,:] ens_p_np = np.asarray(<double[:dim_p[0]:1,:dim_ens[0]]> ens_p, order="F")

    state_p_np,uinv_np,ens_p_np,flag[0] = (<object>init_ens_pdaf)(
                                                                  filtertype[0],
                                                                  dim_p[0],
                                                                  dim_ens[0],
                                                                  state_p_np.base,
                                                                  uinv_np.base,
                                                                  ens_p_np.base,
                                                                  flag[0])

    cdef double[::1] state_p_new
    if state_p != &state_p_np[0]:
        state_p_new = np.asarray(<double[:dim_p[0]:1]> state_p, order="F")
        state_p_new[...] = state_p_np
        warnings.warn("The memory address of state_p is changed in c__add_obs_err_pdaf."
            "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)
    cdef double[::1,:] uinv_new
    if uinv != &uinv_np[0,0]:
        uinv_new = np.asarray(<double[:uinv_len:1,:uinv_len]> uinv, order="F")
        uinv_new[...] = uinv_np
        warnings.warn("The memory address of uinv is changed in c__add_obs_err_pdaf."
            "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)
    cdef double[::1,:] ens_p_new
    if ens_p != &ens_p_np[0,0]:
        ens_p_new = np.asarray(<double[:dim_p[0]:1,:dim_ens[0]]> ens_p, order="F")
        ens_p_new[...] = ens_p_np
        warnings.warn("The memory address of ens_p is changed in c__add_obs_err_pdaf."
            "The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)
```
where `init_ens_pdaf` is defined in `src/pyPDAF/pdaf_c_cb_interface.pxd` as a pointer:
`cdef void*  init_ens_pdaf = NULL;`.
The pointer is associated in pyPDAF functions, for example, in `src/pyPDAF/PDAF3/_pdaf3_c.pyx`:
```cython
from pyPDAF cimport pdaf_c_cb_interface as pdaf_cb
pdaf_cb.init_ens_pdaf = <void*>py__init_ens_pdaf
```
where `py__init_ens_pdaf` is the Python call-back function. The C function,
`pdaf_cb.c__init_ens_pdaf`, is
used for calling Fortran subroutines:
```cython
    with nogil:
        c__pdaf3_init(&filtertype, &subtype, &stepnull, &param_int[0],
                     &dim_pint, &param_real[0], &dim_preal,
                     pdaf_cb.c__init_ens_pdaf, &in_screen, &outflag)
```

#### Caveats
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

#### Exposing the function to pyPDAF and Mypy
To expose your function to pyPDAF or subpackages of pyPDAF, you need to import
it in `__init__.py` in corresponding directories. Further, as pyPDAF supports
type checking and other Python support features. You can add typing and docstring
to stub files ending with `.pyi`.
