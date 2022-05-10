# Design Details

The interface makes use of the `iso_c_binding` introduced in Fortran 2003 standard. This allows interoperability between C and Fortran. Meanwhile, Cython can call C functions based on given C declarations. Hence, the Fortran subroutines can be called as a C function based on its declarations in Cython.

One of the biggest benefit of using Cython instead of f2py is that C declarations conform to C standard and allow us to avoid undefined behaviors.

## pyPDAF wrapper to Fortran call
In `pyPDAF.PDAF`, `.pxd` files define Fortran subroutines as external C functions. For example, the following Fortran subroutine interface:
```Fortran
subroutine c__PDAF_get_state(steps, time, doexit, &
                             U_next_observation, &
                             U_distribute_state, &
                             U_prepoststep, flag) bind(c)
      ! Flag and number of time steps
      INTEGER(c_int), INTENT(inout) :: steps
      ! current model time
      REAL(c_double), INTENT(out) :: time
      ! Whether to exit from forecasts
      INTEGER(c_int), INTENT(inout) :: doexit
      ! Status flag
      INTEGER(c_int), INTENT(inout) :: flag
      ! Provide time step and time of next observation
      procedure(c__next_observation_pdaf) :: U_next_observation 
      ! Routine to distribute a state vector
      procedure(c__distribute_state_pdaf) :: U_distribute_state 
      ! User supplied pre/poststep routine
      procedure(c__prepoststep_pdaf) :: U_prepoststep
end subroutine c__PDAF_get_state
```
can be translated into:
```C
cdef extern void c__pdaf_get_state (int* steps, double* time, int* doexit,
                                    void (*c__next_observation_pdaf)(int*,   int*,
                                                                     int*, double*),
                                    void (*c__distribute_state_pdaf)(int*, double*),
                                    void (*c__prepoststep_pdaf)(int*, int*, int*,
                                                                int*, int*, double*,
                                                                double*, double*, int*),
                                    int* flag);
```
Here, all arguments in the C declaration are pointers because Fortran is by default pass by reference, and user-defined procedures are provided with explicit interface in the C declaration. This is a good benefits compared to the f2py approach, where more tuning is required in the signature file to enable user-supplied external functions.

To provide a more Pythonic interface, the above C call is wrapped by a Python definition:
```Python
get_state (int steps, int doexit,
           py__next_observation_pdaf,
           py__distribute_state_pdaf,
           py__prepoststep_pdaf,
           int flag
          )
```
such that the user can have:
```Python
import pyPDAF.PDAF
PDAF.get_state(...)
```

## Treatment of arguments
There are input and output arguments in Fortran subroutines. In `pyPDAF`, input arguments with `intent(in)` or `intent(inout)` are declared in Python function arguments, and arguments with `intent(inout)` or `intent(out)` are treated as returned values. For example, `double* time` is not present in the `get_state` Python function argument as it is an `intent(out)` argument. They are defined as `cdef double time` before calling `c__pdaf_get_state`.

For some detailed explanation of the following sections, it is recommended to read Cython documentation for [ memory view](https://cython.readthedocs.io/en/stable/src/userguide/memoryviews.html), and [syntax](https://cython.readthedocs.io/en/stable/src/userguide/language_basics.html#differences-between-c-and-cython-expressions).

### Treatment of scalar arguments
In Cython, `int steps` declares the argument as integer value, which can be directly passed to C function calls by accesing its reference using `&`.

### Treatment of array arguments
In `pyPDAF`, all input arrays are numpy arrays. These arrays have to be convert to memoryview in Cython and then passed their reference to Fortran regardless of the shape of the array. To avoid problems with the different treatment of C ordering and Fortran ordering of multidimensional arrays, they're all raveled using Fortran order, e.g. `np.ravel(a, order='F')`.

The same applies to the returned array, where we reshape the memoryview back to its original shape using `np.reshape(a, order='F')`.

### Treatment of procedure arguments (user-supplied routines)
One of the most important aspects of PDAF is the user-supplied routines. This provides the flexibility to different models and observations. In `pyPDAF.UserFunc`, C functions are defined and declared with dummy Python functions. For example, 
```C
cdef void c__init_ens_pdaf (int* filtertype, int* dim_p,
                            int* dim_ens, double* state_p,
                            double* uinv, double* ens_p,
                            int* flag
                           );
``` 
the above C function calls the dummy Python function defined as:
```Python
def py__init_ens_pdaf(filtertype, dim_p, dim_ens, state_p, uinv, ens_p, flag):
    raise RuntimeError('...Wrong py__init_ens_pdaf is called!!!...')
```
Hence, `pyPDAF` can raise an error if required user-supplied function is not given. In `pyPDAF.PDAF`, the user-supplied function has to be given to `pyPDAF.UserFunc`. Taking the previous example for `c__PDAF_get_state`, we have the following treatment:
```Python
import pyPDAF.UserFunc as PDAFcython
cimport pyPDAF.UserFunc as c__PDAFcython
# giving user-supplied functions to UserFunc
PDAFcython.py__next_observation_pdaf = py__next_observation_pdaf
PDAFcython.py__distribute_state_pdaf = py__distribute_state_pdaf
PDAFcython.py__prepoststep_pdaf = py__prepoststep_pdaf
# call the actual functions
c__pdaf_get_state (&steps,
                   &time,
                   &doexit,
                   c__PDAFcython.c__next_observation_pdaf,
                   c__PDAFcython.c__distribute_state_pdaf,
                   c__PDAFcython.c__prepoststep_pdaf,
                   &flag
                  )
```
#### Treatment of array input in user-supplied routines
Because user-supplied routines are given as C pointers, we need to translate them into numpy array for Python function in a similar manner as [Treatment of array arguments](implementation.md#treatment-of-array-arguments), e.g. `state_p_np = np.asarray(<double[:np.prod((dim_p[0]))]> state_p).reshape((dim_p[0]), order='F')`.

##### Caveat
1. Fortran subroutines use pass by reference  by default
2. Python always pass/assignment by reference
The above two points can lead to undefined behavior. For example, in the user-supplied routines, we have Fortran interface:
```Fortran
subroutine c__collect_state_pdaf(dim_p, state_p) bind(c)
  use iso_c_binding, only: c_double, c_int
  implicit none
  ! pe-local state dimension
  integer(c_int), intent(in) :: dim_p
  ! local state vector
  real(c_double), intent(inout) :: state_p(dim_p)
end subroutine c__collect_state_pdaf
```
If we provide the following Python routine:
```Python
def collect_state_pdaf(model, assim_dim, dim_p, state_p):
    state_p = model.field_p.reshape(assim_dim.dim_state_p, order='F')
    return state_p
```
The Python function changes the reference of `state_p` as the `reshape` method generates a new numpy array instance. Hence, it is necessary to add a layer of safety in `pyPDAF.UserFunc`, where we pass the return value to the original input numpy array. We further safeguard this behavor with an assertion.
```Python
state_p_np_tmp = py__collect_state_pdaf(dim_p[0], state_p_np)
state_p_np[:] = state_p_np_tmp[:]
cdef double[::1] state_p_view = state_p_np.ravel(order='F')
assert state_p == &state_p_view[0], 'reference (memory address) of state_p has changed in c__collect_state_pdaf.'
```



## Note
In `pyPDAF/fortran`, `pyPDAF` provides `bind(C)` wrapper for common PDAF public subroutines with `iso_c_binding` arguments. Some subroutine uses assumed size arrays. The interoperability with C is recently supported in Fortran 2018 standard via `ISO_Fortran_binding.h` in C. To avoid unexpected issues for these new features, the dimension sizes are still given explicitly in the Fortran wrapper. 

The same issue applies to allocatable array in derived types in PDAFomi feature, where `obs_f` and `obs_l` are used. Hence, an array of the derived type are defined in the fortran wrapper with setters for many class members. These can be changed later when Fortran 2018 standard is better supported.