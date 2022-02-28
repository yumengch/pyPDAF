# Declaration of user-defined functions

Here, user-defined functions are declared. All user-defined function should have the same input and return arguments as defined here. 

## Note:
For multidimensional arrays, Cython assumes C-like (row-major) ordering even though the arrays are passed by Fortran (column-major ordering). Hence, transpose
is used for 2D arrays if each dimension has the same size, otherwise the dimension is swapped.

For example，the dimension of ens_p in fortran is `ens_p(dim_p, dim_ens)` as defined [here](https://github.com/yumengch/pyPDAF/blob/9848473114abc4e8ae5d7f22339f22f520e06102/pyPDAF/fortran/interface_pdaf.F90#L58)
Yet, in Cython, it is written as:
```python

if (dim_ens[0] != dim_p[0]):
    ens_p_numpy = np.asarray(
                    <double[:dim_ens[0], :dim_p[0]]> ens_p)
else:
    ens_p_numpy = np.asarray(
                    <double[:dim_ens[0], :dim_ens[0]]> ens_p).T
```
