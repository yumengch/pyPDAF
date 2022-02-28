# Python wrapper for C-type function calls

The wrapper convert Python object to memoryview for numpy arrays (see [here](https://cython.readthedocs.io/en/latest/src/userguide/memoryviews.html) 
for an introduction to memory view). 

For PDAFomi, setter functions are used since there are difficulties in interoperability for iso_c_binding in fortran with derived types using allocatable arrays.
