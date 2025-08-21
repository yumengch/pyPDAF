
pypdaf_sources = files('assim.pyx',
                       'callback.pyx',
                       'diag.pyx',
                       'get.pyx',
                       'iau_internal.pyx',
                       'iau.pyx',
                       'internal.pyx',
                       '_pdaf_c.pyx',
                       'put.pyx',
                       'setter.pyx',
                     )

link_args = [
  '-Wl,-rpath,$ORIGIN/../',
]
foreach file :pypdaf_sources
  cython_ext = python.extension_module(
    fs.stem(file),
    file,
    link_with: pdafc_lib,
    dependencies : [mpi, dep_py, blas_dep, fruntime_dep],
    include_directories: include_directories(incdir_numpy),
    c_args: c_cython_args,
    fortran_args: fortran_args,
    link_args: link_args,
    cython_args : cython_args,
    link_language: 'fortran',
    subdir: 'pyPDAF' / 'PDAF' ,
    install: true,
  )
endforeach