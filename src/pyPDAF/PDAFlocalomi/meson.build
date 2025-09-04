pypdaf_sources = files('assim.pyx',
                       'put.pyx',
                     )

foreach file :pypdaf_sources

  cython_ext = python.extension_module(
    fs.stem(file),
    file,
    link_with: [pdafc_lib],
    dependencies : [mpi, dep_py, blas_dep, fruntime_dep],
    include_directories: include_directories(incdir_numpy),
    c_args: c_cython_args,
    link_args: link_args,
    fortran_args: fortran_args,
    cython_args : cython_args,
    link_language: 'fortran',
    subdir: 'pyPDAF' / 'PDAFlocalomi',
    install: true,
    install_rpath: '$ORIGIN'
  )
endforeach

python.install_sources(['__init__.py'], subdir: 'pyPDAF/PDAFlocalomi/')