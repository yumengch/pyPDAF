pypdaf_sources = files('assim.pyx',
                       'put.pyx',
                     )

link_args = [
  '-Wl,-rpath,$ORIGIN/../',
]

foreach file :pypdaf_sources

  cython_ext = python.extension_module(
    fs.stem(file),
    file,
    link_with: [pdafc_lib],
    dependencies : [mpi, dep_py, blas_dep],
    include_directories: include_directories(incdir_numpy),
    c_args: ['-O3',],
    fortran_args: fortran_args,
    link_args: link_args,
    cython_args : ['-Xboundscheck=False', '-Xbinding=False', '-Xwraparound=False', '-Xinitializedcheck=False'],
    link_language: 'fortran',
    subdir: 'pyPDAF' / 'PDAF3',
    install: true,
    install_rpath: '$ORIGIN'
  )
endforeach