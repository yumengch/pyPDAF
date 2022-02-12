from setuptools import setup, Extension, find_packages
from setuptools.command.develop import develop
from setuptools.command.install import install

from Cython.Distutils import build_ext
from Cython.Build import cythonize
from Cython.Compiler import Options

import numpy
import glob
import os

pwd = os.getcwd()
# assuming the directory of PDAF directory
PDAFdir=f'{pwd}/../PDAF-D_V1.16'

# mpi linker requierments
mpi_lib_opnempi = ['mpi_usempif08', 'mpi_mpifh', 'mpi']
mpi_lib_mpich = ['mpifort', 'mpi']
# compiler 
os.environ["CC"] = "gcc-7"
compilier_options = ['-fPIC']
# include directory
inc_dirs = [numpy.get_include(), f'{pwd}/pyPDAF/PDAF/', f'{pwd}/build/']
# linking options
lib_dirs = [f'{PDAFdir}/lib', '/lib/x86_64-linux-gnu/']
libs = ['pdaf-d', ':libgfortran.so.4', 'm', 'lapack', 'blas'] + mpi_lib_opnempi
extra_link_args = []
objs = []
f90_files = glob.glob(os.path.join('pyPDAF', 'fortran', '*.F90'))
for f in f90_files:
    objs.append(f'build/{os.path.basename(f[:-4])}.o')


# fortran compiler options
FC = 'mpif90'
# compiler options
OPT = '-O3 -fdefault-real-8 -fPIC -Jbuild/'
# include directory
INC = f'-I/{PDAFdir}/include'
# CPP
CPP_DEFS = '-DUSE_PDAF'

class PreDevelopCommand(develop):
    """Pre-installation for development mode.
        compiling PDAF iso_c_binding interface
    """
    def run(self):
        os.makedirs('build', exist_ok=True)
        for src, obj in zip(f90_files, objs):
            cmd = f'{FC} {OPT} {CPP_DEFS} -I{PDAFdir}/include -c {src} -o {obj}'
            os.system(cmd)
        develop.run(self)


ext_modules = [Extension('*',
                         [f'{pwd}/pyPDAF/PDAF/*.pyx'],
                         extra_compile_args=compilier_options,
                         library_dirs = lib_dirs,
                         libraries=libs,
                         extra_objects=objs,
                         extra_link_args=extra_link_args
                         ), 
                Extension('*',
                         [f'{pwd}/pyPDAF/Cython/*.pyx'])
                ]

setup(
    name = 'pyPDAF',
    version='0.0.1',
    ext_modules = cythonize(ext_modules),
    include_dirs=inc_dirs,
    cmdclass={
        'develop': PreDevelopCommand
    },
)