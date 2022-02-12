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
PDAFdir='/home/yumengch/PDAF-D_V1.16'

class PreDevelopCommand(develop):
    """Pre-installation for development mode.
        compiling interface_pdaf.F90
    """
    def run(self):
        os.makedirs('build', exist_ok=True)
        # compile interface
        FC = 'mpif90'
        OPT = '-O3 -g -fdefault-real-8 -fPIC -Jbuild/'
        INC = f'-I/{PDAFdir}/include'
        CPP_DEFS = '-DUSE_PDAF'
        f90_files = glob.glob(os.path.join('pyPDAF', 'fortran', '*.F90'))
        for f in f90_files:
            cmd = f'{FC} {OPT} {CPP_DEFS} -I{PDAFdir}/include -c {f} -o build/{os.path.basename(f[:-4])}.o'
            print(cmd)
            os.system(cmd)
        develop.run(self)


mpi_lib_opnempi = ['mpi_usempif08', 'mpi_mpifh', 'mpi']
mpi_lib_mpich = ['mpifort', 'mpi']
os.environ["CC"] = "gcc-7"

ext_modules = [Extension('*',
                         [f'{pwd}/pyPDAF/PDAF/*.pyx'],
                         extra_compile_args=['-g', '-fPIC'],
                         # libraries=["gfortran", 'm'],
                         library_dirs = [f'{PDAFdir}/lib', '/lib/x86_64-linux-gnu/'],
                         # '/home/yumengch/spack/opt/spack/linux-ubuntu20.04-skylake/gcc-7.5.0/netlib-scalapack-2.1.0-y36qb4uycjwhcjz5d6demraj55iwthhk/lib',
                         # '/home/yumengch/spack/opt/spack/linux-ubuntu20.04-skylake/gcc-7.5.0/openblas-0.3.18-7irzwmgibpdjihbmqwhrwu4e3kctfemo'],
                         libraries=['pdaf-d', ':libgfortran.so.4', 'm', 'lapack', 'blas'] + mpi_lib_opnempi,
                         extra_objects=[f'build/interface_pdaf.o'],
                         extra_link_args=['-fexceptions', '-pthread']
                         ), 
                Extension('*',
                         [f'{pwd}/pyPDAF/Cython/*.pyx'])
                ]

setup(
    name = 'pyPDAF',
    version='0.0.1',
    ext_modules = cythonize(ext_modules),
    include_dirs=[numpy.get_include(), f'{pwd}/pyPDAF/PDAF/', f'{pwd}/pyPDAF/fortran/', f'{pwd}/build/'],
    cmdclass={
        'develop': PreDevelopCommand
    },
)