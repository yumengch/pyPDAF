from setuptools import setup, find_packages
# from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
from Cython.Compiler import Options

import numpy
import os

pwd = os.getcwd()

ext_modules = [Extension('*',
                         [f'{pwd}/pyPDAF/PDAF/*.pyx'],
                         extra_compile_args=['-g', '-fPIC'],
                         libraries=['PDAFc', "gfortran", 'm'],
                         library_dirs = [f'{pwd}/lib'],
                         ), 
                Extension('*',
                         [f'{pwd}/pyPDAF/Cython/*.pyx'])
                ]

setup(
    name = 'pyPDAF',
    version='0.0.1',
    ext_modules = cythonize(ext_modules),
    include_dirs=[numpy.get_include(), f'{pwd}/pyPDAF/PDAF/', f'{pwd}/pyPDAF/fortran/'],
    packages=find_packages()
)