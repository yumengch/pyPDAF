"""This file is part of pyPDAF

Copyright (C) 2022 University of Reading and
National Centre for Earth Observation

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
from setuptools import setup, Extension
from setuptools.command.develop import develop
from setuptools.command.install import install

from Cython.Build import cythonize

import numpy
import glob
import os

pwd = os.getcwd()
# assuming the directory of PDAF directory
PDAFdir = f'{pwd}/../PDAF-D_V1.16'

# compiler
os.environ["CC"] = "gcc-7"
compilier_options = ['-fPIC']
# include directory
inc_dirs = [numpy.get_include(), f'{pwd}/pyPDAF/PDAF/']
# linking options
lib_dirs = ['lib']
# Cython set-up will automatically add -l as a prefix
# For example, 'PDAFc' becomes -lPDAFc in final compilation
libs = ['PDAFc']
extra_link_args = []
objs = []


def compile_interface():
    # fortran compiler options for
    FC = 'mpif90'
    # compiler options
    OPT = '-O3 -fdefault-real-8 -fPIC'
    # include directory
    INC = f'-I/{PDAFdir}/include'
    LINK_LIBS = '-llapack -lblas'
    # CPP
    CPP_DEFS = '-DUSE_PDAF'

    f90_files = glob.glob(os.path.join('pyPDAF', 'fortran', '*.F90'))
    objs = []
    for src in f90_files:
        objs.append(f'{os.path.basename(src[:-4])}.o')
        cmd = f'{FC} {OPT} {CPP_DEFS} {INC} -c {src} -o {objs[-1]}'
        print(cmd)
        os.system(cmd)
    objs = ' '.join(objs)
    cmd = f'{FC} {objs} -shared -L{PDAFdir}/lib -lpdaf-d' \
          f' {LINK_LIBS} -o lib/libPDAFc.so'
    print(cmd)
    os.system(cmd)


class PreDevelopCommand(develop):
    """Pre-installation for development mode.
        compiling PDAF iso_c_binding interface
    """
    def run(self):
        compile_interface()
        develop.run(self)


class PreInstallCommand(install):
    """Pre-installation for development mode.
        compiling PDAF iso_c_binding interface
    """
    def run(self):
        compile_interface()
        install.run(self)


ext_modules = [Extension('*',
                         [f'{pwd}/pyPDAF/PDAF/*.pyx'],
                         extra_compile_args=compilier_options,
                         library_dirs=lib_dirs,
                         libraries=libs,
                         extra_objects=objs,
                         extra_link_args=extra_link_args),
               Extension('*',
                         [f'{pwd}/pyPDAF/Cython/*.pyx'])
               ]

setup(
    name='pyPDAF',
    version='0.0.1',
    ext_modules=cythonize(ext_modules),
    include_dirs=inc_dirs,
    cmdclass={
        'develop': PreDevelopCommand,
        'install': PreInstallCommand
    },
)
