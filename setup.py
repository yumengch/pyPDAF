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
from setuptools.dist import Distribution
from setuptools.command.develop import develop
from setuptools.command.install import install

from Cython.Build import cythonize

import configparser

import numpy
import glob
import os

pwd = os.getcwd()
# compiler
os.environ["CC"] = "gcc"
compilier_options = ['-fPIC']
# include directory
inc_dirs = [numpy.get_include(), f'{pwd}/pyPDAF/PDAF/']
# linking options
lib_dirs = ['lib']
# Cython set-up will automatically add -l as a prefix
# For example, 'PDAFc' becomes -lPDAFc in final compilation
libs = ['PDAFc']
extra_link_args = ['-Llib', '-Wl,-rpath=lib']
objs = []


def compile_interface():
    # Get our own instance of Distribution
    dist = Distribution()
    dist.parse_config_files()
    dist.parse_command_line()

    # Get prefix from either config file or command line
    PDAFdir = dist.get_option_dict('PDAF')['directory'][1]
    # fortran compiler options for
    FC = dist.get_option_dict('PDAF')['FC'][1]
    # compiler options
    OPT = dist.get_option_dict('PDAF')['OPT'][1]
    # include directory
    INC = dist.get_option_dict('PDAF')['INC'][1]
    LINK_LIBS = dist.get_option_dict('PDAF')['LINK_LIBS'][1]
    # CPP
    CPP_DEFS = dist.get_option_dict('PDAF')['CPP_DEFS'][1]

    PDAFsrc = os.path.join(PDAFdir, 'src')
    cmd = f'rm -f {PDAFsrc}/*.o {PDAFsrc}/*.mod {PDAFdir}/lib/libpdaf-d.a {PDAFdir}/include/*.mod'
    os.system(cmd)
    f90_files = glob.glob(os.path.join(PDAFdir, 'src', '*.F90'))
    f90_files += glob.glob(os.path.join('pyPDAF', 'fortran', '*.F90'))

    objs = []
    for src in f90_files:
        objs.append(f'{os.path.basename(src[:-4])}.o')
        cmd = f'{FC} {OPT} {CPP_DEFS} {INC} -c {src} -o {objs[-1]}'
        print(cmd)
        os.system(cmd)
    objs = ' '.join(objs)
    cmd = f'{FC} {objs} -shared {LINK_LIBS} -o lib/libPDAFc.so'
    print(cmd)
    os.makedirs(name, exist_ok=True)
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
                         [f'{pwd}/pyPDAF/UserFunc/*.pyx'])
               ]

setup(
    ext_modules=cythonize(ext_modules),
    include_dirs=inc_dirs,
    cmdclass={
        'develop': PreDevelopCommand,
        'install': PreInstallCommand
    },
)
