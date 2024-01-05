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
from setuptools.command.build_ext import build_ext as build_ext_orig

from Cython.Build import cythonize

import configparser

import numpy
import glob
import sys
import os
import subprocess

# Get our own instance of Distribution
dist = Distribution()
dist.parse_config_files()
dist.parse_command_line()
# get the path to current directory
pwd = dist.get_option_dict('pyPDAF')['pwd'][1]
print ('pwd', pwd)
# get path to PDAF directory
PDAFdir = dist.get_option_dict('PDAF')['directory'][1]
if not os.path.isabs(PDAFdir):
    PDAFdir = os.path.join(pwd, PDAFdir)
    print ('input PDAF directory is not absolute path, changing to: ', PDAFdir)
# set up C compiler for cython and Python
os.environ["CC"] = dist.get_option_dict('pyPDAF')['CC'][1]
result = subprocess.run([os.environ["CC"], '--version'], stdout=subprocess.PIPE)
if result.stdout[:3] == b'icc':
    print ('....using Intel compiler....')
    os.environ["LDSHARED"] = "mpiicc -shared"
else:
    print ('....using GNU compiler....')
# compiler options for cython
extra_compile_args=['-Wno-unreachable-code-fallthrough']
# linking static PDAF library and interface objects
extra_objects=['-Wl,--whole-archive', f'{PDAFdir}/lib/libpdaf-var.a',
               f'lib/libPDAFc.a', '-Wl,--no-whole-archive']
# PDAF library contains multiple same .o files
# multiple-definition is thus necessary 
extra_link_args=['-Wl,--allow-multiple-definition']
# setup library to MPI-fortran 
LAPACK_PATH=dist.get_option_dict('pyPDAF')['LAPACK_PATH'][1]
print ('LAPACK_PATH', LAPACK_PATH)
library_dirs=['/usr/lib', ]
if LAPACK_PATH != '': library_dirs += LAPACK_PATH.split(',')
result = subprocess.run(['mpifort', '-show'], stdout=subprocess.PIPE)
result = result.stdout.decode()[:-1].split(' ')
s = [l[2:] for l in result if l[:2] == '-L']
if len(s) > 0: library_dirs += s
print ('library_dirs', library_dirs)
LAPACK_Flag=dist.get_option_dict('pyPDAF')['LAPACK_Flag'][1]
print ('LAPACK_Flag', LAPACK_Flag)
libraries=['gfortran', 'm', *LAPACK_Flag.split(',')]
s = [l[2:] for l in result if l[:2] == '-l']
if len(s) > 0: libraries += s
print ('libraries', libraries)

def compilePDAFLibraryInterface():
    """This function is used to compile PDAF library and its C interface
    """
    os.chdir(pwd)
    os.system('rm -rf pyPDAF.egg-info lib/*')

    options = {}
    # Get compiler options
    for key in dist.get_option_dict('PDAF'):
        options[key] = dist.get_option_dict('PDAF')[key][1]
    # generate configuration file for pyPDAF
    with open(f'{PDAFdir}/make.arch/pyPDAF.h', 'w') as the_file:
        for key in options:
            if key == 'directory':
                continue
            else:
                the_file.write(f'{key}={options[key]}\n')

    os.chdir(f'{PDAFdir}/src')
    status = os.system('make clean PDAF_ARCH=pyPDAF')
    if status:
        raise RuntimeError('failed to clean old PDAF installation')
    # modify the Makefile
    status = os.system('mv -v Makefile  Makefile.tmp')
    # remove _si files as they're not part of pyPDAF
    status = os.system('grep -v "_si.o" Makefile.tmp > Makefile')
    # this is for the .f files
    status = os.system('sed -i \'s/$(FC) -O3 -o/$(FC) -O3 -fPIC -o/g\' Makefile')
    # compile PDAF
    status = os.system('make pdaf-var PDAF_ARCH=pyPDAF')
    if status:
        raise RuntimeError('failed to install PDAF')
    # restore the original Makefile
    status = os.system('mv -v Makefile.tmp  Makefile')
    os.chdir(pwd)

    # compile the C interface to PDAF
    f90_files = ['pyPDAF/fortran/U_PDAF_interface_c_binding.F90',
                 'pyPDAF/fortran/PDAF_c_binding.F90',
                 'pyPDAF/fortran/PDAFomi_obs_c_binding.F90']
    # compile
    objs = []
    for src in f90_files:
        objs.append(f'{os.path.basename(src[:-4])}.o')
        cmd = f'{options["FC"]} {options["OPT"]} {options["CPP_DEFS"]} {options["INC"]} -I. -c {src} -o {objs[-1]}'
        print(cmd)
        os.system(cmd)
    objs = ' '.join(objs)
    # generate static library
    os.makedirs('lib', exist_ok=True)
    cmd = f'{options["AR"]} rc lib/libPDAFc.a {objs}'
    status = os.system(cmd)
    cmd = f'{options["RANLIB"]} lib/libPDAFc.a'
    status = os.system(cmd)


class build_ext(build_ext_orig):
    """Pre-installation for pre build command.
        compiling PDAF iso_c_binding interface
    """
    def run(self):
        for ext in self.extensions:
            if ext.name == 'PDAFc':
                compilePDAFLibraryInterface()
        super().run()


ext_modules = [Extension('PDAFc',
                          ['pyPDAF/fortran/PDAFc.pyx']),
               Extension('*',
                         ['pyPDAF/UserFunc.pyx']),
               Extension('*',
                         ['pyPDAF/PDAF.pyx'],
                         extra_compile_args=extra_compile_args,
                         extra_objects=extra_objects,
                         extra_link_args=extra_link_args,
                         library_dirs=library_dirs,
                         libraries=libraries
                         ),
               ]

setup(name='pyPDAF',
    ext_modules=cythonize(ext_modules,
                          compiler_directives={'language_level': "3"}),
    include_dirs= [numpy.get_include(),],
    cmdclass={
        'build_ext': build_ext,
    },
)
