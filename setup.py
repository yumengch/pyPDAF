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
pwd = dist.get_option_dict('pyPDAF')['pwd'][1]
print ('pwd', pwd)
# compiler
os.environ["CC"] = dist.get_option_dict('pyPDAF')['CC'][1]
result = subprocess.run([os.environ["CC"], '--version'], stdout=subprocess.PIPE)
if result.stdout[:3] == b'icc':
    print ('....using Intel compiler....')
    os.environ["LDSHARED"] = "mpiicc -shared"
else:
    print ('....using GNU compiler....')

compilier_options = ['-fPIC', '-Wno-unreachable-code-fallthrough']
# include directory
inc_dirs = [numpy.get_include(), f'{pwd}/pyPDAF/PDAF/']
# linking options
lib_dirs = [f'{pwd}/lib']
# Cython set-up will automatically add -l as a prefix
# For example, 'PDAFc' becomes -lPDAFc in final compilation
libs = ['PDAFc']
extra_link_args = [f'-L{pwd}/lib',]
if sys.platform == 'darwin':
    pass
#    extra_link_args += [f'-rpath {pwd}/lib',]
else:
    extra_link_args += [f'-Wl,-rpath={pwd}/lib', ]
objs = []


def compile_interface():
    os.chdir(pwd)
    os.system('rm -rf pyPDAF.egg-info lib/*')

    # Get prefix from either config file or command line
    PDAFdir = dist.get_option_dict('PDAF')['directory'][1]
    if not os.path.isabs(PDAFdir):
        PDAFdir = os.path.join(pwd, PDAFdir)
        print ('input PDAF directory is not absolute path, changing to: ', PDAFdir)
    options = {}

    # Get compiler options
    for key in dist.get_option_dict('PDAF'):
        options[key] = dist.get_option_dict('PDAF')[key][1]


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
    status = os.system('sed \'s/$(FC) -O3 -o/$(FC) -O3 -fPIC -o/g\' Makefile > Makefile.tmp')
    status = os.system('mv -v Makefile.tmp  Makefile')
    status = os.system('rm -vf Makefile.tmp')
    status = os.system('make pdaf-var PDAF_ARCH=pyPDAF')
    if status:
        raise RuntimeError('failed to install PDAF')
    status = os.system('sed \'s/$(FC) -O3 -fPIC -o/$(FC) -O3 -o/g\' Makefile > Makefile.tmp')
    status = os.system('mv -v Makefile.tmp  Makefile')
    status = os.system('rm -vf Makefile.tmp')
    os.chdir(pwd)

    f90_files = ['pyPDAF/fortran/U_PDAF_interface_c_binding.F90',
                 'pyPDAF/fortran/PDAF_c_binding.F90',
                 'pyPDAF/fortran/PDAFomi_obs_c_binding.F90']
    objs = []
    for src in f90_files:
        objs.append(f'{os.path.basename(src[:-4])}.o')
        cmd = f'{options["FC"]} {options["OPT"]} {options["CPP_DEFS"]} {options["INC"]} -I. -c {src} -o {objs[-1]}'
        print(cmd)
        os.system(cmd)
    objs = ' '.join(objs)
    if sys.platform == 'darwin':
        cmd = f'{options["FC"]} {objs} -shared  -L{PDAFdir}/lib -lpdaf-var '\
              f'{options["LINK_LIBS"]} -o {pwd}/lib/libPDAFc.dylib'
    else:
        cmd = f'{options["FC"]} {objs} -shared -L{PDAFdir}/lib -lpdaf-var '\
              f'{options["LINK_LIBS"]} -o {pwd}/lib/libPDAFc.so'
    print(cmd)
    os.makedirs('lib', exist_ok=True)
    os.system(cmd)


class build_ext(build_ext_orig):
    """Pre-installation for pre build command.
        compiling PDAF iso_c_binding interface
    """
    def run(self):
        for ext in self.extensions:
            if ext.name == 'PDAFc':
                compile_interface()
        super().run()

ext_modules = [ Extension('PDAFc',
                          [f'{pwd}/pyPDAF/fortran/PDAFc.pyx']),
                Extension('*',
                         [f'{pwd}/pyPDAF/PDAF/*.pyx'],
                         extra_compile_args=compilier_options,
                         library_dirs=lib_dirs,
                         libraries=libs,
                         #extra_objects=objs,
                         extra_link_args=extra_link_args,
                         runtime_library_dirs=lib_dirs),
               Extension('*',
                         [f'{pwd}/pyPDAF/UserFunc/*.pyx'])
               ]

setup(name='pyPDAF',
    ext_modules=cythonize(ext_modules,
                          compiler_directives={'language_level': "3"}),
    include_dirs=inc_dirs,
    cmdclass={
        'build_ext': build_ext,
    },
)
