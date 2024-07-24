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
from setuptools import Distribution
from setuptools.command.build_ext import build_ext as build_ext_orig

from Cython.Build import cythonize

import logging

import numpy
import sys
import os
import subprocess
import shutil

# logging 
logging.getLogger().setLevel(logging.DEBUG)

# Get our own instance of Distribution
dist = Distribution()
dist.parse_config_files()
dist.parse_command_line()
# get the path to current directory
cwd = os.getcwd()
pwd = dist.get_option_dict('pyPDAF')['pwd'][1]
logging.info(f'pwd: {pwd}; getcwd: {os.getcwd()}')
# get path to PDAF directory
PDAFdir = dist.get_option_dict('pyPDAF')['PDAF_dir'][1]
if not os.path.isabs(PDAFdir):
    PDAFdir = os.path.join(pwd, PDAFdir)
    logging.info (f'....input PDAF directory is not absolute path, changing to: {PDAFdir}....')

condaBuild = dist.get_option_dict('pyPDAF')['condaBuild'][1]
cmake_config_path = dist.get_option_dict('pyPDAF')['cmake_config_path'][1]
logging.info (f'....condaBuild: {condaBuild}....')
logging.info (f'....cmake_config_path: {cmake_config_path}....')

# set up C compiler for cython and Python
c_compiler = dist.get_option_dict('pyPDAF')['c_compiler'][1]
assert c_compiler in ['gcc', 'msvc', 'icc', 'clang'], f'{c_compiler} is not a supported C compiler,'\
                                                        ' check c_compiler option in setup.cfg'
logging.info (f'....using {c_compiler} compiler for C language....')
fortran_compiler = dist.get_option_dict('pyPDAF')['fortran_compiler'][1]
assert fortran_compiler in ['gfortran', 'ifort'], f'{fortran_compiler} is not a supported C compiler,'\
                                                   ' check fortran_compiler option in setup.cfg'
logging.info (f'....using {fortran_compiler} compiler for C language....')
if c_compiler == 'icc': os.environ["LDSHARED"] = "mpiicc -shared"

extra_compile_args=[]
extra_link_args = []
extra_objects = []
library_dirs=[]
libraries = []

# linking static PDAF library and interface objects
if sys.platform == 'darwin':
    # Linking static library in mac
    extra_objects+=['-Wl,-force_load', f'{PDAFdir}/lib/libpdaf-var.a', '-Wl,-force_load', f'{pwd}/lib/libPDAFc.a',]
elif os.name == 'nt':
    # linking static library in windows
    library_dirs+=[os.path.join(PDAFdir, 'lib'), os.path.join(pwd, 'lib'),]
    libraries += ['pdaf-var', 'pdafc']
else:
    # linking static library in linux
    extra_objects+=['-Wl,--whole-archive', f'{PDAFdir}/lib/libpdaf-var.a', f'{pwd}/lib/libPDAFc.a', '-Wl,--no-whole-archive']

# add mpi library path
if os.name == 'nt':
    # always use external msmpi as msmpi from conda cannot be linked
    # there seems to be no mpi compiler wrapper in windows
    MPI_LIB_PATH=dist.get_option_dict('pyPDAF')['MPI_LIB_PATH'][1]
    if MPI_LIB_PATH != '': library_dirs += MPI_LIB_PATH.split(',')
    libraries += ['msmpi', 'msmpifec']
else:
    # using mpi compiler wrapper is easier for linux and mac
    # we do not consider cray etc. at the moment
    mpifortran = 'mpifort' if fortran_compiler == 'gfortran' else 'mpiifort'
    result = subprocess.run([mpifortran, '-show'], stdout=subprocess.PIPE)
    result = result.stdout.decode()[:-1].split(' ')
    s = [l[2:].replace('"', '') for l in result if l[:2] == '-L']
    if len(s) > 0: library_dirs += s
    s = [l[2:] for l in result if l[:2] == '-l']
    if len(s) > 0: libraries += s

# linking BLAS/LAPACK/MKL
use_MKL=dist.get_option_dict('pyPDAF')['use_MKL'][1]
if use_MKL == 'True':
    # if condaBuild == 'True':
    #     MKLROOT = os.environ['LIBRARY_LIB'] if os.name == 'nt' else \
    #               os.path.join(os.environ['PREFIX'], 'lib')
    # else:
    MKLROOT = dist.get_option_dict('pyPDAF')['MKLROOT'][1]
    # assert MKLROOT != '', 'MKLROOT must not be empty, check setup.cfg file'

    if os.name == 'nt':
        if MKLROOT != '': library_dirs+=[MKLROOT,]
        libraries += ['mkl_core', 'mkl_sequential', 'mkl_intel_lp64']
    elif sys.platform == "linux" or sys.platform == "linux2":
        if condaBuild == 'True': MKLROOT = os.path.join(os.environ['PREFIX'], 'lib')
        extra_objects+=['-Wl,--start-group', 
                        f'{MKLROOT}/libmkl_intel_lp64.a',
                        f'{MKLROOT}/libmkl_sequential.a',
                        f'{MKLROOT}/libmkl_core.a',
                        '-Wl,--end-group']
    else:
        if condaBuild == 'True': MKLROOT = os.path.join(os.environ['PREFIX'], 'lib')
        extra_objects+=[ 
                        f'{MKLROOT}/libmkl_intel_lp64.a',
                        f'{MKLROOT}/libmkl_sequential.a',
                        f'{MKLROOT}/libmkl_core.a']
else:
    # setup library to LAPACKA and BLAS
    LAPACK_PATH=dist.get_option_dict('pyPDAF')['LAPACK_PATH'][1]
    if LAPACK_PATH != '': library_dirs += LAPACK_PATH.split(',')
    LAPACK_Flag=dist.get_option_dict('pyPDAF')['LAPACK_Flag'][1]
    logging.info (f'LAPACK_Flag: {LAPACK_Flag}')
    if LAPACK_Flag != '': libraries += LAPACK_Flag.split(',')

# add gfortran library to the linking
if os.name != 'nt':
    suffix = 'dylib' if sys.platform == 'darwin' else 'so'
    FC = os.environ['FC'] if condaBuild == 'True' else 'gfortran'
    result = subprocess.run([FC, '--print-file', 'libgfortran.'+suffix], stdout=subprocess.PIPE)
    result = result.stdout.decode()
    result = result[:-18] if sys.platform == 'darwin' else result[:-15]
    library_dirs+=[result,]
    library_dirs+=['/usr/lib', ]
    # somehow gfortran is always necessary
    libraries += ['gfortran', 'm']
    if fortran_compiler == 'ifort': libraries += ['ifcore', 'ifcoremt']

logging.info (f'extra_compile_args: {extra_compile_args}')
logging.info (f'extra_link_args: {extra_link_args}')
logging.info (f'extra_objects: {extra_objects}')
logging.info (f'library_dirs: {library_dirs}')
logging.info (f'libraries: {libraries}')

class build_ext(build_ext_orig):
    """Pre-installation for pre build command.
        compiling PDAF iso_c_binding interface
    """
    def run(self):
        for ext in self.extensions:
            if ext.name == 'pyPDAF.PDAFc':
                cwd = os.getcwd()
                # compile PDAF and PDAFc
                PDAFc_build_dir = os.path.join(pwd, 'pyPDAF', 'fortran', 'build')
                os.makedirs(PDAFc_build_dir, exist_ok=True)
                os.chdir(PDAFc_build_dir)
                shutil.copyfile(os.path.join(pwd, 'PDAFBuild', 'CMakeLists.txt'),
                                os.path.join(PDAFdir, 'src', 'CMakeLists.txt')
                                )
                logging.info(f'cmake -DConfig_PATH={cmake_config_path} -DPDAF_PATH={PDAFdir} ..')
                os.system(f'cmake -DConfig_PATH={cmake_config_path} -DPDAF_PATH={PDAFdir} ..')
                # --config Release argument is used for multi-configuration generator, e.g., Visual Studio
                # This should just be a dummy argument in Linux and Mac where 
                # CMAKE_BUILD_TYPE should be used to set the DEBUG/RELEASE versions
                os.system('cmake --build . --verbose --target install --config Release')
                os.chdir(cwd)
        super().run()


ext_modules = [Extension('pyPDAF.PDAFc',
                          ['pyPDAF/fortran/PDAFc.pyx']),
               Extension('pyPDAF.UserFunc',
                         ['pyPDAF/UserFunc.pyx']),
               Extension('pyPDAF.PDAF',
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
    packages=["pyPDAF"]
)
