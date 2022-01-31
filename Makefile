# $Id: Makefile 746 2009-08-04 12:16:28Z lnerger $

#######################################################
# Generic Makefile for to build PDAF with dummy model #
# To choose the architecture set $PDAF_ARCH           #
#######################################################

######################################################


# User specifications
# 1. Set BASEDIR, the directory where the PDAF package resides
# 2. Set PDAF_ARCH to include compile definitions
#    (See directory BASEDIR/make.arch for files. PDAF_ARCH is filename without .h)

# Root directory of PDAF package
BASEDIR=../PDAF-D_V1.16/
# current directory
PWD=$(shell pwd)
# fortran source code directory
FORTRAN_DIR=$(PWD)/pyPDAF/fortran
CythonDIR=$(PWD)/pyPDAF
LibDir=$(PWD)/lib
# Include machine-specific definitions
# For available include files see directory make.arch
# To choose a file, set PDAF_ARCH either here or by an
# environment variable.
include $(PWD)/linux_gfortran_openmpi.h

# EXT_SUFFIX := $(shell python3-config --extension-suffix)

# End of user specifications
######################################################

######################################################

all: libpdaf-d.a libPDAFc # pybinding

info:
	@echo "Makefile to build PDAF tutorial online implementation";
	@echo "Example: 2D serial model (without parallelization)";

libPDAFc: 
	@cd $(FORTRAN_DIR); $(MAKE) libPDAFc

pybinding:
	python setup.py build_ext --inplace

######################################################

libpdaf-d.a: 
	@echo "++++++ Generate Filter library ++++++"
	@cd $(BASEDIR)/src; make;

clean :
	@cd $(FORTRAN_DIR); $(MAKE) clean
	rm -rf $(CythonDIR)/PDAF/*.so $(CythonDIR)/PDAF/*.c $(CythonDIR)/Cython/*.so $(CythonDIR)/Cython/*.c build/

