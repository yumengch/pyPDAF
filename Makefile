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
pyPDAF_DIR=$(PWD)/pyPDAF
# Include machine-specific definitions
# For available include files see directory make.arch
# To choose a file, set PDAF_ARCH either here or by an
# environment variable.

# End of user specifications
######################################################

######################################################

all: libpdaf-d.a pyPDAF

info:
	@echo "Makefile to build PDAF tutorial online implementation";
	@echo "Example: 2D serial model (without parallelization)";

pyPDAF:
	python setup.py build_ext --inplace

######################################################

libpdaf-d.a: 
	@echo "++++++ Generate Filter library ++++++"
	@cd $(BASEDIR)/src; make;

clean_PDAF :
	@cd $(BASEDIR)/src; make clean;

clean :
	rm -rf build/ pyPDAF.egg-info/
	rm -rf $(pyPDAF_DIR)/PDAF/*.so $(pyPDAF_DIR)/PDAF/*.c \
			$(pyPDAF_DIR)/Cython/*.so $(pyPDAF_DIR)/Cython/*.c

