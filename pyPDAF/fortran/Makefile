# Root directory of PDAF package
PDAFDIR=../../../PDAF-D_V1.16/
# current directory
FORTRAN_DIR=$(shell pwd)
LibDir=$(PWD)/../../lib
# Include machine-specific definitions
# For available include files see directory make.arch
# To choose a file, set PDAF_ARCH either here or by an
# environment variable.
include $(PWD)/../../linux_gfortran_openmpi.h

OBJ_INTERFACE=interface_pdaf.o

libPDAFc : $(OBJ_INTERFACE)
	$(FC) $^ -shared -L$(PDAFDIR)/lib -lpdaf-d $(LINK_LIBS) -o $(LibDir)/libPDAFc.so

%.o : %.F90
	@cd $(FORTRAN_DIR)
	$(FC) $(OPT) -fPIC $(MPI_INC) $(CPP_DEFS) -I$(PDAFDIR)/src -c $^ 

clean:
	rm *.o *.mod  $(LibDir)/*.so