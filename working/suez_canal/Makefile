#-----------BEGIN MAKEFILE---------------------------
FUNWAVE_DIR = ../..
WORK_DIR    = .
COMPILER    = gnu
PARALLEL    = true
EXEC        = funwave
PRECISION   = single

#-----------DEFINE FLAGS-----------------------------
FLAG_1  		= -DVESSEL
FLAG_2			= -DPARALLEL
FLAG_3			= -DDEEP_DRAFT_VESSEL

#----------------uncommon options---------------------
DEF_FC      =
DEF_FC_FLAG =
SPHERICAL   = false
MPI         = openmpi
DEBUG       = false
INCS        = $(IOINCS) $(GOTMINCS)
LIBS        = $(PV3LIB) $(CLIB)  $(PARLIB) $(IOLIBS) $(MPILIB) $(GOTMLIB)
CLIB        =
MDEPFLAGS   = --cpp --fext=f90 --file=-
RANLIB      = ranlib

#----------include the essential makefiles------------
include $(FUNWAVE_DIR)/GNUMake/Essential/Make_Essential
