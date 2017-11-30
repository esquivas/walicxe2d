#---------------------------------------------------
# WALICXE-2D
#---------------------------------------------------
#   Name of the executable
PROGRAM=WALICXE-2D
#
#  Seaerh Path, object files are searched in the order below
VPATH = . : ../src
#
#   Choice of compiler: ifort and gfortran are the only implemented
#   Choice of compiler: ifort and gfortran are tested
#   if MPI enabled, do not use mpif90, but rather the base compiler
#   of the mpi build (see 'mpif90 --version')
#COMPILER= gfortran
COMPILER= ifort
#
#   Compiler flags, make sure they are compatible with 
#   the previous selection
#  -g -traceback
FLAGS= -O3 -cpp
# -i_dynamic -vec-report0
#
#---------------------------------------------------
#   compilation time parameters (Y=on, N=off)
#   watch out all of the following is case sensitive
#---------------------------------------------------
#   MPI paralalelization
MPIP = Y
#   Double Precision (Default is single)
DOUBLEP= Y
#   Riemann Solver for the fluxes (available: HLL, HLL_HLLC, HLLC)
FLUXES = HLLC
#   Type of output (BIN are raw unformatted, VTK are binary)
OUTBIN = Y
OUTVTK = N
#   additional equations (i.e. passive scalars)?
PASSIVES = Y
#   cylindrical geometry source terms (axisymmetric)
CYLINDRICAL = Y
#   adiabatic or type of cooling (choose only one)
#   options are ADIABATIC, H, BBC, CORONAL
#   ADIABATIC: of course turns off the cooling
#   H: single parametrized cooling function
#   BBC: cooling function of Benjamin, Benson and Cox (2003)
#   DMC: Coronal eq. (tabulated) from Dalgarno & Mc Cray (1972)
COOLING   = DMC
#   boundary conditions
#   options are OUTFLOW, REFLECTIVE, PERIODIC , INFLOW, OTHERB
#   INFLOW refers to the outer grid boundary  (e.g. to
#   impose a plane parallel side wind)
#   OTHER is used to add sink/sources at any point in the
#   computational domain, calls 'otherbounds.pro'
LEFTX   = REFLECTIVE
RIGHTX  = OUTFLOW
BOTTOMY = REFLECTIVE
TOPY    = OUTFLOW
OTHERB  = Y
#   choice of slope limiter, available limiters are:
#   limiter =-1: no average, 0 : no limiter, 1: Minmod,
#   2: Van Leer, 3: Van Albada,4: UMIST, 5: Woodward
#   6: Superbee
LIMITER=  2
#####################################################
# There should be no need to modify below this line #
#####################################################
#   list of objects to compile
OBJECTS =\
parameters.o	\
globals.o		\
dmc_module.o	\
sn.o			\
jet.o 			\
user_mod.o 		\
main.o			\
initmain.o		\
basegrid.o		\
refine.o        \
updatelpup.o	\
calcprim.o 		\
uprim.o 		\
loadbalance.o  	\
output.o		\
timestep.o		\
sound.o 		\
tstep.o			\
hllcfluxes.o	\
hllfluxes.o		\
hll_hllcfluxes.o\
swapy.o			\
step.o			\
limiter.o 		\
primf.o 		\
primu.o 		\
locatebounds.o	\
viscosity.o		\
boundaryI.o		\
boundaryII.o	\
coolingh.o		\
coolingdmc.o	\
criteria.o		\
markref.o		\
critneighup.o 	\
critneighdown.o	\
admesh.o		\
coarse.o		\
updatelpdown.o
#---------------------------------------------------
ifeq ($(DOUBLEP),Y)
FLAGS += -DDOUBLEP
ifeq ($(COMPILER),ifort)
FLAGS += -r8
endif
ifeq ($(COMPILER),gfortran)
FLAGS += -fdefault-real-8
endif
endif
ifeq ($(MPIP),Y)
FLAGS += -DMPIP
COMPILER = mpif90
endif
ifeq ($(FLUXES),HLLC)
FLAGS += -DHLLC
endif
ifeq ($(FLUXES),HLL)
FLAGS += -DHLL
endif
ifeq ($(FLUXES),HLL_HLLC)
FLAGS += -DHLL_HLLC
endif
ifeq ($(OUTBIN),Y)
FLAGS += -DOUTBIN
endif
ifeq ($(OUTVTK),Y)
FLAGS += -DOUTVTK
endif
ifeq ($(PASSIVES),Y)
FLAGS += -DPASSIVES
endif
ifeq ($(CYLINDRICAL),Y)
FLAGS += -DCYLINDRICAL
endif
ifeq ($(COOLING),ADIABATIC)
FLAGS += -DADIABATIC
endif
ifeq ($(COOLING),H)
FLAGS += -DCOOLINGH
endif
ifeq ($(COOLING),BBC)
FLAGS += -DCOOLINGBBC
endif
ifeq ($(COOLING),DMC)
FLAGS += -DCOOLINGDMC
endif
ifeq ($(LEFTX),PERIODIC)
FLAGS += -DPERIODX
endif
ifeq ($(BOTTOMY),PERIODIC)
FLAGS += -DPERIODY
endif
ifeq ($(LEFTX),REFLECTIVE)
FLAGS += -DREFXL
endif
ifeq ($(RIGHTX),REFLECTIVE)
FLAGS += -DREFXR
endif
ifeq ($(BOTTOMY),REFLECTIVE)
FLAGS += -DREFYB
endif
ifeq ($(TOPY),REFLECTIVE)
FLAGS += -DREFYT
endif
ifeq ($(LEFTX),OUTFLOW)
FLAGS += -DOUTFXL
endif
ifeq ($(RIGHTX),OUTFLOW)
FLAGS += -DOUTFXR
endif
ifeq ($(TOPY),OUTFLOW)
FLAGS += -DOUTFYT
endif
ifeq ($(BOTTOMY),OUTFLOW)
FLAGS += -DOUTFYB
endif
ifeq ($(LEFTX),INFLOW)
FLAGS += -DINFXL
endif
ifeq ($(RIGHTX),INFLOW)
FLAGS += -DINFXR
endif
ifeq ($(TOPY),INFLOW)
FLAGS += -DINFYT
endif
ifeq ($(BOTTOMY),INFLOW)
FLAGS += -DINFYB
endif
ifeq ($(OTHERB),Y)
FLAGS += -DOTHERB
endif
FLAGS += -DLIMITER=$(LIMITER)
#---------------------------------------------------
$(PROGRAM)  : ${OBJECTS}
	$(COMPILER) $(FLAGS) -o $@ $(OBJECTS)

%.o:%.f90
	$(COMPILER) $(FLAGS) -c $^ -o $@

.PHONY : clean
clean :
	rm -f *.o *.mod src/*.o src/*.mod
	rm -f $(PROGRAM)
#---------------------------------------------------
