# Makefile for Capreole.
#
# Author: Garrelt Mellema

# This Makefile can make different versions of Capreole.
# These versions differ in their parallelization and/or
# different initial conditions (called "problems").
#
# Note 1: Parallelization
# The parallelization intended is specified in the name
# of the executable: _omp means OpenMP (shared memory), 
# _mpi means MPI (distributed memory). Both can also be
# used at the same time (if your architecture supports
# it.
#
# Note 2: Initial conditions
# Different initial conditions are specified by modules
# called "geometry"-"problem name", where geometry for
# 3D Capreole is almost always "cart" for cartesian.
# So, cart-densityfield is a problem module (initial
# conditions) for a triaxial cloud, possibly hit by
# a shock wave.
#
# Note 3: Compiler & Flags
# The compiler is specified by the FC variable (MPIFC for the MPI
# compiler). We have only extensively used the Intel F90 compiler. 
# Support for other compilers will have to be added.
# Parts of the code need to know about the compiler, this is
# done through preprocessor statements. So when compiling with
# intel compiler, -DIFORT needs to be specified. Support for
# new compilers thus needs to be added in the code too.
#
# Note 4: Recompiling
# Some dependencies are through module parameters, and thus
# not recognized by make. Best practise is to run "make clean"
# before running "make".
#-------------------------------------------------------

# Compiler: gfortran
#FC = gfortran # GNU compiler
#MPIFC = mpif90 # MPI compiler

# F90 options (gfortran)
#GFORTFLAGS = -O3 -DGFORT -DMPILOG
# Processor dependent optimization
#F90FLAGS1 = $(GFORTFLAGS) 

# These flags should be added to the F90FLAGS1 depending on the executable
# made. Specify this below on a per executable basis.
#MPI_FLAGS = -I/usr/include/lam -DMPI # For LAM mpi (Stockholm)
#MPI_FLAGS = -DMPI # 
#MPI_FLAGS = -DMPI -DMPILOG # Add more (MPI node) diagnostic output
#OPENMP_FLAGS = -openmp # For gfortran compiler

#-------------------------------------------------------
# Compiler: ifort (Intel) best tested
FC = ifort # Intel compiler
MPIFC = mpif90 # MPI compiler

# F90 options (ifort)
#IFORTFLAGS = -O0 -g -DIFORT
IFORTFLAGS = -O3 -vec_report -u -fpe0 -ipo -DIFORT -shared-intel #-check all -traceback
#IFORTFLAGS = -O3 -vec_report -u -fpe0 -ipo -mcmodel=medium -shared-intel -DIFORT #-check all -traceback
# Processor dependent optimization
#F90FLAGS1 = $(IFORTFLAGS) 
F90FLAGS1 = -xW $(IFORTFLAGS) 
#F90FLAGS1 = -xO $(IFORTFLAGS) 
#F90FLAGS1 = -xT $(IFORTFLAGS) # Laptop 
#F90FLAGS1 = -xB $(IFORTFLAGS)

# These flags should be added to the F90FLAGS1 depending on the executable
# made. Specify this below on a per executable basis.
#MPI_FLAGS = -I/usr/include/lam -DMPI # For LAM mpi (Stockholm)
MPI_FLAGS = -DMPI # 
#MPI_FLAGS = -DMPI -DMPILOG # Add more (MPI node) diagnostic output
OPENMP_FLAGS = -openmp # For Intel compiler

#-------------------------------------------------------

# Compiler: Sun
#(possible problems with constant definition. Cannot have sqrt in constant
# definition)
#FC = f95 # Sun compiler
#MPIFC = mpif90 # MPI compiler

# F90 options (ifort)
#SUNFLAGS = -O3 -DSUN
# Processor dependent optimization
#F90FLAGS1 = $(SUNFLAGS) 
#F90FLAGS1 = -xW $(SUNFLAGS) 

# These flags should be added to the F90FLAGS1 depending on the executable
# made. Specify this below on a per executable basis.
#MPI_FLAGS = -I/usr/include/lam -DMPI # For LAM mpi (Stockholm)
#MPI_FLAGS = -DMPI # 
#MPI_FLAGS = $(MPI_FLAGS) -DMPILOG # Add more (MPI node) diagnostic output
#OPENMP_FLAGS = -openmp # For Sun compiler

#-------------------------------------------------------

# PGI compiler (not recently used/tested)
#FC = pf90
#MPIFC = mpif77
#MPIFC = mpif90

# F90 options (pgi)
#PGIFLAGS = -O3 -fast -DPGI
#F90FLAGS1 = -tp barcelona-64  $(PGIFLAGS) # ranger processors

# These flags should be added to the F90FLAGS1 depending on the executable
# made. Specify this below on a per executable basis.
#MPI_FLAGS = -DMPI 
#MPI_FLAGS = $(MPI_FLAGS) -DMPILOG # Add more (MPI node) diagnostic output
#OPENMP_FLAGS = -mp 

#-------------------------------------------------------

# Absoft compilers (not recently used/tested)
#FC = f90
#ABSOFTF90FLAGS = -O3 -cpu:opteron #-DMPI
#ABSOFTF90FLAGS = -keep -Wv,-Pprocs,1 -O3 -cpu:opteron #-DMPI
LIBS = #-lU77

#-------------------------------------------------------

#LDR     = $(F90)

OPTIONS = $(F90FLAGS)

LDFLAGS = $(OPTIONS) #-L/afs/astro.su.se/pkg/intel/Compiler/11.1/056/lib/intel64/
LIBS = -lirc

#-------------------------------------------------------

# list of objects we're using

CONSTANTS = mathconstants.o cgsconstants.o cgsphotoconstants.o cgsastroconstants.o

COSMOPARMS = cosmoparms.o

PROT = protection2.o

LOF = lof-2.2.o

ROESOL = roesol-adv.o

INTEGRATE_TVD = integrate_tvd.o

TVDSOL = tvdsolver.o

INTEGRATE_VLFVS = integrate_vlfvs.o

# CART problems

CART-CONSTANT = cart-constant.o

CART-ELLIPSECLUMP = cart-ellipseclump.o

CART-DENSITYFIELD = cart-densityfield.o

CART-ENRIQUE = cart-enrique2.o

CART-HALO = cart-halo.o

CART-MINIHALO = cart-minihalo.o

# RT-3D
COOLING = cooling.o
HCOOLING = cooling_h.o
HMCOOLING = cooling_hm.o
RT-3D_BASIC = romberg.o tped.o radiation.o clumping.o doric.o thermal.o 
RT-3D_BASIC_ME = romberg.o tped.o radiation_me.o clumping.o doric.o thermal.o 
RT-3D = $(RT-3D_BASIC) sourceprops.o ionic.o
RT-3D_OMP = $(RT-3D_BASIC) sourceprops.o ionic_openmp.o
RT-3D-PP = $(RT-3D_BASIC) sourceprops_pp.o ionic_pp.o
RT-3D-PP_ME = $(RT-3D_BASIC_ME) sourceprops_pp.o ionic_pp.o
RT-3D-PP-NA = sourceprops_pp.o rad_evolve_planeparallel_noav.o $(RT-3D_BASIC)

# amuse -------------------------------------------------------------------

cart-amuse: F90=$(FC)
cart-amuse: F90FLAGS = $(F90FLAGS1)
cart-amuse : precision.o file_admin.o string.o sizes.o noscaling.o $(CONSTANTS) abundances.o atomic.o no_mpi.o clocks.o mesh.o cart-coords.o hydro.o time.o cart-routines.o $(PROT) boundary.o $(LOF) $(ROESOL) no_ionic.o cart-amuse.o integrate-strang.o  ah3out.o evolve.o capreole.o
	$(F90) $(OPTIONS) -o $@ file_admin.o string.o no_mpi.o clocks.o mesh.o cart-coords.o hydro.o time.o cart-routines.o $(PROT) boundary.o $(LOF) $(ROESOL) no_ionic.o cart-amuse.o integrate-strang.o ah3out.o evolve.o capreole.o $(LIBS)

mpi_cart-amuse: F90=$(MPIFC)
mpi_cart-amuse: F90FLAGS = $(F90FLAGS1) $(MPI_FLAGS)
mpi_cart-amuse : precision.o file_admin.o string.o sizes.o noscaling.o $(CONSTANTS) abundances.o atomic.o mpi.o clocks.o mesh.o cart-coords.o hydro.o time.o cart-routines.o $(PROT) boundary.o $(LOF) $(ROESOL) no_ionic.o cart-amuse.o integrate-strang.o  ah3out.o evolve.o capreole.o
	$(F90) $(OPTIONS) -o $@ file_admin.o string.o mpi.o clocks.o mesh.o cart-coords.o hydro.o time.o cart-routines.o $(PROT) boundary.o $(LOF) $(ROESOL) no_ionic.o cart-amuse.o integrate-strang.o ah3out.o evolve.o capreole.o $(LIBS)

# ellipseclump------------------------------------------------------------------

cart-ellipseclump : F90=$(FC)
cart-ellipseclump : F90FLAGS = $(F90FLAGS1)
cart-ellipseclump : precision.o file_admin.o string.o sizes.o noscaling.o $(CONSTANTS) abundances.o atomic.o no_mpi.o clocks.o mesh.o cart-coords.o hydro.o time.o cart-routines.o $(PROT) boundary.o $(LOF) $(ROESOL) no_ionic.o cart-ellipseclump.o integrate-strang.o  ah3out.o evolve.o capreole.o
	$(F90) $(OPTIONS) -o $@ file_admin.o string.o no_mpi.o clocks.o mesh.o cart-coords.o hydro.o time.o cart-routines.o $(PROT) boundary.o $(LOF) $(ROESOL) no_ionic.o cart-ellipseclump.o integrate-strang.o ah3out.o evolve.o capreole.o $(LIBS)

mpi_cart-ellipseclump : F90=$(MPIFC)
mpi_cart-ellipseclump : F90FLAGS = $(F90FLAGS1) $(MPI_FLAGS)
mpi_cart-ellipseclump : precision.o file_admin.o string.o sizes.o noscaling.o $(CONSTANTS) abundances.o atomic.o mpi.o clocks.o mesh.o cart-coords.o hydro.o time.o cart-routines.o $(PROT) boundary.o $(LOF) $(ROESOL) no_ionic.o cart-ellipseclump.o integrate-strang.o  ah3out.o evolve.o capreole.o
	$(F90) $(OPTIONS) -o $@ file_admin.o string.o mpi.o clocks.o mesh.o cart-coords.o hydro.o time.o cart-routines.o $(PROT) boundary.o $(LOF) $(ROESOL) no_ionic.o cart-ellipseclump.o integrate-strang.o ah3out.o evolve.o capreole.o $(LIBS)

cart-ellipseclump_hc_c2ray : $(PRECISION) $(FILES) $(STRINGS) $(SIZES) $(SCALING) $(CONSTANTS) $(CGSPHOTOCONS) $(ABUNDANCES) $(ATOMIC) $(NOMPI) $(MESH) $(CART-COORDS) $(HYDRO) $(HCOOLING) $(RT-3D-PP) $(TIME) $(CART-ROUTINES) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-ELLIPSECLUMP) $(INTEGRATE)  $(OUTPUT) $(EVOLVE) $(CAPREOLE) 
	$(F90) $(OPTIONS) -o $@ $(STRINGS) $(NOMPI) $(MESH) $(CART-COORDS) $(HYDRO) $(TIME) $(CART-ROUTINES) $(HCOOLING) $(RT-3D-PP) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-ELLIPSECLUMP) $(INTEGRATE) $(OUTPUT) $(EVOLVE) $(CAPREOLE) $(LIBS)

cart-ellipseclump_c2ray : $(PRECISION) $(FILES) $(STRINGS) $(SIZES) $(SCALING) $(CONSTANTS) $(CGSPHOTOCONS) $(ABUNDANCES) $(ATOMIC) $(NOMPI) $(MESH) $(CART-COORDS) $(HYDRO) $(COOLING) $(RT-3D-PP) $(TIME) $(CART-ROUTINES) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-ELLIPSECLUMP) $(INTEGRATE)  $(OUTPUT) $(EVOLVE) $(CAPREOLE) 
	$(F90) $(OPTIONS) -o $@ $(STRINGS) $(NOMPI) $(MESH) $(CART-COORDS) $(HYDRO) $(TIME) $(CART-ROUTINES) $(COOLING) $(RT-3D-PP) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-ELLIPSECLUMP) $(INTEGRATE) $(OUTPUT) $(EVOLVE) $(CAPREOLE) $(LIBS)

mpi_cart-ellipseclump_hc_c2ray_me : $(PRECISION) $(FILES) $(STRINGS) $(SIZES) $(SCALING) $(CONSTANTS) $(CGSPHOTOCONS) $(ABUNDANCES) $(ATOMIC) $(MPI) $(MESH) $(CART-COORDS) $(HYDRO) $(HCOOLING) $(RT-3D-PP_ME) $(TIME) $(CART-ROUTINES) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-ELLIPSECLUMP) $(INTEGRATE)  $(OUTPUT) $(EVOLVE) $(CAPREOLE) 
	$(F90) $(OPTIONS) -o $@ $(STRINGS) $(MPI) $(MESH) $(CART-COORDS) $(HYDRO) $(TIME) $(CART-ROUTINES) $(HCOOLING) $(RT-3D-PP_ME) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-ELLIPSECLUMP) $(INTEGRATE) $(OUTPUT) $(EVOLVE) $(CAPREOLE) $(LIBS)

mpi_cart-ellipseclump_hc_c2ray : $(PRECISION) $(FILES) $(STRINGS) $(SIZES) $(SCALING) $(CONSTANTS) $(CGSPHOTOCONS) $(ABUNDANCES) $(ATOMIC) $(MPI) $(MESH) $(CART-COORDS) $(HYDRO) $(HCOOLING) $(RT-3D-PP) $(TIME) $(CART-ROUTINES) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-ELLIPSECLUMP) $(INTEGRATE)  $(OUTPUT) $(EVOLVE) $(CAPREOLE) 
	$(F90) $(OPTIONS) -o $@ $(STRINGS) $(MPI) $(MESH) $(CART-COORDS) $(HYDRO) $(TIME) $(CART-ROUTINES) $(HCOOLING) $(RT-3D-PP) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-ELLIPSECLUMP) $(INTEGRATE) $(OUTPUT) $(EVOLVE) $(CAPREOLE) $(LIBS)

mpi_cart-ellipseclump_c2ray : $(PRECISION) $(FILES) $(STRINGS) $(SIZES) $(SCALING) $(CONSTANTS) $(CGSPHOTOCONS) $(ABUNDANCES) $(ATOMIC) $(MPI) $(MESH) $(CART-COORDS) $(HYDRO) $(COOLING) $(RT-3D-PP) $(TIME) $(CART-ROUTINES) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-ELLIPSECLUMP) $(INTEGRATE)  $(OUTPUT) $(EVOLVE) $(CAPREOLE) 
	$(F90) $(OPTIONS) -o $@ $(STRINGS) $(MPI) $(MESH) $(CART-COORDS) $(HYDRO) $(TIME) $(CART-ROUTINES) $(COOLING) $(RT-3D-PP) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-ELLIPSECLUMP) $(INTEGRATE) $(OUTPUT) $(EVOLVE) $(CAPREOLE) $(LIBS)

cart-halo_hc_c2ray : $(PRECISION) $(FILES) $(STRINGS) $(SIZES) $(SCALING) $(CONSTANTS) $(CGSPHOTOCONS) $(ABUNDANCES) $(ATOMIC) $(NOMPI) $(MESH) $(CART-COORDS) $(HYDRO) $(HCOOLING) $(RT-3D) $(TIME) $(CART-ROUTINES) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-HALO) $(INTEGRATE)  $(OUTPUT) $(EVOLVE) $(CAPREOLE) 
	$(F90) $(OPTIONS) -o $@ $(STRINGS) $(NOMPI) $(MESH) $(CART-COORDS) $(HYDRO) $(TIME) $(CART-ROUTINES) $(HCOOLING) $(RT-3D) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-HALO) $(INTEGRATE) $(OUTPUT) $(EVOLVE) $(CAPREOLE) $(LIBS)

cart-halo : $(PRECISION) $(FILES) $(STRINGS) $(SIZES) $(SCALING) $(CONSTANTS) $(CGSPHOTOCONS) $(ABUNDANCES) $(ATOMIC) $(NOMPI) $(MESH) $(CART-COORDS) $(HYDRO) $(HCOOLING) $(NO_IONIC) $(TIME) $(CART-ROUTINES) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-HALO) $(INTEGRATE)  $(OUTPUT) $(EVOLVE) $(CAPREOLE) 
	$(F90) $(OPTIONS) -o $@ $(STRINGS) $(NOMPI) $(MESH) $(CART-COORDS) $(HYDRO) $(TIME) $(CART-ROUTINES) $(HCOOLING) $(NO_IONIC) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-HALO) $(INTEGRATE) $(OUTPUT) $(EVOLVE) $(CAPREOLE) $(LIBS)

cart-densityfield_hc_c2ray : $(PRECISION) $(FILES) $(STRINGS) $(SIZES) $(SCALING) $(CONSTANTS) $(CGSPHOTOCONS) $(ABUNDANCES) $(ATOMIC) $(NOMPI) $(MESH) $(CART-COORDS) $(HYDRO) $(HCOOLING) $(RT-3D-PP) $(TIME) $(CART-ROUTINES) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-DENSITYFIELD) $(INTEGRATE)  $(OUTPUT) $(EVOLVE) $(CAPREOLE) 
	$(F90) $(OPTIONS) -o $@ $(STRINGS) $(NOMPI) $(MESH) $(CART-COORDS) $(HYDRO) $(TIME) $(CART-ROUTINES) $(HCOOLING) $(RT-3D-PP) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-DENSITYFIELD) $(INTEGRATE) $(OUTPUT) $(EVOLVE) $(CAPREOLE) $(LIBS)

cart-densityfield_c2ray : $(PRECISION) $(FILES) $(STRINGS) $(SIZES) $(SCALING) $(CONSTANTS) $(CGSPHOTOCONS) $(ABUNDANCES) $(ATOMIC) $(NOMPI) $(MESH) $(CART-COORDS) $(HYDRO) $(COOLING) $(RT-3D-PP) $(TIME) $(CART-ROUTINES) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-DENSITYFIELD) $(INTEGRATE)  $(OUTPUT) $(EVOLVE) $(CAPREOLE) 
	$(F90) $(OPTIONS) -o $@ $(STRINGS) $(NOMPI) $(MESH) $(CART-COORDS) $(HYDRO) $(TIME) $(CART-ROUTINES) $(COOLING) $(RT-3D-PP) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-DENSITYFIELD) $(INTEGRATE) $(OUTPUT) $(EVOLVE) $(CAPREOLE) $(LIBS)

mpi_cart-densityfield_hc_c2ray : $(PRECISION) $(FILES) $(STRINGS) $(SIZES) $(SCALING) $(CONSTANTS) $(CGSPHOTOCONS) $(ABUNDANCES) $(ATOMIC) $(MPI) $(MESH) $(CART-COORDS) $(HYDRO) $(HCOOLING) $(RT-3D-PP) $(TIME) $(CART-ROUTINES) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-DENSITYFIELD) $(INTEGRATE)  $(OUTPUT) $(EVOLVE) $(CAPREOLE) 
	$(F90) $(OPTIONS) -o $@ $(STRINGS) $(MPI) $(MESH) $(CART-COORDS) $(HYDRO) $(TIME) $(CART-ROUTINES) $(HCOOLING) $(RT-3D-PP) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-DENSITYFIELD) $(INTEGRATE) $(OUTPUT) $(EVOLVE) $(CAPREOLE) $(LIBS)

mpi_cart-densityfield_c2ray : $(PRECISION) $(FILES) $(STRINGS) $(SIZES) $(SCALING) $(CONSTANTS) $(CGSPHOTOCONS) $(ABUNDANCES) $(ATOMIC) $(MPI) $(MESH) $(CART-COORDS) $(HYDRO) $(COOLING) $(RT-3D-PP) $(TIME) $(CART-ROUTINES) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-DENSITYFIELD) $(INTEGRATE)  $(OUTPUT) $(EVOLVE) $(CAPREOLE) 
	$(F90) $(OPTIONS) -o $@ $(STRINGS) $(MPI) $(MESH) $(CART-COORDS) $(HYDRO) $(TIME) $(CART-ROUTINES) $(COOLING) $(RT-3D-PP) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-DENSITYFIELD) $(INTEGRATE) $(OUTPUT) $(EVOLVE) $(CAPREOLE) $(LIBS)

cart-minihalo_hc_c2ray : $(PRECISION) $(FILES) $(STRINGS) $(SIZES) $(SCALING) $(CONSTANTS) $(CGSPHOTOCONS) $(COSMOPARMS) $(ABUNDANCES) $(ATOMIC) $(NOMPI) $(MESH) $(CART-COORDS) $(HYDRO) $(HCOOLING) $(RT-3D-PP) $(TIME) $(CART-ROUTINES) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-MINIHALO) $(INTEGRATE)  $(OUTPUT) $(EVOLVE) $(CAPREOLE) 
	$(F90) $(OPTIONS) -o $@ $(STRINGS) $(NOMPI) $(MESH) $(CART-COORDS) $(HYDRO) $(TIME) $(CART-ROUTINES) $(HCOOLING) $(RT-3D-PP) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-MINIHALO) $(INTEGRATE) $(OUTPUT) $(EVOLVE) $(CAPREOLE) $(LIBS)

cart-minihalo_c2ray : $(PRECISION) $(FILES) $(STRINGS) $(SIZES) $(SCALING) $(CONSTANTS) $(CGSPHOTOCONS) $(COSMOPARMS) $(ABUNDANCES) $(ATOMIC) $(NOMPI) $(MESH) $(CART-COORDS) $(HYDRO) $(COOLING) $(RT-3D-PP) $(TIME) $(CART-ROUTINES) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-MINIHALO) $(INTEGRATE)  $(OUTPUT) $(EVOLVE) $(CAPREOLE) 
	$(F90) $(OPTIONS) -o $@ $(STRINGS) $(NOMPI) $(MESH) $(CART-COORDS) $(HYDRO) $(TIME) $(CART-ROUTINES) $(COOLING) $(RT-3D-PP) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-MINIHALO) $(INTEGRATE) $(OUTPUT) $(EVOLVE) $(CAPREOLE) $(LIBS)

mpi_cart-minihalo_hc_c2ray : $(PRECISION) $(FILES) $(STRINGS) $(SIZES) $(SCALING) $(CONSTANTS) $(CGSPHOTOCONS) $(COSMOPARMS) $(ABUNDANCES) $(ATOMIC) $(MPI) $(MESH) $(CART-COORDS) $(HYDRO) $(HCOOLING) $(RT-3D-PP) $(TIME) $(CART-ROUTINES) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-MINIHALO) $(INTEGRATE)  $(OUTPUT) $(EVOLVE) $(CAPREOLE) 
	$(F90) $(OPTIONS) -o $@ $(STRINGS) $(MPI) $(MESH) $(CART-COORDS) $(HYDRO) $(TIME) $(CART-ROUTINES) $(HCOOLING) $(RT-3D-PP) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-MINIHALO) $(INTEGRATE) $(OUTPUT) $(EVOLVE) $(CAPREOLE) $(LIBS)

mpi_cart-minihalo_c2ray : $(PRECISION) $(FILES) $(STRINGS) $(SIZES) $(SCALING) $(CONSTANTS) $(CGSPHOTOCONS) $(COSMOPARMS) $(ABUNDANCES) $(ATOMIC) $(MPI) $(MESH) $(CART-COORDS) $(HYDRO) $(COOLING) $(RT-3D-PP) $(TIME) $(CART-ROUTINES) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-MINIHALO) $(INTEGRATE)  $(OUTPUT) $(EVOLVE) $(CAPREOLE) 
	$(F90) $(OPTIONS) -o $@ $(STRINGS) $(MPI) $(MESH) $(CART-COORDS) $(HYDRO) $(TIME) $(CART-ROUTINES) $(COOLING) $(RT-3D-PP) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-MINIHALO) $(INTEGRATE) $(OUTPUT) $(EVOLVE) $(CAPREOLE) $(LIBS)

cart-enrique_hc_c2ray : $(PRECISION) $(FILES) $(STRINGS) $(SIZES) $(SCALING) $(MATHCONS) $(CGSCONS) $(CGSASTROCONS) $(CGSPHOTOCONS) $(ABUNDANCES) $(ATOMIC) $(NOMPI) $(MESH) $(CART-COORDS) $(HYDRO) $(HMCOOLING) $(RT-3D) $(TIME) $(CART-ROUTINES) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-ENRIQUE) $(INTEGRATE)  $(OUTPUT) $(EVOLVE) $(CAPREOLE) 
	$(F90) $(OPTIONS) -o $@ $(STRINGS) $(NOMPI) $(MESH) $(CART-COORDS) $(HYDRO) $(TIME) $(CART-ROUTINES) $(HMCOOLING) $(RT-3D) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-ENRIQUE) $(INTEGRATE) $(OUTPUT) $(EVOLVE) $(CAPREOLE) $(LIBS)

cart-enrique_hc_c2ray_openmp : $(PRECISION) $(FILES) $(STRINGS) $(SIZES) $(SCALING) $(MATHCONS) $(CGSCONS) $(CGSASTROCONS) $(CGSPHOTOCONS) $(ABUNDANCES) $(ATOMIC) $(NOMPI) $(MESH) $(CART-COORDS) $(HYDRO) $(HMCOOLING) $(RT-3D_OMP) $(TIME) $(CART-ROUTINES) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-ENRIQUE) $(INTEGRATE)  $(OUTPUT) $(EVOLVE) $(CAPREOLE) 
	$(F90) $(OPTIONS) -o $@ $(STRINGS) $(NOMPI) $(MESH) $(CART-COORDS) $(HYDRO) $(TIME) $(CART-ROUTINES) $(HMCOOLING) $(RT-3D_OMP) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-ENRIQUE) $(INTEGRATE) $(OUTPUT) $(EVOLVE) $(CAPREOLE) $(LIBS)

cart-enrique_c2ray : $(PRECISION) $(FILES) $(STRINGS) $(SIZES) $(SCALING) $(MATHCONS) $(CGSCONS) $(CGSASTROCONS) $(CGSPHOTOCONS) $(ABUNDANCES) $(ATOMIC) $(NOMPI) $(INPUT) $(MESH) $(CART-COORDS) $(HYDRO) $(HMCOOLING) $(RT-3D) $(TIME) $(CART-ROUTINES) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-ENRIQUE) $(INTEGRATE)  $(OUTPUT) $(EVOLVE) $(CAPREOLE) 
	$(F90) $(OPTIONS) -o $@ $(STRINGS) $(NOMPI) $(INPUT) $(MESH) $(CART-COORDS) $(HYDRO) $(TIME) $(CART-ROUTINES) $(HMCOOLING) $(RT-3D) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-ENRIQUE) $(INTEGRATE) $(OUTPUT) $(EVOLVE) $(CAPREOLE) $(LIBS)

cart-constant_hc_c2ray : $(PRECISION) $(FILES) $(STRINGS) $(SIZES) $(SCALING) $(MATHCONS) $(CGSCONS) $(CGSASTROCONS) $(CGSPHOTOCONS) $(ABUNDANCES) $(ATOMIC) $(NOMPI) $(MESH) $(CART-COORDS) $(HYDRO) $(HMCOOLING) $(RT-3D) $(TIME) $(CART-ROUTINES) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-CONSTANT) $(INTEGRATE)  $(OUTPUT) $(EVOLVE) $(CAPREOLE) 
	$(F90) $(OPTIONS) -o $@ $(STRINGS) $(NOMPI) $(MESH) $(CART-COORDS) $(HYDRO) $(TIME) $(CART-ROUTINES) $(HMCOOLING) $(RT-3D) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-CONSTANT) $(INTEGRATE) $(OUTPUT) $(EVOLVE) $(CAPREOLE) $(LIBS)

cart-ellipseclump_tvd : $(PRECISION) $(FILES) $(STRINGS) $(SIZES) $(SCALING) $(CONSTANTS) $(ABUNDANCES) $(ATOMIC) $(NOMPI) $(MESH) $(CART-COORDS) $(HYDRO) $(TIME) $(CART-ROUTINES) $(PROT) $(BOUNDARY) $(LOF) $(TVDSOL) $(NO_IONIC) $(CART-ELLIPSECLUMP) $(INTEGRATE_TVD)  $(OUTPUT) $(EVOLVE) $(CAPREOLE) 
	$(F90) $(OPTIONS) -o $@ $(STRINGS) $(NOMPI) $(MESH) $(CART-COORDS) $(HYDRO) $(TIME) $(CART-ROUTINES) $(PROT) $(BOUNDARY) $(LOF) $(TVDSOL) $(NO_IONIC) $(CART-ELLIPSECLUMP) $(INTEGRATE_TVD) $(OUTPUT) $(EVOLVE) $(CAPREOLE) $(LIBS)

cart-ellipseclump_vlfvs : $(PRECISION) $(FILES) $(STRINGS) $(SIZES) $(SCALING) $(CONSTANTS) $(ABUNDANCES) $(ATOMIC) $(NOMPI) $(MESH) $(CART-COORDS) $(HYDRO) $(TIME) $(CART-ROUTINES) $(PROT) $(BOUNDARY) $(LOF) $(NO_IONIC) $(CART-ELLIPSECLUMP) $(INTEGRATE_VLFVS)  $(OUTPUT) $(EVOLVE) $(CAPREOLE) 
	$(F90) $(OPTIONS) -o $@ $(STRINGS) $(NOMPI) $(MESH) $(CART-COORDS) $(HYDRO) $(TIME) $(CART-ROUTINES) $(PROT) $(BOUNDARY) $(LOF) $(NO_IONIC) $(CART-ELLIPSECLUMP) $(INTEGRATE_VLFVS) $(OUTPUT) $(EVOLVE) $(CAPREOLE) $(LIBS)

clean:
	rm -f *.o *.mod *.l *.il *.vo

.f.o:
	$(F90) -c $(OPTIONS) $<

.f90.o:
	$(F90) -c $(OPTIONS) $<

.F90.o:
	$(F90) -c $(OPTIONS) $<

f.mod:
	$(F90) -c $(OPTIONS) $<

.SUFFIXES: .f90 .F90 .mod .o




