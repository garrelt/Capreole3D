# Makefile for MPI hydro
#F90 = pf90
F90 = f95
#F90 = mpif77
#F90 = ifort
#F90 = /disk/sn-12/garrelt//mpich/bin/mpif90
#F90 = ifort
LDR     = $(F90)
#PP = cpp -P

# F90 options
F90FLAGS = -O3 -cpu:opteron #-DMPI
#F90FLAGS = -keep -Wv,-Pprocs,1 -O3 -cpu:opteron #-DMPI
#F90FLAGS = -assume 2underscores -xW -O3 -vec_report -u -ipo #-DMPI
#F90FLAGS = -O3 -ipo

OPTIONS = $(F90FLAGS)

LDFLAGS = $(OPTIONS) 
LIBS = -lU77

# list of objects we're using

PRECISION = precision.o

FILES = file_admin.o

MESH = mesh.o

CART-COORDS = cart-coords.o

CART-ROUTINES = cart-routines.o

PROT = protection.o

LOF = lof-2.2.o

BOUNDARY = boundary.o

ROESOL = roesol-adv.o

# CART problems

CART-ELLIPSECLUMP = cart-ellipseclump.o

HYDRO = hydro.o

TIME = time.o

INTEGRATE = integrate-strang.o

EVOLVE = evolve.o

OUTPUT = output.o

CAPREOLE = capreole.o

MPI = mpi.o

NOMPI = no_mpi.o

SIZES = sizes.o

SCALING = scaling.o

NOSCALING = noscaling.o

ATOMIC = atomic.o

CGSCONS = cgsconstants.o

CGSASTROCONS = cgsastroconstants.o

STRINGS = string.o

.f90.o:
	$(F90) -c $(OPTIONS) $<

.F90.o:
	$(F90) -c $(OPTIONS) $<

#.F90.f90:
#	$(PP) $< $*.f90

.f.mod:
	$(F90) -c $(OPTIONS) $<

.SUFFIXES: .f90 .F90 .mod .o

# ISW-------------------------------------------------------------------

cart-ellipseclump : $(PRECISION) $(FILES) $(STRINGS) $(SIZES) $(SCALING) $(CGSCONS) $(CGSASTROCONS) $(ATOMIC) $(NOMPI) $(MESH) $(CART-COORDS) $(HYDRO) $(TIME) $(CART-ROUTINES) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(NO_IONIC) $(CART-ELLIPSECLUMP) $(INTEGRATE)  $(OUTPUT) $(EVOLVE) $(CAPREOLE) 
	$(F90) $(OPTIONS) -o $@ $(STRINGS) $(NOMPI) $(MESH) $(CART-COORDS) $(HYDRO) $(TIME) $(CART-ROUTINES) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(NO_IONIC) $(CART-ELLIPSECLUMP) $(INTEGRATE) $(OUTPUT) $(EVOLVE) $(CAPREOLE) $(LIBS)

mpi_cart-ellipseclump : $(PRECISION) $(FILES) $(STRINGS) $(SIZES) $(SCALING) $(CGSCONS) $(CGSASTROCONS) $(ATOMIC) $(MPI) $(MESH) $(CART-COORDS) $(HYDRO) $(TIME) $(CART-ROUTINES) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(NO_IONIC) $(CART-ELLIPSECLUMP) $(INTEGRATE) $(OUTPUT) $(EVOLVE) $(CAPREOLE) 
	$(F90) $(OPTIONS) -o $@ $(MPI) $(STRINGS) $(MESH) $(CART-COORDS) $(HYDRO) $(TIME) $(CART-ROUTINES) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(NO_IONIC) $(CART-ELLIPSECLUMP) $(INTEGRATE)  $(OUTPUT) $(EVOLVE) $(CAPREOLE) $(LIBS)

clean:
	rm -f *.o *.mod *.l *.il *.vo


