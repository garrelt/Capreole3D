# Makefile for MPI hydro
F90 = ifort
LDR     = $(F90)

# F90 options
F90FLAGS = -O3 -u -fpe0
#F90FLAGS = -xW -O3 -vec_report -u -ipo
#F90FLAGS = -O3 -ipo

OPTIONS = $(F90FLAGS)

LDFLAGS = $(OPTIONS)
LIBS =

# list of objects we're using

PRECISION = precision.o

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

.f90.o:
	$(F90) -c $(OPTIONS) $<
.F90.o:
	$(F90) -c $(OPTIONS) $<

.f.mod:
	$(F90) -c $(OPTIONS) $<

.SUFFIXES: .f90 .F90 .mod .o

# ISW-------------------------------------------------------------------

cart-ellipseclump : $(PRECISION) $(SIZES) $(SCALING) $(CGSCONS) $(ATOMIC) $(NOMPI) $(MESH) $(CART-COORDS) $(HYDRO) $(TIME) $(CART-ROUTINES) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(NO_IONIC) $(CART-ELLIPSECLUMP) $(INTEGRATE)  $(OUTPUT) $(EVOLVE) $(CAPREOLE) 
	$(F90) $(OPTIONS) -o $@ $(NOMPI) $(MESH) $(CART-COORDS) $(HYDRO) $(TIME) $(CART-ROUTINES) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(NO_IONIC) $(CART-ELLIPSECLUMP) $(INTEGRATE) $(OUTPUT) $(EVOLVE) $(CAPREOLE) $(LIBS)

mpi_cart-ellipseclump : $(PRECISION) $(SIZES) $(SCALING) $(ATOMIC) $(MPI) $(MESH) $(CART-COORDS) $(HYDRO)  $(CART-ROUTINES) $(TIME) $(PROT) $(BOUNDARY) $(ROESOL) $(CGSCONS) $(NO_IONIC) $(CART-ELLIPSECLUMP) $(INTEGRATE)  $(OUTPUT) $(HYDRO) $(MPI_HYDRO) 
	$(F90) $(OPTIONS) -o $@ $(MPI) $(GRID) $(CART-COORDS) $(CART-ROUTINES) $(NO_IONIC) $(PROT) $(BOUNDARY) $(ROESOL) $(CART-ELLIPSECLUMP) $(INTEGRATE)  $(OUTPUT) $(HYDRO) $(MPI_HYDRO) $(LIBS)

clean:
	rm -f *.o *.mod *.l *.il


