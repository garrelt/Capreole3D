# Makefile for MPI hydro
#F90 = pf90
#F90 = f95
#F90 = mpif77
F90 = ifort
#F90 = /disk/sn-12/garrelt//mpich/bin/mpif90
#F90 = ifort
LDR     = $(F90)
#PP = cpp -P

# F90 options
#F90FLAGS = -O3 -cpu:opteron #-DMPI
#F90FLAGS = -keep -Wv,-Pprocs,1 -O3 -cpu:opteron #-DMPI
F90FLAGS = -assume 2underscores -xW -O3 -vec_report -u -ipo -fpe0 #-DMPI
#F90FLAGS = -O3 -ipo

OPTIONS = $(F90FLAGS)

LDFLAGS = $(OPTIONS) 
LIBS = #-lU77

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

ABUNDANCES = abundances.o

CGSCONS = cgsconstants.o

CGSASTROCONS = cgsastroconstants.o

CGSPHOTOCONS = cgsphotoconstants.o

STRINGS = string.o

# RT-3D
COOLING = cooling.o
HCOOLING = cooling_h.o
RT-3D_BASIC = romberg.o tped.o radiation.o clumping.o doric.o thermal.o 
RT-3D = $(RT-3D_BASIC) sourceprops.o ionic.o
RT-3D-PP = $(RT-3D_BASIC) sourceprops_pp.o ionic_pp.o
RT-3D-PP-NA = sourceprops_pp.o rad_evolve_planeparallel_noav.o $(RT-3D_BASIC)

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

cart-ellipseclump_hc_c2ray : $(PRECISION) $(FILES) $(STRINGS) $(SIZES) $(SCALING) $(CGSCONS) $(CGSASTROCONS) $(CGSPHOTOCONS) $(ABUNDANCES) $(ATOMIC) $(NOMPI) $(MESH) $(CART-COORDS) $(HYDRO) $(HCOOLING) $(RT-3D-PP) $(TIME) $(CART-ROUTINES) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-ELLIPSECLUMP) $(INTEGRATE)  $(OUTPUT) $(EVOLVE) $(CAPREOLE) 
	$(F90) $(OPTIONS) -o $@ $(STRINGS) $(NOMPI) $(MESH) $(CART-COORDS) $(HYDRO) $(TIME) $(CART-ROUTINES) $(HCOOLING) $(RT-3D-PP) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(NO_IONIC) $(CART-ELLIPSECLUMP) $(INTEGRATE) $(OUTPUT) $(EVOLVE) $(CAPREOLE) $(LIBS)

cart-ellipseclump_c2ray : $(PRECISION) $(FILES) $(STRINGS) $(SIZES) $(SCALING) $(CGSCONS) $(CGSASTROCONS) $(CGSPHOTOCONS) $(ABUNDANCES) $(ATOMIC) $(NOMPI) $(MESH) $(CART-COORDS) $(HYDRO) $(COOLING) $(RT-3D-PP) $(TIME) $(CART-ROUTINES) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-ELLIPSECLUMP) $(INTEGRATE)  $(OUTPUT) $(EVOLVE) $(CAPREOLE) 
	$(F90) $(OPTIONS) -o $@ $(STRINGS) $(NOMPI) $(MESH) $(CART-COORDS) $(HYDRO) $(TIME) $(CART-ROUTINES) $(COOLING) $(RT-3D-PP) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(NO_IONIC) $(CART-ELLIPSECLUMP) $(INTEGRATE) $(OUTPUT) $(EVOLVE) $(CAPREOLE) $(LIBS)

clean:
	rm -f *.o *.mod *.l *.il *.vo


