# Makefile for MPI hydro
#F90 = pf90
#F90 = f95
#F90 = mpif77
F90 = ifort
#F90 = mpif90

LDR     = $(F90)
#PP = cpp -P

# F90 options
IFORTFLAGS = -xW -O3 -vec_report -u -ipo -fpe0
#IFORTFLAGS = -assume 2underscores -xW -O3 -vec_report -u -ipo -fpe0
#F90FLAGS = -O3 -cpu:opteron #-DMPI
#F90FLAGS = -keep -Wv,-Pprocs,1 -O3 -cpu:opteron #-DMPI
F90FLAGS = $(IFORTFLAGS)
#F90FLAGS = $(IFORTFLAGS) -DMPI 

OPTIONS = $(F90FLAGS)

LDFLAGS = $(OPTIONS) 
LIBS = #-lU77

# list of objects we're using

PRECISION = precision.o

FILES = file_admin.o

MESH = mesh.o

CART-COORDS = cart-coords.o

CART-ROUTINES = cart-routines.o

PROT = protection2.o

LOF = lof-2.2.o

BOUNDARY = boundary.o

ROESOL = roesol-adv.o

# CART problems

CART-CONSTANT = cart-constant.o

CART-ELLIPSECLUMP = cart-ellipseclump.o

CART-DENSITYFIELD = cart-densityfield.o

CART-ENRIQUE = cart-enrique_wind.o

CART-HALO = cart-halo.o

CART-MINIHALO = cart-minihalo.o

HYDRO = hydro.o

TIME = time.o

INTEGRATE = integrate-strang.o

EVOLVE = evolve.o

INPUT = ah3in.o

OUTPUT = ah3out.o

CAPREOLE = capreole.o

MPI = mpi.o

NOMPI = no_mpi.o

SIZES = sizes.o

SCALING = scaling.o

NOSCALING = noscaling.o

ATOMIC = atomic.o

ABUNDANCES = abundances.o

MATHCONS = mathconstants.o

CGSCONS = cgsconstants.o

CGSASTROCONS = cgsastroconstants.o

CGSPHOTOCONS = cgsphotoconstants.o

CONSTANTS = $(MATHCONS) $(CGSCONS) $(CGSASTROCONS)

COSMOPARMS = cosmoparms.o

STRINGS = string.o

NO_IONIC = no_ionic.o

# RT-3D
COOLING = cooling.o
HCOOLING = cooling_h.o
HMCOOLING = cooling_hm.o
RT-3D_BASIC = romberg.o tped.o radiation.o clumping.o doric.o thermal.o 
RT-3D_BASIC_ME = romberg.o tped.o radiation_me.o clumping.o doric.o thermal.o 
RT-3D = $(RT-3D_BASIC) sourceprops.o ionic.o
RT-3D-PP = $(RT-3D_BASIC) sourceprops_pp.o ionic_pp.o
RT-3D-PP_ME = $(RT-3D_BASIC_ME) sourceprops_pp.o ionic_pp.o
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

cart-ellipseclump : $(PRECISION) $(FILES) $(STRINGS) $(SIZES) $(SCALING) $(CONSTANTS) $(ATOMIC) $(NOMPI) $(MESH) $(CART-COORDS) $(HYDRO) $(TIME) $(CART-ROUTINES) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(NO_IONIC) $(CART-ELLIPSECLUMP) $(INTEGRATE)  $(OUTPUT) $(EVOLVE) $(CAPREOLE) 
	$(F90) $(OPTIONS) -o $@ $(STRINGS) $(NOMPI) $(MESH) $(CART-COORDS) $(HYDRO) $(TIME) $(CART-ROUTINES) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(NO_IONIC) $(CART-ELLIPSECLUMP) $(INTEGRATE) $(OUTPUT) $(EVOLVE) $(CAPREOLE) $(LIBS)

mpi_cart-ellipseclump : $(PRECISION) $(FILES) $(STRINGS) $(SIZES) $(SCALING) $(CONSTANTS) $(ABUNDANCES) $(ATOMIC) $(MPI) $(MESH) $(CART-COORDS) $(HYDRO) $(TIME) $(CART-ROUTINES) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(NO_IONIC) $(CART-ELLIPSECLUMP) $(INTEGRATE) $(OUTPUT) $(EVOLVE) $(CAPREOLE) 
	$(F90) $(OPTIONS) -o $@ $(MPI) $(STRINGS) $(MESH) $(CART-COORDS) $(HYDRO) $(TIME) $(CART-ROUTINES) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(NO_IONIC) $(CART-ELLIPSECLUMP) $(INTEGRATE)  $(OUTPUT) $(EVOLVE) $(CAPREOLE) $(LIBS)

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

cart-enrique_c2ray : $(PRECISION) $(FILES) $(STRINGS) $(SIZES) $(SCALING) $(MATHCONS) $(CGSCONS) $(CGSASTROCONS) $(CGSPHOTOCONS) $(ABUNDANCES) $(ATOMIC) $(NOMPI) $(INPUT) $(MESH) $(CART-COORDS) $(HYDRO) $(HMCOOLING) $(RT-3D) $(TIME) $(CART-ROUTINES) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-ENRIQUE) $(INTEGRATE)  $(OUTPUT) $(EVOLVE) $(CAPREOLE) 
	$(F90) $(OPTIONS) -o $@ $(STRINGS) $(NOMPI) $(INPUT) $(MESH) $(CART-COORDS) $(HYDRO) $(TIME) $(CART-ROUTINES) $(HMCOOLING) $(RT-3D) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-ENRIQUE) $(INTEGRATE) $(OUTPUT) $(EVOLVE) $(CAPREOLE) $(LIBS)

cart-constant_hc_c2ray : $(PRECISION) $(FILES) $(STRINGS) $(SIZES) $(SCALING) $(MATHCONS) $(CGSCONS) $(CGSASTROCONS) $(CGSPHOTOCONS) $(ABUNDANCES) $(ATOMIC) $(NOMPI) $(MESH) $(CART-COORDS) $(HYDRO) $(HMCOOLING) $(RT-3D) $(TIME) $(CART-ROUTINES) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-CONSTANT) $(INTEGRATE)  $(OUTPUT) $(EVOLVE) $(CAPREOLE) 
	$(F90) $(OPTIONS) -o $@ $(STRINGS) $(NOMPI) $(MESH) $(CART-COORDS) $(HYDRO) $(TIME) $(CART-ROUTINES) $(HMCOOLING) $(RT-3D) $(PROT) $(BOUNDARY) $(LOF) $(ROESOL) $(CART-CONSTANT) $(INTEGRATE) $(OUTPUT) $(EVOLVE) $(CAPREOLE) $(LIBS)

clean:
	rm -f *.o *.mod *.l *.il *.vo


