Program Capreole
  
  ! Author: Garrelt Mellema
  ! version 2004.05.11 (F90)
  ! (previous version 2003.06.01) 
  ! (previous version 2001.10.01)
  ! (descends from F77 version 2001.03.21)

  ! This is general surpose hydrodynamics programme in Fortran 90/95.
  ! It can be compiled for use in an MPI environment (MPICH).
  ! For this use compiler option -DMPI, and make sure you use the
  ! routines in mpi.f90. Without MPI use no_mpi.f90.
  !
  ! All initial/boundary conditions should be supplied in separate routines, 
  ! as well as some other functions, noted in the programme.
  ! Modules needed:
  ! atomic
  ! boundary
  ! evolution
  ! geometry
  ! grid
  ! hydro
  ! hydrosolver
  ! integrator
  ! mesh
  ! my_mpi
  ! output
  ! precision
  ! problem
  ! protection
  ! scaling
  ! sizes
  ! times

  ! This version avoids passing large dynamic arrays as subroutine arguments
  ! since this tends to fill up the stack space (and on some machines
  ! 'limit stacksize unlimited' does not help). Instead most arrays are 
  ! defined as allocatable in modules (e.g. hydro) and allocated once the 
  ! mesh size is known. In order to use generic names in subroutines, 
  ! pointers are used.
  !----------------------------------------------------------------------------

  use my_mpi
  use mesh
  use output
  use grid
  use hydro
  use times
  use problem
  use evolution

  implicit none

  ! Start and end time for CPU report
  real :: tstart,tend

  ! Restart (disabled)
  logical :: restart=.false.

  ! Error flag
  integer :: ierror

  !----------------------------------------------------------------------------

  call cpu_time(tstart)

  call mpi_setup()

  call init_mesh()

  ! Initialize output routines
  call init_output(restart,ierror)

  ! Initialize the coordinate system
  call init_coords(restart)
  
  ! Initialize the hydrodynamic variables
  call init_hydro(restart) ! allocates arrays
  call init_problem(restart) ! sets initial conditions

  ! Initialiaze time variables
  call init_time(restart)

  ! Evolve the problem
  call evolve()

  ! End the run
  call mpi_end()

  ! Find out CPU time
  call cpu_time(tend)

  write(*,*) 'CPU time: ',tend-tstart,' s'

end Program Capreole







