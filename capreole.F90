Program Capreole
  
  ! Author: Garrelt Mellema
  ! version 2007.10.22
  ! (previous 2004.05.11)
  ! (previous version 2003.06.01) 
  ! (previous version 2001.10.01)
  ! (descends from F77 version 2001.03.21)

  ! This is general surpose hydrodynamics programme in Fortran 90/95.
  ! It can be compiled for use in an MPI environment (MPICH).
  ! For this use compiler option -DMPI, and make sure you use the
  ! routines in mpi.f90. For compilation without MPI do not use
  ! the -DMPI option and use no_mpi.f90.
  !
  ! All initial/boundary conditions should be supplied in separate routines, 
  ! as well as some other functions, noted in the programme.

  ! This version avoids passing large dynamic arrays as subroutine arguments
  ! since this tends to fill up the stack space (and on some machines
  ! 'limit stacksize unlimited' does not help). Instead most arrays are 
  ! defined as allocatable in modules (e.g. hydro) and allocated once the 
  ! mesh size is known. In order to use generic names in subroutines, 
  ! pointers are used.
  !----------------------------------------------------------------------------

  use my_mpi ! too many variables inherited from the MPI library
  !            to use 'only' here.
  use file_admin, only: stdinput, log_unit
  use mesh, only: init_mesh
  use output, only: init_output
  use grid, only: init_coords
  use hydro, only: init_hydro
  use times, only: init_time
  use problem, only: init_problem
  use evolution, only: evolve

  implicit none

  ! Start and end time for CPU and wall clock report
  real :: tstart,tend
  integer :: cntr1,cntr2,countspersec

  ! Restart variables
  logical :: restart!=.false.
  character(len=19) :: restartfile
  character(len=1) :: answer

  ! Error flag
  integer :: ierror

  ! iargc library function (number of command line arguments)
  integer :: iargc

  ! File with input values
  character(len=512) :: inputfile

  !----------------------------------------------------------------------------

  ! Initialize cpu timer
  call cpu_time(tstart)

  ! Initialize wall cock timer
  call system_clock(cntr1)

  ! Set up MPI structure
  call mpi_setup()

  ! Set up input stream (either standard input or from file given
  ! by first argument)
  if (iargc() > 0) then
     call getarg(1,inputfile)
     if (rank == 0) then
        write(30,*) 'reading input from ',trim(adjustl(inputfile))
        open(unit=stdinput,file=inputfile)
     endif
  endif

  ! Ask if this is a restart
  if (rank == 0) then
     write (*,'(A,$)') 'Restart of old run: (y/n): '
     read (unit=stdinput,fmt='(A1)') answer
     if (answer == 'y' .or. answer == 'Y') then
        restart=.true.
        write (*,'(A,$)') 'Output file from which to restart: '
        read (unit=stdinput,fmt=*) restartfile
     else
        restart=.false.
     endif
  endif
#ifdef MPI
    ! Tell the other nodes
    call MPI_BCAST(restart,1,MPI_LOGICAL,0,MPI_COMM_NEW,ierror)
    call MPI_BCAST(restartfile,19,MPI_CHARACTER,0,MPI_COMM_NEW,ierror)
#endif

  ! Initialize computational mesh
  call init_mesh(restart,restartfile)

  ! Initialize output routines
  call init_output(restart,restartfile,ierror)

  ! Initialize the coordinate system (grid)
  call init_coords(restart,restartfile)
  
  ! Initialize the hydrodynamic variables
  call init_hydro( ) ! allocates arrays
  call init_problem(restart,restartfile) ! sets initial conditions

  ! Initialiaze time variables
  call init_time(restart,restartfile)

  ! Evolve the problem
  call evolve()

  ! Find out CPU and wall clock time
  call cpu_time(tend)
  call system_clock(cntr2,countspersec)

  ! Report times
  write(log_unit,*) 'CPU time: ',tend-tstart,' s'
  write(log_unit,*) 'Wall clock time: ',(cntr2-cntr1)/countspersec,' s'

  ! End the run
  call mpi_end()

end Program Capreole







