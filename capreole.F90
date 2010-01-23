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

  use precision, only: dp
  use my_mpi ! too many variables inherited from the MPI library
  !            to use 'only' here.
  use file_admin, only: stdinput, log_unit, file_input, flag_for_file_input
  use mesh, only: init_mesh
  use output, only: init_output
  use grid, only: init_coords
  use hydro, only: init_hydro
  use times, only: init_time
  use problem, only: init_problem
  use evolution, only: evolve

  implicit none

  ! Start and end time for CPU report
  real :: cputime1 !< Start time for CPU report
  real :: cputime2 !< End time for CPU report
  real(kind=dp) :: cpu_seconds=0.0
  integer :: cpu_hours=0
  integer :: cpu_minutes=0

  ! Wall clock time variables
  integer :: cntr1 !< Start time wall clock
  integer :: cntr2 !< End time wall clock
  integer :: countspersec !< counts per second (for wall clock time)
  real(kind=dp) :: clock_seconds=0.0
  integer :: clock_hours=0
  integer :: clock_minutes=0

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
  call cpu_time(cputime1)

  ! Initialize wall cock timer
  call system_clock(cntr1)

  ! Set up MPI structure
  call mpi_setup()

  ! Set up input stream (either standard input or from file given
  ! by first argument)
  if (rank == 0) then
     write(log_unit,*) "input or input?"
     flush(log_unit)
     if (COMMAND_ARGUMENT_COUNT () > 0) then
        call GET_COMMAND_ARGUMENT(1,inputfile)
        write(log_unit,*) "reading input from ",trim(adjustl(inputfile))
        open(unit=stdinput,file=inputfile)
        call flag_for_file_input(.true.)
     else
        write(log_unit,*) "reading input from command line"
     endif
     flush(log_unit)
  endif

  ! Ask if this is a restart
  if (rank == 0) then
     if (.not.file_input) write (*,"(A,$)") "Restart of old run: (y/n): "
     read (unit=stdinput,fmt="(A1)") answer
     if (answer == "y" .or. answer == "Y") then
        restart=.true.
        if (.not.file_input) &
             write (*,"(A,$)") "Output file from which to restart: "
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

  ! Find out CPU time
  call cpu_time(cputime2)
  cpu_seconds=cpu_seconds+real(cputime2-cputime1,dp)
  cpu_minutes = cpu_minutes + int(cpu_seconds) / 60
  cpu_seconds = MOD ( cpu_seconds , 60.0 )
  cpu_hours = cpu_hours + cpu_minutes / 60
  cpu_minutes = MOD ( cpu_minutes , 60 )

  ! Find out wall clock time
  call system_clock(cntr2,countspersec)
  clock_seconds=clock_seconds+real(cntr2-cntr1,dp)/real(countspersec,dp)
  clock_minutes = clock_minutes + int(clock_seconds) / 60
  clock_seconds = MOD ( clock_seconds , 60.0 )
  clock_hours = clock_hours + clock_minutes / 60
  clock_minutes = MOD ( clock_minutes , 60 )
  
  if (rank == 0) then
     write(log_unit,*) "CPU time: ",cpu_hours," hours",cpu_minutes," minutes", &
          cpu_seconds," seconds."
     write(log_unit,*) "Wall clock time: ",clock_hours," hours", &
          clock_minutes," minutes",clock_seconds," seconds."
  endif

  ! End the run
  call mpi_end()

end Program Capreole







