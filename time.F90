module times

  ! Module for Capreole (f90)
  ! Author: Garrelt Mellema
  ! Date: 2004-05-11
  ! This module needs to be checked for the F compiler

  ! This module contains the routine which allocates the hydro
  ! arrays. This way they do not end up on stack
  
  use precision
  use sizes
  use mesh
  use my_mpi
  use scaling

  implicit none
  private

  real(kind=dp),public :: frametime
  integer,public       :: LastFrame
  real(kind=dp),public :: time
  real(kind=dp),public :: dt

  public :: init_time

contains

  subroutine init_time (restart)
    
    ! This routine initializes all time variables
    
    ! This may be a fresh start or a restart of a saved run
    
    logical,intent(in) :: restart
    
    if (.not.restart) then ! Fresh start
       
       ! Ask for the input if you are processor 0.
       
       if (rank == 0) then
          write (*,'(//,A,/)') '----- Output -----'
          write (*,'(A,$)') '1) Every how many simulation seconds: '
          read (*,*) frametime
          write (*,'(A,$)') '2) How many output frames: '
          read (*,*) LastFrame
       endif
       
       ! report input parameters
       if (rank == 0) then
          write (30,'(//,A,/)') '----- Output -----'
          write (30,'(A,1PE10.3)') '1) Every how many simulation seconds: ', & 
               frametime
          write (30,'(A,I4)') '2) Number of output frames: ',LastFrame
       endif
       
#ifdef MPI       
       ! Distribute the input parameters to the other nodes
       call MPI_BCAST(frametime,1,MPI_DOUBLE_PRECISION,0,&
            MPI_COMM_NEW,ierror)
       call MPI_BCAST(LastFrame,1,MPI_INTEGER,0,MPI_COMM_NEW,ierror)
#endif
       
       ! Scale time
       frametime=frametime/sctime
       
       ! Initialize the time and time step to zero
       time=0.0d0
       dt=0.0d0
    endif

  end subroutine init_time

end module times
