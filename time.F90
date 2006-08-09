module times

  ! Module for Capreole (f90)
  ! Author: Garrelt Mellema
  ! Date: 2004-05-11
  ! This module needs to be checked for the F compiler

  ! This module contains the routine which allocates the hydro
  ! arrays. This way they do not end up on stack
  
  use file_admin, only: stdinput
  use precision
  use sizes
  use mesh
  use my_mpi
  use scaling
  use string_manipulation
  use astroconstants

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
    character(len=10) :: str_time_unit
#ifdef MPI       
    integer :: ierror
#endif

    if (.not.restart) then ! Fresh start
       
       ! Ask for the input if you are processor 0.
       
       if (rank == 0) then
          write (*,'(//,A,/)') '----- Output -----'
          write (*,'(A,$)') '1) Time between outputs (specify unit): '
          read (stdinput,*) frametime,str_time_unit
          write (*,'(A,$)') '2) How many output frames: '
          read (stdinput,*) LastFrame
       endif
       
       ! report input parameters
       if (rank == 0) then
          write (30,'(//,A,/)') '----- Output -----'
          write (30,'(A,1PE10.3,A)') '1) Time between outputs: ', & 
               frametime,str_time_unit
          write (30,'(A,I4)') '2) Number of output frames: ',LastFrame
          ! Convert to seconds
          call convert_case(str_time_unit,0) ! conversion to lower case
          select case (trim(adjustl(str_time_unit)))
          case ('s','sec','second','secs','seconds')
          case ('years','yrs','year','yr')
             frametime=frametime*YEAR
          case ('myears','myrs','myear','myr')
             frametime=frametime*1e6*YEAR
          case ('gyears','gyrs','gyear','gyr')
             frametime=frametime*1e9*YEAR
          case default
             write(*,*) 'Time unit not recognized, assuming seconds'
          end select
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
