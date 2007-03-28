module evolution

  ! Module for Capreole
  ! Author: Garrelt Mellema
  ! Date: 2004-05-11

  ! This module contains the routines for the framework of the
  ! hydrodynamics calculation.

  ! Contents:
  ! evolve - basic integration frame work

  use precision
  use scaling
  use sizes
  use mesh
  use grid
  use hydro
  use times
  use output
  use geometry
  use integrator

  implicit none

  real(kind=dp),parameter :: cfl=0.4d0  ! CFL number for time step (<1.0)

contains

  subroutine evolve ()
    
    ! This routine handles the basic integration frame work

    integer :: nstep,inegative,istop ! various control
                                                    !   integers
    real(kind=dp)    :: dtlocal    ! processor local time step
    real(kind=dp)    :: nexttime   ! timer for output
    integer :: nframe              ! integers for output
#ifdef MPI
    integer :: ierror
#endif

    !--------------------------------------------------------------------------

    ! Report initial conditions
    nframe=int(time/frametime) ! initialize output counter
    if (nframe == 0) call make_output(nframe)

    ! Set time for next output
    nexttime=real(nframe+1,dp)*frametime
    nframe=nframe+1    ! update output counter

    nstep=0            ! timestep counter
    istop=0            ! initialize stop flag
    inegative=0        ! initialize negative density/energy flag
    

    ! Integration loop
    do
       nstep=nstep+1 ! count the integration loops
       
       if (mod(nstep,2) == 0) then
          stold => state2   ! copy new state to old state
          stnew => state1   ! copy new state to old state
       else
          stold => state1   ! copy new state to old state
          stnew => state2   ! copy new state to old state
       endif
       !stold(:,:,:,:)=stnew(:,:,:,:)   ! copy new state to old state
       
       ! Determine the time step on the local grid;
       ! the timestep function has to be supplied
       state => stold
       dtlocal=timestep(cfl,OLD)
       dt=dtlocal

#ifdef MPI
       ! communicate with all other processors to find
       ! the global minimal time step
       call MPI_ALLREDUCE(dtlocal,dt,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
            MPI_COMM_NEW,ierror)
#endif
       ! Make sure that dt does not take us beyond nexttime
       dt=min(nexttime-time,dt)

       ! Integrate one time step
       call integrate(nstep,istop)

       time=time+dt          ! update time

       ! Report time step info
       write(30,'(A,I6,3(X,1PE10.3))') 'Time info: ',&
            nstep,time*sctime,dt*sctime,nexttime*sctime
       call flush(30)

       ! Check the need for new output
       if (time >= nexttime) then
          state => stnew
          call make_output(nframe)
          nframe=nframe+1
          nexttime=nexttime+frametime
       endif

       if (nframe.gt.LastFrame.or.istop.ne.0) exit ! end the integration loop
    enddo

    if (istop.ne.0) then
       ! Record conditions where error occured
       write(30,*) 'Stop triggered by correction routine'
       state => stnew
       call make_output(nframe)
    else
       write(30,*) 'Maximum number of frames reached; nframe = ',nframe
    endif

  end subroutine evolve

end module evolution
