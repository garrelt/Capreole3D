module problem

  ! Module for Capreole2D (F90)
  ! Author: Garrelt Mellema
  ! Date: 2004-05-11

  ! This is the problem module. It contains the routines which define the
  ! problem being solved:

  ! Contents:
  ! init_problem - sets up the hydro variables according to the specified
  !                  problem
  ! inflow - sets the inflow boundary conditions

  ! This version: read in density field (cartesian coordinates).

  use file_admin, only: stdinput
  use precision, only: dp
  use my_mpi
  use sizes
  use scaling
  use atomic
  use cgsconstants
  use abundances, only: mu
  use mesh
  use grid
  use hydro
  use boundary
  use ionic

  implicit none

  real(kind=dp),private :: sedensity,sevelocity,sepressr

contains

  subroutine init_problem (restart)
    
    ! This routine initializes all hydro variables
    
    ! This may be a fresh start or a restart of a saved run
    
    ! Case: steady shock interacting with elliptical cloud; (x,y)
    
    ! smoothing parameter for making soft-edged clump
    real(kind=dp),parameter :: eta=0.1_dp
    
    logical,intent(in) :: restart ! tells you whether it's a new run or a restart
    ! e variables are environment
    real(kind=dp)    :: etemperature
    real(kind=dp)    :: epressr

    integer :: mx1,my1,mz1,temperature_setup
    character(len=512) :: densityfield_file
    real,allocatable,dimension(:,:,:) :: tempdens
    real(kind=dp) :: maxdens
    integer :: i,j,k,nitt,ieq
    real(kind=dp)    :: r_interface ! dummy needed for calling init_ionic

#ifdef MPI       
    integer :: status(MPI_STATUS_SIZE)
    integer :: request,nextproc
    integer,parameter :: outputcircle=601
    integer :: ierror
#endif

    if (.not.restart) then ! Fresh start
       
       ! Ask for the input if you are processor 0.
       if (rank == 0) then
          write (*,'(A,$)') '1) File with density (cm^-3): '
          read (stdinput,*) densityfield_file
          write (*,'(A,$)') '2) Temperature setup: '
          write (*,'(A,$)') '   1: isobaric: '
          write (*,'(A,$)') '   2: isothermal: '
          read (stdinput,*) temperature_setup
          if (temperature_setup == 1) then
             write (*,'(A,$)') '3) Minimum temperature: '
             read (stdinput,*) etemperature
          else
             write (*,'(A,$)') '3) Temperature: '
             read (stdinput,*) etemperature
          endif
       endif

       ! report input parameters
       if (rank == 0) then
          write(30,'(A)') & 
               'Problem: cart-density field (cartesian: density field)'
          write (30,'(2A)') '1) Density file (cm^-3): ', &
               densityfield_file
          write (30,'(A)') '2) Temperature setup:'
          if (temperature_setup == 1) then
             write (30,'(A)') 'isobaric'
             write (30,'(A,1PE10.3)') '3) Minimum Temperature: ',etemperature
          else
             write (30,'(A)') '  isothermal'
             write (30,'(A,1PE10.3)') '3) Temperature: ',etemperature
          endif
       endif
          
#ifdef MPI       
       ! Distribute the input parameters to the other nodes
       ! (densityfield_file is distributed in the MPI ring below)
       call MPI_BCAST(temperature_setup,1,MPI_INTEGER,0, & 
            MPI_COMM_NEW,ierror)
       call MPI_BCAST(etemperature,1,MPI_DOUBLE_PRECISION,0, & 
            MPI_COMM_NEW,ierror)
#endif
       
       ! Read in density field (MPI ring)
       if (rank.eq.0) then
          
          open(unit=40,file=densityfield_file,form='UNFORMATTED',status='OLD')
          read(40) mx1,my1,mz1
          write(30,*)  mx1,my1,mz1
          if (mx1 .ne. meshx .or. my1.ne.meshy .or. mz1.ne.meshz) then
             write(30,*) 'Error: ', &
                  'density field does not match specified mesh size'
             stop
          endif
          allocate(tempdens(1-mbc:mx1+mbc,1-mbc:my1+mbc,1-mbc:mz1+mbc))
          read(40) tempdens(1:mx1,1:my1,1:mz1)
          close(40)
       
#ifdef MPI
          ! Send filename to next processor
          if (npr > 1) then
             call MPI_ISSEND(densityfield_file,512,MPI_CHARACTER,rank+1, &
                  outputcircle,MPI_COMM_NEW,request,ierror)
             write(30,*) rank,'sent filename to ',rank+1
             ! Wait for the circle to complete
             call MPI_RECV(densityfield_file,512,MPI_CHARACTER,npr-1, &
                  outputcircle,MPI_COMM_NEW,status,ierror)
             write(30,*) rank,'received filename from ',npr-1
             call flush(30)
          endif

       else
       
          ! Receive filename from previous processor
          call MPI_RECV(densityfield_file,512,MPI_CHARACTER,rank-1, &
               outputcircle,MPI_COMM_NEW,status,ierror)
          write(30,*) rank,'received filename from ',rank-1
          call flush(30)
          
          if (ierror == 0) then ! if ok
             open(unit=40,file=densityfield_file,form='UNFORMATTED', &
                  status='OLD')
             read(40) mx1,my1,mz1
             write(30,*)  mx1,my1,mz1
             if (mx1 .ne. meshx .or. my1.ne.meshy .or. mz1.ne.meshz) then
                write(30,*) 'Error: ', &
                     'density field does not match specified mesh size'
                stop
             endif
             allocate(tempdens(1-mbc:mx1+mbc,1-mbc:my1+mbc,1-mbc:mz1+mbc))
             read(40) tempdens(1:mx1,1:my1,1:mz1)
             close(40)

             ! Determine next processor (npr -> 0)
             nextproc=mod(rank+1,npr)
             ! Send filename along
             call MPI_ISSEND(densityfield_file,512,MPI_CHARACTER,nextproc, &
                  outputcircle,MPI_COMM_NEW,request,ierror)
             write(30,*) rank,'sent filename to ',nextproc
          else
             write(30,*) 'MPI error in output routine'
             ierror=ierror+1
          endif
#endif
       endif

       call flush(30)
       ! Assume outflow boundaries
       do k=1-mbc,mz1+mbc
          do j=1-mbc,my1+mbc
             do i=1-mbc,0
                tempdens(i,j,k)=tempdens(1,j,k)
             enddo
          enddo
       enddo
       do k=1-mbc,mz1+mbc
          do j=1-mbc,0
             do i=1-mbc,mx1+mbc
                tempdens(i,j,k)=tempdens(i,1,k)
             enddo
          enddo
       enddo
       do k=1-mbc,0
          do j=1-mbc,my1+mbc
             do i=1-mbc,mx1+mbc
                tempdens(i,j,k)=tempdens(i,j,1)
             enddo
          enddo
       enddo
       
       do k=1-mbc,mz1+mbc
          do j=1-mbc,my1+mbc
             do i=1-mbc,mx1+mbc
                if (tempdens(i,j,k) == 0.0) tempdens(i,j,k)=tempdens(1,1,1)
             enddo
          enddo
       enddo
       ! Find maximum density (for isobaric initial conditions)
       write(30,*) maxval(tempdens),minval(tempdens),scdens
       call flush(30)
       maxdens=maxval(tempdens)*mu*m_p/scdens 

       ! Set density in state array
       ! Scale to program scaling
       state(sx-mbc:ex+mbc,sy-mbc:ey+mbc,sz-mbc:ez+mbc,RHO)= &
            tempdens(sx-mbc:ex+mbc,sy-mbc:ey+mbc,sz-mbc:ez+mbc)* &
            mu*m_p/scdens ! densities are in cm^-3, 
       deallocate(tempdens)
                
       ! Set the initial momenta
       do k=sz-mbc,ez+mbc
          do j=sy-mbc,ey+mbc
             do i=sx-mbc,ex+mbc
                state(i,j,k,RHVX)=0.0d0
                state(i,j,k,RHVY)=0.0d0
                state(i,j,k,RHVZ)=0.0d0
             enddo
          enddo
       enddo
       
       write(30,*) 'setting pressure'
       write(30,*) minval(state(:,:,:,RHO))
       call flush(30)
       ! Set the initial pressure and energy density
       ! (depends on the temperature setup, isobaric or isothermal)
       if (temperature_setup == 1) then
          epressr=maxdens*boltzm*etemperature/mu
          do k=sz-mbc,ez+mbc
             do j=sy-mbc,ey+mbc
                do i=sx-mbc,ex+mbc
                   pressr(i,j,k)=epressr
                   state(i,j,k,EN)=pressr(i,j,k)/gamma1+ &
                        0.5d0*(state(i,j,k,RHVX)*state(i,j,k,RHVX)+ &
                        state(i,j,k,RHVY)*state(i,j,k,RHVY)+ &
                        state(i,j,k,RHVZ)*state(i,j,k,RHVZ))/state(i,j,k,RHO)
                   !state(i,j,k,TRACER1)=1.0d0
                enddo
             enddo
          enddo
       else
          do k=sz-mbc,ez+mbc
             do j=sy-mbc,ey+mbc
                do i=sx-mbc,ex+mbc
                   pressr(i,j,k)=state(i,j,k,RHO)*boltzm*etemperature/mu
                   state(i,j,k,EN)=pressr(i,j,k)/gamma1+ &
                        0.5d0*(state(i,j,k,RHVX)*state(i,j,k,RHVX)+ &
                        state(i,j,k,RHVY)*state(i,j,k,RHVY)+ &
                        state(i,j,k,RHVZ)*state(i,j,k,RHVZ))/state(i,j,k,RHO)
                   !state(i,j,k,TRACER1)=1.0d0
                enddo
             enddo
          enddo
       endif
          
       ! Initialize the ionic concentrations
       call init_ionic(restart,r_interface)

       ! Record variables which remain constant during a run in a file
       ! runparams to be used at restarts
       
       if (rank == 0) then
          open(unit=80,file='runparams', &
               status='unknown',form='unformatted')
          write(80) 0,0
              
          ! Disabled until I've figured out how to do this
          ! distributed
          
          !!write(80) dx,dy,frametime
          !!write(80) sedensity,sevelocity,sepressr
          !!write(80) (x(i),i=1-mbc,meshx+mbc)
          !!write(80) (y(j),j=1-mbc,meshy+mbc)
          !!write(80) ((volx(i,j),i=1-mbc,meshx+mbc), &
          !!            j=1-mbc,meshy+mbc)
              !!write(80) ((voly(i,j),i=1-mbc,meshx+mbc), &
          !!            j=1-mbc,meshy+mbc)
       endif
    else
       write(30,*) 'No restart implemented yet!'
    endif

  end subroutine init_problem

  !==========================================================================

  subroutine inflow (newold)
    
    ! This routine resets the inner boundary to the inflow condition
    
    integer,intent(in) :: newold

    integer :: i,j,k

    ! Point state to appropriate array
    state => set_state_pointer(newold)

    if (sx == 1 .and. sevelocity > 0.0) then
       do k=sz-1,ez+1
          do j=sy-1,ey+1
             do i=1-mbc,1
                state(i,j,k,RHO)=sedensity
                state(i,j,k,RHVX)=sedensity*sevelocity
                state(i,j,k,RHVY)=0.0d0
                state(i,j,k,RHVZ)=0.0d0
                pressr(i,j,k)=sepressr
                   state(i,j,k,EN)=pressr(i,j,k)/gamma1+ &
                        0.5d0*(state(i,j,k,RHVX)*state(i,j,k,RHVX)+ &
                        state(i,j,k,RHVY)*state(i,j,k,RHVY)+ &
                        state(i,j,k,RHVZ)*state(i,j,k,RHVZ))/state(i,j,k,RHO)
                !state(i,j,k,TRACER1)=-1.0d0
             enddo
          enddo
       enddo
    endif
    
  end subroutine inflow
  
end module problem
