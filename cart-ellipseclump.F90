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

  ! This version: steady shock interacting with elliptical cloud
  !  (cartesian coordinates).

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
    ! w variables are clump 
    ! se variables are shocked environment
    real(kind=dp)    :: edensity,etemperature,vblast
    real(kind=dp)    :: wdensity,wtemperature,x0,y0,z0,axis
    real(kind=dp)    :: axis1,axis2,axis3,rotangle1,rotangle2
    real(kind=dp)    :: epressr,wpressr,vs1,xm1,fp,shckdist
    real(kind=dp)    :: xc,yc,zc,xr1,yr1,xr2,zr2

    integer :: i,j,k,nitt,ieq
    real(kind=dp)    :: r_interface ! dummy needed for calling init_ionic

#ifdef MPI       
    integer :: ierror
#endif

    if (.not.restart) then ! Fresh start
       
       ! Ask for the input if you are processor 0.
       if (rank == 0) then
          write (*,'(//,A,/)') '----- Environment -----'
          write (*,'(A,$)') '1) Density (cm^-3): '
          read (stdinput,*) edensity
          write (*,'(A,$)') '2) Temperature: '
          read (stdinput,*) etemperature
          write (*,'(A,$)') '3) Blast velocity (km/s): '
          read (stdinput,*) vblast
          
          write (*,'(//,A,/)') '----- Knot -----'
          write (*,'(A,$)') '1) Density (cm^-3): '
          read (stdinput,*) wdensity
          write (*,'(A,$)') '2) Temperature: '
          read (stdinput,*) wtemperature
          write (*,'(A,$)') '3) Position of centre x,y,z (cm): '
          read (stdinput,*) x0,y0,z0
          write (*,'(A,$)') '4) Axis 1: '
          read (stdinput,*) axis1
          write (*,'(A,$)') '5) Axis 2: '
          read (stdinput,*) axis2
          write (*,'(A,$)') '6) Axis 3: '
          read (stdinput,*) axis3
          write (*,'(A,$)') '7) Rotation angle1: '
          read (stdinput,*) rotangle1
          write (*,'(A,$)') '8) Rotation angle2: '
          read (stdinput,*) rotangle2
       endif

       ! report input parameters
       if (rank == 0) then
          write(30,'(A)') & 
               'Problem: cart-ellipseclump (cartesian shock - elliptical cloud interaction)'
          write (30,'(//,A,/)') '----- Environment -----'
          write (30,'(A,1PE10.3)') '1) Density (cm^-3): ',edensity
          write (30,'(A,1PE10.3)') '2) Temperature: ',etemperature
          write (30,'(A,F8.3)') '3) Blast velocity (km/s): ',vblast
          
          write (30,'(//,A,/)') '----- Knot -----'
          write (30,'(A,1PE10.3)') '1) Density (cm^-3): ',wdensity
          write (30,'(A,1PE10.3)') '2) Temperature: ',wtemperature
          write (30,'(A,3(1PE10.3,X))') '3) Position of centre x,y (cm): ',&
               x0,y0,z0
          write (30,'(A,1PE10.3)') '4) Axis 1: ',axis1
          write (30,'(A,1PE10.3)') '5) Axis 2: ',axis2
          write (30,'(A,1PE10.3)') '6) Axis 3: ',axis3
          write (30,'(A,F8.3)') '7) Rotation angle1: ',rotangle1
          write (30,'(A,F8.3)') '8) Rotation angle2: ',rotangle2
       endif
          
#ifdef MPI       
       ! Distribute the input parameters to the other nodes
       call MPI_BCAST(edensity,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
            ierror)
       call MPI_BCAST(etemperature,1,MPI_DOUBLE_PRECISION,0, & 
            MPI_COMM_NEW,ierror)
       call MPI_BCAST(vblast,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
            ierror)
       call MPI_BCAST(wdensity,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
            ierror)
       call MPI_BCAST(wtemperature,1,MPI_DOUBLE_PRECISION,0, &
            MPI_COMM_NEW,ierror)
       call MPI_BCAST(x0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
            ierror)
       call MPI_BCAST(y0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
            ierror) 
       call MPI_BCAST(z0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
            ierror) 
       call MPI_BCAST(axis1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
            ierror)
       call MPI_BCAST(axis2,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
            ierror)
       call MPI_BCAST(axis3,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
            ierror)
       call MPI_BCAST(rotangle1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
            ierror)
       call MPI_BCAST(rotangle2,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
            ierror)
#endif
       
       ! Scale physical parameters

       vblast=vblast*1d5/scvelo ! vblast is in km/s, not cm/s
       edensity=mu*m_p*edensity/scdens ! densities are in cm^-3, 
       wdensity=mu*m_p*wdensity/scdens ! 
       
       ! Calculate the pressures
       
       epressr=edensity*boltzm*etemperature/xmu
       wpressr=wdensity*boltzm*wtemperature/xmu
       wpressr=epressr       ! Pressure equilibrium!!!!
       
       ! Calculate the properties of the shock wave 
       vs1=sqrt(gamma*epressr/edensity) ! Sound speed in unshocked gas
       xm1=vblast/vs1                   ! Mach numnber of shock
       
       write(30,*) 'Pressure= ',epressr*scener
       write(30,*) 'Density= ',edensity*scdens
       write(30,*) 'Mach number= ',xm1
       
       ! Hugoniot conditions for density and pressure jumps
       
       if (xm1 <= 1.0) then
          sedensity=edensity
          sevelocity=0.0
          sepressr=epressr
       else
          sedensity=(gamma+1.0d0)*xm1*xm1/((gamma-1.0d0)*xm1*xm1+2.0d0)* &
               edensity
          sevelocity=vblast*(1.0d0-edensity/sedensity) ! post shock velocity
          sepressr=epressr*(2.0d0*gamma*xm1*xm1/(gamma+1.0d0)-(gamma-1.0d0)/ & 
               (gamma+1.0d0))
       endif

       ! Scale the geometrical parameters of the clump
       x0=x0/scleng           ! centre position along z-axis 
       y0=y0/scleng           ! centre position along z-axis 
       z0=z0/scleng           ! centre position along z-axis 
       axis1=axis1/scleng     ! radius; both are in cm's 
       axis2=axis2/scleng     ! radius; both are in cm's 
       axis3=axis3/scleng     ! radius; both are in cm's 
       axis=max(axis1,axis2,axis3)
       rotangle1=rotangle1/180.0*acos(-1.0)
       rotangle2=rotangle2/180.0*acos(-1.0)

       ! Calculate the initial state
       ! The clump is elliptical, we calculate the ellipse in fp and set
       ! the clump properies if we are inside and the unshocked environment
       ! if outdside
       
       do k=sz-mbc,ez+mbc
          do j=sy-mbc,ey+mbc
             do i=sx-mbc,ex+mbc
                xc=x(i)-x0
                yc=y(j)-y0
                zc=z(k)-z0

                xr1 =  cos(rotangle1)*xc+sin(rotangle1)*yc
                yr1 = -sin(rotangle1)*xc+cos(rotangle1)*yc
                 
                xr2 = cos(rotangle2)*xr1+sin(rotangle2)*zc
                zr2 = -sin(rotangle2)*xr1+cos(rotangle2)*zc

                fp=sqrt(xr2*xr2/(axis1*axis1)+yr1*yr1/(axis2*axis2)+ &
                     zr2*zr2/(axis3*axis3))

                if (fp <= 1.0) then
                   state(i,j,k,RHO)=wdensity
                   state(i,j,k,RHVX)=0.0d0
                   state(i,j,k,RHVY)=0.0d0
                   state(i,j,k,RHVZ)=0.0d0
                   pressr(i,j,k)=wpressr
                   state(i,j,k,EN)=pressr(i,j,k)/gamma1+ &
                        0.5d0*(state(i,j,k,RHVX)*state(i,j,k,RHVX)+ &
                        state(i,j,k,RHVY)*state(i,j,k,RHVY)+ &
                        state(i,j,k,RHVZ)*state(i,j,k,RHVZ))/state(i,j,k,RHO)
                   state(i,j,k,TRACER1)=1.0d0
                else
                   state(i,j,k,RHO)=edensity
                   state(i,j,k,RHVX)=0.0d0
                   state(i,j,k,RHVY)=0.0d0
                   state(i,j,k,RHVZ)=0.0d0
                   pressr(i,j,k)=epressr
                   state(i,j,k,EN)=pressr(i,j,k)/gamma1+ &
                        0.5d0*(state(i,j,k,RHVX)*state(i,j,k,RHVX)+ &
                        state(i,j,k,RHVY)*state(i,j,k,RHVY)+ &
                        state(i,j,k,RHVZ)*state(i,j,k,RHVZ))/state(i,j,k,RHO)
                   state(i,j,k,TRACER1)=-1.0d0
                endif
             enddo
          enddo
       enddo
          
       ! To avoid numerical problems the edge of the clump needs to be
       ! smoothed. We achieve this by diffusing the state variables a number
       ! of time with a high diffuse parameter eta. 
       
       tmpstate => stold ! Use stold for temporary storage of state
       if (eta > 0.0d0) then
          do nitt=1,10
             do ieq=1,neq
                do k=sz,ez
                   do j=sy,ey
                      do i=sx,ex
                         tmpstate(i,j,k,ieq)=state(i,j,k,ieq)+ & 
                              eta*( &
                              state(i-1,j,k,ieq)+ &
                              state(i+1,j,k,ieq)+ &
                              state(i,j-1,k,ieq)+ &
                              state(i,j+1,k,ieq)+ &
                              state(i,j,k-1,ieq)+ &
                              state(i,j,k+1,ieq)- &
                              6.0d0*state(i,j,k,ieq))
                      enddo
                   enddo
                enddo
             enddo
             do ieq=1,neq
                do k=sz,ez
                   do j=sy,ey
                      do i=sx,ex
                         state(i,j,k,ieq)=tmpstate(i,j,k,ieq)
                      enddo
                   enddo
                enddo
             enddo
             
             ! In order for the grid boundaries to be consistent, we need
             ! to reset them after each smoothing loop.
             ! exchange boundaries with neighbours
             call exchngxy(NEW)
          enddo
       endif

       ! Introduce the shockwave. This happens outside the smoothing
       ! The shock wave is initially located at a distance shckdist from
       ! the edge of the clump, but not further than the left edge of the
       ! the grid (x(2))
       
       shckdist=10.0d0*dx
       do k=sz-mbc,ez+mbc
          do j=sy-mbc,ey+mbc
             do i=sx-mbc,ex+mbc
                if (x(i) < max(1.5d0*dx,x0-axis-shckdist)) then
                   state(i,j,k,RHO)=sedensity
                   state(i,j,k,RHVX)=sedensity*sevelocity
                   state(i,j,k,RHVY)=0.0d0
                   state(i,j,k,RHVZ)=0.0d0
                   pressr(i,j,k)=sepressr
                   state(i,j,k,EN)=pressr(i,j,k)/gamma1+ &
                        0.5d0*(state(i,j,k,RHVX)*state(i,j,k,RHVX)+ &
                        state(i,j,k,RHVY)*state(i,j,k,RHVY)+ &
                        state(i,j,k,RHVZ)*state(i,j,k,RHVZ))/state(i,j,k,RHO)
                   state(i,j,k,TRACER1)=-1.0d0
                endif
             enddo
          enddo
       enddo
       
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

    if (sx == 1) then
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
                state(i,j,k,TRACER1)=-1.0d0
             enddo
          enddo
       enddo
    endif
    
  end subroutine inflow
  
end module problem
