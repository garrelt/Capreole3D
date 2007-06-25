module problem

  ! Module for Capreole3D (F90)
  ! Author: Garrelt Mellema
  ! Date: 2006-08-21

  ! This is the problem module. It contains the routines which define the
  ! problem being solved:

  ! Contents:
  ! init_problem - sets up the hydro variables according to the specified
  !                  problem
  ! inflow - sets the inflow boundary conditions

  ! This version: 1/r^2 density profile with flat core
  !    (cartesian coordinates).

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
    real(kind=dp) :: r_core,dens_core,etemperature
    real(kind=dp) :: dens_val,xs,ys,zs,dist

    integer :: i,j,k,nitt,ieq
    real(kind=dp)    :: r_interface ! dummy needed for calling init_ionic

#ifdef MPI       
    integer :: ierror
#endif

    if (.not.restart) then ! Fresh start
       
       ! Ask for the input if you are processor 0.
       if (rank == 0) then
          write (*,'(//,A,/)') '----- Halo -----'
          write(*,'(A,$)') '1) Reference (core) radius (pc): '
          read(stdinput,*) r_core
          write(*,'(A,$)') '2) Density at reference (core) radius (cm^-3): '
          read(stdinput,*) dens_core
          write (*,'(A,$)') '3) Initial temperature (K): '
          read (stdinput,*) etemperature
       endif

       ! report input parameters
       if (rank == 0) then
          write(30,'(A)') & 
               'Problem: cart-halo (halo with constant density core', &
               ' and 1/r^2 outside of the core)'
          write (30,'(//,A,/)') '----- Halo -----'
          write (30,'(A,1PE10.3)') '1) Reference (core) radius (pc): ',r_core
          write (30,'(A,1PE10.3)') '2) Core density (cm^-3): ',dens_core
          write (30,'(A,1PE10.3)') '3) Temperature: ',etemperature
          
       endif
          
#ifdef MPI       
       ! Distribute the input parameters to the other nodes
       call MPI_BCAST(r_core,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
            ierror)
       call MPI_BCAST(dens_core,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
            ierror)
       call MPI_BCAST(etemperature,1,MPI_DOUBLE_PRECISION,0, & 
            MPI_COMM_NEW,ierror)
#endif
       
       ! Scale physical parameters (density scaled below)
       r_core=r_core*pc/scleng ! core radius is in pc
       
       ! Initialize the ionic concentrations and all radiation
       ! quantities (need to do this before initializing the
       ! state since this needs the source position
       call init_ionic(restart,r_interface)

       ! Calculate the initial state
       ! Source position determines the core location
       do k=sz-mbc,ez+mbc
          do j=sy-mbc,ey+mbc
             do i=sx-mbc,ex+mbc
                xs=x(i)-rsrcpos(1)
                ys=y(j)-rsrcpos(2)
                zs=z(k)-rsrcpos(3)
                dist=sqrt(xs*xs+ys*ys+zs*zs)
                if (dist.le.r_core) then
                   ! Flat core
                   dens_val=dens_core
                else
                   dens_val=dens_core*(dist/r_core)**(-2.0)
                endif
                state(i,j,k,RHO)=dens_val*m_p/scdens
                state(i,j,k,RHVX)=0.0d0
                state(i,j,k,RHVY)=0.0d0
                state(i,j,k,RHVZ)=0.0d0
                pressr(i,j,k)=dens_val*kb*etemperature/scener
                state(i,j,k,EN)=pressr(i,j,k)/gamma1+ &
                     0.5d0*(state(i,j,k,RHVX)*state(i,j,k,RHVX)+ &
                     state(i,j,k,RHVY)*state(i,j,k,RHVY)+ &
                     state(i,j,k,RHVZ)*state(i,j,k,RHVZ))/state(i,j,k,RHO)
             enddo
          enddo
       enddo
       
    endif

  end subroutine init_problem

  !==========================================================================

  subroutine inflow (newold)
    
    ! This routine resets the inner boundary to the inflow condition
    
    integer,intent(in) :: newold

  end subroutine inflow
  
end module problem
