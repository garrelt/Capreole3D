module problem

  ! Module for Capreole3D (F90)
  ! Author: Garrelt Mellema
  ! Date: 2007-10-02

  ! This is the problem module. It contains the routines which define the
  ! problem being solved:

  ! Contents:
  ! init_problem - sets up the hydro variables according to the specified
  !                  problem
  ! inflow - sets the inflow boundary conditions

  ! This version: constant density

  use file_admin, only: stdinput
  use precision, only: dp
  use my_mpi
  use sizes, only: mbc,neq,RHO,RHVX,RHVY,RHVZ,EN
  use scaling, only: SCDENS, SCVELO, SCLENG, SCMOME, SCENER
  use atomic, only: gamma,gamma1
  use cgsconstants, only: m_p
  use astroconstants, only: pc
  use abundances, only: mu
  use mesh, only: sx,ex,sy,ey,sz,ez,meshx,meshy,meshz
  use grid, only: x,y,z,dx,dy,dz
  use hydro, only: state,pressr,set_state_pointer,NEW,OLD,restart_state
  use ionic, only: init_ionic
  use tped, only: n2rho,rho2n,pressr2temper,temper2pressr
  use protection, only: presprot
  use geometry, only: presfunc
  use boundary, only: exchngxy

  implicit none

contains

  subroutine init_problem (restart,restartfile)
    
    ! This routine initializes all hydro variables
    
    ! This may be a fresh start or a restart of a saved run
    
    logical,intent(in) :: restart ! tells you whether it's a new run or a restart
    character(len=19),intent(in) :: restartfile

    ! Local variables
    real(kind=dp) :: r_interface ! dummy needed for calling init_ionic
    integer :: ierror

    if (.not.restart) then ! Fresh start

       call fresh_start_state( )
       
    else

       call restart_state(restartfile,ierror)
       state(:,:,:,RHO)=state(:,:,:,RHO)/scdens
       state(:,:,:,RHVX)=state(:,:,:,RHVX)/scmome
       state(:,:,:,RHVY)=state(:,:,:,RHVY)/scmome
       state(:,:,:,RHVZ)=state(:,:,:,RHVZ)/scmome
       state(:,:,:,EN)=state(:,:,:,EN)/scener
       call exchngxy(OLD)

    endif
       
    ! Initialize the ionic concentrations
    call init_ionic(restart,r_interface)

  end subroutine init_problem

  !==========================================================================

  subroutine fresh_start_state ( )
    
    ! This routine initializes all hydro variables for a fresh start
    
    ! Case: Enrique's density and velocity fields for a turbulent
    !  molecular cloud.
    
    ! Local variables
    real(kind=dp) :: numdens0,temper0
    real(kind=dp) :: rhoconversion,pressrconversion
    integer :: i,j,k,ieq,ierror,nitt,idim, inegative

    ! Ask for the input if you are processor 0.
    if (rank == 0) then
       write (*,'(//,A,/)') '----- Environment -----'
       write (*,'(A,$)') '1) Number density: '
       read (*,*) numdens0
       write (*,'(A,$)') '2) Temperature: '
       read (*,*) temper0
    endif
    
    ! report input parameters
    if (rank == 0) then
       write(30,'(A)') & 
            'Problem: cart-constant'
       write (30,'(//,A,/)') '----- Environment -----'
       write (30,'(A,F10.3)') '1) Number density: ',numdens0
       write (30,'(A,F10.3)') '2) Temperature: ',temper0
    endif
    
    rhoconversion=n2rho(numdens0)
    pressrconversion=temper2pressr(temper0,numdens0,0.0d0)

    ! Calculate the initial state: pressure and energy density
    do k=sz-mbc,ez+mbc
       do j=sy-mbc,ey+mbc
          do i=sx-mbc,ex+mbc
             state(i,j,k,RHO)=rhoconversion/SCDENS
             state(i,j,k,RHVX)=0.0
             state(i,j,k,RHVY)=0.0
             state(i,j,k,RHVZ)=0.0
             pressr(i,j,k)=pressrconversion/SCENER
             state(i,j,k,EN)=pressr(i,j,k)/gamma1+ &
                  0.5d0*(state(i,j,k,RHVX)*state(i,j,k,RHVX)+ &
                  state(i,j,k,RHVY)*state(i,j,k,RHVY)+ &
                  state(i,j,k,RHVZ)*state(i,j,k,RHVZ))/state(i,j,k,RHO)
          enddo
       enddo
    enddo
    
  end subroutine fresh_start_state

  !==========================================================================

  subroutine inflow (newold)
    
    ! This routine resets the inner boundary to the inflow condition
    ! Dummy version

    integer,intent(in) :: newold

  end subroutine inflow
  
end module problem
