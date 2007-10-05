module problem

  ! Module for Capreole3D (F90)
  ! Author: Garrelt Mellema
  ! Date: 2007-10-03

  ! This is the problem module. It contains the routines which define the
  ! problem being solved:

  ! Contents:
  ! init_problem - sets up the hydro variables according to the specified
  !                  problem
  ! inflow - sets the inflow boundary conditions

  ! This version: density, velocity, temperature fields (to be read in) for 
  ! HII region evolution (cartesian coordinates). The fields are supposed
  ! to come from a GADGET turbulence calculation of 
  ! Enrique Vazquez-Semadeni et al. 

  use file_admin, only: stdinput
  use precision, only: dp
  use my_mpi
  use sizes, only: mbc,neq,RHO,RHVX,RHVY,RHVZ,EN
  use scaling, only: SCDENS, SCVELO, SCLENG, SCMOME, SCENER
  use atomic, only: gamma,gamma1
  use cgsconstants, only: m_p
  use astroconstants, only: M_SOLAR,pc
  use abundances, only: mu
  use mesh, only: sx,ex,sy,ey,sz,ez,meshx,meshy,meshz
  use grid, only: x,y,z,dx,dy,dz
  use hydro, only: state,pressr,set_state_pointer,NEW,OLD,restart_state
  use ionic, only: init_ionic
  use tped, only: n2rho,rho2n,temper2pressr
  use protection, only: presprot
  use geometry, only: presfunc
  use boundary, only: exchngxy
  use sourceprops, only: srcpos

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

    ! Initialize the radiative bits
    ! (ionization fractions, sources, radiation tables)
    call init_ionic(restart,r_interface)

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
       
    ! Set the boundary conditions
    call exchngxy(OLD)

  end subroutine init_problem

  !==========================================================================

  subroutine fresh_start_state ( )
    
    ! This routine initializes all hydro variables for a fresh start
    
    ! Case: Density, temperature and velocity fields for a turbulent
    !  molecular cloud (GADGET simulation)
    
    ! Scaling parameters of GADGET
    real(kind=dp),parameter :: density_scaling=M_SOLAR/pc**3
    real(kind=dp),parameter :: velocity_scaling=736200.0

    ! Local variables
    real :: r1,r2,r3,r4
    character(len=512) :: densityfield,temperaturefield
    character(len=512) :: velxfield,velyfield,velzfield
    real,allocatable,dimension(:,:,:) :: temparray
    integer :: i,j,k,ieq,ierror,nitt,idim, inegative

    ! Ask for the input if you are processor 0.
    if (rank == 0) then
       write (*,'(//,A,/)') '----- Fields -----'
       write (*,'(A,$)') '1) Density File: '
       read (*,*) densityfield
       write (*,'(A,$)') '2) Temperature File: '
       read (*,*) temperaturefield
       write (*,'(A,$)') '3) X-velocity File: '
       read (*,*) velxfield
       write (*,'(A,$)') '4) Y-velocity File: '
       read (*,*) velyfield
       write (*,'(A,$)') '5) Z-velocity File: '
       read (*,*) velzfield
    endif
    
    ! report input parameters
    if (rank == 0) then
       write(30,'(A)') & 
            'Problem: cart-gadget'
       write (30,'(//,A,/)') '----- Fields -----'
       write (30,'(2A)') '1) Density File: ', densityfield
       write (30,'(2A)') '2) Temperature File: ', temperaturefield
       write (30,'(2A)') '3) X-velocity File: ', velxfield
       write (30,'(2A)') '4) Y-velocity File: ', velyfield
       write (30,'(2A)') '5) Z-velocity File: ', velzfield
    endif
    
    ! Read in fields
    open(unit=21,file=densityfield,form='unformatted',status='old')
    read(21) r1,r2,r3
    allocate(temparray(meshx,meshy,meshz))
    read(21) temparray
    close(21)
    
    write(30,*) 'Maximum density: ',maxval(temparray), &
         ' at ',maxloc(temparray)
    
    state(sx:ex,sy:ey,sz:ez,RHO)=temparray*density_scaling/SCDENS

    open(unit=21,file=temperaturefield,form='unformatted',status='old')
    read(21) r1,r2,r3,r4
    allocate(temparray(meshx,meshy,meshz))
    read(21) temparray
    close(21)
    
    write(30,*) 'Maximum temperature: ',maxval(temparray), &
         ' at ',maxloc(temparray)
    do k=sz,ez
       do j=sy,ey
          do i=sx,ex
             pressr(i,j,k)=temper2pressr(real(temparray(i,j,k),kind=dp), &
                  rho2n(state(i,j,k,RHO)),0.0d0)/SCENER
          enddo
       enddo
    enddo

    open(unit=21,file=velxfield,form='unformatted',status='old')
    read(21) r1,r2,r3,r4
    allocate(temparray(meshx,meshy,meshz))
    read(21) temparray
    close(21)

    ! Subtract the velocity of the source point
    temparray(:,:,:)=temparray(:,:,:)-temparray(srcpos(1),srcpos(2),srcpos(3))
    state(sx:ex,sy:ey,sz:ez,RHVX)=state(sx:ex,sy:ey,sz:ez,RHO)*temparray/SCVELO

    open(unit=21,file=velyfield,form='unformatted',status='old')
    read(21) r1,r2,r3,r4
    allocate(temparray(meshx,meshy,meshz))
    read(21) temparray
    close(21)

    ! Subtract the velocity of the source point
    temparray(:,:,:)=temparray(:,:,:)-temparray(srcpos(1),srcpos(2),srcpos(3))
    state(sx:ex,sy:ey,sz:ez,RHVY)=state(sx:ex,sy:ey,sz:ez,RHO)*temparray/SCVELO

    open(unit=21,file=velzfield,form='unformatted',status='old')
    read(21) r1,r2,r3,r4
    allocate(temparray(meshx,meshy,meshz))
    read(21) temparray
    close(21)

    ! Subtract the velocity of the source point
    temparray(:,:,:)=temparray(:,:,:)-temparray(srcpos(1),srcpos(2),srcpos(3))
    state(sx:ex,sy:ey,sz:ez,RHVZ)=state(sx:ex,sy:ey,sz:ez,RHO)*temparray/SCVELO

    ! Deallocate the temporary arrays
    deallocate(temparray)
    
    ! Calculate the initial state: energy density
    do k=sz-mbc,ez+mbc
       do j=sy-mbc,ey+mbc
          do i=sx-mbc,ex+mbc
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
  
  !==========================================================================

  subroutine apply_grav_force(dt,newold)

    ! Dummy routine

    real(kind=dp),intent(in) :: dt
    integer,intent(in) :: newold

  end subroutine apply_grav_force
    
end module problem
