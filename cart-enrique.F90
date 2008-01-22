module problem

  ! Module for Capreole3D (F90)
  ! Author: Garrelt Mellema
  ! Date: 2007-10-01

  ! This is the problem module. It contains the routines which define the
  ! problem being solved:

  ! Contents:
  ! init_problem - sets up the hydro variables according to the specified
  !                  problem
  ! inflow - sets the inflow boundary conditions

  ! This version: density and velocity fields (to be read in) for 
  ! HII region evolution (cartesian coordinates). The density and
  ! and velocity fields are supposed to come from a turbulence
  ! calculation of Enrique Vazquez-Semadeni et al. 

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
  use tped, only: n2rho,rho2n,pressr2temper
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
    
    ! Parameters of turbulence simulation
    real(kind=dp),parameter :: Jeans=4.0 ! Jeans number
    real(kind=dp),parameter :: csound=0.2e5 ! sound speed
    real(kind=dp),parameter :: mu_simulation=2.4 ! mean weight per particle

    ! Local variables
    real(kind=dp) :: rhomean,vx_c,vy_c,vz_c
    integer :: m1,m2,m3
    character(len=4) :: file_id
    character(len=512) :: densityfield,momxfield,momyfield,momzfield
    character(len=3) :: center_answer
    real,allocatable,dimension(:,:,:) :: tempdens,tempmomx,tempmomy,tempmomz
    integer,dimension(3) :: center_here
    integer :: shft
    integer :: i,j,k,ieq,ierror,nitt,idim, inegative

    ! Ask for the input if you are processor 0.
    if (rank == 0) then
       write (*,'(//,A,/)') '----- Environment -----'
       write (*,'(A,$)') '1) File identifier (e.g. 0020): '
       read (stdinput,*) file_id
       ! center on maximum density
       center_answer='max'
    endif
    
    ! report input parameters
    if (rank == 0) then
       write(30,'(A)') & 
            'Problem: cart-enrique'
       write (30,'(//,A,/)') '----- Environment -----'
       write (30,'(2A)') '1) File identifier: ',file_id
       write (30,'(2A)') '3) Grid centered on: ',center_answer
    endif
    
    ! Construct file names from file identifier
    densityfield='ro'//file_id//'.bin'
    momxfield='rovx'//file_id//'.bin'
    momyfield='rovy'//file_id//'.bin'
    momzfield='rovz'//file_id//'.bin'
    
    ! Read in density fields
    open(unit=21,file=densityfield,form='unformatted',status='old')
    read(21) m1,m2,m3
    if (m1 .ne. meshx .or. m2.ne.meshy .or. m3.ne.meshz) then
       write(*,*) 'Error: density field does not match specified mesh size'
       stop
    endif
    allocate(tempdens(m1,m2,m3))
    read(21) tempdens
    close(21)
    
    ! Read in momenta fields
    allocate(tempmomx(m1,m2,m3))
    allocate(tempmomy(m1,m2,m3))
    allocate(tempmomz(m1,m2,m3))
    open(unit=21,file=momxfield,form='unformatted',status='old')
    read(21) m1,m2,m3
    read(21) tempmomx
    close(21)
    open(unit=21,file=momyfield,form='unformatted',status='old')
    read(21) m1,m2,m3
    read(21) tempmomy
    close(21)
    open(unit=21,file=momzfield,form='unformatted',status='old')
    read(21) m1,m2,m3
    read(21) tempmomz
    close(21)
    
    ! Take care of centering
    select case (center_answer)
    case ('max') 
       center_here=maxloc(tempdens)
    case ('min') 
       center_here=minloc(tempdens)
    case default 
       center_here=(/ m1/2, m2/2, m3/2 /)
    end select
    
    do idim=1,3
       shft=center_here(idim)-m1/2
       tempdens=cshift(tempdens,shft,idim)
       tempmomx=cshift(tempmomx,shft,idim)
       tempmomy=cshift(tempmomy,shft,idim)
       tempmomz=cshift(tempmomz,shft,idim)
    enddo
    
    write(*,*) 'Maximum density: ',maxval(tempdens), &
         ' at ',maxloc(tempdens)
    
    ! Subtract velocity of source point
    vx_c=tempmomx(m1/2,m2/2,m3/2)/tempdens(m1/2,m2/2,m3/2)
    vy_c=tempmomy(m1/2,m2/2,m3/2)/tempdens(m1/2,m2/2,m3/2)
    vz_c=tempmomz(m1/2,m2/2,m3/2)/tempdens(m1/2,m2/2,m3/2)
    tempmomx=(tempmomx/tempdens-vx_c)*tempdens
    tempmomy=(tempmomy/tempdens-vy_c)*tempdens
    tempmomz=(tempmomz/tempdens-vz_c)*tempdens
    write(*,*) minval(tempmomx/tempdens),maxval(tempmomx/tempdens)
    write(*,*) minval(tempmomy/tempdens),maxval(tempmomy/tempdens)
    write(*,*) minval(tempmomz/tempdens),maxval(tempmomz/tempdens)
    
    ! Scaling: the mean number density is given by 
    ! 500*(J/L(pc))^2, see Vazquez-Semanami et al. (2005).
    rhomean=n2rho(500.0*(Jeans/(x(meshx)*scleng)*PC)**2)
    write(*,*) minval(tempdens)*(500.0*(Jeans/(x(meshx)*scleng)*PC)**2), &
         maxval(tempdens)*(500.0*(Jeans/(x(meshx)*scleng)*PC)**2)
    ! The turbulence simulation uses molecular hydrogen
    ! with a mean weight per particle of 2.4. We assume
    ! atomic hydrogen with a mean weight per particle of
    ! 1.3 (including the helium).
    rhomean=mu_simulation/mu*rhomean
    
    state(sx:ex,sy:ey,sz:ez,RHO)=tempdens*rhomean
    state(sx:ex,sy:ey,sz:ez,RHVX)=tempmomx*csound*rhomean
    state(sx:ex,sy:ey,sz:ez,RHVY)=tempmomy*csound*rhomean
    state(sx:ex,sy:ey,sz:ez,RHVZ)=tempmomz*csound*rhomean
    
    ! Deallocate the temporary arrays
    deallocate(tempdens)
    deallocate(tempmomx)
    deallocate(tempmomy)
    deallocate(tempmomz)
    
    ! periodic box
    do ieq=RHO,RHVZ
       state(sx-mbc:sx-1,:,:,ieq)=state(ex-mbc:ex-1,:,:,ieq)
       state(ex+1:ex+mbc,:,:,ieq)=state(sx:sx+mbc,:,:,ieq)
       state(:,sy-mbc:sy-1,:,ieq)=state(:,ey-mbc:ey-1,:,ieq)
       state(:,ey+1:ey+mbc,:,ieq)=state(:,sy:sy+mbc,:,ieq)
       state(:,:,sz-mbc:sz-1,ieq)=state(:,:,ez-mbc:ez-1,ieq)
       state(:,:,ez+1:ez+mbc,ieq)=state(:,:,sz:sz+mbc,ieq)
    enddo
    
    ! Scale physical parameters
    state(:,:,:,RHO)=state(:,:,:,RHO)/scdens
    state(:,:,:,RHVX)=state(:,:,:,RHVX)/scmome
    state(:,:,:,RHVY)=state(:,:,:,RHVY)/scmome
    state(:,:,:,RHVZ)=state(:,:,:,RHVZ)/scmome
    
    ! Calculate the initial state: pressure and energy density
    do k=sz-mbc,ez+mbc
       do j=sy-mbc,ey+mbc
          do i=sx-mbc,ex+mbc
             pressr(i,j,k)=state(i,j,k,RHO)*csound*csound/gamma/ &
                  scvelo/scvelo
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
