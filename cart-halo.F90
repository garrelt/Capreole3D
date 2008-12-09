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

  ! This version: 1/r^2 density profile with flat core (cartesian coordinates).
  ! This is Test 6 from the Cosmological Radiative Transfer Comparison Project

  use file_admin, only: stdinput
  use precision, only: dp
  use string_manipulation, only: convert_case
  use my_mpi
  use sizes, only: mbc,neq,RHO,RHVX,RHVY,RHVZ,EN,XHI,XHII
  use scaling, only: SCDENS, SCVELO, SCLENG, SCENER, SCMOME
  use cgsconstants, only: m_p, kb
  use astroconstants, only: pc, kpc, Mpc
  use abundances, only: mu
  use atomic, only: gamma, gamma1, boltzm
  use mesh, only: sx,ex,sy,ey,sz,ez,meshx,meshy,meshz
  use grid, only: x,y,z,dx,dy,dz
  use hydro, only: state,pressr,set_state_pointer,NEW,OLD,stnew,tmpstate,restart_state
  use boundary, only: exchngxy
  use ionic, only: init_ionic
  use sourceprops, only: rsrcpos
  use tped, only: electrondens,temper2pressr,n2rho

  implicit none

  real(kind=dp),private :: sedensity,sevelocity,sepressr

contains

  subroutine init_problem (restart,restartfile)
    
    ! This routine initializes all hydro variables
    
    ! This may be a fresh start or a restart of a saved run
    
    logical,intent(in) :: restart ! tells you whether it's a new run or a restart
    character(len=19),intent(in) :: restartfile

    ! Local variables
    real(kind=dp) :: r_interface ! dummy needed for calling init_ionic
    integer :: ierror

    ! Initialize the ionic concentrations and all radiation
    ! quantities (need to do this before initializing the
    ! state since this needs the source position
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
       
  end subroutine init_problem

  !==========================================================================

  subroutine fresh_start_state ( )
    
    ! This routine initializes all hydro variables for a fresh start
    
    ! Case: 1/r^2 halo
    
    real(kind=dp) :: r_core,dens_core,etemperature
    real(kind=dp) :: dens_val,xs,ys,zs,dist
    real(kind=dp) :: mindist,maxdist
    real(kind=dp) :: hdx,hdy,hdz,r_core2
    real(kind=dp),dimension(0:8) :: dist2
    real(kind=dp),dimension(10,10,10) :: dist3d
    
    integer :: i,j,k,nitt,ieq,ii,jj,kk
    character(len=10) :: str_length_unit

#ifdef MPI       
    integer :: ierror
#endif

    ! Ask for the input if you are processor 0.
    if (rank == 0) then
       write (*,'(//,A,/)') '----- Halo -----'
       write(*,'(A,$)') '1) Reference (core) radius (specify unit): '
       read(stdinput,*) r_core,str_length_unit
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
       write (30,'(A,1PE10.3,A)') '1) Reference (core) radius: ',r_core, &
            str_length_unit
       write (30,'(A,1PE10.3)') '2) Core density (cm^-3): ',dens_core
       write (30,'(A,1PE10.3)') '3) Temperature: ',etemperature
       
       ! Convert to cm
       r_core=r_core*unit_conversion(str_length_unit)
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
    r_core=r_core/scleng
    
    ! Calculate the initial state
    ! Source position determines the core location
    hdx=0.5*dx
    hdy=0.5*dy
    hdz=0.5*dz
    r_core2=r_core*r_core
    do k=sz-mbc,ez+mbc
       do j=sy-mbc,ey+mbc
          do i=sx-mbc,ex+mbc
             xs=x(i)-rsrcpos(1)
             ys=y(j)-rsrcpos(2)
             zs=z(k)-rsrcpos(3)
             dist2(0)=((xs)*(xs)+(ys)*(ys)+(zs)*(zs))
             dist2(1)=((xs+hdx)*(xs+hdx)+(ys+hdy)*(ys+hdy)+(zs+hdz)*(zs+hdz))
             dist2(2)=((xs-hdx)*(xs-hdx)+(ys+hdy)*(ys+hdy)+(zs+hdz)*(zs+hdz))
             dist2(3)=((xs+hdx)*(xs+hdx)+(ys-hdy)*(ys-hdy)+(zs+hdz)*(zs+hdz))
             dist2(4)=((xs+hdx)*(xs+hdx)+(ys+hdy)*(ys+hdy)+(zs-hdz)*(zs-hdz))
             dist2(5)=((xs-hdx)*(xs-hdx)+(ys-hdy)*(ys-hdy)+(zs+hdz)*(zs+hdz))
             dist2(6)=((xs+hdx)*(xs+hdx)+(ys-hdy)*(ys-hdy)+(zs-hdz)*(zs-hdz))
             dist2(7)=((xs-hdx)*(xs-hdx)+(ys+hdy)*(ys+hdy)+(zs-hdz)*(zs-hdz))
             dist2(8)=((xs-hdx)*(xs-hdx)+(ys-hdy)*(ys-hdy)+(zs-hdz)*(zs-hdz))
             mindist=minval(dist2)
             maxdist=maxval(dist2)
             if (maxdist <= r_core2) then
                ! Flat core
                dens_val=dens_core
             elseif (mindist > r_core2) then
                ! r^-2 environment
                dens_val=dens_core*(dist2(0)/r_core2)
             else
                ! do weighting
                dens_val=0.0
                do kk=1,10
                   do jj=1,10
                      do ii=1,10
                         dist3d(ii,jj,kk)= &
                              (xs-hdx +(real(ii-1)+0.5)*dx*0.1)**2 +  &
                              (ys-hdy +(real(jj-1)+0.5)*dy*0.1)**2 +  & 
                              (zs-hdz +(real(kk-1)+0.5)*dz*0.1)**2
                         if (dist3d(ii,jj,kk) <= r_core2) then
                            dens_val=dens_val + dens_core
                         else
                            dens_val=dens_val + &
                                 dens_core*(dist3d(ii,jj,kk)/r_core2)
                         endif
                      enddo
                   enddo
                enddo
                dens_val=dens_val*1e-3
             endif
             state(i,j,k,RHO)=n2rho(dens_val)/scdens
             state(i,j,k,RHVX)=0.0d0
             state(i,j,k,RHVY)=0.0d0
             state(i,j,k,RHVZ)=0.0d0
             pressr(i,j,k)=temper2pressr(etemperature,dens_val,&
                  electrondens(dens_val,state(i,j,k,XHI:XHII)))/scener
             state(i,j,k,EN)=pressr(i,j,k)/gamma1+ &
                  0.5d0*(state(i,j,k,RHVX)*state(i,j,k,RHVX)+ &
                  state(i,j,k,RHVY)*state(i,j,k,RHVY)+ &
                  state(i,j,k,RHVZ)*state(i,j,k,RHVZ))/state(i,j,k,RHO)
          enddo
       enddo
    enddo

  end subroutine fresh_start_state

  !==========================================================================
  
  subroutine apply_grav_force(dt,newold)
    
    ! Dummy routine
    
    real(kind=dp),intent(in) :: dt
    integer,intent(in) :: newold
    
  end subroutine apply_grav_force
  
  !==========================================================================
  
  subroutine inflow (newold)
    
    ! This routine resets the inner boundary to the inflow condition
    ! Version: dummy
    
    integer,intent(in) :: newold
    
  end subroutine inflow
  
  !==========================================================================
  
  function unit_conversion(in_str_unit)
    
    real(kind=dp) :: unit_conversion
    
    real(kind=dp) :: conversion_factor
    
    character(len=10),intent(in) :: in_str_unit
    character(len=10) :: str_unit
    
    str_unit=in_str_unit
    call convert_case(str_unit,0) ! conversion to lower case
    select case (trim(adjustl(str_unit)))
    case ('cm','centimeter','cms','centimeters')
       conversion_factor=1.0
    case ('m','meter','ms','meters')
       conversion_factor=100.0
    case ('km','kilometer','kms','kilometers','clicks')
       conversion_factor=1000.0
    case ('pc','parsec','parsecs')
       conversion_factor=pc
    case ('kpc','kiloparsec','kiloparsecs')
       conversion_factor=kpc
    case ('mpc','megaparsec','megaparsecs')
       conversion_factor=Mpc
    case default
       write(*,*) 'Length unit not recognized, assuming cm'
       conversion_factor=1.0
    end select
    
    unit_conversion=conversion_factor
    
  end function unit_conversion

end module problem
