module problem

  ! Module for Capreole3D (F90)
  ! Author: Garrelt Mellema
  ! Date: 2010-01-25

  ! This is the problem module. It contains the routines which define the
  ! problem being solved:

  ! Contents:
  ! init_problem - sets up the hydro variables according to the specified
  !                  problem
  ! inflow - sets the inflow boundary conditions

  ! This version: steady shock interacting with elliptical cloud
  !  (cartesian coordinates).
  ! Can be used for Test 7 of the Cosmological Radiative Transfer Comparison
  ! Project

  use file_admin, only: stdinput, log_unit, file_input
  use precision, only: dp
  use string_manipulation, only: convert_case
  use my_mpi
  use sizes, only: mbc,neq,RHO,RHVX,RHVY,RHVZ,EN,nrofDim
  use scaling, only: SCDENS, SCVELO, SCLENG, SCENER, SCMOME
  use cgsconstants, only: m_p, kb
  use astroconstants, only: pc, kpc, Mpc
  use abundances, only: mu
  use atomic, only: gamma, gamma1, boltzm
  use mesh, only: sx,ex,sy,ey,sz,ez,meshx,meshy,meshz
  use grid, only: x,y,z,dx,dy,dz
  use hydro, only: state,pressr,set_state_pointer,NEW,OLD,stnew,tmpstate,restart_state
  use boundary, only: exchngxy,REFLECTIVE,OUTFLOW,PROBLEM_DEF,X_IN,X_OUT,Y_IN, &
       Y_OUT,Z_IN,Z_OUT
  use ionic, only: init_ionic

  implicit none

  integer, dimension(nrofDim,2) :: domainboundaryconditions

  real(kind=dp),private :: sedensity,sevelocity,sepressr

contains

  subroutine init_problem (restart,restartfile)
    
    ! This routine initializes all hydro variables
    
    ! This may be a fresh start or a restart of a saved run
    
    !> tells you whether it's a new run or a restart
    logical,intent(in) :: restart 
    character(len=19),intent(in) :: restartfile !< file from which to restart

    ! Local variables
    real(kind=dp) :: r_interface !< dummy needed for calling init_ionic
    integer :: ierror !< error flag

    ! Set domain boundary conditions
    domainboundaryconditions(2:3,:)=OUTFLOW
    domainboundaryconditions(1,2)=OUTFLOW
    domainboundaryconditions(1,1)=PROBLEM_DEF

    ! Fill the state variable
    if (.not.restart) then 

       ! Fresh start
       call fresh_start_state( )
       
    else

       ! Read state from restartfile
       call restart_state(restartfile,ierror)

       ! Scale to code scaling
       state(:,:,:,RHO)=state(:,:,:,RHO)/scdens
       state(:,:,:,RHVX)=state(:,:,:,RHVX)/scmome
       state(:,:,:,RHVY)=state(:,:,:,RHVY)/scmome
       state(:,:,:,RHVZ)=state(:,:,:,RHVZ)/scmome
       state(:,:,:,EN)=state(:,:,:,EN)/scener

       call exchngxy(OLD,domainboundaryconditions,problemboundary) ! Fill boundary conditions

    endif
       
    ! Initialize the ionic concentrations
    call init_ionic(restart,r_interface)

  end subroutine init_problem

  !==========================================================================

  subroutine fresh_start_state ( )
    
    ! This routine initializes all hydro variables for a fresh start
    
    ! Case: steady shock interacting with elliptical cloud; (x,y,z)
    
    ! smoothing parameter for making soft-edged clump
    real(kind=dp),parameter :: eta=0.1_dp
    
    ! e variables are environment
    ! w variables are clump 
    ! se variables are shocked environment
    real(kind=dp)    :: edensity,etemperature,vblast
    real(kind=dp)    :: wdensity,wtemperature,x0,y0,z0,axis
    real(kind=dp)    :: axis1,axis2,axis3,rotangle1,rotangle2
    real(kind=dp)    :: epressr,wpressr,vs1,xm1,shckdist
    real(kind=dp)    :: dens_val,pres_val
    real(kind=dp),dimension(8) :: fp
    real(kind=dp),dimension(10,10,10) :: fp3d
    real(kind=dp)    :: hdx,hdy,hdz

    integer :: i,j,k,nitt,ieq,ii,jj,kk

    character(len=10) :: str_length_unit,str_a1_unit,str_a2_unit,str_a3_unit
#ifdef MPI       
    integer :: ierror
#endif

    ! Ask for the input if you are processor 0.
    if (rank == 0) then
       if (.not.file_input) then
          write (*,'(//,A,/)') '----- Environment -----'
          write (*,'(A,$)') '1) Density (cm^-3): '
       endif
       read (stdinput,*) edensity
       if (.not.file_input) write (*,'(A,$)') '2) Temperature: '
       read (stdinput,*) etemperature
       if (.not.file_input) write (*,'(A,$)') '3) Blast velocity (km/s): '
       read (stdinput,*) vblast
       
       if (.not.file_input) then
          write (*,'(//,A,/)') '----- Knot -----'
          write (*,'(A,$)') '1) Density (cm^-3): '
       endif
       read (stdinput,*) wdensity
       if (.not.file_input) write (*,'(A,$)') '2) Temperature: '
       read (stdinput,*) wtemperature
       if (.not.file_input) write (*,'(A,$)') '3) Position of centre x,y,z (specify units): '
       read (stdinput,*) x0,y0,z0,str_length_unit
       if (.not.file_input) write (*,'(A,$)') '4) Axis 1: '
       read (stdinput,*) axis1,str_a1_unit
       if (.not.file_input) write (*,'(A,$)') '5) Axis 2: '
       read (stdinput,*) axis2,str_a2_unit
       if (.not.file_input) write (*,'(A,$)') '6) Axis 3: '
       read (stdinput,*) axis3,str_a3_unit
       if (.not.file_input) write (*,'(A,$)') '7) Rotation angle1: '
       read (stdinput,*) rotangle1
       if (.not.file_input) write (*,'(A,$)') '8) Rotation angle2: '
       read (stdinput,*) rotangle2
    endif

    ! report input parameters
    if (rank == 0) then
       write(log_unit,'(A)') & 
            'Problem: cart-ellipseclump (cartesian shock - elliptical cloud interaction)'
       write (log_unit,'(//,A,/)') '----- Environment -----'
       write (log_unit,'(A,1PE10.3)') '1) Density (cm^-3): ',edensity
       write (log_unit,'(A,1PE10.3)') '2) Temperature: ',etemperature
       write (log_unit,'(A,F8.3)') '3) Blast velocity (km/s): ',vblast
       
       write (log_unit,'(//,A,/)') '----- Knot -----'
       write (log_unit,'(A,1PE10.3)') '1) Density (cm^-3): ',wdensity
       write (log_unit,'(A,1PE10.3)') '2) Temperature: ',wtemperature
       write (log_unit,'(A,3(1PE10.3,X),A)') '3) Position of centre x,y,z: ',&
            x0,y0,z0,str_length_unit
       write (log_unit,'(A,1PE10.3,A)') '4) Axis 1: ',axis1,str_a1_unit
       write (log_unit,'(A,1PE10.3,A)') '5) Axis 2: ',axis2,str_a2_unit
       write (log_unit,'(A,1PE10.3,A)') '6) Axis 3: ',axis3,str_a3_unit
       write (log_unit,'(A,F8.3)') '7) Rotation angle1: ',rotangle1
       write (log_unit,'(A,F8.3)') '8) Rotation angle2: ',rotangle2
       
       ! Convert to cm
       x0=x0*unit_conversion(str_length_unit)
       y0=y0*unit_conversion(str_length_unit)
       z0=z0*unit_conversion(str_length_unit)
       axis1=axis1*unit_conversion(str_a1_unit)
       axis2=axis2*unit_conversion(str_a2_unit)
       axis3=axis3*unit_conversion(str_a3_unit)
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
    
    epressr=edensity*boltzm*etemperature/mu
    wpressr=wdensity*boltzm*wtemperature/mu
    !wpressr=epressr       ! Pressure equilibrium!!!!
    
    ! Calculate the properties of the shock wave 
    vs1=sqrt(gamma*epressr/edensity) ! Sound speed in unshocked gas
    xm1=vblast/vs1                   ! Mach numnber of shock
    
    write(log_unit,*) 'Pressure= ',epressr*scener
    write(log_unit,*) 'Density= ',edensity*scdens
    write(log_unit,*) 'Mach number= ',xm1
    
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
    
    hdx=0.5*dx
    hdy=0.5*dy
    hdz=0.5*dz
    do k=sz-mbc,ez+mbc
       do j=sy-mbc,ey+mbc
          do i=sx-mbc,ex+mbc
             fp(1)=ellipse(x(i)+hdx,y(j)+hdy,z(k)+hdz, &
                  x0,y0,z0,rotangle1,rotangle2,axis1, &
                  axis2,axis3)
             fp(2)=ellipse(x(i)-hdx,y(j)+hdy,z(k)+hdz, &
                  x0,y0,z0,rotangle1,rotangle2,axis1, &
                  axis2,axis3)
             fp(3)=ellipse(x(i)+hdx,y(j)-hdy,z(k)+hdz, &
                  x0,y0,z0,rotangle1,rotangle2,axis1, &
                  axis2,axis3)
             fp(4)=ellipse(x(i)-hdx,y(j)-hdy,z(k)+hdz, &
                  x0,y0,z0,rotangle1,rotangle2,axis1, &
                  axis2,axis3)
             fp(5)=ellipse(x(i)+hdx,y(j)+hdy,z(k)-hdz, &
                  x0,y0,z0,rotangle1,rotangle2,axis1, &
                  axis2,axis3)
             fp(6)=ellipse(x(i)-hdx,y(j)+hdy,z(k)-hdz, &
                  x0,y0,z0,rotangle1,rotangle2,axis1, &
                  axis2,axis3)
             fp(7)=ellipse(x(i)+hdx,y(j)-hdy,z(k)-hdz, &
                  x0,y0,z0,rotangle1,rotangle2,axis1, &
                  axis2,axis3)
             fp(8)=ellipse(x(i)-hdx,y(j)-hdy,z(k)-hdz, &
                  x0,y0,z0,rotangle1,rotangle2,axis1, &
                  axis2,axis3)

             if (maxval(fp) <= 1.0) then
                state(i,j,k,RHO)=wdensity
                state(i,j,k,RHVX)=0.0d0
                state(i,j,k,RHVY)=0.0d0
                state(i,j,k,RHVZ)=0.0d0
                pressr(i,j,k)=wpressr
                !state(i,j,k,TRACER1)=1.0d0
             elseif (minval(fp) >= 1.0) then
                state(i,j,k,RHO)=edensity
                state(i,j,k,RHVX)=0.0d0
                state(i,j,k,RHVY)=0.0d0
                state(i,j,k,RHVZ)=0.0d0
                pressr(i,j,k)=epressr
                !state(i,j,k,TRACER1)=-1.0d0
             else
                ! do weighting
                dens_val=0.0
                pres_val=0.0
                do kk=1,10
                   do jj=1,10
                      do ii=1,10
                         fp3d(ii,jj,kk)=ellipse( &
                              x(i)-hdx+(real(ii-1)+0.5)*dx*0.1, &
                              y(j)-hdy+(real(jj-1)+0.5)*dy*0.1, &
                              z(k)-hdz+(real(kk-1)+0.5)*dz*0.1, &
                              x0,y0,z0,rotangle1,rotangle2,axis1, &
                              axis2,axis3)
                         if (fp3d(ii,jj,kk) <= 1.0) then
                            dens_val=dens_val + wdensity
                            pres_val=pres_val + wpressr
                         else
                            dens_val=dens_val + edensity
                            pres_val=pres_val + epressr
                         endif
                      enddo
                   enddo
                enddo
                dens_val=dens_val*1e-3
                pres_val=pres_val*1e-3
                state(i,j,k,RHO)=dens_val
                state(i,j,k,RHVX)=0.0d0
                state(i,j,k,RHVY)=0.0d0
                state(i,j,k,RHVZ)=0.0d0
                pressr(i,j,k)=pres_val
             endif
             state(i,j,k,EN)=pressr(i,j,k)/gamma1+ &
                  0.5d0*(state(i,j,k,RHVX)*state(i,j,k,RHVX)+ &
                  state(i,j,k,RHVY)*state(i,j,k,RHVY)+ &
                  state(i,j,k,RHVZ)*state(i,j,k,RHVZ))/state(i,j,k,RHO)
          enddo
       enddo
    enddo
    
    ! To avoid numerical problems the edge of the clump needs to be
    ! smoothed. We achieve this by diffusing the state variables a number
    ! of time with a high diffuse parameter eta. 
    
    tmpstate => stnew ! Use stnew for temporary storage of state
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
       call exchngxy(OLD,domainboundaryconditions,problemboundary) ! Fill boundary conditions
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
                !state(i,j,k,TRACER1)=-1.0d0
             endif
          enddo
       enddo
    enddo

  end subroutine fresh_start_state

  !==========================================================================

  subroutine problemboundary (boundary,newold)
    
    ! This routine resets the inner boundary to the inflow condition
    
    integer,intent(in) :: boundary
    integer,intent(in) :: newold

    integer :: i,j,k

    ! Point state to appropriate array
    state => set_state_pointer(newold)

    select case (boundary)
    case (X_IN)
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
    case (X_OUT)
    case (Y_IN)
    case (Y_OUT)
    case (Z_IN)
    case (Z_OUT)
    end select
    
  end subroutine problemboundary

  !==========================================================================

  subroutine apply_grav_force(dt,newold)

    ! Dummy routine

    real(kind=dp),intent(in) :: dt
    integer,intent(in) :: newold

  end subroutine apply_grav_force
    
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
       write(log_unit,*) 'Length unit not recognized, assuming cm'
       conversion_factor=1.0
    end select
    
    unit_conversion=conversion_factor
    
  end function unit_conversion

  !==========================================================================

  function ellipse (x,y,z,x0,y0,z0,angle1,angle2,axis1,axis2,axis3)
    
    real(kind=dp) :: ellipse

    real(kind=dp),intent(in) :: x
    real(kind=dp),intent(in) :: y
    real(kind=dp),intent(in) :: z
    real(kind=dp),intent(in) :: x0
    real(kind=dp),intent(in) :: y0
    real(kind=dp),intent(in) :: z0
    real(kind=dp),intent(in) :: angle1
    real(kind=dp),intent(in) :: angle2
    real(kind=dp),intent(in) :: axis1
    real(kind=dp),intent(in) :: axis2
    real(kind=dp),intent(in) :: axis3
    
    real(kind=dp) :: xc,yc,zc
    real(kind=dp) :: xr1,yr1,xr2,zr2

    xc=x-x0
    yc=y-y0
    zc=z-z0
    
    xr1 =  cos(angle1)*xc+sin(angle1)*yc
    yr1 = -sin(angle1)*xc+cos(angle1)*yc
    
    xr2 = cos(angle2)*xr1+sin(angle2)*zc
    zr2 = -sin(angle2)*xr1+cos(angle2)*zc
    
    ellipse=sqrt(xr2*xr2/(axis1*axis1)+yr1*yr1/(axis2*axis2)+ &
         zr2*zr2/(axis3*axis3))
    
  end function ellipse

end module problem
