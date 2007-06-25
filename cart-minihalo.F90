module problem

  ! Module for Capreole3D (F90)
  ! Author: Garrelt Mellema
  ! Date: 2007-06-22

  ! This is the problem module. It contains the routines which define the
  ! problem being solved:

  ! Contents:
  ! init_problem - sets up the hydro variables according to the specified
  !                  problem
  ! inflow - sets the inflow boundary conditions

  ! This version: Mini-halo (cartesian coordinates).

  use file_admin, only: stdinput
  use precision, only: dp
  use my_mpi
  use sizes, only: mbc,neq,RHO,RHVX,RHVY,RHVZ,EN
  use scaling, only: SCDENS, SCVELO, SCLENG
  use mathconstants, only: pi
  use cgsconstants, only: G_grav, m_p, kb
  use cosmology_parameters, only: rho_crit_0, Omega0, Omega_b, H0
  use astroconstants, only: M_solar, kpc
  use atomic, only: gamma1, boltzm
  use abundances, only: mu
  use mesh, only: sx,ex,sy,ey,sz,ez,meshx,meshy,meshz
  use grid, only: x,y,z,dx,dy,dz
  use hydro, only: state,pressr,gforce,set_state_pointer,NEW,OLD
  !use boundary
  use ionic, only: init_ionic

  implicit none

  ! Variables associated with the 1D Bertschinger infall solution
  integer,private :: ng
  real(kind=dp),dimension(:),allocatable,private :: lambdal
  real(kind=dp),dimension(:),allocatable,private :: dens1d
  real(kind=dp),dimension(:),allocatable,private :: vel1d
  real(kind=dp),dimension(:),allocatable,private :: mass1d

  ! TIS parameters
  ! Truncated Isothermal Sphere (TIS) parameters
  real(kind=dp),parameter,private :: zeta_t = 29.400d0
  real(kind=dp),parameter,private :: Mtil   = 61.48d0
  real(kind=dp),parameter,private :: etaTIS = 0.5544d0
  real(kind=dp),parameter,private :: etaSUS = 0.5d0  
  ! TIS radius in lambda=r/r_ta coordinates (i.e. in units of the 
  ! turnaround radius for the self-similar solution of Bertschinger 1985)
  real(kind=dp),parameter,private :: lambdat=0.3327339 
  ! TIS core radius in lambda units
  real(kind=dp),parameter,private :: lambda0=lambdat/zeta_t
  ! Overdensity for uniform tophat
  ! it is approximately EdS at high z; 
  ! replace with more precise value, if needed
  ! TIS mean overdensity
  real(kind=dp),parameter,private :: Delta_c = (etaSUS/etaTIS)**3*18.*pi**2
  
  ! Initialized parameters needed when calculating the (time dependent)
  ! gravitational force
  real(kind=dp),private :: xc,yc,zc
  real(kind=dp),private :: tcross, zcross,lambdacross,tinit
  real(kind=dp),private :: r0,densnorm

contains
  
  subroutine init_problem (restart)
    
    ! This routine initializes all hydro variables
    
    ! This may be a fresh start or a restart of a saved run
    
    ! Case: Mini-halo (cartesian coordinates)
    
    logical,intent(in) :: restart ! new run or a restart

    real(kind=dp) :: mass, zcoll, temper_out
    real(kind=dp) :: rho_crit_z, rt, rho0, sigma_V, tv
    real(kind=dp) :: tt,zz,lambdafactor,densnorm
    real(kind=dp) :: xl,yl,zl,lambda,xi,d,g, rhotot_b

    real(kind=dp)    :: r_interface ! dummy needed for calling init_ionic
    real(kind=dp),dimension(:,:,:),allocatable :: velr,temp

    integer :: lpos,lpos1
    real(kind=dp) :: dlpos
    real(kind=dp) :: rad
    integer :: i,j,k
    
#ifdef MPI       
    integer :: ierror
#endif

    if (.not.restart) then ! Fresh start
       
       ! Ask for the input if you are processor 0.
       if (rank == 0) then
          write (*,'(//,A,/)') '----- Mini-halo -----'
          write (*,'(A,$)') '1) Mass (solar masses): '
          read (stdinput,*) mass
          write (*,'(A,$)') '2) Redshift: '
          read (stdinput,*) zcoll
          write (*,'(A,$)') '3) Temperature outside (K): '
          read (stdinput,*) temper_out
          write (*,'(A,$)') '4) Position of centre x,y,z (cm): '
          read (stdinput,*) xc,yc,zc
       endif

       ! report input parameters
       if (rank == 0) then
          write(30,'(A)') & 
               'Problem: Mini-halo (cartesian)'
          write (30,'(//,A,/)') '----- Mini-halo -----'
          write (30,'(A,1PE10.3)') '1) Mass (solar masses): ',mass
          write (30,'(A,1PE10.3)') '2) Redshift: ',zcoll
          write (30,'(A,F8.3)')    '3) Temperature outside (K): ',temper_out
          write (30,'(A,3(1PE10.3))') '4) Position of centre x,y,z (cm): ', &
               xc,yc,zc
       endif
          
#ifdef MPI       
       ! Distribute the input parameters to the other nodes
       call MPI_BCAST(mass,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
            ierror)
       call MPI_BCAST(zcoll,1,MPI_DOUBLE_PRECISION,0, & 
            MPI_COMM_NEW,ierror)
       call MPI_BCAST(temper_out,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
            ierror)
       call MPI_BCAST(xc,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
            ierror)
       call MPI_BCAST(yc,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
            ierror)
       call MPI_BCAST(zc,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
            ierror)
#endif
       ! Scale
       xc=xc/SCLENG
       yc=yc/SCLENG
       zc=zc/SCLENG

       call infall_1d (meshx)
       
       ! critical density at redshift z
       rho_crit_z = rho_crit_0*(1.-Omega0+Omega0*(zcoll+1.)**3)

       ! TIS radius
       rt = (3.*mass*M_solar/(4.*pi*Delta_c*rho_crit_z))**(1./3.)
       ! TIS core radius
       r0  = rt/zeta_t
       ! TIS central density
       rho0 = zeta_t**3/3./Mtil*Delta_c*rho_crit_z 
       ! TIS velocity dispersion 
       sigma_V = sqrt(4.*pi*G_grav*rho0*r0**2)
       ! TIS virial temperature
       tv   = (mu*m_p/kb)*sigma_V**2

       write(*,*) 'check TIS',r0/kpc,rt/kpc,sigma_V/1.d5,tv

       ! redshift at which shock is crossed 
       zcross=1.0797*(zcoll+1.)-1.
       ! Cosmological time corresponding to redshift zcross
       tcross = 2.*(1.+zcross)**(-1.5)/(3.*H0*sqrt(Omega0))
       
       ! Initial time
       tinit=tcross

       ! Current time and redshift
       tt=tinit!+time
       zz=-1.0d0+(1.0d0+zcross)*(tcross/tt)**(2./3.)

       ! Effective radius of the mini-halo
       lambdacross=lambdat
       
       ! Factor for converting physical length to dimensionless length
       ! lambda (time dependent (through tt))
       lambdafactor=(lambda0/r0)*(tinit/tt)**(8./9.)

       ! normalization of the TIS profile
       densnorm=((1.+zcoll)/(1.+zcross))**3*rho0/rho_crit_z

       write(*,*) 'check',rho0/rho_crit_z,densnorm
    
       allocate(velr(sx-mbc:ex+mbc,sy-mbc:ey+mbc,sz-mbc:ez+mbc))
       allocate(temp(sx-mbc:ex+mbc,sy-mbc:ey+mbc,sz-mbc:ez+mbc))

        write(*,*) 'check TIS size',lambdacross,lambda0,lambdacross/lambda0
       ! Set the 3D solution (dimensionless)
       do k=sz-mbc,ez+mbc
          zl=(z(k)-zc)*lambdafactor*SCLENG
          do j=sy-mbc,ey+mbc
             yl=(y(j)-yc)*lambdafactor*SCLENG
             do i=sx-mbc,ex+mbc
                xl=(x(i)-xc)*lambdafactor*SCLENG
                lambda=sqrt(xl*xl+yl*yl+zl*zl)
                if (lambda <= lambdacross) then
                   ! Solution for the hydrostatic TIS minihalo 
                   xi=lambda/lambda0
                   call tis(xi,g,d)
                   state(i,j,k,RHO)=d*densnorm
                   velr(i,j,k)=0.
                   temp(i,j,k)=tv
                else 
                   ! infall solution: interpolate onto 3D uniform grid
                   call interpolate(lambda,lpos,lpos1,dlpos)
                   state(i,j,k,RHO)=dens1d(lpos)+ &
                        (dens1d(lpos1)-dens1d(lpos))*dlpos
                   velr(i,j,k)=vel1d(lpos)+ &
                        (vel1d(lpos1)-vel1d(lpos))*dlpos
                   temp(i,j,k)=temper_out
                end if
             end do
          end do
       end do

       ! baryon density at the current redshift (zz)
       ! Time dependent (via zz)
       rhotot_b=rho_crit_0*Omega_b*(1.0d0+zz)**3
       
       ! Set the 3D solution (dimensional)
       do k=sz-mbc,ez+mbc
          do j=sy-mbc,ey+mbc
             do i=sx-mbc,ex+mbc
                rad=sqrt((x(i)-xc)**2+(y(j)-yc)**2+(z(k)-zc)**2)
                state(i,j,k,RHO)=state(i,j,k,RHO)*rhotot_b/SCDENS
                velr(i,j,k)=velr(i,j,k)*r0/(lambda0*tcross)/SCVELO
                state(i,j,k,RHVX)=state(i,j,k,RHO)*velr(i,j,k)* &
                     (xc-x(i))/max(rad,1.0)
                state(i,j,k,RHVY)=state(i,j,k,RHO)*velr(i,j,k)* &
                     (yc-y(j))/max(rad,1.0)
                state(i,j,k,RHVZ)=state(i,j,k,RHO)*velr(i,j,k)* &
                     (zc-z(k))/max(rad,1.0)
                pressr(i,j,k)=state(i,j,k,RHO)*boltzm*temp(i,j,k)/mu
                state(i,j,k,EN)=pressr(i,j,k)/gamma1+ &
                        0.5d0*(state(i,j,k,RHVX)*state(i,j,k,RHVX)+ &
                        state(i,j,k,RHVY)*state(i,j,k,RHVY)+ &
                        state(i,j,k,RHVZ)*state(i,j,k,RHVZ))/state(i,j,k,RHO)
             end do
          end do
       end do

       ! Deallocate temporary arrays
       deallocate(velr)
       deallocate(temp)

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

  subroutine dm_grav_force (time)

    ! This subroutine calculates the gravitational force due to
    ! the dark matter profile.
    ! It is time dependent.

    real(kind=dp),intent(in) :: time

    real(kind=dp) :: tt,zz,lambdafactor
    real(kind=dp) :: zl,yl,xl,lambda,xi,massin,d,g,rhotot_b,rad,force
    integer :: lpos,lpos1
    real(kind=dp) :: dlpos
    integer :: i,j,k

    ! Time and redshift
    tt=tinit+time
    zz=-1.0d0+(1.0d0+zcross)*(tcross/tt)**(2./3.)

    ! Factor for converting physical length to dimensionless length
    ! lambda (time dependent (through tt))
    lambdafactor=(lambda0/r0)*(tinit/tt)**(8./9.)
    
    ! baryon density at the current redshift (zz)
    rhotot_b=rho_crit_0*Omega_b*(1.0d0+zz)**3

    ! Set the gravitational force.
    ! Two step procedure: dimensionless, followed by dimensional
    do k=sz-mbc,ez+mbc
       zl=(z(k)-zc)*lambdafactor*SCLENG
       do j=sy-mbc,ey+mbc
          yl=(y(j)-yc)*lambdafactor*SCLENG
          do i=sx-mbc,ex+mbc
             xl=(x(i)-xc)*lambdafactor*SCLENG
             ! Dimensionless interior mass
             lambda=sqrt(xl*xl+yl*yl+zl*zl)
             if (lambda <= lambdacross) then
                ! Solution for the hydrostatic TIS minihalo 
                xi=lambda/lambda0
                call tis(xi,g,d)
                massin=3.*lambda0*densnorm*lambda**2*g
             else 
                ! infall solution: interpolate onto 3D uniform grid
                call interpolate(lambda,lpos,lpos1,dlpos)
                massin=mass1d(lpos)+ &
                     (mass1d(lpos1)-mass1d(lpos))*dlpos
             endif

             ! Convert to dimensional mass and force
             rad=sqrt((1.d-10*(x(i)-xc))**2+(1.d-10*(y(j)-yc))**2+ &
                  (1.d-10*(z(k)-zc))**2)*1.d10
             !              print*,'check radius', rad
             massin=massin*(4.*pi/3.)*(rhotot_b*1d30)* &
                  (r0/1.d10/lambda0)**3!/M_solar
             ! Scaling:
             ! G is in cm^3 g^-1 s^-2, massin is in g
             ! G*massin is thus cm (cm/s)^2, hence the scaling.
             ! rad is already scaled.
             force=G_grav*massin/(SCVELO*SCVELO*SCLENG)
             force=force/max(rad,0.1*dx)**2
             gforce(i,j,k,1)=force*(xc-x(i))/max(rad,0.1*dx)
             gforce(i,j,k,2)=force*(yc-y(j))/max(rad,0.1*dx)
             gforce(i,j,k,3)=force*(zc-z(k))/max(rad,0.1*dx)

          end do
       end do
    end do

  end subroutine dm_grav_force

  !==========================================================================

  subroutine apply_grav_force(dt,newold)
    
    real(kind=dp),intent(in) :: dt
    integer,intent(in) :: newold

    integer :: i,j,k

    ! Point state to appropriate array
    state => set_state_pointer(newold)

    do k=sz-mbc,ez+mbc
       do j=sy-mbc,ey+mbc
          do i=sx-mbc,ex+mbc
             state(i,j,k,RHVX)=state(i,j,k,RHVX)+ &
                  dt*state(i,j,k,RHO)*gforce(i,j,k,1)
             state(i,j,k,RHVY)=state(i,j,k,RHVY)+ &
                  dt*state(i,j,k,RHO)*gforce(i,j,k,2)
             state(i,j,k,RHVZ)=state(i,j,k,RHVZ)+ &
                  dt*state(i,j,k,RHO)*gforce(i,j,k,3)
             state(i,j,k,EN)=state(i,j,k,EN)+ &
                  dt*(state(i,j,k,RHVX)*gforce(i,j,k,1)+ &
                  state(i,j,k,RHVY)*gforce(i,j,k,2)+ &
                  state(i,j,k,RHVZ)*gforce(i,j,k,3))
          enddo
       enddo
    enddo

  end subroutine apply_grav_force
    
  !==========================================================================

  subroutine infall_1d (size1d)

    ! One dimensional infall solution (Bertschinger).
    ! The solution is part of the module data

    integer,intent(in) :: size1d
    
    real(kind=dp) :: dtheta,theta,stheta,ctheta,s2theta
    real(kind=dp) :: lambda,v,ds,xi,beta,d,m
    integer :: i

    ! 1D Spherical solution of infall solution
  
    ng=size1d
    allocate(dens1d(ng))
    allocate(vel1d(ng))
    allocate(mass1d(ng))
    allocate(lambdal(ng))

    dtheta=2.*0.999*pi/real(ng,kind=dp) !radial grid size
    do i=1,ng
       theta=real(i,kind=dp)*dtheta
       stheta=sin(theta)
       ctheta=cos(theta)
       s2theta=(sin(0.5*theta))**2
       !note that this diverges for theta->0, but 
       !for small theta we would substitute the TIS solution anyway
       lambda=s2theta*(pi/(theta-stheta))**(8./9.)
       v=lambda*stheta*(theta-stheta)/(1.-ctheta)**2
       ds=(3./4.)*(theta-stheta)
       xi=1.-1.5*v/lambda
       beta=s2theta  
       d=ds**2/beta**3/(1.+3.*xi)
       m=4.5*lambda**2*d*((8./9.)*lambda-v)
       lambdal(i)=lambda !radius
       dens1d(i)=d       !density
       vel1d(i)=v        !velocity
       mass1d(i)=m       !enclosed mass 
       !       print*,'check',real(lambda),real(v),real(d),real(m)
       !       write(2,*) real(theta),real(lambda),real(v),real(d),real(m)
    enddo
    
  end subroutine infall_1d

  !==========================================================================

  subroutine interpolate(lambda,lpos,lpos1,dlpos)

    ! Find the interpolation parameters 

    implicit none

    real(kind=dp),intent(in) :: lambda
    integer,intent(out) :: lpos,lpos1
    real(kind=dp),intent(out) :: dlpos

    integer :: l,lp

    do l=1,ng-1
       lp=l+1
       if (lambdal(l) >= lambda .and. lambdal(lp) <= lambda) then
          lpos=l
          lpos1=lp
          dlpos=(lambda-lambdal(l))/(lambdal(lp)-lambdal(l))
       end if
    end do

  end subroutine interpolate

  !==========================================================================

  subroutine tis(xi,g,d)
     
    ! Subroutine providing an approximation to the TIS solution
    ! (dimensionless)

    ! xi = r/r0
    ! g  = gravity force
    ! d  = density
    real(kind=dp),intent(in) :: xi
    real(kind=dp),intent(out) :: g, d

    real(kind=dp), parameter :: aa=21.38,bb=19.81,a2=9.08,b2=14.62

    real(kind=dp) :: xi2
    
    !xi=r/r0 is the TIS radius in units of the core radius
    xi2=xi**2 

    !TIS density profile normalized to the central density 
    ! (Shapiro, Iliev and Raga 99, Appendix)
    d=aa/(a2+xi2)-bb/(b2+xi2)

    !gravity force (GM/r^2, in dimensionless units)
    g=2.*xi*(aa/(a2+xi2)**2-bb/(b2+xi2)**2)/d

  end subroutine tis

  !==========================================================================

  subroutine inflow (newold)
    
    ! This routine resets the inner boundary to the inflow condition
    
    integer,intent(in) :: newold

  end subroutine inflow
  
end module problem
