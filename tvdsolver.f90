module hydrosolver
  
  ! Module for Capreole 3D (f90)
  ! Author: Garrelt Mellema
  ! Date: 2007-10-09

  ! This module contains the relaxing TVD solver for solving 
  ! the Euler equations in three dimensions, including advected quantities.
  ! 
  ! The main routine is split into many subroutines for readability
  ! and to aid the compiler.

  ! This version has routines for constructing and destructing the arrays 
  ! needed for the solver, since otherwise they would be placed on the
  ! stack, which can cause problems.
  
  use precision
  use scaling
  use sizes
  use mesh
  use grid
  use atomic
  use geometry

  !use cgsconstants

  implicit none

  private
  
  real(kind=dp),parameter,private    :: HALF=0.5_dp ! 1/2
  real(kind=dp),parameter,private    :: ONE=1.0_dp  ! 1
  real(kind=dp),parameter,private    :: ZERO=0.0_dp ! 0
  real(kind=dp),parameter,private    :: TWO=2.0_dp ! 0

  ! Limiters
  integer,parameter,private :: NO_LIMITER=0
  integer,parameter,private :: MON_CEN=1
  integer,parameter,private :: VAN_LEER=2
  integer,parameter,private :: SUPERBEE=3
  integer,parameter,private :: MINMOD=4


  ! The following variables are public so that they can be filled
  ! or used outside of the solver. They correspond to the one
  ! spatial dimension
  real(kind=dp),dimension(:,:),allocatable,public :: state1d
  ! pressure
  real(kind=dp),dimension(:),allocatable,public :: wp 
  ! dstate will contain the state changes:
  real(kind=dp),dimension(:,:),allocatable,public :: dstate
  !$OMP THREADPRIVATE(state1d,wp,dstate)  

  ! The following variables are private since they are only
  ! needed inside the solver.
  real(kind=dp),dimension(:),allocatable :: c 
  real(kind=dp),dimension(:,:),allocatable :: w
  real(kind=dp),dimension(:,:),allocatable :: fu
  real(kind=dp),dimension(:,:),allocatable :: fl
  real(kind=dp),dimension(:,:),allocatable :: fr
  real(kind=dp),dimension(:,:),allocatable :: dfl
  real(kind=dp),dimension(:,:),allocatable :: dfr
  !$OMP THREADPRIVATE(c,fu,fl,fr,dfl,dfr)

  public :: constr_solver, destr_solver, solver
  
contains

  !============================================================================
  subroutine constr_solver (mesh)
    
    integer,intent(in) :: mesh

    ! Initializes solver variables

    allocate(state1d(1-mbc:mesh+mbc,neq))
    allocate(wp(1-mbc:mesh+mbc))
    allocate(dstate(1-mbc:mesh+mbc,neq))

    allocate(c(1-mbc:mesh+mbc))
    allocate(w(1-mbc:mesh+mbc,neq))
    allocate(fu(1-mbc:mesh+mbc,neq))
    allocate(fl(2-mbc:mesh+mbc,neq))
    allocate(fr(2-mbc:mesh+mbc,neq))
    allocate(dfl(2-mbc:mesh+mbc,neq))
    allocate(dfr(2-mbc:mesh+mbc,neq))

  end subroutine constr_solver

!-----------------------------------------------------------------------------

  subroutine destr_solver

    ! Destructs solver variables

    deallocate(state1d)
    deallocate(wp)
    deallocate(dstate)

    deallocate(c)
    deallocate(w)
    deallocate(fu)
    deallocate(fl)
    deallocate(fr)
    deallocate(dfl)
    deallocate(dfr)

  end subroutine destr_solver

!-----------------------------------------------------------------------------

  subroutine solver (mesh,dt,dx,dy,dz,V1,V2,V3,ij,ik,ierror)
    
    ! The TVD solver routine, supplies dstate
    
    integer,intent(in) :: mesh
    real(kind=dp),intent(in) :: dt,dx,dy,dz
    integer,intent(in) :: ij,ik
    integer,intent(in) :: V1,V2,V3 ! indices of the the two velocities:
                                        ! V1: integration direction
                                        ! V2,V3: perpendicular direction
    integer,intent(out) :: ierror ! control integer

    integer :: i
    real(kind=dp) :: dtdx!,dtdy

    !------------------------------------------------------------------------

    ierror=0 ! initialise the control variable to 0

    dtdx=dt/dx

    call calculate_fluxes (state1d,V1,V2,V3)

    !! Do half step using first-order upwind scheme
    fr=(state1d*spread(c,2,neuler)+w)*HALF
    fl=cshift(state1d*spread(c,2,neuler)-w,1,1)*HALF

    fu=(fr-fl)
    dstate=-(fu-cshift(fu,-1,1))*HALF*dtdx
    
    !! Do full step using second-order TVD scheme
    call calculate_fluxes (state1d+dstate,V1,V2,V3)

    !! Right-moving waves
    fr=((state1d+dstate)*spread(c,1,5)+w)*HALF
    dfl=(fr-cshift(fr,-1,2))*HALF
    dfr=cshift(dfl,1,2)
    fr=fr+limiter(dfl,dfr,VAN_LEER)

    !! Left-moving waves
    fl=cshift((state1d+dstate)*spread(c,1,5)-w,1,2)*HALF
    dfl=(cshift(fl,-1,2)-fl)*HALF
    dfr=cshift(dfl,1,2)
    fl=fl+limiter(dfl,dfr,VAN_LEER)

    fu=(fr-fl)
    dstate=-(fu-cshift(fu,-1,2))*dtdx

    !ipres_error=pressure_fix ()
    
    !--------------------------------------------------------------------------

  contains
    
    !--------------------------------------------------------------------------

    subroutine calculate_fluxes (u,V1,V2,V3)

      real(kind=dp),dimension(1-mbc:mesh+mbc,neuler),intent(in) :: u
      integer,intent(in) :: V1,V2,V3

      integer :: i,ieq

    ! ---------------------------------------------------------------------
    ! calculate the fluxes at the cell centre
    ! ---------------------------------------------------------------------
    do i=1-mbc,mesh+mbc
       w(i,RHO)=u(i,V1)
       w(i,V1)=u(i,V1)*u(i,V1)/u(i,RHO) + wp(i)
       w(i,V2)=u(i,V1)*u(i,V2)/u(i,RHO)
       w(i,V3)=u(i,V1)*u(i,V3)/u(i,RHO)
       w(i,EN)=u(i,V1)*(wp(i)+u(i,EN))/u(i,RHO)
       c(i)=abs(u(i,V1)/u(i,RHO))+ &
            sqrt(gamma*wp(i)/u(i,RHO))
    enddo
    do ieq=neuler+1,neq
       do i=1-mbc,mesh+mbc
          w(i,ieq)=u(i,ieq)*u(i,V1)/u(i,RHO)
       enddo
    enddo

    end subroutine calculate_fluxes

    !--------------------------------------------------------------------------
    
    function limiter(a,b,limfunc)
      
      ! Limiter function used in flux limiting.
      ! It limits a ratio using the function limfunc.
      
      real(kind=dp),dimension(1-mbc:mesh+mbc,neq) :: limiter
      
      ! superbee parameter (1 to 2)
      real(kind=dp),parameter :: sbpar=1.8
      
      real(kind=dp),dimension(1-mbc:mesh+mbc,neq),intent(in) :: a,b
      integer,intent(in) :: limfunc
      
      select case (limfunc)
      case (NO_LIMITER)
         limiter = ZERO
      case(VAN_LEER)
         ! van Leer
         limiter=TWO*max(ZERO,a*b)/(a+b)
      case(MINMOD)
         ! MinMod
         limiter=(sign(HALF,a)+sign(HALF,b))*min(abs(a),abs(b))
      case(SUPERBEE)
         ! Superbee (sbpar)
         where (abs(a) > abs(b))
            limiter=(sign(HALF,a)+sign(HALF,b))*min(abs(a),abs(sbpar*b))
         elsewhere
            limiter=(sign(HALF,a)+sign(HALF,b))*min(abs(sbpar*a),abs(b))
         endwhere
      ! This one needs to be fixed
      !case(MON_CEN)
      !   ! monotinized centered 
      !   c = (1.d0 + ratio)/2.d0
      !   limiter = max(0.d0, min(c, 2.d0, 2.d0*ratio))
      end select
      
    end function limiter

  end subroutine solver
    
end module hydrosolver





