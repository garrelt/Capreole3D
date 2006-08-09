module geometry

  ! Module for Capreole (f90)
  ! Author: Garrelt Mellema
  ! Date: 2004-05-12 (previous 2003-06-01)
  ! This module is also accepted by the F compiler (Dec 9, 2003)
  !
  ! This module contains the routines related to or dependent on
  ! the choice of the 2D coordinate system.
  !
  ! Version: cartesian coordinates (x,y).
  ! 
  ! Contents:
  ! timestep : calculates the CFL time step
  ! xg1g2    : sets the geometric factors for the 1D hydro solver (x-direction)
  ! yg1g2    : sets the geometric factors for the 1D hydro solver (y-direction)
  ! xgeosource : sets the geometric source terms (x-direction)
  ! ygeosource : sets the geometric source terms (y-direction)
  ! presfunc : calculates the pressure from the state variables
  !
  ! History:
  ! 2004-05-12: adapted for new approach without large (dynamic) arrays
  !             as arguments to timestep and presfunc. state is now a
  !             pointer, pointing to stnew or stold (see hydro.F90).

  use precision
  use scaling
  use sizes
  use mesh
  use grid
  use atomic
  use hydro
  use times

  private

  real(kind=dp),public :: maxdt=1.0e20_dp ! non CFL limit on time step

  public :: timestep,presfunc

contains

  function absvel(st0d) result(absvel_result)

    ! Function to calculate absolute velocity
    ! Version: cartesian 3D

    real(kind=dp) :: absvel_result
    real(kind=dp),dimension(neq) :: st0d

    absvel_result=(st0d(RHVX)*st0d(RHVX)+st0d(RHVY)*st0d(RHVY)+ &
         st0d(RHVZ)*st0d(RHVZ))/(st0d(RHO)*st0d(RHO))

  end function absvel

!--------------------------------------------------------------------

  subroutine presfunc(ist,ifi,jst,jfi,kst,kfi,newold,ierror)
    
    real(kind=dp),parameter :: precision=1.0d-13

    integer,intent(in) :: ist,ifi,jst,jfi,kst,kfi,newold
    integer,intent(out) :: ierror

    integer i,j,k

    ierror=0 ! nullify error variable

    ! Point state to appropriate array
    state => set_state_pointer(newold)

    do k=kst,kfi
       do j=jst,jfi
          do i=ist,ifi
             pressr(i,j,k)=state(i,j,k,EN)-0.5*state(i,j,k,RHO)* &
                  absvel(state(i,j,k,:))
             if (abs(pressr(i,j,k)/state(i,j,k,EN)) < precision) then
                write(30,'(A,2(1PE10.3,X),A,3(I4,X),A,E10.3)') &
                     'PRECISION ISSUE: ',pressr(i,j,k), &
                     state(i,j,k,EN),' at ',i,j,k,' time = ',time
                ierror=1
             endif
             pressr(i,j,k)=gamma1*pressr(i,j,k)
             if (pressr(i,j,k) <= 0.0_dp) then
                !!write(30,'(A,2(1PE10.3,X),A,3(I4,X),A,E10.3)') &
                !!     'Presfunc reports negative pressure: ', &
                !!     pressr(i,j,k),state(i,j,k,EN),' at ',i,j,k, &
                !!     ' time = ',time
                ierror=ierror+1
             endif
          enddo
       enddo
    enddo
    
  end subroutine presfunc

!------------------------------------------------------------------------------

  function timestep(cfl,newold) result(timestep_result)

    ! Function to calculate maximum time step from the CFL condition
    ! Version: spherical coordinates r, theta

    real(kind=dp) :: timestep_result

    real(kind=dp),intent(in) :: cfl
    integer,intent(in) :: newold

    real(kind=dp) :: vs ! sound speed
    real(kind=dp) :: vabs ! absolute velocities

    integer :: i,j,k ! grid loop counters
    real(kind=dp) :: vmax ! maximum signal velocity
    real(kind=dp) :: dcell ! cell size

    ! Point state to appropriate array
    state => set_state_pointer(newold)
    
    vmax=-1.0d0
    dcell=min(dy,dx,dz)
    
    do k=sz,ez
       do j=sy,ey
          do i=sx,ex
             vs=sqrt(max(0.0d0,gamma*pressr(i,j,k)/state(i,j,k,RHO)))
             vabs=sqrt(absvel(state(i,j,k,:)))
             vmax=max(vmax,(vs+vabs)/dcell)
          enddo
       enddo
    enddo

    timestep_result=min(maxdt,cfl/vmax)

  end function timestep

end module geometry
