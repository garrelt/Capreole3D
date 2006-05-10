module protection

  ! Module for Capreole 
  ! Author: Garrelt Mellema
  ! Date: 2004-05-11
  !
  ! This module contains the routines related to the negative
  ! pressure / density / total energy correction

  use precision
  use scaling
  use sizes
  use mesh
  use grid
  use hydro
  use times
  use atomic
  use geometry

  implicit none
  
contains
  !========================================================================

  subroutine presprot(inegative,icall,newold)

    ! This routine protects the pressure, density, and energy density 
    ! from becoming negative.
    ! We do this in two stages, if 1) fails, do 2).
    ! 1) diffuse  
    ! 2) setting a minimal temperature tmin
    
    ! the diffusion coefficient used in stage 1)
    real(kind=dp),parameter :: eta=0.05d0 

    ! the minimum temperature applied in stage 2)
    real(kind=dp),parameter :: tmin=1.0d-10

    integer,intent(out) :: inegative ! control integer; /= 0 if fixes failed
    integer,intent(in)  :: icall ! call id (to distinguish between different
                                 !  calls to this subroutine)
    integer,intent(in) :: newold
    integer :: i,j,k,ieq
    real(kind=dp)    :: pnew
    real(kind=dp),dimension(:,:,:,:),allocatable,save :: fdiff
    real(kind=dp),dimension(:,:,:,:),allocatable,save :: gdiff
    real(kind=dp),dimension(:,:,:,:),allocatable,save :: hdiff
    !real(kind=dp),dimension(sx-1:ex+2,sy-1:ey+2,neq) :: fdiff
    !real(kind=dp),dimension(sx-1:ex+2,sy-1:ey+2,neq) :: gdiff
    logical :: problem
    integer :: ix,jx,imin,jmin,iplus,jplus,kmin,kplus
    integer :: ierror
    
    ! Point state to appropriate array
    state=set_state_pointer(newold)

    inegative=0
    problem=.false.

    if (.not.(allocated(fdiff))) allocate(fdiff(sx-1:ex+2,sy-1:ey+2, &
         sz-1:ez+2,neq))
    if (.not.(allocated(gdiff))) allocate(gdiff(sx-1:ex+2,sy-1:ey+2, &
         sz-1:ez+2,neq))
    if (.not.(allocated(hdiff))) allocate(hdiff(sx-1:ex+2,sy-1:ey+2, &
         sz-1:ez+2,neq))
    fdiff(:,:,:,:)=0.0d0
    gdiff(:,:,:,:)=0.0d0
    hdiff(:,:,:,:)=0.0d0

    do k=sz,ez
       do j=sy,ey
          do i=sx,ex
             ! Check for negative pressure/density/energy
             if (pressr(i,j,k) <= 0.0.or.state(i,j,k,RHO) <= 0.0.or. &
                  state(i,j,k,EN) <= 0.0) then
                ! Report to log file
                write(30,'(A,1PE10.3,2(2X,E10.3),$)') &
                     'Negative pressure/density/energy: ', &
                     pressr(i,j,k),state(i,j,k,RHO),state(i,j,k,EN)
                write(30,'(A,3(I4,X),A,1PE10.3)') ' at ', &
                     i,j,k,' time = ',time
                write(30,*) 'call ',icall
                call flush(30)
                ! Set a control variable
                problem=.true.

                ! Pressure fix 1: diffuse with four neighbours
                ! diffusion parameter is eta (used below)
                ! To keep things conservative, express everything
                ! as fluxes
                imin=max(i-1,1) ! diffusive flux at edge 1 is zero
                iplus=min(i+1,meshx) ! diffusive flux at edge 1 is zero
                jmin=max(j-1,1)
                jplus=min(j+1,meshy)
                kmin=max(k-1,1)
                kplus=min(k+1,meshz)
                ! Only set diffuse flux if the neighbouring cell has
                ! enough positive pressure itself, or if we are
                ! correcting a negative density or energy.
                do ieq=1,neq
                   if (-(pressr(imin,j,k)/pressr(i,j,k)) > eta .or. &
                        pressr(i,j,k) > 0.0d0 )  &
                        fdiff(i,j,k,ieq)=(state(imin,j,k,ieq)-state(i,j,k,ieq))
                   if (-(pressr(iplus,j,k)/pressr(i,j,k)) > eta .or. &
                        pressr(i,j,k) > 0.0d0 ) &
                        fdiff(i+1,j,k,ieq)=(-state(iplus,j,k,ieq)+&
                        state(i,j,k,ieq))
                   if (-(pressr(i,jmin,k)/pressr(i,j,k)) > eta .or. &
                        pressr(i,j,k) > 0.0d0 )  &
                        gdiff(i,j,k,ieq)=(state(i,jmin,k,ieq)-state(i,j,k,ieq))
                   if (-(pressr(i,jplus,k)/pressr(i,j,k)) > eta .or. &
                        pressr(i,j,k) > 0.0d0 ) &
                        gdiff(i,j+1,k,ieq)=(-state(i,jplus,k,ieq)+&
                        state(i,j,k,ieq))
                   if (-(pressr(i,j,kmin)/pressr(i,j,k)) > eta .or.  &
                        pressr(i,j,k) > 0.0d0 )  &
                        hdiff(i,j,k,ieq)=(state(i,j,kmin,ieq)-state(i,j,k,ieq))
                   if (-(pressr(i,j,kplus)/pressr(i,j,k)) > eta .or. &
                        pressr(i,j,k) > 0.0d0 ) &
                        hdiff(i,j,k+1,ieq)=(-state(i,j,kplus,ieq)+ &
                        state(i,j,k,ieq))
                enddo
             endif
          enddo
       enddo
    enddo

    if (problem) then
       ! Apply fluxes
       do ieq=1,neq
          do k=sz,ez
             do j=sy,ey
                do i=sx,ex
                   state(i,j,k,ieq)=state(i,j,k,ieq)+eta* &
                        (fdiff(i,j,k,ieq)-fdiff(i+1,j,k,ieq) + &
                         gdiff(i,j,k,ieq)-gdiff(i,j+1,k,ieq) + &
                         hdiff(i,j,k,ieq)-hdiff(i,j,k+1,ieq))
                enddo
             enddo
          enddo
       enddo

       ! Recalculate the pressure after the diffusion
       ! the presfunc routine needs to be supplied
       call presfunc(sx,ex,sy,ey,sz,ez,newold,ierror)
       
       do k=sz,ez
          do j=sy,ey
             do i=sx,ex
                ! Check if the pressure is still negative
                if (pressr(i,j,k) <= 0.0) then ! check result of fix 1
                   write(30,*) 'Still Negative pressure: ', &
                        pressr(i,j,k),' at ',i,j,k
                   call flush(30)

                   ! Pressure fix 2: set temperature to minimum value
                   pnew=state(i,j,k,RHO)*boltzm*tmin/xmu
                   state(i,j,k,EN)=state(i,j,k,EN)+(pnew-pressr(i,j,k))/gamma1
                   pressr(i,j,k)=pnew
                endif
                
                ! Check for negative densities and energies
                ! These are fatal. inegative is used to communicate
                ! this to the calling program.
                if (state(i,j,k,RHO) <= 0.0) then
                   write(30,*) 'Still negative density: ', &
                        state(i,j,k,RHO),' at ',i,j,k
                   inegative=1
                endif
                if (state(i,j,k,EN) <= 0.0) then
                   write(30,*) 'Still negative energy: ', &
                        state(i,j,k,EN),' at ',i,j,k
                   inegative=2
                endif
             enddo
          enddo
       enddo
    endif

  end subroutine presprot

end module protection
