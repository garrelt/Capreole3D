module oddeven

  ! Module for Capreole (2D)
  ! Author: Garrelt Mellema
  ! Date: 2003-08-26
  !
  ! This module contains the routines related to the correction
  ! of the odd-even decoupling / carbuncle phenomenon
 

  use precision
  use scaling
  use sizes
  use grid
  use coords
  use atomic
  use geometry

  private
  
  ! Shock detection parameter
  real(kind=dp),parameter,private :: shockDetectFactor = 0.5_dp
  ! The diffusion coefficient used
  real(kind=dp),parameter,private :: eta=0.05_dp

  integer,private :: i,j,k,ieq,m
  real(kind=dp),private :: deltaPressure

  public :: odd_even
!  private :: set_odd_even_x,set_odd_even_y,apply_odd_even_x,apply_odd_even_y
!  private :: detect_shock_x,detect_shock_y

contains
  !========================================================================
  
  subroutine odd_even (state,pressr,vol,time,action)

    ! This routine handles the odd-even fix interface with the
    ! integration routine. The actual work is done in the set and
    ! apply routines for the different coordinate directions.

    real(kind=dp),dimension(sx-mbc:ex+mbc,sy-mbc:ey+mbc,neq), &
         intent(inout) :: state
    real(kind=dp),dimension(sx-mbc:ex+mbc,sy-mbc:ey+mbc), &
         intent(inout)     :: pressr
    real(kind=dp),dimension(sx-mbc:ex+mbc,sy-mbc:ey+mbc), &
         intent(in)        :: vol
    real(kind=dp),intent(in)     :: time
    integer :: action ! what to do:
    ! 1: fix x-direction
    ! 2: fix y-direction
    ! 0: fix all directions
    ! other values: do nothing

    ! Flag for shock detection
    ! ( exported via state(neq) )
    integer,dimension(sx-mbc:ex+mbc,sy-mbc:ey+mbc) :: flag

    ! Diffusive fluxes
    real(kind=dp),dimension(sx-mbc:ex+mbc,sy-mbc:ey+mbc,neq) :: fdiff
    real(kind=dp),dimension(sx-mbc:ex+mbc,sy-mbc:ey+mbc,neq) :: gdiff

    ! Counters
    integer :: nrOfXCorrections,nrOfYCorrections

    ! Figure out what to do.
    select case (action)
    case (0)
       call set_odd_even_x ()
       state(sx:ex,sy:ey,neq)=state(sx:ex,sy:ey,neq)+real(flag(sx:ex,sy:ey),dp)
       call set_odd_even_y ()
       state(sx:ex,sy:ey,neq)=state(sx:ex,sy:ey,neq)+real(flag(sx:ex,sy:ey),dp)
       call apply_odd_even_x ()
       call apply_odd_even_y ()
    case (1)
       call set_odd_even_x ()
       call apply_odd_even_x ()
       state(sx:ex,sy:ey,neq)=state(sx:ex,sy:ey,neq)+real(flag(sx:ex,sy:ey),dp)
    case (2)
       call set_odd_even_y ()
       call apply_odd_even_y ()
       state(sx:ex,sy:ey,neq)=state(sx:ex,sy:ey,neq)+real(flag(sx:ex,sy:ey),dp)
    case default
       ! Do nothing
    end select

  contains

    subroutine set_odd_even_x ()
      
      ! This routine find shocks in the perpendicular direction,
      ! then looks for density oscillations parallel to those shocks,
      ! and calculates a diffusive flux if found.
      
      integer :: nrOfFlags,cellCount
      integer :: imin,iplus
      
      ! Reset flags
      flag(:,:) = 0
      cellCount=0
      nrOfXCorrections=0
      
      ! Reset diffusive flux
      fdiff(:,:,:)=0.0_dp

      ! Find the shocks in y direction, set flag to 1.
      nrOfFlags=detect_shock_y (1)

      ! Test for odd-even pattern.
      ! Use 4-point up and down pattern for this
      do j=sy,ey
         do i=sx,ex
            cellCount = 0
            do k=-2,1
               if(flag(i+k,j) == 1 .or. flag(i+k,j) == 2) & 
                    cellCount = cellCount + 1
            enddo
            if(cellCount == 4) then
               if(state(i-2,j,1) > state(i-1,j,1).and. &
                    state(i-1,j,1) < state(i  ,j,1).and. &
                    state(i  ,j,1) > state(i+1,j,1)) & 
                    flag(i,j) = 2
               if(state(i-2,j,1) < state(i-1,j,1).and. & 
                    state(i-1,j,1) > state(i  ,j,1).and.  &
                    state(i  ,j,1) < state(i+1,j,1)) & 
                    flag(i,j) = 2
            endif
         enddo
      enddo
      
      do j=sy,ey
         do i=sx,ex
            if( flag(i,j) == 2 ) then
               nrOfXCorrections = nrOfXCorrections + 1
               imin=max(i-1,1) ! diffusive flux at edge 1 is zero
               iplus=min(i+1,meshx) ! diffusive flux at edge 1 is zero
               do ieq=1,neq
                  fdiff(i,j,ieq)=(state(imin,j,ieq)-state(i,j,ieq))
                  fdiff(i+1,j,ieq)=(-state(iplus,j,ieq)+state(i,j,ieq))
               enddo 
            end if 
         end do 
      enddo 
       
    end subroutine set_odd_even_x
    
    !--------------------------------------------------------------------------

    subroutine apply_odd_even_x ()

      ! Apply the diffusive flux
      
      if ( nrOfXCorrections > 0 ) then
         do ieq=1,neq
            do j=sy,ey
               do i=sx,ex
                  state(i,j,ieq)=state(i,j,ieq) + &
                       eta*(fdiff(i,j,ieq)-fdiff(i+1,j,ieq))
               enddo
            enddo
         enddo
      end if
      
    end subroutine apply_odd_even_x
    
    !--------------------------------------------------------------------------
    
    subroutine set_odd_even_y ()
      
      ! This routine find shocks in the perpendicular direction,
      ! then looks for density oscillations parallel to those shocks,
      ! and calculates a diffusive flux if found.
      
      integer :: nrOfFlags,cellCount
      integer :: jmin,jplus
      
      ! Reset flags
      flag(:,:) = 0
      cellCount=0
      nrOfYCorrections=0
      
      ! Reset diffusive flux
      gdiff(:,:,:)=0.0_dp
      
      ! Find the shocks in x direction, set flag to 3
      nrOfFlags=detect_shock_x (3)

      ! Test for odd-even pattern.
      ! Use 4-point up and down pattern for this
      do j=sy,ey
         do i=sx,ex
            cellCount = 0
            do k=-2,1
               if(flag(i,j+k) == 3 .or. flag(i,j+k) == 4) & 
                    cellCount = cellCount + 1
            enddo
            if(cellCount == 4) then
               if(state(i,j-2,1) > state(i,j-1,1).and. &
                    state(i,j-1,1) < state(i,j  ,1).and. &
                    state(i,j  ,1) > state(i,j+1,1)) &
                    flag(i,j) = 4
               if(state(i,j-2,1) < state(i,j-1,1).and. &
                    state(i,j-1,1) > state(i,j  ,1).and.  & 
                    state(i,j  ,1) < state(i,j+1,1)) & 
                    flag(i,j) = 4
            endif
         enddo
      enddo

      do j=sy,ey
         do i=sx,ex
            if(flag(i,j) == 4) then
               nrOfYCorrections = nrOfYCorrections + 1
               jmin=max(j-1,1)      ! diffusive flux at edge 1 is zero
               jplus=min(j+1,meshy) ! diffusive flux at edge 1 is zero
               !if (j == 169) write(*,*) i,j, nrOfYCorrections
               do ieq=1,neq
                  gdiff(i,j,ieq)=(state(i,jmin,ieq)-state(i,j,ieq))
                  gdiff(i,j+1,ieq)=(-state(i,jplus,ieq)+state(i,j,ieq))
               end do
            end if
         end do
      end do
    
    end subroutine set_odd_even_y

    !--------------------------------------------------------------------------
    
    subroutine apply_odd_even_y ()
      
      ! Apply the diffusive flux
      
      if (nrOfYCorrections > 0) then
         do ieq=1,neq
            do j=sy,ey
               do i=sx,ex
                  state(i,j,ieq)=state(i,j,ieq) + & 
                       eta*(gdiff(i,j,ieq)-gdiff(i,j+1,ieq))
               enddo
            enddo
         enddo
         
      end if
      
    end subroutine apply_odd_even_y
    
    !--------------------------------------------------------------------------
    
    function detect_shock_x (marker) result(shock_x)
      
      ! Detects shocks in x-direction and flags post-shock region with marker
      ! This is done by looking at the pressure jumps.

      ! Output: number of cells flagged
      
      integer :: shock_x
      
      integer,intent(in) :: marker
      
      integer :: i1,i2,i3,i4,i5

      shock_x=0

      do j=sy,ey
         do i=sx,ex
            deltaPressure = (pressr(i+1,j)-pressr(i,j)) / &
                 max(pressr(i+1,j),pressr(i,j))
            ! Right going shock
            if(-(deltaPressure) > shockDetectFactor) then
               i4=max(i-4,sx-mbc)
               i3=max(i-3,sx-mbc)
               i2=max(i-2,sx-mbc)
               i1=max(i-1,sx-mbc)
               flag(i4,j) = marker
               flag(i3,j) = marker
               flag(i2,j) = marker
               flag(i1,j) = marker
               flag(i,  j) = marker
               shock_x = shock_x + 5
            endif
            ! Left going shock
            if((deltaPressure) > shockDetectFactor) then
               i5=min(i+5,ex+mbc)
               i4=min(i+4,ex+mbc)
               i3=min(i+3,ex+mbc)
               i2=min(i+2,ex+mbc)
               i1=min(i+1,ex+mbc)
               flag(i5,j) = marker
               flag(i4,j) = marker
               flag(i3,j) = marker
               flag(i2,j) = marker
               flag(i1,j) = marker
               shock_x = shock_x + 5
            endif
         enddo
      enddo

    end function detect_shock_x
    
    !--------------------------------------------------------------------------
    
    function detect_shock_y (marker) result(shock_y)
      
      ! Detects shocks in y-direction and flags post-shock region with marker
      ! This is done by looking at the pressure jumps.
      
      ! Output: number of cells flagged
      
      integer :: shock_y
      
      integer,intent(in) :: marker

      integer :: j1,j2,j3,j4,j5
      shock_y=0
      do j=sy,ey
         do i=sx,ex
            deltaPressure = (pressr(i,j+1)-pressr(i,j)) / &
                 max(pressr(i,j+1),pressr(i,j))
            ! Right going shock
            if(-(deltaPressure) > shockDetectFactor) then
               j4=max(j-4,sy-mbc)
               j3=max(j-3,sy-mbc)
               j2=max(j-2,sy-mbc)
               j1=max(j-1,sy-mbc)
               flag(i,j4) = marker
               flag(i,j3) = marker
               flag(i,j2) = marker
               flag(i,j1) = marker
               flag(i,  j) = marker
               shock_y = shock_y + 5
            endif
            ! Left going shock
            if((deltaPressure) > shockDetectFactor) then
               j5=min(j+5,ey+mbc)
               j4=min(j+4,ey+mbc)
               j3=min(j+3,ey+mbc)
               j2=min(j+2,ey+mbc)
               j1=min(j+1,ey+mbc)
               flag(i,j5) = marker
               flag(i,j4) = marker
               flag(i,j3) = marker
               flag(i,j2) = marker
               flag(i,j1) = marker
               shock_y = shock_y + 5
            endif
         enddo
      enddo
      
    end function detect_shock_y
    
  end subroutine odd_even
    
end module oddeven
