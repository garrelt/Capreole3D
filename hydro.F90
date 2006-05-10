module hydro

  ! Module for Capreole (3D)
  ! Author: Garrelt Mellema
  ! Date: 2005-04-20 (previous 2004-05-10)
  ! This module needs to be checked for the F compiler

  ! This module contains the routine which allocates the hydro
  ! arrays. This way they do not end up on stack
  
  use precision
  use sizes
  use mesh

  implicit none
  private

  real(kind=dp),dimension(:,:,:,:),allocatable,target,public :: stold
  real(kind=dp),dimension(:,:,:,:),allocatable,target,public :: stnew
  real(kind=dp),dimension(:,:,:),allocatable,public   :: pressr

  ! These pointers are either used as a generic name (state)
  ! in some routines
  ! or as temporary storage during initialization.
  real(kind=dp),pointer,dimension(:,:,:,:),public :: state
  real(kind=dp),pointer,dimension(:,:,:,:),public :: tmpstate

  integer,parameter,public :: NEW=1
  integer,parameter,public :: OLD=0

  public :: init_hydro,set_state_pointer

contains

  subroutine init_hydro (restart)
    
    ! This routine allocates the hydrodynamic variables
    
    logical,intent(in) :: restart

    ! Allocate the arrays
    allocate(stold(sx-mbc:ex+mbc,sy-mbc:ey+mbc,sz-mbc:ez+mbc,neq))
    allocate(stnew(sx-mbc:ex+mbc,sy-mbc:ey+mbc,sz-mbc:ez+mbc,neq))
    allocate(pressr(sx-mbc:ex+mbc,sy-mbc:ey+mbc,sz-mbc:ez+mbc))

    ! point generic state variable to stnew
    state => stnew

  end subroutine init_hydro

  function set_state_pointer (whichone) result(state_pointer_result)

    integer,intent(in) :: whichone
    real(kind=dp),pointer,dimension(:,:,:,:) :: state_pointer_result

    ! Point state to appropriate array
    if (whichone == NEW) then
       state_pointer_result => stnew
    else
       state_pointer_result => stold
    endif

  end function set_state_pointer

end module hydro
