module sizes

  ! Module for Capreole (f90)
  ! Author: Garrelt Mellema
  ! Date: 2003-06-01
  ! This module is also accepted by the F compiler (Dec 9, 2003)

  ! This module contains basic parameter definitions

  use precision
  private

  integer,public,parameter :: nrOfDim=3        ! number of dimensions
  integer,public,parameter :: mbc=2            ! number of ghost cells
  integer,public,parameter :: neuler=2+NrOfDim ! number of Euler equations
  integer,public,parameter :: neq=neuler+1     ! number of equations (Euler +
                                               !   advected quantities)
  ! Indices of state array
  integer,public,parameter :: RHO=1
  integer,public,parameter :: RHVX=2
  integer,public,parameter :: RHVY=3
  integer,public,parameter :: RHVZ=4
  integer,public,parameter :: EN=5

  ! The following define constants for identifying the coordinate system
  integer,public,parameter :: CART =1
  integer,public,parameter :: CYL  =2
  integer,public,parameter :: SPHER=3

end module sizes
