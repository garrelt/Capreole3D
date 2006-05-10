module grid

  ! Module for Capreole (3D)
  ! Author: Garrelt Mellema
  ! Date: 2005-04-20 (previous 2004-05-11, 2003-06-01)
  ! This module is also accepted by the F compiler (Dec 9, 2003)

  ! This module contains routines dealing with the physical 
  ! coordinate system
  ! - initcoords sets up the coordinate system
  ! - inner/outer x/y/z bound takes care of the boundary conditions

  ! Version: cartesian coordinates (x,y,z)

  ! History:
  ! 2004-05-11 - adapted to the new approach of not passing large
  !              arrays as subroutine arguments. To use the generic
  !              name vol, it is now a pointer which can point to 
  !              volx, voly, volz
  ! 2005-04-20 - adapted for 3D
  
  use precision
  use scaling
  use sizes
  use my_mpi
  use mesh

  private
  !implicit none

  ! Identify type of coordinate system
  integer,parameter,public :: coordsys=CART

  ! dx,dy - cell sizes
  ! xlength,ylength - length of the entire grid
  real(kind=dp),public :: dx,dy,dz,xlength,ylength,zlength

  ! x : x-coordinate
  ! y : y-coordinate
  ! z : z-coordinate
  real(kind=dp),dimension(:),allocatable,public :: x
  real(kind=dp),dimension(:),allocatable,public :: y
  real(kind=dp),dimension(:),allocatable,public :: z

  ! volx : volume factor (for x-integration)
  ! voly : volume factor (for y-integration)
  ! volz : volume factor (for z-integration)
  ! vol  : generic name for volume, can point to volx or voly
  real(kind=dp),dimension(:,:,:),allocatable,target,public :: volx
  real(kind=dp),dimension(:,:,:),allocatable,target,public :: voly
  real(kind=dp),dimension(:,:,:),allocatable,target,public :: volz
  real(kind=dp),pointer,dimension(:,:,:),public :: vol

  real(kind=dp),dimension(2),public :: xedge
  real(kind=dp),dimension(2),public :: yedge
  real(kind=dp),dimension(2),public :: zedge

  public :: init_coords

contains

  subroutine init_grid ( )
    
    ! Allocate the arrays
    allocate(x(sx-mbc:ex+mbc))
    allocate(y(sy-mbc:ey+mbc))
    allocate(z(sz-mbc:ez+mbc))
    allocate(volx(sx-mbc:ex+mbc,sy-mbc:ey+mbc,sz-mbc:ez+mbc))
    allocate(voly(sx-mbc:ex+mbc,sy-mbc:ey+mbc,sz-mbc:ez+mbc))
    allocate(volz(sx-mbc:ex+mbc,sy-mbc:ey+mbc,sz-mbc:ez+mbc))
    
    ! point generic volume variable to volx
    vol => volx

  end subroutine init_grid

  subroutine init_coords (restart)
    
    ! This routine initializes the coordinate variables
    
    ! This may be a fresh start or a restart of a saved run
    
    ! Case: cartesian coordinates (x,y,z)
    
    logical,intent(in) :: restart

    integer :: i,j,k,ierror

    if (.not.restart) then ! Fresh start
       
       ! Ask for the input if you are processor 0.
       
       if (rank == 0) then
          write (unit=*,fmt="(a)",advance="no") "2) Size of grid box (cm): "
          read (unit=*,fmt=*) xlength,ylength,zlength
          write (unit=30,fmt="(a,2(es10.3))") "2) Size of grid box(cm): ", &
               xlength,ylength,zlength
          ! Record variables which remain constant during a run in a file
          ! runparams to be used at restarts
          open(unit=80,file="runparams",status="unknown",form="unformatted",&
               action="write")
          write(unit=80) xlength,ylength,zlength
       endif
    else
       open(unit=80,file="runparams",status="old",form="unformatted", &
            action="write")
       read(unit=80) xlength,ylength,zlength
    endif

#ifdef MPI
    ! Distribute the input parameters to the other nodes
    call MPI_BCAST(xlength,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,&
         ierror)
    call MPI_BCAST(ylength,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,&
         ierror)
    call MPI_BCAST(zlength,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,&
         ierror)
#endif    

    ! Setup the coordinate grid by allocating the coordinate arrays
    call init_grid ()

    ! Scale the physical grid lengths
    xlength=xlength/scleng
    ylength=ylength/scleng
    zlength=zlength/scleng
    
    ! Cell sizes
    dx=xlength/real(max(1,meshx),dp)
    dy=ylength/real(max(1,meshy),dp)
    dz=zlength/real(max(1,meshz),dp)

    ! Fill arrays
    do i=sx-mbc,ex+mbc            
       x(i)=dx*(real(i-1,dp)+0.5_dp)
    enddo
    do j=sy-mbc,ey+mbc
       y(j)=dy*(real(j-1,dp)+0.5_dp)
    enddo
    do k=sz-mbc,ez+mbc
       z(k)=dz*(real(k-1,dp)+0.5_dp)
    enddo

    ! For higher accuracy define separate volumes for x and y
    ! This is not used in the cylindrical version
    
    do k=sz-mbc,ez+mbc
       do j=sy-mbc,ey+mbc
          do i=sx-mbc,ex+mbc
             volx(i,j,k)=1.0_dp
             voly(i,j,k)=1.0_dp
             volz(i,j,k)=1.0_dp
          enddo
       enddo
    enddo

    ! These are edge factors which are fed to the solver
    ! in [xyz]integrate. If the boundary is a singularity,
    ! it is good to set these to zero, otherwise they
    ! should be one.
    ! edge(1) = inner boundary
    ! edge(2) = outer boundary
    xedge(1:2)=1.0
    yedge(1:2)=1.0
    zedge(1:2)=1.0

  end subroutine init_coords

end module grid
