module mesh

  ! Module for Capreole (3D)
  ! Author: Garrelt Mellema
  ! Date: 2005-04-20 (prev 2003-06-01)
  ! This module is also accepted by the F compiler (Dec 9, 2003)
  !
  ! This module contains the routines to set up a computational mesh,
  ! which can be distributed over several processors using
  ! MPI routines.
  !
  ! 3D version

  use precision, only: dp
  use my_mpi
  use file_admin, only: stdinput

  implicit none

  private

  integer,public :: meshx,meshy,meshz ! total grid size in x, y, z
  integer,public :: sx,ex,sy,ey,sz,ez ! local grid start and end coordinates

  public :: init_mesh

contains

  !----------------------------------------------------------------------------

  subroutine init_mesh (restart,restartfile)
    ! Read in the grid dimensions, and distribute them over the
    ! processors

    logical,intent(in) :: restart
    character(len=19),intent(in) :: restartfile

    integer :: ierror=0

    if (.not.restart) then ! Fresh start
       
       if (rank == 0) then
          print "(2/,A,/)", "----- Grid -----"
          write(unit=*,fmt="(a)",advance="no") "1) Number of grid points: "
          read (unit=stdinput,fmt=*) meshx,meshy,meshz
          
          ! Report
          write(unit=30,fmt="(2/,A,/)") "----- Grid -----"
          write(unit=30,fmt="(A,3I5)") "1) Number of grid points: ", &
               meshx,meshy,meshz
       endif

    else
       ! Ask for the input if you are processor 0.
       if (rank == 0) then
          call restart_mesh(restartfile,meshx,meshy,meshz,ierror)
       endif
    endif

#ifdef MPI
    ! Distribute the total grid size over all processors
    call MPI_BCAST(meshx,1,MPI_INTEGER,0,MPI_COMM_NEW,ierror)
    call MPI_BCAST(meshy,1,MPI_INTEGER,0,MPI_COMM_NEW,ierror)
    call MPI_BCAST(meshz,1,MPI_INTEGER,0,MPI_COMM_NEW,ierror)
#endif

    ! Compute the decomposition of the grid over the processors
    ! (find 3D decomposition)
    ! sx = start x coordinate, ex = end x coordinate
    ! sy = start y coordinate, ey = end y coordinate
    ! sz = start z coordinate, ez = end z coordinate
    call fnd3ddecomp ()

    ! Report the grid for the local processor
    write(unit=30,fmt=*) "Grid: ",sx,ex,sy,ey,sz,ez
    
  end subroutine init_mesh

  !========================================================================
  subroutine restart_mesh(filename,xmesh,ymesh,zmesh,ierror)
    
    ! This routine retrieved the mesh size
    ! (xmesh,ymesh,zmesh) from the ah3 file filename.
    ! Should be called from module mesh

    use sizes, only: nrOfDim, neq
    use atomic, only: gamma

    character(len=19),intent(in) :: filename ! name of ah3 file
    integer,intent(out)    :: xmesh,ymesh,zmesh ! 3D size of mesh
    integer,intent(out) :: ierror

    ! AH3D header variables
    character(len=80) :: banner
    integer :: nrOfDim_in ! corresponds to parameter nrOfDim (no. of dimensions)
    integer :: neq_in     ! corresponds to parameter neq (no. of equations)
    integer :: npr_in     ! corresponds to parameter npr (no. of processors)
    integer :: refinementFactor ! not used
    integer :: nframe           ! output counter
    real(kind=dp) :: gamma_in  ! corresponds to parameter gamma (adiab. index)
    real(kind=dp) :: time      ! output time

    ! AH3D grid variables
    integer :: igrid ! counters
    real(kind=dp) :: x_corner,y_corner,z_corner
    real(kind=dp) :: dx_in,dy_in,dz_in
    integer :: level

    ierror=0
    
    ! Read in header
    if (rank.eq.0) then
       open(unit=40,file=filename,form='UNFORMATTED',status='old')
       read(40) banner
       read(40) nrOfDim_in
       read(40) neq_in
       read(40) npr_in
       read(40) refinementFactor
       read(40) nframe
       read(40) gamma_in
       read(40) time
       
       ! Check for consistency
       if (nrOfDim_in /= nrOfDim .or. neq_in /= neq .or. npr_in /= npr .or. &
            gamma_in /= gamma ) then
          ierror=1
          write(*,*) "Error: ah3 file inconsistent with program parameters"
       endif
       
       if (ierror == 0) then
          ! Read in grids
          ! (each processor has its grid, we read in all and find the
          !  largest value to obtain the physical size of the full grid).
          read(40) xmesh,ymesh,zmesh
       endif

       close(40)

    endif
    
  end subroutine restart_mesh

  !----------------------------------------------------------------------------

  subroutine fnd3ddecomp ()
    
    ! This routine makes the decomposition of the 3D grid into 
    ! local processor grids

    call MPE_DECOMP1D( meshx, dims(1), grid_struct(1), sx, ex )
    call MPE_DECOMP1D( meshy, dims(2), grid_struct(2), sy, ey )
    call MPE_DECOMP1D( meshz, dims(3), grid_struct(3), sz, ez )
    
  end subroutine fnd3ddecomp
  
  !----------------------------------------------------------------------------

  subroutine MPE_DECOMP1D( n, numprocs, myid, s, e )

    ! This routine distributes n over the number of processors
    
    ! Input:
    ! n        - number of grid cells
    ! numprocs - number of processors (to distribute over)
    ! myid     - id of current processor
    
    ! Output:
    ! s - start index for current processor
    ! e - end index for current processor

    integer,intent(in)  :: n, numprocs, myid
    integer,intent(out) :: s, e
    integer ::  nlocal
    integer ::  deficit
    
    nlocal  = n / numprocs
    s = myid * nlocal + 1
    deficit = n-int(n/numprocs)*numprocs !mod(n,numprocs)
    s = s + min(myid,deficit)
    if (myid < deficit) then
       nlocal = nlocal + 1
    endif
    e = s + nlocal - 1
    if (e > n .or. myid == numprocs-1) e = n
    
  end subroutine MPE_DECOMP1D
  
end module mesh

