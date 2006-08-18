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

  use my_mpi
  use file_admin, only: stdinput
  private

  integer,public :: meshx,meshy,meshz ! total grid size in x, y, z
  integer,public :: sx,ex,sy,ey,sz,ez ! local grid start and end coordinates

  public :: init_mesh

contains

  !----------------------------------------------------------------------------

  subroutine init_mesh ()
    ! Read in the grid dimensions, and distribute them over the
    ! processors

    integer :: ierror=0

    if (rank == 0) then
       print "(2/,A,/)", "----- Grid -----"
       write(unit=*,fmt="(a)",advance="no") "1) Number of grid points: "
       read (unit=stdinput,fmt=*) meshx,meshy,meshz

    ! Report
       write(unit=30,fmt="(2/,A,/)") "----- Grid -----"
       write(unit=30,fmt="(A,3I5)") "1) Number of grid points: ", &
            meshx,meshy,meshz
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

