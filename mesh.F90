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
  use file_admin, only: stdinput,log_unit,ah3

  implicit none

  private

  integer,public :: meshx,meshy,meshz ! total grid size in x, y, z
  integer,public :: sx,ex,sy,ey,sz,ez ! local grid start and end coordinates

  public :: init_mesh

#ifdef MPI
  integer :: mpi_ierror
#endif

contains

  !----------------------------------------------------------------------------

  subroutine init_mesh (restart,restartfile)
    ! Read in the grid dimensions, and distribute them over the
    ! processors

    logical,intent(in) :: restart
    character(len=*),intent(in) :: restartfile

    integer :: ierror=0

    if (.not.restart) then ! Fresh start
       
       if (rank == 0) then
          print "(2/,A,/)", "----- Grid -----"
          write(unit=*,fmt="(a)",advance="no") "1) Number of grid points: "
          read (unit=stdinput,fmt=*) meshx,meshy,meshz
          
          ! Report
          write(unit=log_unit,fmt="(2/,A,/)") "----- Grid -----"
          write(unit=log_unit,fmt="(A,3I5)") "1) Number of grid points: ", &
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
    call MPI_BCAST(meshx,1,MPI_INTEGER,0,MPI_COMM_NEW,mpi_ierror)
    call MPI_BCAST(meshy,1,MPI_INTEGER,0,MPI_COMM_NEW,mpi_ierror)
    call MPI_BCAST(meshz,1,MPI_INTEGER,0,MPI_COMM_NEW,mpi_ierror)
#endif

    ! Compute the decomposition of the grid over the processors
    ! (find 3D decomposition)
    ! sx = start x coordinate, ex = end x coordinate
    ! sy = start y coordinate, ey = end y coordinate
    ! sz = start z coordinate, ez = end z coordinate
    call fnd3ddecomp ()

    ! Report the grid for the local processor
    write(unit=log_unit,fmt=*) "Grid: ",sx,ex,sy,ey,sz,ez
    
  end subroutine init_mesh

  !========================================================================

  subroutine restart_mesh(filename,xmesh,ymesh,zmesh,ierror)
    
    ! This routine retrieved the mesh size
    ! (xmesh,ymesh,zmesh) from the ah3 file filename.
    ! Should be called from module mesh

    use sizes, only: nrOfDim, neq
    use atomic, only: gamma

    character(len=*),intent(in) :: filename ! name of ah3 file
    integer,intent(out) :: xmesh,ymesh,zmesh ! 3D size of mesh
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

    ierror=0
    
    ! Read in header
    if (rank == 0) then
       open(unit=ah3,file=filename,form="unformatted",status="old", &
            action="read")
       read(unit=ah3) banner
       read(unit=ah3) nrOfDim_in
       read(unit=ah3) neq_in
       read(unit=ah3) npr_in
       read(unit=ah3) refinementFactor
       read(unit=ah3) nframe
       read(unit=ah3) gamma_in
       read(unit=ah3) time
       
       ! Check for consistency
       if (nrOfDim_in /= nrOfDim .or. neq_in /= neq .or. npr_in /= npr .or. &
            gamma_in /= gamma ) then
          ierror=1
          write(unit=log_unit,fmt=*) &
               "Error: ah3 file inconsistent with program parameters"
       endif
       
       if (ierror == 0) then
          ! Read in grids
          ! (each processor has its grid, we read in all and find the
          !  largest value to obtain the physical size of the full grid).
          read(unit=ah3) xmesh,ymesh,zmesh
       endif

       close(unit=ah3)

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

  subroutine MPE_DECOMP1D (nmesh, numprocs, myid, startpnt, endpnt)

    ! This routine distributes n over the number of processors
    
    ! Input:
    ! nmesh    - number of mesh cells
    ! numprocs - number of processors (to distribute over)
    ! myid     - id of current processor
    
    ! Output:
    ! startpnt - start index for current processor
    ! endpnt   - end index for current processor

    integer,intent(in)  :: nmesh, numprocs, myid
    integer,intent(out) :: startpnt, endpnt
    integer ::  nlocal
    integer ::  deficit
    
    nlocal  = nmesh / numprocs
    startpnt = myid * nlocal + 1
    deficit = nmesh-int(nmesh/numprocs)*numprocs !mod(nmesh,numprocs)
    startpnt = startpnt + min(myid,deficit)
    if (myid < deficit) then
       nlocal = nlocal + 1
    endif
    endpnt = startpnt + nlocal - 1
    if (endpnt > nmesh .or. myid == numprocs-1) then 
       endpnt = nmesh
    endif

  end subroutine MPE_DECOMP1D
  
end module mesh

