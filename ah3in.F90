module input
	
  ! Module for Capreole (f90)
  ! Author: Garrelt Mellema
  ! Date: 2005-12-02

  ! This module contains routines dealing with input from
  ! AH3 files.
  ! - restart_grid: obtain full grid size from ah3 file
  ! - restart_vars: obtain state variables from ah3 file

  ! Version: 3D / MPI not fully enabled.

  !--------------------------------------
  ! AH3D file format (24-10-2002)
  
  ! File name: 24102002_a_0000.ah3

    ! Header:
    ! string 80 bytes
    ! int    nrOfDim
    ! int    nrOfVars
    ! int    nrOfGrids
    ! int    refinementFactor
    ! int    frameNr
    ! double gamma
    ! double time

    ! Grids:
    ! int    nrOfCells1 [, nrOfCells2, nrOfCells3]
    ! double corner1 [, corner2, corner3]
    ! double cellSize1 [,cellSize2, cellSize3]
    ! int    level

    ! Cells;
    ! double rho
    ! double rho*v1 [, rho*v2, rho*v3]
    ! double rho*e
    ! [double var(nrOfVars-(2+nrOfDim))]
  !--------------------------------------

  use precision, only: dp
  use sizes
  use my_mpi
  use atomic

  implicit none
  
contains

  !========================================================================
  subroutine restart_mesh(filename,xmesh,ymesh,zmesh,ierror)
    
    ! This routine retrieved the mesh size
    ! (xmesh,ymesh,zmesh) from the ah3 file filename.
    ! Should be called from module mesh

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

  !========================================================================
  subroutine restart_grid(filename,xgrid,ygrid,zgrid,ierror)
    
    ! This routine constructs the physical size of the grid
    ! (xgrid,ygrid,zgrid) from the ah3 file filename.
    ! Should be called from module coords

    character(len=19),intent(in) :: filename ! name of ah3 file
    real(kind=dp),intent(out)    :: xgrid,ygrid,zgrid ! 3D size of grid
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
    integer :: xmesh,ymesh,zmesh 
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
          xgrid=0.0
          ygrid=0.0
          zgrid=0.0
          do igrid=1,npr_in
             read(40) xmesh,ymesh,zmesh
             read(40) x_corner,y_corner,z_corner
             read(40) dx_in,dy_in,dz_in
             read(40) level
             xgrid=max(xgrid,x_corner+dx_in*(real(xmesh)-0.5))
             ygrid=max(ygrid,y_corner+dy_in*(real(ymesh)-0.5))
             zgrid=max(zgrid,z_corner+dz_in*(real(zmesh)-0.5))
          enddo
       endif

       close(40)

    endif
    
  end subroutine restart_grid

  !========================================================================
  subroutine restart_state(filename,ierror)
    
    ! This routine reads in the state variables and the time
    ! from the ah3 file filename.
    ! It should be called from the problem module.
    ! It handles multi-processor input by reading it from rank 0
    ! and sending it on the appropriate processor.
    
    character(len=19),intent(in) :: filename ! name of output file
    real(kind=dp),dimension(sx-mbc:ex+mbc,sy-mbc:ey+mbc,sz-mbc:ez+mbc,neq), & 
         intent(out) :: state
    real(kind=dp),intent(out) :: time
    integer,intent(out) :: ierror

    ! AH3D header variables
    character(len=80) :: banner
    integer :: nrOfDim_in,neq_in,npr_in
    integer :: refinementFactor
    integer :: nframe
    real(kind=dp) :: gamma_in

    ! AH3D grid variables
    integer :: igrid ! counters
    integer :: xmesh,ymesh,zmesh
    real(kind=dp) :: x_corner,y_corner,z_corner
    real(kind=dp) :: dx_in,dy_in,dz_in
    integer :: level

    ierror=0
    
    ! Header
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
    endif

    ! Check for consistency
    if (nrOfDim_in /= nrOfDim .or. neq_in /= neq .or. npr_in /= npr .or. &
         gamma_in /= gamma ) then
       ierror=1
       write(*,*) "Error: ah3 file inconsistent with program parameters"
    endif

    if (ierror == 0 ) then
#ifdef MPI
       ! Distribute runid over nodes
       call MPI_BCAST(time,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,ierror)
#endif     

       ! Grid: read in, but ignore. Grid should be handled by
       ! restart_grid and module coords
       if (rank == 0) then
          do igrid=1,npr_in
             read(40) xmesh,ymesh,zmesh
             read(40) x_corner,y_corner,z_corner
             read(40) dx_in,dy_in,dz_in
             read(40) level
          enddo
       endif
       
       if (rank == 0) then
          ! read in state for processor 0
          ! (this can always be done)
          read(40) state(sx:ex,sy:ey,sz:ez,RHO)
          read(40) state(sx:ex,sy:ey,sz:ez,RHVX)
          read(40) state(sx:ex,sy:ey,sz:ez,RHVY)
          read(40) state(sx:ex,sy:ey,sz:ez,RHVZ)
          read(40) state(sx:ex,sy:ey,sz:ez,EN)
          if (neq > neuler) then
             read(40) state(sx:ex,sy:ey,sz:ez,neuler+1:neq)
          endif
       endif
       
#ifdef MPI
       ! Allocate temporary arrays for reading the state variables
       allocate(temparray(ex-sx+1,ey-sy+1,ez-sz+1)
       allocate(temparray2(ex-sx+1,ey-sy+1,ez-sz+1,neuler+1:neq)
       if (rank == 0) then
          ! read in for other processors
          do igrid=2,npr_in
             do ieq=RHO,EN
                ! read in state variables one by one and send them to
                ! appropriate processor
                read(40) temparray
                call MPI_ISSEND(temparray,(ex-sx+1)*(ey-sy+1)*(ez-sz+1), &
                     MPI_DOUBLE_PRECISION,igrid-1,ieq,MPI_COMM_NEW, &
                     request,ierror)
             enddo
          enddo
          ! The non-euler variables are written in one big array, read it
          ! in and send it on
          read(40) temparray2
          call MPI_ISSEND(temparray,(ex-sx+1)*(ey-sy+1)*(ez-sz+1)*(neq-neuler), &
               MPI_DOUBLE_PRECISION,igrid-1,neuler+1,MPI_COMM_NEW, &
               request,ierror)
       else
          ! if you are a non rank 0 processor, wait for the arrays to arrive
          ! and attach it to the appropriate state variable
          do ieq=RHO,EN
             call MPI_RECV(temparray,(ex-sx+1)*(ey-sy+1)*(ez-sz+1), &
                  MPI_DOUBLE_PRECISION,0,ieq,MPI_COMM_NEW, &
                  status,ierror)
             state(sx:ex,sy:ey,sz:ez,ieq)=temparray
          enddo
          call MPI_RECV(temparray2,(ex-sx+1)*(ey-sy+1)*(ez-sz+1)*(neq-neuler), &
               MPI_DOUBLE_PRECISION,0,neuler+1,MPI_COMM_NEW, &
               status,ierror)
          state(sx:ex,sy:ey,sz:ez,neuler+1:neq)=temparray2
       endif
       deallocate(temparray)
       deallocate(temparray2)
#endif
    endif
    
  end subroutine restart_state
  
  !========================================================================
  subroutine restart_time(filename,simtime,ierror)
    
    ! This routine retrieves the time (simtime)
    ! from the ah3 file filename.
    ! Should be called from module time

    character(len=19),intent(in) :: filename ! name of ah3 file
    integer,intent(out)    :: simtime ! time of restart
    integer,intent(out) :: ierror

    ! AH3D header variables
    character(len=80) :: banner
    integer :: nrOfDim_in ! corresponds to parameter nrOfDim (no. of dimensions)
    integer :: neq_in     ! corresponds to parameter neq (no. of equations)
    integer :: npr_in     ! corresponds to parameter npr (no. of processors)
    integer :: refinementFactor ! not used
    integer :: nframe           ! output counter
    real(kind=dp) :: gamma_in  ! corresponds to parameter gamma (adiab. index)

    ierror=0
    
    ! Read in header
    if (rank == 0) then
       open(unit=40,file=filename,form='UNFORMATTED',status='old')

       read(40) banner
       read(40) nrOfDim_in
       read(40) neq_in
       read(40) npr_in
       read(40) refinementFactor
       read(40) nframe
       read(40) gamma_in
       read(40) simtime
       
       close(40)
    endif
    
  end subroutine restart_time

end module input
   
