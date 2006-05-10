module output
  
  use precision
  use scaling
  use sizes
  use my_mpi
  use mesh
  use grid
  use atomic
  use hydro
  use times

  !implicit none

  private

  character(len=1),private :: runid
  character(len=8),private :: revdate

  public :: init_output, make_output

contains
  
  !========================================================================
  subroutine init_output(restart,ierror)
    
    ! This routine takes care of the output
    
    logical,intent(in) :: restart
    integer,intent(out) :: ierror

    integer  :: i,j,ieq ! counters and error flags
    character(len=19) :: filename ! name of output file
    character(len=6) :: string_nframe,string_rank ! string version of integers
    character(len=8) :: date
    character(len=4) :: string_frame
    logical :: test_exist

    ! Construct date-identifier: d-m-y
    call date_and_time(date)
    revdate(1:2)=date(7:8)
    revdate(3:4)=date(5:6)
    revdate(5:8)=date(1:4)
    ! Run ID
    if (rank.eq.0) then
       write(*,'(A,$)') 'Run ID (one letter): '
       read(*,*) runid
       ! Inquire if it exists
       if (.not.restart) then
          do
             filename=revdate//'_'//runid//'_0000.ah3'
             inquire(file=filename,exist=test_exist)
             if (.not.test_exist) exit
             if (runid < 'z') then
                runid=achar(iachar(runid)+1)
                write(30,*) 'RunID already in use, trying ',runid
             else
                write(30,*) 'All runids are in use?'
                ierror=1
             endif
          enddo
       endif
    endif

#ifdef MPI
    ! Distribute runid over nodes
    call MPI_BCAST(runid,1,MPI_CHARACTER,0,MPI_COMM_NEW,ierror)
#endif     

  end subroutine init_output

  !========================================================================
  subroutine make_output(nframe)
    
    ! This routine takes care of the output
    
    integer,intent(in) :: nframe
    
    ! AH3D Output variables
    integer,parameter :: refinementFactor=1
    integer,parameter :: level=1
    character(len=80),parameter :: banner='Capreole (F90) 3D Hydrodynamics'

    integer   :: ierror,i,j,k,ieq ! counters and error flags
    character(len=19) :: filename ! name of output file
    character(len=6)  :: string_nframe,string_rank ! string version of integers
    character(len=4)  :: string_frame
    integer :: space

    integer,parameter :: outputcircle=601
    integer :: request,nextproc
#ifdef MPI
    integer :: status(MPI_STATUS_SIZE)
#endif

    ierror=0
    write(30,*) 'Writing Frame: ',nframe
      
    ! Construct file name and open the file
    write(string_frame,'(i4)') nframe
    do
       space=index(string_frame,' ')
       if (space == 0) exit
       string_frame(space:space)='0'
    enddo
    filename=revdate//'_'//runid//'_'//string_frame//'.ah3'
    
    ! AH3D output format (24-10-2002)

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
    ! double rho*v1 [, rho*v1, rho*v3]
    ! double rho*e
    ! [double var(nrOfVars-(2+nrOfDim))]

    if (rank.eq.0) then
       ! Header
       open(unit=40,file=filename,form='UNFORMATTED')
       write(40) banner
       write(40) nrOfDim
       write(40) neq
       write(40) npr
       write(40) refinementFactor
       write(40) nframe
       write(40) gamma
       write(40) time*sctime
       close(40)
    endif

    if (rank.eq.0) then

       ! Grid 
       !!open(unit=40,file=filename,form='UNFORMATTED',status='OLD', &
       !!        position='APPEND')
       open(unit=40,file=filename,form='UNFORMATTED',status='OLD', &
               position='APPEND')
       write(40) ex-sx+1,ey-sy+1,ez-sz+1
       write(40) x(sx)*scleng,y(sy)*scleng,z(sz)*scleng
       write(40) dx*scleng,dy*scleng,dz*scleng
       write(40) level
       close(40)
       
#ifdef MPI
       ! Send filename to next processor
       if (npr > 1) then
          call MPI_ISSEND(filename,19,MPI_CHARACTER,rank+1,outputcircle, &
               MPI_COMM_NEW,request,ierror)
          write(30,*) rank,'sent filename to ',rank+1
          ! Wait for the circle to complete
          call MPI_RECV(filename,19,MPI_CHARACTER,npr-1,outputcircle, &
               MPI_COMM_NEW,status,ierror)
          write(30,*) rank,'received filename from ',npr-1
       endif
#endif

       ! Cells
       !!open(unit=40,file=filename,form='UNFORMATTED',status='OLD', &
       !!        position='APPEND')
       open(unit=40,file=filename,form='UNFORMATTED',status='OLD', &
            position='APPEND')
       write(40) (((state(i,j,k,RHO)*scdens,i=sx,ex),j=sy,ey),k=sz,ez)
       write(40) (((state(i,j,k,RHVX)*scmome,i=sx,ex),j=sy,ey),k=sz,ez)
       write(40) (((state(i,j,k,RHVY)*scmome,i=sx,ex),j=sy,ey),k=sz,ez)
       write(40) (((state(i,j,k,RHVZ)*scmome,i=sx,ex),j=sy,ey),k=sz,ez)
       write(40) (((state(i,j,k,EN)*scener,i=sx,ex),j=sy,ey),k=sz,ez)
       if (neq > neuler) then
          write(40) ((((state(i,j,k,ieq),i=sx,ex),j=sy,ey),k=sz,ez), &
               ieq=neuler+1,neq)
       endif
       close(40)
       ! Send filename to next processor
#ifdef MPI
       if (npr > 1) then
          call MPI_ISSEND(filename,19,MPI_CHARACTER,rank+1,outputcircle, &
               MPI_COMM_NEW,request,ierror)
          write(30,*) rank,'sent filename to ',rank+1
          ! Wait for the circle to complete
          call MPI_RECV(filename,19,MPI_CHARACTER,npr-1,outputcircle, &
               MPI_COMM_NEW,status,ierror)
          write(30,*) rank,'received filename from ',npr-1,' end of loop.'
       endif
#endif
#ifdef MPI
    else
       
       ! Grids
       ! Receive filename from previous processor
       call MPI_RECV(filename,19,MPI_CHARACTER,rank-1,outputcircle, &
            MPI_COMM_NEW,status,ierror)
       write(30,*) rank,'received filename from ',rank-1
       
       if (ierror == 0) then ! if ok
          open(unit=40,file=filename,form='UNFORMATTED',status='OLD', &
               position='APPEND')
          ! Grid 
          write(40) ex-sx+1,ey-sy+1,ez-sz+1
          write(40) x(sx)*scleng,y(sy)*scleng,z(sz)*scleng
          write(40) dx*scleng,dy*scleng,dz*scleng
          write(40) level
          close(40)
          
          ! Determine next processor (npr -> 0)
          nextproc=mod(rank+1,npr)
          ! Send filename along
          call MPI_ISSEND(filename,19,MPI_CHARACTER,nextproc, &
               outputcircle,MPI_COMM_NEW,request,ierror)
          write(30,*) rank,'sent filename to ',nextproc
       endif

       ! Cells
       ! Receive filename from previous processor
       call MPI_RECV(filename,19,MPI_CHARACTER,rank-1,outputcircle, &
            MPI_COMM_NEW,status,ierror)
       write(30,*) rank,'received filename from ',rank-1
          
       if (ierror == 0) then ! if ok
          open(unit=40,file=filename,form='UNFORMATTED',status='OLD', &
               position='APPEND')
          ! Cells
          write(40) (((state(i,j,k,RHO)*scdens,i=sx,ex),j=sy,ey),k=sz,ez)
          write(40) (((state(i,j,k,RHVX)*scmome,i=sx,ex),j=sy,ey),k=sz,ez)
          write(40) (((state(i,j,k,RHVY)*scmome,i=sx,ex),j=sy,ey),k=sz,ez)
          write(40) (((state(i,j,k,RHVZ)*scmome,i=sx,ex),j=sy,ey),k=sz,ez)
          write(40) (((state(i,j,k,EN)*scener,i=sx,ex),j=sy,ey),k=sz,ez)
          if (neq > neuler) then
             write(40) ((((state(i,j,k,ieq),i=sx,ex),j=sy,ey),k=sz,ez), &
                  ieq=neuler+1,neq)
          endif
          close(40)
          
          ! Determine next processor (npr -> 0)
          nextproc=mod(rank+1,npr)
          ! Send filename along
          call MPI_ISSEND(filename,19,MPI_CHARACTER,nextproc, &
               outputcircle,MPI_COMM_NEW,request,ierror)
          write(30,*) rank,'sent filename to ',nextproc
       else
          write(30,*) 'MPI error in output routine'
          ierror=ierror+1
       endif
#endif
    endif

  end subroutine make_output
  
end module output
