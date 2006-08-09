module integrator

  ! Module for Capreole (3D)
  ! Author: Garrelt Mellema
  ! Date: 2003-08-14
 
  ! Integrates the set of equations over one time step (dt).

  ! Version: Strang splitting
 
  use precision
  use sizes
  use mesh
  use hydro
  use times
  use problem
  use protection
  use oddeven
  use boundary
  use hydrosolver

  implicit none
  
  integer,parameter,private :: transverse=1
  integer,parameter,private :: notransverse=0

contains
  
  subroutine integrate(nstep,istop)

    ! This routine integrates the grid for one time step

    integer,intent(in)  :: nstep  ! time step counter
    integer,intent(out) :: istop  ! control integer

    ! Local variables
    integer :: istop1,istop2
    integer :: inegative
    integer :: ierror,num,mpierror

    integer,dimension(nrOfDim) :: instep
    
    ! Set flags to zero
    inegative=0
    ierror=0
    mpierror=0
    istop1=0
    istop2=0
    istop=0

    ! reset this to store flag for oe-fix
    !stold(:,:,:,neq)=0.0d0
    !stnew(:,:,:,neq)=0.0d0

    ! Take one time step Strang splitting
    ! Alternate the order between time steps
    instep=ran123()

    do num=1,nrofDim
       istoptest: if (istop == 0) then
          
          ! Point generic volume array to volx
          vol => volx

          select case (instep(num))
          case (1)
             if (meshx > 1) call xintegr(ierror)
          case (2)
             if (meshy > 1) call yintegr(ierror)
          case (3)
             if (meshz > 1) call zintegr(ierror)
          end select

          ! check for hydro errors (ierror <> 0)
          ! (find the maximum of ierror and store this in istop1)
          istop1=ierror
#ifdef MPI	
          call MPI_ALLREDUCE(ierror,istop1,1,MPI_INTEGER,MPI_MAX, &
               MPI_COMM_NEW,mpierror)
#endif

          ! Point generic state array to stnew
          state => stnew
          
          ! set the inflow condition
          ! the inflow routine has to be supplied
          call inflow(NEW)

          ! exchange boundaries with neighbours
          ! This routine also calculates the new pressure
          call exchngxy(NEW)

          ! Odd-even fix
          !call odd_even(-1,NEW)

          ! exchange boundaries with neighbours
          ! This routine also calculates the new pressure
          !call exchngxy(NEW)

          ! Protect against negative pressures
          call presprot(inegative,2,NEW)
          istop2=inegative
#ifdef MPI	
          call MPI_ALLREDUCE(inegative,istop2,1,MPI_INTEGER,MPI_MAX, &
               MPI_COMM_NEW,mpierror)
#endif	
          istop=istop1+istop2

          ! exchange boundaries with neighbours
          ! This routine also calculates the new pressure
          call exchngxy(NEW)

          ! Copy new state to old state
          stold=stnew
       endif istoptest
    enddo

    if (istop == 0) then ! otherwise serious error occurred
       call rad_evolve3D( )
       ! exchange boundaries with neighbours
       ! This routine also calculates the new pressure
       call exchngxy(NEW)
    endif

  end subroutine integrate
  
  !=======================================================================
     
  subroutine xintegr(itoterror)
    
    ! This routine integrates the grid in the x direction
    
    integer,intent(out) :: itoterror
    
    ! Local variables
    integer ::  i,j,k,ieq,ij,ik,mesh,ioff,ierror
    !real(kind=dp),dimension(2) :: edge

    ! Initialize error status to zero
    itoterror=0

    mesh=ex-sx+1     ! the size of the y row
    ioff=sx-1        ! offset between grid and y rows
    call constr_solver(mesh)

    ! remap to 1D state variable running from 1-mbc to ex-sx+1+mbc
    kloop: do k=sz,ez
       jloop :do j=sy,ey  ! integrate over all rows 
          
          ij=j             ! to export the current j position
          ik=k             ! to export the current k position
          
          ! Calculate volume weighted state and pressure
          do i=1-mbc,mesh+mbc
             state1d(i,1:neq)=stold(i+ioff,j,k,1:neq)
             wp(i)=pressr(i+ioff,j,k)
          enddo
          
          ! Call Roe solver
          call solver(mesh,dt,dx,dy,dz,2,3,4,ij,ik,ierror)
          
          ! Report problems
          if (ierror /= 0) then
             write(30,*) 'Roe solver error (x), j,k= ',j,k
             call flush(30)
             itoterror=itoterror+ierror
          endif

          ! Update the stnew variable
          do ieq=1,neq
             do i=sx,ex
                stnew(i,j,k,ieq)=stold(i,j,k,ieq)+dstate(i-ioff,ieq)
             enddo
          enddo

       enddo jloop
    enddo kloop
    
    ! Destroy solver variables
    call destr_solver()

  end subroutine xintegr

  !========================================================================

  subroutine yintegr(itoterror)

    ! This routine integrates the grid in the y direction
    
    integer,intent(out) :: itoterror
    
    integer :: i,j,k,ieq,ij,ik,mesh,joff,ierror

    ! Initialize status to zero
    itoterror=0

    ! remap to 1D state variable running from 1-mbc to ey-sy+1+mbc      
    mesh=ey-sy+1    ! size of the y row
    joff=sy-1       ! offset between grid and row
    call constr_solver(mesh)
    
    kloop: do k=sz,ez
       iloop: do i=sx,ex      ! integrate over all rows 
          
          ij=i            ! to export the current i position
          ik=k            ! to export the current k position
          
          ! Calculate volume weighted state and pressure
          do j=1-mbc,mesh+mbc
             state1d(j,1:neq)=stold(i,j+joff,k,1:neq)
             wp(j)=pressr(i,j+joff,k)
          enddo

          ! Call the Roe solver
          call solver(mesh,dt,dy,dz,dx,3,4,2,ij,ik,ierror)

          ! Report problems
          if (ierror /= 0) then
             write(30,*) 'Roe solver error (y), i,k= ',i,k
             call flush(30)
             itoterror=itoterror+ierror
          endif

          ! Update the stnew variable
          do ieq=1,neq
             do j=sy,ey
                stnew(i,j,k,ieq)=stold(i,j,k,ieq)+dstate(j-joff,ieq)
             enddo
          enddo

       enddo iloop
    enddo kloop
    
    ! Destroy solver variables
    call destr_solver()
    
  end subroutine yintegr

  !========================================================================

  subroutine zintegr(itoterror)

    ! This routine integrates the grid in the z direction
    
    integer,intent(out) :: itoterror
    
    integer :: i,j,k,ieq,ij,ik,mesh,koff,ierror

    ! Initialize status to zero
    itoterror=0

    ! remap to 1D state variable running from 1-mbc to ez-sz+1+mbc      
    mesh=ez-sz+1    ! size of the x row
    koff=sz-1       ! offset between grid and row
    call constr_solver(mesh)

    jloop: do j=sy,ey
       iloop: do i=sx,ex      ! integrate over all rows 
          
          ij=i            ! to export the current i position
          ik=j            ! to export the current j position
          
          ! Calculate volume weighted state and pressure
          do k=1-mbc,mesh+mbc
             state1d(k,1:neq)=stold(i,j,k+koff,1:neq)
             wp(k)=pressr(i,j,k+koff)
          enddo

          ! Call the Roe solver
          call solver(mesh,dt,dz,dx,dy,4,2,3,ij,ik,ierror)

          ! Report problems
          if (ierror /= 0) then
             write(30,*) 'Roe solver error (z), i,j= ',i,j
             call flush(30)
             itoterror=itoterror+ierror
          endif

          ! Update the stnew variable
          do ieq=1,neq
             do k=sz,ez
                stnew(i,j,k,ieq)=stold(i,j,k,ieq)+dstate(k-koff,ieq)
             enddo
          enddo

       enddo iloop
    enddo jloop
    
    ! Destroy solver variables
    call destr_solver()

  end subroutine zintegr

  !========================================================================

  function ran123 ()

    ! Find a random order of the numbers 1, 2, 3

    integer,dimension(3) :: ran123

    real :: ran_num,order

    call random_number(ran_num)
    order=6.0*ran_num
    if (order.lt.2.0) then
       ran123(1)=1
       if (order.lt.1.0) then
          ran123(2)=2
          ran123(3)=3
       else
          ran123(2)=3
          ran123(3)=2
       endif
    elseif (order.lt.4.0) then
       ran123(1)=2
       if (order.lt.3.0) then
          ran123(2)=1
          ran123(3)=3
       else
          ran123(2)=3
          ran123(3)=1
       endif
    else
       ran123(1)=3
       if (order.lt.5.0) then
          ran123(2)=1
          ran123(3)=2
       else
          ran123(2)=2
          ran123(3)=1
       endif
    endif

  end function ran123

end module integrator




