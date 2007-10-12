module integrator

  ! Module for Capreole (3D)
  ! Author: Garrelt Mellema
  ! Date: 2007-10-09 (2D: 2007-02-28 (2006-02-11, 2003-08-14)
 
  ! Integrates the set of equations over one time step (dt).

  ! Version: van Leer Flux Vector Splitting
 
  use precision, only: dp
  use sizes, only: neuler, neq, mbc, RHO, RHVX, RHVY, RHVZ, EN
  use mesh, only: sx,ex,sy,ey,sz,ez
  use grid, only: dx,dy,dz
  use atomic, only: gamma,gamma1
  use hydro, only: state,stnew,stold,NEW,OLD,pressr
  use times, only: dt
  use problem, only: inflow
  use protection, only: presprot
  use boundary, only: exchngxy
  use ionic, only: rad_evolve3d

  implicit none
  
contains
  
  subroutine integrate(nstep,istop)
    
    !! This routine integrates the hydrodynamic quantities (state)
    !! in three dimensions, using the van Leer Flux Vector Splitting 
    !! (VLFVS) method.
    
    !! Author: Garrelt Mellema (mellema@strw.leidenuniv.nl)
    !! Date:   2007-10-09 (2007-02-28 (14 August 2001)
    
    integer,intent(in)  :: nstep  ! time step counter
    integer,intent(out) :: istop  ! control integer
    
    real(kind=dp) :: dt2 ! half a time time step
    
    real(kind=dp) :: st1(1-mbc:1+mbc,1-mbc:1+mbc,1-mbc:1+mbc,neq) !! stencil for state
    real(kind=dp) :: p1(1-mbc:1+mbc,1-mbc:1+mbc,1-mbc:1+mbc)      !! stencil for pressure
    real(kind=dp) :: source(neq)
    real(kind=dp) :: df(neq),dg(neq),dh(neq)  !! state changes 
    integer :: inegative,istop1,istop2
    integer :: i,j,k,ii,jj,kk,ieq,ij
    
    ! Set flags to zero
    inegative=0

    istoptest: if (istop == 0) then
       
       ! To achieve second order in space and time, 
       ! do half a time step first.
       dt2=0.5d0*dt
       
       do k=sz,ez
          do j=sy,ey
             do i=sx,ex
                ! make (2*mbc+1)^3 stencil for state and pressure
                do kk=1-mbc,1+mbc
                   do jj=1-mbc,1+mbc
                      do ii=1-mbc,1+mbc
                         do ieq=1,neuler
                            st1(ii,jj,kk,ieq)=stold(i+ii-1,j+jj-1,k+kk-1,ieq)
                         enddo
                         do ieq=neuler+1,neq
                            st1(ii,jj,kk,ieq)=stold(i+ii-1,j+jj-1,k+kk-1,ieq)* &
                                 stold(i+ii-1,j+jj-1,k+kk-1,RHO)
                         enddo
                         p1(ii,jj,kk)=pressr(i+ii-1,j+jj-1,k+kk-1)
                      enddo
                   enddo
                enddo
                ! Calculate the VLFVS fluxes (1st order)
                call fluxes(st1,p1,df,dg,dh)
                source=0.0
                ! Update the state
                do ieq=1,neuler
                   stnew(i,j,k,ieq)=stold(i,j,k,ieq)- &
                        dt2*(df(ieq)/dx+dg(ieq)/dy+dh(ieq)/dz-source(ieq))
                enddo
                do ieq=neuler+1,neq
                   stnew(i,j,k,ieq)=(st1(1,1,1,ieq)- &
                        dt2*(df(ieq)/dx+dg(ieq)/dy+dh(ieq)/dz-source(ieq))) &
                        /stnew(i,j,k,RHO)
                enddo
                
             enddo
          enddo
       enddo

       ! exchange boundaries with neighbours
       ! This routine also calculates the new pressure
       call exchngxy(NEW)
       ! Protect against negative pressures
       call presprot(inegative,2,NEW)
       
       ! check for hydro errors (inegative <> 0)
       ! (find the maximum of inegative and store this in istop)
       istop2=inegative
#ifdef MPI	
       call MPI_ALLREDUCE(inegative,istop2,1,MPI_INTEGER,MPI_MAX, &
            MPI_COMM_NEW,mpierror)
#endif	
       istop=istop1+istop2
    endif istoptest
    
    if (istop.eq.0) then
       
       ! Integrate a full time step
       do k=sz,ez
          do j=sy,ey
             do i=sx,ex
                ! make (2*mbc+1)^3 stencil for state and pressure
                do kk=1-mbc,1+mbc
                   do jj=1-mbc,1+mbc
                      do ii=1-mbc,1+mbc
                         do ieq=1,neuler
                            st1(ii,jj,kk,ieq)=stnew(i+ii-1,j+jj-1,k+kk-1,ieq)
                         enddo
                         do ieq=neuler+1,neq
                            st1(ii,jj,kk,ieq)=stnew(i+ii-1,j+jj-1,k+kk-1,ieq)* &
                                 stnew(i+ii-1,j+jj-1,k+kk-1,RHO)
                         enddo
                         p1(ii,jj,kk)=pressr(i+ii-1,j+jj-1,k+kk-1)
                      enddo
                   enddo
                enddo
                ! Calculate the VLFVS fluxes (1st order)
                call fluxes2(st1,p1,df,dg,dh)
                source=0.0
                ! Update the state
                ! Update stold
                ! Note: do not update stnew here, since we are stepping
                !        through the grid one would get `feed-back' effects
                !        if stnew changes
                ! GM (2007-10-09): this comment seems out of date
                do ieq=1,neuler
                   stnew(i,j,k,ieq)=stold(i,j,k,ieq)- &
                        dt*(df(ieq)/dx+dg(ieq)/dy+dh(ieq)/dz-source(ieq))
                enddo
                do ieq=neuler+1,neq
                   stnew(i,j,k,ieq)=(stold(i,j,k,ieq)*stold(i,j,k,RHO)- &
                        dt*(df(ieq)/dx+dg(ieq)/dy+dh(ieq)/dz-source(ieq))) &
                        /stnew(i,j,k,RHO)
                enddo
                
             enddo
          enddo
       enddo
       
       ! set the inflow condition
       ! the inflow routine has to be supplied
       call inflow(NEW) 
       ! exchange boundaries with neighbours
       ! This routine also calculates the new pressure
       call exchngxy(NEW)
       ! Protect against negative pressures
       call presprot(inegative,2,NEW)
       istop2=inegative
#ifdef MPI	
       call MPI_ALLREDUCE(inegative,istop2,1,MPI_INTEGER,MPI_MAX, &
            MPI_COMM_NEW,mpierror)
#endif	
       istop=istop1+istop2
    endif
    
    if (istop == 0) then ! otherwise serious error occurred
       ! Point generic state array to stnew (the newest at this point)
       state => stnew
       ! Apply radiative processes.
       ! rad_evolve changes state in place (so no changing from stold
       ! to stnew is involved).
       call rad_evolve3D(dt)
       ! exchange boundaries with neighbours
       ! This routine also calculates the new pressure
       call exchngxy(NEW)
    endif
    
  end subroutine integrate

  !=======================================================================
  
  subroutine fluxes(st1,p1,df,dg,dh)
    
    !! This routine calculates the first order fluxes for 
    !! Van Leer's Flux Splitting method.
    
    !! Author: Garrelt Mellema (mellema@strw.leidenuniv.nl)
    !! Date:   14 August 2001
    
    !! Input:
    !! st1 : array containing the state in a (mbc+1)*(mbc+1) 
    !!         stencil around the grid point being updates
    !! p1  : the same for the pressure
    
    !! Output:
    !! df  : state change due to fluxes in the x-direction
    !! dg  : state change due to fluxes in the y-direction
    !! dh  : state change due to fluxes in the z-direction
    
    integer,parameter :: nodir=0
    integer,parameter :: null=0

    real(kind=dp),intent(in) :: st1(1-mbc:1+mbc,1-mbc:1+mbc,1-mbc:1+mbc,neq) !! stencil for state
    real(kind=dp),intent(in) :: p1(1-mbc:1+mbc,1-mbc:1+mbc,1-mbc:1+mbc)      !! stencil for pressure
    real(kind=dp),intent(out) :: df(neq),dg(neq),dh(neq)          !! state changes 

    real(kind=dp) :: fl(neq),fr(neq),fpl(neq),fmr(neq)
    real(kind=dp) :: gl(neq),gr(neq),gpl(neq),gmr(neq)
    real(kind=dp) :: hl(neq),hr(neq),hpl(neq),hmr(neq)
    
    real(kind=dp) :: d,u,v,w,cs,dex(neuler+1:neq)
    real(kind=dp) :: dum(neq) ! dummy
    integer :: ieq
      
    ! calculation of the cell centred fluxes
    
    ! First find the primitives (d,u,v,cs,dex)
    call primitives(st1,p1,d,u,v,w,cs,dex,1,1,1,nodir,null)
    ! Find the left and right going fluxes (fl and fr) in the x-direction
    call VanLeerFluxVectorSplit(d,u,v,w,cs,dex,RHVX,RHVY,RHVZ,fl,fr)
    ! Find the down and up going fluxes (gl and gr) in the y-direction
    call VanLeerFluxVectorSplit(d,v,u,w,cs,dex,RHVY,RHVX,RHVZ,gl,gr)
    ! Find the down and up going fluxes (hl and hr) in the z-direction
    call VanLeerFluxVectorSplit(d,w,u,v,cs,dex,RHVZ,RHVX,RHVY,hl,hr)

    ! calculation of the rightgoing flux from the left cell (fmr)
    call primitives(st1,p1,d,u,v,w,cs,dex,0,1,1,nodir,null)
    call VanLeerFluxVectorSplit(d,u,v,w,cs,dex,RHVX,RHVY,RHVZ,dum,fmr)
    
    ! calculation of the left going flux from the right cell (fpl)
    call primitives(st1,p1,d,u,v,w,cs,dex,2,1,1,nodir,null)
    call VanLeerFluxVectorSplit(d,u,v,w,cs,dex,RHVX,RHVY,RHVZ,fpl,dum)

    ! calculation of the up going flux from the lower cell (gmr)
    call primitives(st1,p1,d,u,v,w,cs,dex,1,0,1,nodir,null)
    call VanLeerFluxVectorSplit(d,v,u,w,cs,dex,RHVY,RHVX,RHVZ,dum,gmr)

    ! calculation of the down going flux from the upper cell (gpl)
    call primitives(st1,p1,d,u,v,w,cs,dex,1,2,1,nodir,null)
    call VanLeerFluxVectorSplit(d,v,u,w,cs,dex,RHVY,RHVX,RHVZ,gpl,dum)

    ! calculation of the up going flux from the lower cell (hmr)
    call primitives(st1,p1,d,u,v,w,cs,dex,1,1,0,nodir,null)
    call VanLeerFluxVectorSplit(d,w,u,v,cs,dex,RHVZ,RHVX,RHVY,dum,hmr)

    ! calculation of the down going flux from the upper cell (hpl)
    call primitives(st1,p1,d,u,v,w,cs,dex,1,1,2,nodir,null)
    call VanLeerFluxVectorSplit(d,w,u,v,cs,dex,RHVZ,RHVX,RHVY,hpl,dum)

    ! Assemble all fluxes into one df and dg
    do ieq=1,neq
       df(ieq)=(fr(ieq)+fpl(ieq)-fl(ieq)-fmr(ieq))
       dg(ieq)=(gr(ieq)+gpl(ieq)-gl(ieq)-gmr(ieq))
       dh(ieq)=(hr(ieq)+hpl(ieq)-hl(ieq)-hmr(ieq))
    enddo

  end subroutine fluxes

  !=======================================================================
  
  subroutine fluxes2(st1,p1,df,dg,dh)

    !! This routine calculates the second order fluxes for 
    !! Van Leer's Flux Splitting method.
    !! Second order accuracy is achieved by finding the values
    !! of the primitive variables at the cell interfaces. This
    !! is done through (limited) interpolation.

    !! Author: Garrelt Mellema (mellema@strw.leidenuniv.nl)
    !! Date:   14 August 2001
    
    !! Input:
    !! st1 : array containing the state in a (mbc+1)*(mbc+1) 
    !!         stencil around the grid point being updates
    !! p1  : the same for the pressure
    
    !! Output:
    !! df  : state change due to fluxes in the x-direction
    !! dg  : state change due to fluxes in the y-direction
    
    ! define the directional pointers xdir and ydir
    integer,parameter :: xdir=1  ! (for primitives)
    integer,parameter :: ydir=2 ! (for primitives)
    integer,parameter :: zdir=3 ! (for primitives)
    ! define plus and minus (for passing to primitives)
    integer,parameter :: minus=-1
    integer,parameter :: plus=1
    
    real(kind=dp),intent(in) :: st1(1-mbc:1+mbc,1-mbc:1+mbc,1-mbc:1+mbc,neq) !! stencil for state
    real(kind=dp),intent(in) :: p1(1-mbc:1+mbc,1-mbc:1+mbc,1-mbc:1+mbc)      !! stencil for pressure
    real(kind=dp),intent(out) :: df(neq),dg(neq),dh(neq)                  !! state changes 

    real(kind=dp) :: fl(neq),fr(neq),fpl(neq),fmr(neq)
    real(kind=dp) :: gl(neq),gr(neq),gpl(neq),gmr(neq)
    real(kind=dp) :: hl(neq),hr(neq),hpl(neq),hmr(neq)
    
    real(kind=dp) :: d,u,v,w,cs,dex(neuler+1:neq)
    real(kind=dp) :: dum(neq) ! dummy
    integer :: ieq

    ! calculation of the left going flux from the cell (fl)
    call primitives(st1,p1,d,u,v,w,cs,dex,1,1,1,xdir,minus) 
    call VanLeerFluxVectorSplit(d,u,v,w,cs,dex,RHVX,RHVY,RHVZ,fl,dum)
    ! calculation of the right going flux from the cell (fr)
    call primitives(st1,p1,d,u,v,w,cs,dex,1,1,1,xdir,plus)
    call VanLeerFluxVectorSplit(d,u,v,w,cs,dex,RHVX,RHVY,RHVZ,dum,fr)
    
    ! calculation of the down going flux from the cell (gl)
    call primitives(st1,p1,d,u,v,w,cs,dex,1,1,1,ydir,minus)
    call VanLeerFluxVectorSplit(d,v,u,w,cs,dex,RHVY,RHVX,RHVZ,gl,dum)
    ! calculation of the up going flux from the cell (gr)
    call primitives(st1,p1,d,u,v,w,cs,dex,1,1,1,ydir,plus)
    call VanLeerFluxVectorSplit(d,v,u,w,cs,dex,RHVY,RHVX,RHVZ,dum,gr)
    
    ! calculation of the down going flux from the cell (hl)
    call primitives(st1,p1,d,u,v,w,cs,dex,1,1,1,zdir,minus)
    call VanLeerFluxVectorSplit(d,w,u,v,cs,dex,RHVZ,RHVX,RHVY,hl,dum)
    ! calculation of the up going flux from the cell (hr)
    call primitives(st1,p1,d,u,v,w,cs,dex,1,1,1,zdir,plus)
    call VanLeerFluxVectorSplit(d,w,u,v,cs,dex,RHVZ,RHVX,RHVY,dum,gr)
    

    ! calculation of the rightgoing flux from the left cell (fmr)
    call primitives(st1,p1,d,u,v,w,cs,dex,0,1,1,xdir,plus)
    call VanLeerFluxVectorSplit(d,u,v,w,cs,dex,RHVX,RHVY,RHVZ,dum,fmr)
    
    ! calculation of the left going flux from the right cell (fpl)
    call primitives(st1,p1,d,u,v,w,cs,dex,2,1,1,xdir,minus)
    call VanLeerFluxVectorSplit(d,u,v,w,cs,dex,RHVX,RHVY,RHVZ,fpl,dum)
    
    ! calculation of the up going flux from the lower cell (gmr)
    call primitives(st1,p1,d,u,v,w,cs,dex,1,0,1,ydir,plus)
    call VanLeerFluxVectorSplit(d,v,u,w,cs,dex,RHVY,RHVX,RHVZ,dum,gmr)
    
    ! calculation of the down going flux from the upper cell (gpl)
    call primitives(st1,p1,d,u,v,w,cs,dex,1,2,1,ydir,minus)
    call VanLeerFluxVectorSplit(d,v,u,w,cs,dex,RHVY,RHVX,RHVZ,gpl,dum)
    
    ! calculation of the up going flux from the lower cell (hmr)
    call primitives(st1,p1,d,u,v,w,cs,dex,1,1,0,ydir,plus)
    call VanLeerFluxVectorSplit(d,w,u,v,cs,dex,RHVZ,RHVX,RHVY,dum,hmr)
    
    ! calculation of the down going flux from the upper cell (hpl)
    call primitives(st1,p1,d,u,v,w,cs,dex,1,1,2,ydir,minus)
    call VanLeerFluxVectorSplit(d,w,u,v,cs,dex,RHVZ,RHVX,RHVY,hpl,dum)
    
    ! Assemble all fluxes into one df and dg
    do ieq=1,neq
       df(ieq)=(fr(ieq)+fpl(ieq)-fl(ieq)-fmr(ieq))
       dg(ieq)=(gr(ieq)+gpl(ieq)-gl(ieq)-gmr(ieq))
       dh(ieq)=(hr(ieq)+hpl(ieq)-hl(ieq)-hmr(ieq))
    enddo

  end subroutine fluxes2

  !=======================================================================

  subroutine primitives(st1,p1,d,u,v,w,cs,dex,i,j,k,idir,ikind)

    !! This routine calculates the primtive variables left or right of 
    !! the cell interfaces for use in calculating the fluxes in 
    !! the van Leer Flux Vector Splitting algorithm.
    !! The calculation of the primitive can be zeroth order
    !! (value at cell interface = value at cell centre)
    !! leading to first order fluxes,
    !! or first order (linear interpolation), leading to second
    !! order fluxes. The first order primitives need to be limited
    !! to avoid oscillations in the final result. This is done by
    !! the function LimitedAverage.
    
    !! Input:
    !! st1  - state variables on a (1+2*mbc)*(1+2*mbc) stencil
    !! p1   - pressure on a (1+2*mbc)*(1+2*mbc) stencil
    !! i,j  - indicates the position of the cell interface
    !! idir - indicates the direction of the cell interface (x or y)
    !! ikind - indicates the left or right interface (minus/plus)
    !! Note: there appears to be some redundancy in these last three
    !!       input items. Consider: interface, left or right, and order 
    !!       as input variables instead
    
    !! Output:
    !! d   - density at the interface
    !! u,v - x and y velocity at the interface
    !! cs  - sound speed at the interface
    !! dex - density of tracer at the interface
    
    !! Author: Garrelt Mellema (mellema@strw.leidenuniv.nl)
    !! Date:   14 August 2001
    
    integer,intent(in) :: i,j,k,idir,ikind
    
    real(kind=dp),intent(in) :: st1(1-mbc:1+mbc,1-mbc:1+mbc,1-mbc:1+mbc,neq) !! stencil for state
    real(kind=dp),intent(in) :: p1(1-mbc:1+mbc,1-mbc:1+mbc,1-mbc:1+mbc)      !! stencil for pressure
    real(kind=dp),intent(out) :: d,u,v,w,cs
    real(kind=dp),dimension(neuler+1:neq),intent(out) :: dex

    real(kind=dp) :: gradl,gradr,grad(neq) ! gradients 
    real(kind=dp) :: p  ! pressure (not exported)
    real(kind=dp) :: pc(neq), pr(neq), pl(neq) ! primitives
    integer ieq
    
    real(kind=dp) :: LimitedAverage ! limiting function
    
    ! find primitives in the central cell
    call Conserved2Primitives(st1,p1,i,j,k,pc)
    
    ! find primitives in neighbouring cells (if needed for second order)
    if (ikind /= 0) then
       select case(idir)
       case(1)
          call Conserved2Primitives(st1,p1,i-1,j,k,pl)
          call Conserved2Primitives(st1,p1,i+1,j,k,pr)
       case(2)
          call Conserved2Primitives(st1,p1,i,j-1,k,pl)
          call Conserved2Primitives(st1,p1,i,j+1,k,pr)
       case(3)
          call Conserved2Primitives(st1,p1,i,j,k-1,pl)
          call Conserved2Primitives(st1,p1,i,j,k+1,pr)
       end select
       ! find the gradients, and take a limited average
       do ieq = 1, neq
          gradl = pc(ieq) - pl(ieq)
          gradr = pr(ieq) - pc(ieq)
          grad(ieq) = LimitedAverage(gradl,gradr)
       enddo
       ! find the interpolated primitives
       d = pc(RHO) + ikind*0.5*grad(RHO)
       u = pc(RHVX) + ikind*0.5*grad(RHVX)
       v = pc(RHVY) + ikind*0.5*grad(RHVY)
       w = pc(RHVZ) + ikind*0.5*grad(RHVZ)
       p = pc(EN) + ikind*0.5*grad(EN)
       do ieq=neuler+1,neq
          dex(ieq) = d*(pc(ieq) + ikind*0.5*grad(ieq))
       enddo
    else
       ! if only zeroth order primitives are needed, do not
       ! calculate the gradients
       d = pc(RHO)
       u = pc(RHVX)
       v = pc(RHVY)
       w = pc(RHVZ)
       p = pc(EN)
       do ieq=neuler+1,neq
          dex(ieq) = d*pc(ieq)
       enddo
    endif
    
    ! Find the sound speed at the interface.
    cs=sqrt(gamma*p/d)
    
  end subroutine primitives
  
  !=======================================================================

  subroutine Conserved2Primitives(st1,p1,i,j,k,prim)

    !! This routine finds the primtive variables from the state 
    !! at a given position.
    
    !! Input:
    !! st1 - state variables on a (1+2*mbc)*(1+2*mbc) stencil
    !! p1  - pressure on a (1+2*mbc)*(1+2*mbc) stencil
    !! i,j,k - position (on stencil) for which to calculate primitives
    
    !! Output:
    !! prim - vector of primitive variables (density, velocities, 
    !!         pressure,tracer)
    
    real(kind=dp),intent(in) :: st1(1-mbc:1+mbc,1-mbc:1+mbc,1-mbc:1+mbc,neq) !! stencil for state
    real(kind=dp),intent(in) :: p1(1-mbc:1+mbc,1-mbc:1+mbc,1-mbc:1+mbc)      !! stencil for pressure
    integer,intent(in) :: i,j,k
    real(kind=dp),intent(out) :: prim(neq)

    integer ieq
    
    ! changes conserved variables (st1) to primitive variables (prim)
    prim(RHO)  = st1(i,j,k,RHO)         ! density
    prim(RHVX) = st1(i,j,k,RHVX)/prim(RHO) ! x-velocity
    prim(RHVY) = st1(i,j,k,RHVY)/prim(RHO) ! y-velocity
    prim(RHVZ) = st1(i,j,k,RHVZ)/prim(RHO) ! y-velocity
    prim(EN)   = p1(i,j,k)            ! pressure 
    do ieq=neuler+1,neq
       prim(ieq) = st1(i,j,k,ieq)/prim(RHO) ! tracers
    enddo
    
  end subroutine Conserved2Primitives

  !=======================================================================

  function LimitedAverage(a,b)

    ! This function takes a limited average of two quantities
      
    real(kind=dp) :: LimitedAverage
    real(kind=dp),intent(in) :: a,b
    
    real(kind=dp) :: ab
    
    ab=a*b
    if (ab.le.0.0d0) then
       LimitedAverage=0.0d0
    else
       LimitedAverage=(ab*b+ab*a)/(a*a+b*b)
    endif
    
  end function LimitedAverage

  !=======================================================================

  subroutine VanLeerFluxVectorSplit(d,u,v,w,cs,dex,nu,nv,nw,fm,fp)
    
    !! This routine calculates the fluxes according to van Leer's 
    !! Flux Vector Splitting algorithm. The calculation is done using the
    !! primitive variables. Account has to be taken to know which
    !! direction one is looking (x or y).
    
    !! Input:
    !! d,u,v,cs dex - density, velocity 1, velocity 2, sound speed, tracer(s)
    !! nu,nv        - indexes for momentum fluxes:
    !!                  if velocity 1 = x-velocity nu=2, nv=3 should be used
    !!                  if velocity 1 = y-velocity nu=3, nv=2 should be used
    
    !! Author: Garrelt Mellema (mellema@strw.leidenuniv.nl)
    !! Date:   14 August 2001
    
    real(kind=dp),intent(in) :: d,u,v,w,cs,dex(neuler+1:neq)
    integer,intent(in) :: nu,nv,nw
    real(kind=dp),intent(out) ::  fm(neq),fp(neq) ! fluxes

    real(kind=dp) ::  ft(neq) ! fluxes
    real(kind=dp) :: am ! mach number
    integer :: ieq

    ! calculate mach number
    am=u/cs
    
    ! Find the `supersonic' flux
    ft(RHO)=d*u
    ft(nu)=d*(u*u+cs*cs/gamma)
    ft(nv)=d*u*v
    ft(nw)=d*u*w
    ft(EN)=d*u*(0.5*(u*u+v*v+w*w)+cs*cs/gamma1)
    do ieq=neuler+1,neq
       ft(ieq)=dex(ieq)*u
    enddo
    
    ! Consider the three cases: Mach > 1, -1 < Mach < 1, Mach < -1
    if (am.ge.1.0) then
       ! Mach > 1
       ! assign supersonic flux to `plus' flux
       do ieq=1,neq
          fp(ieq)=ft(ieq)
          fm(ieq)=0.0d0
       enddo
    elseif (abs(am).lt.1.0d0) then
       ! -1 < Mach < 1
       ! Use van Leer's recipe to split the `supersonic' flux vector
       ! between plus and minus fluxes
       fp(RHO)=d*cs*0.25*(am+1.0)*(am+1.0)
       fp(nu)=fp(1)*cs*(gamma1*am+2.0d0)/gamma
       !! fp(nu)=fp(1)*(gamma1*u+2.*cs)/gamma
       fp(nv)=fp(1)*v
       fp(nw)=fp(1)*w
       fp(EN)=fp(1)*(0.5d0*(v*v+w*w)+ &
            (gamma1*u+2.0d0*cs)*(gamma1*u+2.0d0*cs)/ &
            (2.0d0*(gamma*gamma-1.0d0)))
       !**          fp(4)=fp(1)*(cs*cs*(gamma1*am+2.0d0)*(gamma1*am+2.0d0)/
       !**     $        (2.0d0*(gamma*gamma-1.0d0))+0.5d0*v*v)
       do ieq=neuler+1,neq
          fp(ieq)=(dex(ieq)/d)*fp(1)
       enddo
       do ieq=1,neq
          fm(ieq)=ft(ieq)-fp(ieq)
       enddo
    else
       ! Mach < -1
       ! assign supersonic flux to `minus' flux
       do ieq=1,neq
          fp(ieq)=0.0d0
          fm(ieq)=ft(ieq)
       end do
    endif
    
  end subroutine VanLeerFluxVectorSplit

end module integrator
