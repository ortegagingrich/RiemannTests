! =====================================================================
subroutine rp1(maxmx,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,fwave,s,amdq,apdq)
! =====================================================================
!
!solve Riemann problems for the 1D shallow water equations
!    with source-term resulting from variable topography b(x,t)
!    (h)_t + (u h)_x = 0
!    (uh)_t + (uuh + 0.5*gh**2)_x = -g*h*b_x
!
!
!
!     On input,
!     ql contains the state vector at the left edge of each cell
!     qr contains the state vector at the right edge of each cell
!
!     On output, wave contains the fwaves/s,
!                s the speeds,
!                amdq the  left-going flux difference  A**- \Delta q
!               apdq the right-going flux difference  A**+ \Delta q
!
!     Note that the i'th Riemann problem has left state qr(i-1,:)
!     #                                    and right state ql(i,:)
!     From the basic clawpack routine step1, rp1 is called with ql=qr=q.
!
!      This is for use with the Riemann solver(s) used in GeoClaw
!        but for 1D problems. That is, it calls the same solver for left and
!        right states that is used for 2d problems.
!
!        This routine deals with dry-state problems over topography
!        by testing a wall boundary condition, like that done in GeoClaw 2d
!
!        to call other point-wise Riemann solvers alter the call on line 177

    use geoclaw_module, only: dry_tolerance, grav

    implicit none

    ! Input arguments
    integer, intent(in) :: maxmx,meqn,mwaves,mbc,mx,maux

    double precision, intent(in), dimension(meqn, 1-mbc:maxmx+mbc) :: ql,qr
    double precision, intent(in), dimension(maux, 1-mbc:maxmx+mbc) :: auxl
    
    ! Output for debugging purposes
    double precision, intent(out), dimension(maux, 1-mbc:maxmx+mbc) :: auxr

    ! Output arguments
    double precision, intent(out) :: s(mwaves, 1-mbc:maxmx+mbc)
    double precision, intent(out) :: fwave(meqn, mwaves, 1-mbc:maxmx+mbc)
    double precision, intent(out), dimension(meqn, 1-mbc:maxmx+mbc) :: amdq,apdq

    !Local
    integer :: m,i,mw,maxiter
    double precision :: hR,hL,huR,huL,uR,uL,hvR,hvL,vR,vL,phiR,phiL
    double precision :: bR,bL,sL,sR,sRoe1,sRoe2,sE1,sE2,uhat,chat
    double precision :: hstartest,hstarHLL,sLtest,sRtest
    double precision :: wall(2), fw(3,3), sw(3)
    double precision :: g,drytol
    
    !Local: mod
    double precision :: hRr,hLr,huRr,huLr,uRr,uLr
    double precision :: deeptol
    double precision :: a1,a2,deltah,deltahu
    double precision :: fwroe(2,2), swroe(2),fwfull(2,2),swfull(2)
    double precision :: st,ft,roetime,fulltime !time required at time step for each
    double precision :: contime
    !for feeding into full
    double precision :: info

    g=grav
    drytol=dry_tolerance
    deeptol=1.d-3
    
    roetime=0.d0
    fulltime=0.d0
    contime=0.d0
    
    !print *,'meqn=',meqn,'   mwaves=',mwaves

    !loop through Riemann problems at each grid cell
    
    !initialize Riemann problems (should be inside loop; separated for timing)
    do i=2-mbc,mx+mbc
       !inform of a bad riemann problem from the start
       if((qr(1,i-1).lt.0.d0).or.(ql(1,i) .lt. 0.d0)) then
          write(*,*) 'Negative input: hl,hr,i=',qr(1,i-1),ql(1,i),i
       endif
       
       do mw=1,mwaves
          s(mw,i)=0.d0
          do m=1,meqn
             fwave(m,mw,i)=0.d0
          enddo
       enddo
    enddo
    
    
    !roe loop
    call CPU_TIME(st)
      do i=2-mbc,mx+mbc
         !skip problem if in a completely dry area
         if (qr(3,i-1).le.drytol.and.ql(3,i).le.drytol) then
            go to 29
         endif
         
         !If both states exceed depth tolerance, just use a Roe solver:
         !Modularize as a subroutine later
         !simultaneous roe values
         hLr=qr(3,i-1)
         hRr=ql(3,i)
         huLr=qr(4,i-1)
         huRr=ql(4,i)
         bL = auxr(1,i-1)
         bR = auxl(1,i)
         
         
         
         !if water is sufficiently deep, use the roe solver, otherwise the full
         if((hRr.ge.deeptol.and.hLr.ge.deeptol)&
         &.and.(abs(hLr-hRr).lt.(1.d3*(hLr+hRr)))) then
         	!call the roe solver
         	call shallow_roe_jog(hLr,hRr,huLr,huRr,bL,bR,swroe,fwroe)
         	!mark as having used the roe solver
         	auxr(2,i)=1.d0
         else! otherwise, just use the full solver
            call shallow_full(i,info,hLr,hRr,huLr,huRr,bL,bR,swroe,fwroe)
            !mark as having used the full solver
            auxr(2,i)=2.d0
         endif
         
         !set wave speeds
       	 s(3,i)=swroe(1)
       	 s(4,i)=swroe(2)
       	 !compute fwaves
         fwave(3,3,i)=fwroe(1,1)
     	 fwave(4,3,i)=fwroe(2,1)
       	 fwave(3,4,i)=fwroe(1,2)
       	 fwave(4,4,i)=fwroe(2,2)
       	 fwave(1,3,i)=0.d0
       	 fwave(1,4,i)=0.d0
       	 fwave(2,3,i)=0.d0
         fwave(2,4,i)=0.d0
         
         29 continue
      enddo
      call CPU_TIME(ft)
      roetime=roetime+ft-st
      
      
      call CPU_TIME(st)
      do i=2-mbc,mx+mbc
         
         !skip problem if in a completely dry area
         if (qr(1,i-1).le.drytol.and.ql(1,i).le.drytol) then
            go to 30
         endif
         
         !Riemann problem variables
         hL = qr(1,i-1)
         hR = ql(1,i)
         huL = qr(2,i-1)
         huR = ql(2,i)
         bL = auxr(1,i-1)
         bR = auxl(1,i)
         
         !call subroutine to carry out David George's solver
         call shallow_full(i,info,hL,hR,huL,huR,bL,bR,swfull,fwfull)
         !put info in auxilary array
         !print *,i,",",info
         auxr(3,i)=info
         
         
         !set wave speeds
         s(1,i)=swfull(1)
         s(2,i)=swfull(2)
         !set fwaves
         fwave(1,1,i)=fwfull(1,1)
         fwave(2,1,i)=fwfull(2,1)
         fwave(1,2,i)=fwfull(1,2)
         fwave(2,2,i)=fwfull(2,2)
         fwave(3,1,i)=0.d0
         fwave(4,1,i)=0.d0
         fwave(3,2,i)=0.d0
         fwave(4,2,i)=0.d0
         
         
         
         
         30 continue
      enddo
      call CPU_TIME(ft)
      fulltime=fulltime+ft-st
      
      
      call CPU_TIME(st)
      call CPU_TIME(ft)
      contime=contime+ft-st
      
      !print how much time was required by each solver
      !print *,"Time required by full solver: ",fulltime
      !print *,"Time required by Roe solver: ",roetime
      !print *,"Time required for CPU_TIME call: ",contime
      
      !print performance improvement ratio of roe solver
      print *,"Hybrid Solver performance improvement: ",fulltime/roetime
      

      do i=2-mbc,mx+mbc
         do m=1,meqn
            amdq(m,i) = 0.d0
            apdq(m,i) = 0.d0
            do  mw=1,mwaves
               if (s(mw,i).lt.0.d0) then
                  amdq(m,i) = amdq(m,i) + fwave(m,mw,i)
               elseif (s(mw,i).gt.0.d0) then
                  apdq(m,i) = apdq(m,i) + fwave(m,mw,i)
               else
                  amdq(m,i) = amdq(m,i) + .5d0*fwave(m,mw,i)
                  apdq(m,i) = apdq(m,i) + .5d0*fwave(m,mw,i)
               endif
            enddo
         enddo
      enddo

      return
      end subroutine rp1
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
