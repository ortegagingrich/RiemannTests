!Subroutine containing modularized version of the David George solver above
!Note: this is for a single grid cell
subroutine shallow_full(i,hLin,hRin,huLin,huRin,bLin,bRin,swfull,fwfull)
	
	use geoclaw_module, only: dry_tolerance,grav
	implicit none
	
	!input
	integer, intent(in) :: i
	double precision, intent(in) :: hLin,hRin,huLin,huRin,bLin,bRin
	
	!output
	double precision, intent(out) :: swfull(2),fwfull(2,2) !component,wave
	
	!Local
	integer :: m,mw,maxiter
	double precision :: hL,hR,huL,huR,hvL,hvR,bL,bR,uL,uR,vR,vL,phiR,phiL
	double precision :: sL,sR,sRoe1,sRoe2,sE1,sE2,uhat,chat
	double precision :: hstartest,hstarHLL,sLtest,SRtest
	double precision :: wall(2), fw(3,3),sw(3)
	double precision :: g,drytol
	
	g=grav
	drytol=dry_tolerance
	
	!Riemann problem variables
	hL=hLin
	hR=hRin
	huL=huLin
	huR=huRin
	hvL=0.d0
	hvR=0.d0
	bL=bLin
	bR=bRin
	
	
	
	!inform of a bad riemann problem from the start
	if((hL.lt.0.d0).or.(hR.lt.0.d0)) then
		write(*,*) 'Negative input: hL,hR,i=',hL,hR,i
	endif
	
	do mw=1,2
		swfull(mw)=0.d0
		do m=1,2
			fwfull(m,mw)=0.d0
		enddo
	enddo
	
	!if completely dry, do nothing
	if((hL.lt.drytol).and.(hR.lt.drytol)) then
		return
	endif
	
	
	!check for wet/dry boundary
	sE1=1.d99
	sE2=-1.d99!these are overwritten at the earliest opportunity
	wall(2)=1.d0 !used to indicate if the wall solution is used for wet/dry problem
	wall(1)=1.d0
	if(hR.le.drytol) then !perform wall test
		!Einfeldt speed for wall problem:
		sLtest=min(-sqrt(g*hL),huL/hL-sqrt(g*hL)) 
		hstartest=hL-(huL/sLtest) !middle state of wall problem
		if(hstartest+bL.lt.bR) then !use wall test results
			wall(2)=0.d0  !indicate that right cell is a wall
			hR=hL
			huR=-huL
			bR=bL
		endif
	elseif(hL.le.drytol) then !other wall test
		sRtest=max(sqrt(g*hR),huR/hR+sqrt(g*hR))
		hstartest=hR-(huR/sRtest)
		if(hstartest+bR.lt.bl) then
			wall(1)=0.d0
			hL=hR
			huL=-huR
			bL=bR
		endif
	endif
	
	!deal with general dry states; note wall problems are no longer considered dry
	!states at this point
	if(hR.gt.drytol) then
		uR=huR/hR
		vR=hvR/hR
		phiR=0.5d0*g*hR**2+huR**2/hR
	else
		hR=0.d0
		huR=0.d0
		uR=0.d0
		vR=0.d0
		phiR=0.d0
		sE2=max(sE2,huL/hL+2.d0*sqrt(g*hL)) !speed of wet-dry interface
	endif
	
	if(hL.gt.drytol) then
		uL=huL/hL
		vL=hvL/hL
		phiL=0.5d0*g*hL**2+huL**2/hL
	else
		hL=0.d0
		huL=0.d0
		uL=0.d0
		vL=0.d0
		phiL=0.d0
		sE1=min(sE1,huR/hR-2.d0*sqrt(g*hR))
	endif
	
	!determine wave speeds
	sL=uL-sqrt(g*hL) ! 1 wave speed of left state
	sR=uR+sqrt(g*hR) ! 2 wave speed of right state
	
	!calculate Roe speeds
	uhat=(sqrt(g*hL)*uL+sqrt(g*hR)*uR)/(sqrt(g*hR)+sqrt(g*hL)) ! Roe average
	chat=sqrt(0.5d0*g*(hL+hR)) ! Roe average
	sRoe1=uhat-chat
	sRoe2=uhat+chat
	
	!calculate Einfeldt speeds 
	sE1=min(sE1,min(sL,sRoe1))
	sE2=max(sE2,max(sR,sRoe2))
	
	!Now, solve the Riemann problem
	maxiter=1
	
	!We'll deal with this later ...
	call riemann_aug_JCP(maxiter,3,3,hL,hR,huL,huR,hvL,hvR,bL,bR,uL,uR,vL,vR,&
	& phiL,phiR,sE1,sE2,drytol,g,sw,fw)
	
	!don't update using waves going into walls
	swfull(1)=sw(1)*wall(1)
	swfull(2)=sw(3)*wall(2)
	!corrector wave speed is not recorded; only the einfeldt waves
	!fluxes are included with the other waves
	
	do m=1,2
		!einfeldt-type waves
		fwfull(m,1)=fw(m,1)*wall(1)
		fwfull(m,2)=fw(m,3)*wall(2)
		!deal with the corrector wave; update the appropriate fwave
		if(sw(2)>0.0) then
			fwfull(m,2)=fwfull(m,2)+fw(m,2)*wall(2)
		else
			fwfull(m,1)=fwfull(m,1)+fw(m,2)*wall(1)
		endif
	enddo
	
	
	return
end subroutine shallow_full



