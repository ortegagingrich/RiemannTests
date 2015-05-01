subroutine shallow_roe_jog(hL,hR,huL,huR,bL,bR,swroe,fwroe)
	!note: this subroutine assumes that the depth of both cells is sufficient
	!If hL or hR are too close to zero, bad stuff will happen.
	
	use geoclaw_module, only: grav
	
	implicit none
	
	!input
	double precision, intent(in) :: hL,hR,huL,huR,bL,bR
	
	!output
	double precision, intent(out) :: swroe(2),fwroe(2,2)!arg1: component, arg2:wave
	
	
	!local
	double precision :: uL,uR,deltah,deltahu
	double precision :: uhat,chat,sswave,s1,s2,a1,a2
	double precision :: shl,shr,o2ch
	double precision :: g
	g=grav
	
	!To Do: subtract off the steady state wave first so that we can deal
	!with varying bathymetry
	
	!compute velocities
	uL=huL/hL
	uR=huR/hR
	
	
	!roe averaged quantities
	shl=sqrt(hL)
	shr=sqrt(hR)
	uhat=(shl*uL+shr*uR)/(shl+shr)
	chat=sqrt(g*0.5d0*(hL+hR))
	
	!wave speeds
	s1=uhat-chat
	s2=uhat+chat
	
	!decompose jump vector
	deltah=hR-hL
	deltahu=huR-huL
	o2ch=0.5d0/chat
	a1=o2ch*(-deltahu+s2*deltah)
	a2=o2ch*( deltahu-s1*deltah)
	
	!fwaves
	fwroe(1,1)=s1*a1
	fwroe(2,1)=fwroe(1,1)*s1
	fwroe(1,2)=s2*a2
	fwroe(2,2)=fwroe(1,2)*s2
	swroe(1)=s1
	swroe(2)=s2
	
	
	
	return
end subroutine shallow_roe_jog
