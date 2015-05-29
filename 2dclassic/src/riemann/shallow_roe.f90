subroutine shallow_roe_jog(hL,hR,huL,huR,hvL,hvR,bL,bR,swroe,fwroe)
	!note: this subroutine assumes that the depth of both cells is sufficient
	!If hL or hR are too close to zero, bad stuff will happen.
	
	!also, note that v is the transverse velocity, regardless of which direction
	
	!use geoclaw_module, only: grav
	
	!common block
	common /cparam/ g
	
	implicit none
	
	!input
	double precision, intent(in) :: hL,hR,huL,huR,hvL,hvR,bL,bR
	
	!output
	double precision, intent(out) :: swroe(3),fwroe(3,3)!arg1: component, arg2:wave
	
	
	!local
	integer m
	double precision :: uL,uR,vL,vR,deltah,deltahu,deltahv
	double precision :: ghbar,lamlam,ssjump
	double precision :: uhat,vhat,chat,s1,s2,s3,a1,a2,a3
	double precision :: shl,shr,o2ch
	double precision :: lamL1,lamR2,lamM1,lamM2,huM,hM,beta,trans
	double precision :: g
	g=grav
	
	
	!compute velocities
	uL=huL/hL
	uR=huR/hR
	vL=hvL/hL
	vR=hvR/hR
	
	!compute the steady-state wave
	ghbar=g*0.5d0*(hL+hR)
	lamlam=0.25d0*(uL+uR)**2-ghbar
	ssjump=(bR-bL)*ghbar/lamlam
	!roe averaged quantities
	shl=sqrt(hL)
	shr=sqrt(hR)
	uhat=(shl*uL+shr*uR)/(shl+shr)
	vhat=(shl*vL+shr*vR)/(shl+shr)
	chat=sqrt(ghbar)
	
	!wave speeds
	s1=uhat-chat
	s2=uhat
	s3=uhat+chat
	
	!decompose jump vector
	deltah=hR-hL-ssjump
	deltahu=huR-huL
	deltahv=hvR-hvL
	o2ch=0.5d0/chat
	a1=o2ch*(-deltahu+s2*deltah)
	a2=deltahv-vhat*deltah
	a3=o2ch*( deltahu-s1*deltah)
	
	!fwaves
	fwroe(1,1)=s1*a1
	fwroe(2,1)=fwroe(1,1)*s1
	fwroe(3,1)=fwroe(1,1)*vhat
	
	fwroe(1,3)=s3*a3
	fwroe(2,3)=fwroe(1,3)*s3
	fwroe(3,3)=fwroe(1,3)*vhat
	
	fwroe(1,2)=0.d0
	fwroe(2,2)=0.d0
	fwroe(3,2)=s2*a2
	!fwroe(3,2)=hR*uR*vR-hL*uL*vL-fwroe(3,1)-fwroe(3,3)
	
	swroe(1)=s1
	swroe(2)=s2
	swroe(3)=s3
	
	
	
	!now, try the entropy fix
	!return !toggle efix
	
	
	!First, rule out cases with waves in the same direction
	
	!Case: All waves to the right
	lamL1=uL-sqrt(g*hL)!1-wave speed at the left
	if(lamL1>=0.d0 .and. s1>0.d0) then
		!all waves to the right; no need to do anything
		return
	endif
	
	!Case: All waves to the left
	lamR2=uR+sqrt(g*hR)!wave speed at the right
	if(lamR2<=0.d0 .and. s2<0.d0) then
		!all waves to the left; no need to do anything
		return
	endif
	
	!Get middle state (left state+left wave)
	hM=hL+a1
	huM=huL+a1*s1
	!middle state wave speed
	lamM1=huM/hM-sqrt(g*hM)
	lamM2=huM/hM+sqrt(g*hM)
	
	!check which wave (if any is transonic)
	if(lamM1>0.d0) then
		!1-wave is transonic
		
		!compute fraction to update the left
		beta=(lamM1-s1)/(lamM1-lamL1)
		
		!modify the fwaves
		do m=1,2
			!compute flux to be transferred from 1-wave to 2-wave
			trans=(lamM1/s1)*(1-beta)*fwroe(m,1)
			!apply the change to both waves
			fwroe(m,1)=fwroe(m,1)-trans
			fwroe(m,2)=fwroe(m,2)+trans
		enddo
		
	else if(lamM2<0.d0) then
		!2-wave is transonic
		
		!compute fraction to update the left
		beta=(lamR2-s2)/(lamR2-lamM2)
		
		!modify the fwaves
		do m=1,2
			!compute flux to be transferred from 2-wave to 1-wave
			trans=(lamM2/s2)*beta*fwroe(m,2)
			!apply the change to both waves
			fwroe(m,1)=fwroe(m,1)+trans
			fwroe(m,2)=fwroe(m,2)-trans
		enddo
		
	else
		!neither wave is transonic; do nothing
		return
	endif
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	return
end subroutine shallow_roe_jog
