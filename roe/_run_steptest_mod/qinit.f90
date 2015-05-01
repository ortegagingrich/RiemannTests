subroutine qinit(meqn,mbc,mx,xlower,dx,q,maux,aux)

    ! Set initial conditions for the q array.
    ! This default version simply sets eta = max(h + b,0)

    ! For more specific initial conditions
    !  copy this to an application directory and
    !  loop over all grid cells to set values of q(1:meqn, 1:mx).

    !use geoclaw_module, only: dry_tolerance !uncomment if needed
    !use geoclaw_module, only: grav  !uncomment if needed

    implicit none

    integer, intent(in) :: meqn,mbc,mx,maux
    real(kind=8), intent(in) :: xlower,dx
    real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc)
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc)

    !locals
    integer :: i
    real(kind=8) :: xcell,xp

    real(kind=8) :: eta0,eta1,eta2,eta3,damloc1,damloc2,damloc3

       eta0 = 1.0d0
       eta1 = 0.50d0
       eta2 = 0.25d0
       eta3= 0.0d0
       damloc1=-30.0d0!.3d2
       damloc2=0.0d0
       damloc3=30.0d0

    do i=1,mx
      xcell = xlower + (i-0.5d0)*dx
      !call mapc2p(xcell,xp)
      !xp = xcell
      q(2,i)=0.d0
      if (xcell<damloc1) then
         q(1,i) = max(0.0,eta0-aux(1,i))
      else if (xcell<damloc2) then
         q(1,i)=max(0.0,eta1-aux(1,i))
      else if (xcell<damloc3) then
         q(1,i)=max(0.0,eta2-aux(1,i))
      else
         q(1,i) = max(0.0,eta3-aux(1,i))
      endif
      !set values for simulataneous roe test
      q(3,i)=q(1,i)
      q(4,i)=q(2,i)

   enddo


end subroutine qinit
