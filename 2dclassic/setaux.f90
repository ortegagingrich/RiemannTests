subroutine setaux(mbc,mx,my,xlower,ylower,dx,dy,maux,aux)

    ! Called at start of computation before calling qinit, and
    ! when AMR is used, also called every time a new grid patch is created.
    ! Use to set auxiliary arrays aux(1:maux, 1-mbc:mx+mbc, 1-mbc:my+mbc).
    ! Note that ghost cell values may need to be set if the aux arrays
    ! are used by the Riemann solver(s).
    !
    ! This default version does nothing. 
 
    implicit none
    integer, intent(in) :: mbc,mx,my,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy
    real(kind=8), intent(out) ::  aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    
    !set bathymetry
    do i=1-mbc:mx+mbc
    	do j=1-mbc:my+mbc
    		aux(1,i,j)=0.d0
    	enddo
    enddo

end subroutine setaux
