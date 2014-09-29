!======================================================================
!   Adds viscosity in all domain (eta is defined in parameters.f90)
!   nb is the block ID
!   up(n)=up(n)+eta\nabla^2(u(n))
!======================================================================
subroutine viscosity(nb,nbmin,nbmaxtot,u,up)
  use parameters
  integer, intent(in) :: nb, nbmin, nbmaxtot
  real,intent(in), dimension(nbmin:nbmaxtot,neq,nxmin:nxmax,nymin:nymax) & 
                      :: u
  real,intent(inout),dimension(nbmin:nbmaxtot,neq,nxmin:nxmax,nymin:nymax) & 
                      :: up
  integer :: i, j
  !
  !--------------------------------------------------------------------
  do i=1,nx
     do j=1,ny
        !
        up(nb,:,i,j)=up(nb,:,i,j)+eta*( u(nb,:,i+1,j)+u(nb,:,i-1,j)   &
                                       +u(nb,:,i,j+1)+u(nb,:,i,j-1)   &
                                       -4.*u(nb,:,i,j) )
           !
     end do
  end do
end subroutine viscosity
!======================================================================
