!=======================================================================
!  refinement criteria based on grad(P), w/am2  also grad(rho)
!=======================================================================
subroutine criteria(nb,irefup,irefdown)
  use parameters
  use globals
  implicit none
  integer, intent(in) :: nb
  logical, intent(out),   dimension(nblocks*np) :: irefup,irefdown
  integer :: i,j,ip,im,jp,jm
  real, parameter :: floor=1e-15,tol =0.05,tol2=1.0
  real            :: grad,am, am2
  !
  grad=0.
  !
   do i=1,nx
     do j=1,ny
        !
        ip=i+1
        jp=j+1
        im=i-1
        jm=j-1
        !
        am = max(floor, abs( primit(nb,4,i,j)) )
        !
        am2= max(floor, abs( primit(nb,1,i,j)) )
        !
        grad=max( grad, abs(primit(nb,4,ip,j)-primit(nb,4,i,j))/am,    &
                        abs(primit(nb,4,im,j)-primit(nb,4,i,j))/am,    &
                        abs(primit(nb,4,i,jp)-primit(nb,4,i,j))/am,    &
                        abs(primit(nb,4,i,jm)-primit(nb,4,i,j))/am     )
        !
        grad=max( grad, abs(primit(nb,1,ip,j)-primit(nb,1,i,j))/am2,   &
                        abs(primit(nb,1,im,j)-primit(nb,1,i,j))/am2,   &
                        abs(primit(nb,1,i,jp)-primit(nb,1,i,j))/am2,   &
                        abs(primit(nb,1,i,jm)-primit(nb,1,i,j))/am2    )
        !
        if (grad .ge. tol ) then
           !   only mark if nl <=nlevs-1
           if (lp(nb,1).lt.nlevs) then
              !print*,nb,grad, i,j,am
              irefup(nb)=.true.
              irefdown(nb)=.false.
           end if
              return
        end if
        !
     end do
  end do
  if (grad.le.tol2) irefdown(nb)=.true.
  !
  return
end subroutine criteria
!=======================================================================
