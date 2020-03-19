!======================================================================
!   upwind time step u->up
!   nb is the block ID
!======================================================================
subroutine step(nb,dt,f,g)
  use parameters
  use globals
  implicit none
  integer, intent(in) :: nb
  real,    intent(in) :: dt
  real,    intent(in),  dimension(neq,nxmin:nxmax,nymin:nymax) :: f, g
  integer :: i, j, im, jm
  real    :: dtdx, dtdy
  real, dimension (neq) :: S
  !
  dtdx=dt/dx(lp(nb,1))
  dtdy=dt/dy(lp(nb,1))
  !
  !--------------------------------------------------------------------
  do i=1,nx
     do j=1,ny
        !
#ifdef CYLINDRICAL
        !   calculate geometrical source terms
        call prims(nb,i,j,S)
#endif
        !
        im=i-1
        jm=j-1
        !
        up(nb,:,i,j)=u(nb,:,i,j)-dtdx*(f(:,i,j)-f(:,im,j))  &
                                -dtdy*(g(:,i,j)-g(:,i,jm))
        !
#ifdef CYLINDRICAL
        !   add the source terms
        up(nb,:,i,j)= up(nb,:,i,j)+dt*S(:)
#endif
        !
     end do
  end do
  !--------------------------------------------------------------------
#ifdef CYLINDRICAL
  !
contains
  !====================================================================
  !    calculates the geometrical source terms for cylindrical coords
  !    y-->r, x--> z
  !====================================================================
  subroutine prims(nb,i,j,S)
    use parameters
    use globals
    implicit none
    integer, intent(in)                  :: nb, i,j
    real,    intent(out),dimension (neq) :: S
    real :: radius,rho,H
    !
    radius = ( float(j+icoord(nb,2)) - 0.5 ) * dy( lp(nb,1) )
    !if (icoord(nb,2).eq.0) radius=radius+2.*dy(nlevs)
    !
    rho=primit(nb,1,i,j)
    H=0.5*rho*(primit(nb,2,i,j)*primit(nb,2,i,j)                     &
              +primit(nb,3,i,j)*primit(nb,3,i,j) )                   &
              +(cv+1.)*primit(nb,4,i,j)
    !
    S(1) = -      rho            *primit(nb,3,i,j)/radius
    S(2) = - rho*primit(nb,2,i,j)*primit(nb,3,i,j)/radius
    S(3) = - rho*primit(nb,3,i,j)*primit(nb,3,i,j)/radius
    S(4) = -       H             *primit(nb,3,i,j)/radius
    !
#ifdef PASSIVES
    S(5:neqpas) = - primit(nb,5:neqpas,i,j)*primit(nb,3,i,j)/radius
#endif
    !
    return
    !
  end subroutine prims
  !====================================================================
#endif

end subroutine step
!======================================================================
