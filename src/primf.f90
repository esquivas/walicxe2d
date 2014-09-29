!=======================================================================
!   calculates the F fluxes from the primitive variables 
!   at one point
!=======================================================================
subroutine primf(prim,ff)
  use parameters
  implicit none
  real, dimension(neq), intent(in)  :: prim
  real, dimension(neq), intent(out) :: ff
  real :: etot
  !
  etot= 0.5*prim(1)*(prim(2)**2+prim(3)**2)+cv*prim(4)
  !
  ff(1) = prim(1)*prim(2)
  ff(2) = prim(1)*prim(2)*prim(2)+prim(4)
  ff(3) = prim(1)*prim(2)*prim(3)
  ff(4) = prim(2)*(etot+prim(4))
  !
#ifdef PASSIVES
  ff(5:neqpas) = prim(5:neqpas)*prim(2)
#endif
  !
end subroutine primf
!=======================================================================
