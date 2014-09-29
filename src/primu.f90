!=======================================================================
!   calculates conserved variables from the primitives at one point
!=======================================================================
subroutine primu(prim,uu)
  use parameters
  implicit none
  real, dimension(neq), intent(in)  :: prim
  real, dimension(neq), intent(out) :: uu
  !
  uu(1) = prim(1)
  uu(2) = prim(1)*prim(2)
  uu(3) = prim(1)*prim(3)
  uu(4) = 0.5*prim(1)*(prim(2)**2+prim(3)**2)+cv*prim(4)
  !
#ifdef PASSIVES
  uu(5:neqpas) = prim(5:neqpas)
#endif
  !
end subroutine primu
!=======================================================================
