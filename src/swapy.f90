!=======================================================================
!   swaoping routine (in y-direction)
!=======================================================================
subroutine swapy(prim,neq)
  real, intent(inout), dimension(neq) :: prim
  integer, intent(in) :: neq
  real :: aux
  !
  aux=prim(2)
  prim(2)=prim(3)
  prim(3)=aux
  !
end subroutine swapy
!=======================================================================
