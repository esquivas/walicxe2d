!========================================================================
!   calculates the sound speed
!========================================================================
subroutine sound(p,d,cs)
  use parameters
  implicit none
  real, intent(in)  :: p, d
  real, intent(out) :: cs
  !
  cs=sqrt(gamma*p/d)
  !
  return
end subroutine sound
!========================================================================
