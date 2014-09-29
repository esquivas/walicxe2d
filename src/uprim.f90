!========================================================================
!   calculates the primitive variables (rho,u,v,pgas) from the
!   integration variables (rho,u rho, v rho, etot)
!========================================================================
subroutine uprim(prim,uu,t)
  use parameters
  implicit none
  real,    intent(out), dimension(neq)  :: prim
  real,    intent(in),  dimension(neq)  :: uu  
  real,    intent(out)                  :: T
  real :: r
#ifdef COOLINGH
  real :: dentot
#endif
  !
  r=max(uu(1),1e-15)
  prim(1)=r
  prim(2)=uu(2)/r 
  prim(3)=uu(3)/r
  !
#ifdef PASSIVES
  prim(5:neqpas) = uu(5:neqpas)
#endif
  !
  prim(4)=( uu(4)-0.5*r*(prim(2)**2+prim(3)**2) ) /cv 
  prim(4)=max(prim(4),1e-15)
  !
  !-----------------------------------------
#ifdef ADIABATIC
  !
  T=prim(4)/r
  !
#endif
  !-----------------------------------------
#ifdef COOLINGH
  !
  dentot=(2.*r-prim(5))
  dentot=max(dentot,1e-15)
  !
  T=max(1.,(prim(4)/dentot)*Tempsc )
  prim(4)=dentot*T/Tempsc
  !
#endif
  !-----------------------------------------
#ifdef COOLINGDMC
  !
  ! assumes it is fully ionized
  r=max(r,1e-15)
  !
  T=max(1.,(prim(4)/(r))*Tempsc )
  prim(4)=r*T/Tempsc
  !
#endif
  !-----------------------------------------
#ifdef COOLINGBBC
  !
  dentot= prim(5) + prim(6) + prim(7) + prim(8) + prim(9) + prim(11)
  !
  T = (prim(4)/dentot/Rg)*vsc*vsc
  !
#endif
  !-----------------------------------------
end subroutine uprim  
!========================================================================
