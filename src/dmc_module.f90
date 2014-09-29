!========================================================================
!   Dalgarno & Mc Cray coronal equilibrium cooling module
!   The table is filled at initmain.f90
!   cooldmc(T) interpolates the cooling table and returns the coefficient
!========================================================================
module dmc_module
#ifdef COOLINGDMC
  use globals
  implicit none
  real (kind=8), dimension(2,41) :: cooltab
  !
contains
  
  double precision function cooldmc(T)
    !use globals
    !implicit none
    real , intent(in) :: T
    integer           :: if1
    real (kind=8)     :: T0, T1, C0, C1
    !
    if(T.gt.1e8) then
       cooldmc=0.21D-26*Sqrt(dble(T))
    else
       if1=int(log10(T)*10)-39
       T0=cooltab(1,if1)
       c0=cooltab(2,if1)
       T1=cooltab(1,if1+1)
       c1=cooltab(2,if1+1)
       cooldmc=(c1-c0)*(dble(T)-T0)/(T1-T0)+c0
    end if
    !
  end function cooldmc
  !
#endif
end module dmc_module
!========================================================================
