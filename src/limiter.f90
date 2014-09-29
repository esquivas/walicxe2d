!=======================================================================
!   reconstruction of the primitives to the edge of the cell for the
!   second half of the timestep
!   it uses slope limiters, the limiters avaiable are hardwired to
!   ensure performance, but they can be set from the makefile
!=======================================================================
subroutine limiter(PLL,PL,PR,PRR,neq)
  implicit none
  real, dimension(neq), intent(inout) :: pl,  pr
  real, dimension(neq), intent(in)    :: pll, prr
  integer, intent (in)  :: neq
  real :: dl, dm, dr, al, ar
  integer :: ieq
  !
  do ieq=1,neq
     dl=pl(ieq)-pll(ieq)
     dm=pr(ieq)-pl(ieq)
     dr=prr(ieq)-pr(ieq)
     al=average(dl,dm)
     ar=average(dm,dr)
     pl(ieq)=pl(ieq)+al*0.5
     pr(ieq)=pr(ieq)-ar*0.5
  end do
  !
contains
  real function average(a,b)
    implicit none
    real, intent(in)    :: a, b
    !
#if LIMITER==-1
    !   no average (reduces to 1st order)
    average=0.
#endif
    !
#if LIMITER==0
    !   no limiter
    average=0.5*(a+b)
#endif
    !
#if LIMITER==1
    !   Minmod - most diffusive
    real :: s
    s=sign(1.,a)
    average=s*max(0.,min(abs(a),s*b))
#endif
    !
#if LIMITER==2
    !   Falle Limiter (Van Leer)
    if(a*b.le.0.) then
       average=0.
    else
       average=a*b*(a+b)/(a**2+b**2)
    end if
#endif
    !
#if LIMITER==3
    !   Van Albada
    real, parameter :: delta=1.e-7
    average=(a*(b*b+delta)+b*(a*a+delta))/(a*a+b*b+delta)
#endif
#if LIMITER==4
    !   UMIST Limiter - less diffusive
    real :: s, c, d
    s=sign(1.,a)
    c=0.25*a+0.75*b
    d=0.75*a+0.25*b
    average=min(2.*abS(a),2.*s*b,s*c,s*d)
    average=s*max(0.,average)
#endif
#if LIMITER==5
    !    Woodward Limiter (o MC-limiter; monotonized central difference)
    real :: s, c
    s=sign(1.,a)
    c=0.5*(a+b)
    average=min(2.*abs(a),2.*s*b,s*c)
    average=s*max(0.,average)
#endif
#if LIMITER==6
    !   superbee Limiter (tends to flatten circular waves)
    real :: s, av1, av2
    s=sign(1.,b)
    av1=min(2.*abs(b),s*a)
    av2=min(abs(b),2.*s*a)
    average=s*max(0.,av1,av2)
#endif
    !
  end function average
  !
end subroutine limiter
!=======================================================================
