!======================================================================
!   parametrized cooling curve, uses the ionization state of hydrogen
!   and ties the O I and II to it
!======================================================================
subroutine coolingh(dt)
#ifdef COOLINGH
  use parameters
  use globals
  implicit none
  real,    intent(in)  :: dt
  real, dimension(neq) :: uu
  integer :: i, j, nb
  !
  !   *** takes dt in physical units (seconds) ***
  !
  do nb=nbmin,nbmax
     !
     if ( leaf(nb) ) then
        do i=0,nx+1
           do j=0,ny+1
              !
              uu(:)=u(nb,:,i,j)
              call atomic(dt,uu,1.)
              u(nb,:,i,j)=uu(:)
              !
           end do
        end do
     end if
     !
  end do
  !--------------------------------------------------------------------
  !
  !contains
#endif
end subroutine coolingh
!======================================================================
!
#ifdef COOLINGH
!
!--------------------------------------------------------------------
double precision function alpha(T)
  implicit none
  !    calculates the recombination rate (case B)
  !
  real (kind=8), intent(in) :: T
  alpha=2.55d-13*(1.d4/T)**0.79
  !
end function alpha
!--------------------------------------------------------------------
double precision function alpha1(T)
  implicit none
  !    calculates the recombination rate to level 1
  !
  real (kind=8), intent(in) :: T
  alpha1=1.57d-13*(1.d4/T)**0.52
  !
end function alpha1
!--------------------------------------------------------------------
double precision function colf(T)
  implicit none
  !    calculates the collisional ionization rate
  !
  real (kind=8), intent(in) :: T
  colf=5.83d-11*sqrt(T)*exp(-157828./T)
  !
end function colf
!--------------------------------------------------------------------
double precision function betah(T)
  implicit none
  real (kind=8), intent(in) ::  T
  real (kind=8)             ::  a
  !
  a=157890./T
  betah=1.133D-24/sqrt(a)*(-0.0713+0.5*log(a)+0.640*a**(-0.33333))
  !
end function betah
!--------------------------------------------------------------------
!   Non-equilibrium energy loss for low temperatures
!   considering the collisional excitation of [O I] and
!   [O II] lines and radiative recombination of H. This
!   cooling rate is multiplied by a factor of 7.033 so
!   that it has the same value as the "coronal equilibrium"
!   cooling rate at a temperature of 44770 K (at temperatures
!   higher than this value, the equilibrium cooling rate is
!   used). The collisional ionization of H and excitation
!   of Lyman-alpha are computed separately, and added to
!   the cooling rate.
double precision FUNCTION ALOSS(X1,X2,DT,DEN,DH0,TE0)
  implicit none
  !
  real, intent(in)          :: DT
  real (kind=8), intent(in) :: X1,X2,DEN,DH0,TE0
  real, parameter :: XION=2.179e-11, XH=0.9,XO=1.e-3
  real, parameter :: C0=0.5732,C1=1.8288e-5,C2=-1.15822e-10,C3=9.4288e-16
  real, parameter :: D0=0.5856,D1=1.55083e-5,D2=-9.669e-12, D3=5.716e-19
  real, parameter :: ENK=118409.,EN=1.634E-11
  real (kind=8) :: TE, DH,DHP,DE,DOI,DOII,OMEGA,OMEGAL,OMEGAH,FRAC,QLA
  real (kind=8) :: ECOLL,CION,EION,EREC,TM,T2,EOI,EOII,EQUIL,FR,EX2,TANH
  real (kind=8) :: BETAH, BETAF, HIICOOL
  !
  Te=MAX(Te0,1000.)
  DH=DEN
  DHP=(1.-X1)*DH
  DE=DHP+1.E-4*DH
  DOI=XO*DH0
  DOII=XO*DHP
  !
  !     SHP=-DH*(X1-X2)/DT
  !
  !   Collisionally excited Lyman alpha
  !
  IF(TE.LE.55000.) OMEGA=C0+TE*(C1+TE*(C2+TE*C3))
  IF(TE.GE.72000.) OMEGA=D0+TE*(D1+TE*(D2+TE*D3))
  IF(TE.GT.55000..AND.TE.LT.72000.) THEN
     OMEGAL=C0+TE*(C1+TE*(C2+TE*C3))
     OMEGAH=D0+TE*(D1+TE*(D2+TE*D3))
     FRAC=(TE-55000.)/17000.
     OMEGA=(1.-FRAC)*OMEGAL+FRAC*OMEGAH
  END IF
  QLA=8.6287E-6/(2.*SQRT(TE))*OMEGA*EXP(-ENK/TE)
  ECOLL=DE*DH0*QLA*EN
  ECOLL=MAX(ECOLL,0.)
  !
  !   Hydrogen recombination and collisional ionization
  !
  CION=5.834E-11*SQRT(TE)*EXP(-1.579E5/TE)
  !     AREC=2.61E-10/TE**0.7
  !     EREC=DE*DHP*(BETAH(TE)+AREC*XION)
  !     EION=SHP*XION
  EION=DE*DH0*CION*XION
  !     EION=0.
  EREC=DE*DHP*(BETAH(TE))
  EREC=MAX(EREC,0.)
  !
  !   [O I] and [O II] coll. excited lines
  !
  TM=1./TE
  T2=TM*TM
  EOI=DE*DOI*10.**(1381465*T2-12328.69*TM-19.82621)
  EOII=DE*DOII*10.**(-2061075.*T2-14596.24*TM-19.01402)
  EOI=MAX(EOI,0.)
  EOII=MAX(EOII,0.)
  !
  !   free-free cooling    
  !
  BETAF=1.3*1.42E-27*TE**0.5
  HIICOOL=DE*DHP*BETAF
  !
  !   Equilibrium cooling (for high Te)
  !
  EQUIL=(1.0455E-18/TE**0.63)*(1.-EXP(-(TE*1.E-5)**1.63))*DE*DEN+HIICOOL
  !
  !   switch between non-equil. and equil. ionization cooling
  !
  IF(TE.LE.44770.) FR=0.
  IF(TE.GE.54770.) FR=1.
  IF(TE.GT.44770..AND.TE.LT.54770.) THEN
     EX2=EXP(-2.*(TE-49770.)/500.)
     TANH=(1.-EX2)/(1.+EX2)
     FR=0.5*(1.+TANH)
  END IF
  !
  ALOSS=ECOLL+EION+(EREC+7.033*(EOI+EOII))*(1.-FR)+EQUIL*FR
  !
END FUNCTION ALOSS
!--------------------------------------------------------------------
!    calculates the new ionization state and energy density
!      using a time dependent ionization calculation and an
!      approximate time dependent cooling calculation
!--------------------------------------------------------------------
subroutine atomic(dt,uu,tau)
  use parameters
  implicit none
  real, intent(in)                 :: dt, tau
  real, intent(out),dimension(neq) :: uu
  real, dimension(neq) :: prim
  real                 :: T, phi0, psi0

  real (kind=8) :: etau, psi, phi, dh, y0, fpn, g0, e, y1, t1,dh0, al
  real (kind=8) :: gain, tprime, ce, ALOSS
  !COMMON/PHOTO/PHI0,PSI0
  !COMMON/JETENV/DE,TE,DC,TC,XC,RC,AMP,AL,CTH,STH,RAD0
  !
  !    these need to be double precision in order for
  !      the ionization calculation to work
  !
  real (kind=8) :: te, col,rec,a,b,c,d,colf,alpha
  !
  !    parameters
  !      xi - neutral carbon abundance (for non-zero electron density
  !      boltzm - Boltzmann's constant
  !
  real (kind=8), parameter ::  xi=1.d-4,boltzm=1.3807d-16,AMASS=2.158d-24
  !
  !!!=   AE - 18/10/07
  !!!  phi0=0.
  !!!  psi0=0.
  !
  !   atenuate photoionization with optical depth
  !
  !!!   etau=exp(-tau)
  !!!   phi=phi0*etau
  !!!   psi=psi0*etau
  !
  !   solve for the ionization fraction and the internal energy
  !
  call uprim(prim,uu,T)     !# temperature
  col=colf(dble(t))         !# collisional ionization rate 
  rec=alpha(dble(t))        !# rad. recombination rate
  y0=dble( uu(5)/uu(1) )    !# neutral H fraction  
  dh=dble( uu(1) )     !# H density
  !!!   fpn=phi/dh               !# ionizing flux per nucleus
  !
  !print*,T
  !    solve for the new neutral fraction using the analytical
  !    solution (see notes)
  !
  a=rec+col
  !!!b=-((2.+xi)*rec+(1.+xi)*col+fpn)
  b=   -((2.+xi)*rec+(1.+xi)*col    )
  c=(1.+xi)*rec
  d=sqrt(b**2-4.*a*c)
  g0=(2.*a*y0+b+d)/(2.*a*y0+b-d)
  e=exp( -d*dh*dble(dt) )
  !   
  y1=(-b-d*(1.+g0*e)/(1.-g0*e))/(2.*a) !# the new neutral fraction
  y1=min(y1,1.)
  y1=max(y1,0.0001)

  uu(5)=sngl(y1)*uu(1)            !# update the uu array
  !uu(5)=0.0001*uu(1)   
  !y0=0.0001
  !y1=0.0001
  !
  !    find the new total energy using the cooling and heating rates
  !    and assuming that the cooling goes approximately linear with 
  !    temperature
  !    
  dh0=dble( uu(5) )
  al=ALOSS(y0,y1,dt,dh,dh0,dble(t))/dh**2
  !if(al.lt.0.) write(*,*) 'que paso !'
  !if(al.lt.0.) al=0.
  !  al=al*(1.-(0.5e4/max(1.e4,t))**4)
  !f(t.le.1.e4) al=al*dble((t/1.e4)**4)
  !!!  gain=dble(psi)*dh0
  !!!  tprime=gain*dble(t)/(dh**2*al)+100.
  tprime=10.
  !
  ce=(2.*dh*al)/(3.*boltzm*t)
  t1=tprime+(t-tprime)*exp(-ce*dt) !# new temperature
  !
  t1=max(t1,0.1*dble(t) )
  t1=min(t1,10.*dble(t) )
  t1=max(t1,tprime)
  !
  uu(4) = cv*(2.*uu(1)-uu(5))*sngl(t1)/Tempsc        &
       +0.5*prim(1)*(prim(2)**2+prim(3)**2)    
  !
end subroutine atomic
!
#endif
!======================================================================
