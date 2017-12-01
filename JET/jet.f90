!=======================================================================
!  Variable jet module
!=======================================================================
module Jet
    implicit none
  !   Jet parameters
  real, save :: RJ, LJ, denJ, TempJ, vJ0, dvJ, tau
  real, save :: posJ(2)

contains

!--------------------------------------------------------------------
!   Here the parameters of the Jet are initialized, and scaled to
!   code units if needed by the impose_sn subroutine
!--------------------------------------------------------------------
subroutine init_jet()
  use parameters, only : rsc, vsc2, yr, tsc, Tempsc
  implicit none
  
  RJ = 7.5e14/rsc  !  jet radius
  LJ = 7.5e14/rsc  !  jet length

  !  jet position
  posJ(1) = 0.
  posJ(2) = 0.

  denJ= 50.
  TempJ = 1.e3/Tempsc
  vJ0 = 200.e5/sqrt(vsc2)
  dvJ =  15.e5/sqrt(vsc2)

  tau = 40.*yr/tsc

end subroutine init_jet

!--------------------------------------------------------------------
!  imposes the jet at a location in 'pos' variable set above, along 
!  with the parameters of a standard variable jet
!  (not precessing because is 2D ! )
!--------------------------------------------------------------------
subroutine impose_jet(u, nb, time)
  use parameters, only : neq, nxmin, nxmax,nymin,nymax &
                , pi, cv
  use globals, only : dx, dy, icoord, lp
  implicit none
  real,    intent(inout) :: u(neq,nxmin:nxmax,nymin:nymax)
  integer, intent(in)    :: nb
  real,    intent(in)    :: time
  integer :: lev, i, j
  real    :: y, x, vJet

  lev=lp(nb,1)

  !  top hat jet variable velocity
  vJet = vj0 + dvJ*sin(2.*pi*time/tau)

  do i=nxmin,nxmax
     do j=nymin,nymax
        
        x= ( float(i+icoord(nb,1)) -0.5 )* dx(lev) ! normalized units
        y= ( float(j+icoord(nb,2)) -0.5 )* dy(lev) ! normalized units
        
        if( (abs(x-posJ(1)) <= LJ).and.( abs(y-posJ(2)) <= RJ ) )  then
          
          !   inside Jet source
          u(1,i,j)= denJ
          u(2,i,j)= denJ*vJet
          u(3,i,j)= 0.
          u(4,i,j)= 0.5*denJ*vJet**2+cv*denJ*1.0001*TempJ
          
          u(5,i,j)=  denJ*0.9999  !  neutral fraction
          u(6,i,j)=  denJ
          
        endif
     end do
  end do

end subroutine impose_jet
!--------------------------------------------------------------------
end module Jet