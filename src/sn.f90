!=======================================================================
!  SN module
!=======================================================================
module SN
  use parameters, only : pc, Msun
  implicit none
  !   SN parameters
  real, parameter :: Esn = 1.e51
  real, parameter :: Rsn = 1.0*pc
  real, parameter :: Msn = 2.0*Msun
  real, parameter :: chi=0.5  !Fraction of kinetic to total energy
  real, parameter :: nenv = 1.0  !  needed to add swept up mass
contains

!--------------------------------------------------------------------
!  Imposes a SN centered at (xy,yc) for block nb,  with properties 
!  set above
!--------------------------------------------------------------------
subroutine impose_sn(u, nb, xc, yc)
  use parameters, only : neq, nxmin, nxmax,nymin,nymax, rsc, rhosc &
                       , vsc2, Tempsc, psc, nlevs, pi
  use globals, only : dx, dy, icoord, lp
  implicit none
  real, intent(inout) :: u(neq,nxmin:nxmax,nymin:nymax)
  integer, intent(in) :: nb
  real, intent(in)    :: xc, yc ! detonation center
  integer :: lev, i, j
  real    :: ekin, x, y, r, dens, vr, Eth, Mtot

  lev=lp(nb,1)
  
  Mtot=Msn+(4.*pi/3.)*nenv*rhosc*Rsn**3
  !print*,Msn, Mtot, (4.*pi/3.)*nenv*rhosc*Rsn**3
  
  do i=nxmin,nxmax
     do j=nymin,nymax
        
        x= ( float(i+icoord(nb,1)) -0.5 )* dx(lev)*rsc
        y= ( float(j+icoord(nb,2)) -0.5 )* dy(lev)*rsc
        
        R=sqrt((x-xc)**2+(y-yc)**2+0.1*dx(nlevs)*rsc )
        
        if (r.le.Rsn) then
           
           !   inside SN
           dens=(3./4./pi)*Msn/(Rsn**3)*1.3/2.1 /rhosc!  + nenv
           vr=(R/Rsn)*Sqrt(10*chi*Esn/(3*Msn))  /sqrt(vsc2)
           Eth=(3./4./pi)*Esn*(1.-chi)/(Rsn**3) /psc
           
           u(1,i,j)= dens
           u(2,i,j)= dens*vr*(x-xc)/(R)
           u(3,i,j)= dens*vr*(y-yc)/(R)
           u(4,i,j)= 0.5*dens*vr**2+Eth
           
           u(5,i,j)=  dens*0.0  !  neutral fraction
           u(6,i,j)=  dens
           
        endif
     end do
  end do

end subroutine impose_sn
!--------------------------------------------------------------------
end module SN
