!======================================================================
!   Cooling routine foro coronal equilibrium (Dalgarno & Mc Cray 1972)
!   The cooling is turned of bellow 10^4 K
!   And the primitives are updated in this routine as well
!======================================================================
subroutine coolingdmc(dt)
#ifdef COOLINGDMC
  use parameters
  use globals
  use dmc_module
  implicit none
  real,    intent(in)  :: dt
  real, dimension(neq) :: uu, prim
  real                 :: T ,Eth0, dens
  real, parameter :: Tmin=1.e4
  real (kind=8)        :: ALOSS, Ce
  integer :: i, j, nb
  !
  do nb=nbmin,nbmax
     !
     if ( leaf(nb) ) then
        !
        do i=0,nx+1
           do j=0,ny+1
              !
              !   get the primitives (and T)
              call uprim(primit(nb,:,i,j),u(nb,:,i,j),T)
              !
              if(T.gt.Tmin) then
                 !
                 Eth0= cv*primit(nb,4,i,j)
                 !
                 Aloss=cooldmc(T)
                 dens=primit(nb,1,i,j)
                 Ce=(Aloss*dble(dens)**2)/(Eth0*Psc)   !cgs
                 !
                 !  apply cooling to primitive and conserved variables
                 primit(nb,4,i,j)=primit(nb,4,i,j)*exp(-ce*dt)
                 !
                 u(nb,4,i,j)= u(nb,4,i,j)-Eth0+cv*primit(nb,4,i,j)
                 !
              end if
              !call atomic(dt,uu,1.)
              !u(nb,:,i,j)=uu(:)
              !
           end do
        end do
     end if
     !
  end do
  !--------------------------------------------------------------------
  !
#endif
end subroutine coolingdmc
!======================================================================
