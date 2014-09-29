!=======================================================================
!   calculates HLL fluxes from the primitive variables 
!   on all the domain
!   choice --> 1, uses u for the 1st half of timestep (first order)
!          --> 2, uses u and limiter for second order timestep
!=======================================================================
subroutine hllfluxes(nb,choice,f,g)
#ifdef HLL
  use parameters
  use globals
  implicit none
  integer, intent(in) :: nb, choice
  real,    intent(out),dimension(neq,nxmin:nxmax,nymin:nymax) :: f, g
  !
  real, dimension(neq) :: priml, primr, primll, primrr, ff, uu
  real :: csl, csr, sl, sr, sst, slmul, srmur, rholul, rhorur
  integer :: i,j, ip, jp, im, jm ,ip2, jp2
  !
  ! computes the fluxes usig the primitives
  select case(choice)
     !------------------------------------------------------------------
  case(1)        ! 1st half timestep
     !
     do i=0,nx
        do j=0,ny
           !
           ip=i+1
           jp=j+1
           !------------------------------------------------------------
           !   x direction
           priml(:)=primit(nb,:,i ,j )
           primr(:)=primit(nb,:,ip,j )
           !
           call primfhll(priml,primr,ff)
           f(:,i,j)=ff(:)
           !   y direction 
           priml(:)=primit(nb,:,i ,j )
           primr(:)=primit(nb,:,i, jp)
           call swapy(priml,neq)
           call swapy(primr,neq)
           !
           call primfhll(priml,primr,ff)
           call swapy(ff,neq)
           g(:,i,j)=ff(:)
           !------------------------------------------------------------
           !
        end do
     end do
     !------------------------------------------------------------------
  case (2)   !  2nd half timestep
     !------------------------------------------------------------------
     !
     do i=0,nx
        do j=0,ny
         ip=i+1
           ip2=i+2
           im=i-1
           jp=j+1
           jp2=j+2
           jm=j-1
           !------------------------------------------------------------
           !   x direction
           priml (:)=primit(nb,:,i,  j )
           primr (:)=primit(nb,:,ip, j )
           primll(:)=primit(nb,:,im, j )
           primrr(:)=primit(nb,:,ip2,j )
           !
           call limiter(primll,priml,primr,primrr,neq)
           !
           call primfhll(priml,primr,ff)
           f(:,i,j)=ff(:)
           !   y direction
           priml (:)=primit(nb,:,i,j   )
           primr (:)=primit(nb,:,i,jp  )
           primll(:)=primit(nb,:,i,jm  )
           primrr(:)=primit(nb,:,i,jp2 )
           call swapy(priml, neq)
           call swapy(primr, neq)
           call swapy(primll,neq)
           call swapy(primrr,neq)
           call limiter(primll,priml,primr,primrr,neq)
           !
           call primfhll(priml,primr,ff)
           call swapy(ff,neq)
           g(:,i,j)=ff(:)
           !------------------------------------------------------------
        end do
     end do
     !------------------------------------------------------------------  
     !
  end select
  !
contains
  !
!=======================================================================
!   calculates the F HLLC fluxes from the primitive variables
!=======================================================================
  subroutine primfhll(priml,primr,ff)
    use parameters
    implicit none
    real, dimension(neq), intent(in)  :: priml, primr
    real, dimension(neq), intent(inout) :: ff
    real, dimension(neq)                :: fl, fr, ul, ur
    real :: csl, csr, sl, sr, slmul, srmur, rholul, rhorur, sst
    real :: rhost,ek
    !
    call sound(priml(4),priml(1),csl)
    call sound(primr(4),primr(1),csr)
    !
    sr=max(priml(2)+csl,primr(2)+csr)
    sl=min(priml(2)-csl,primr(2)-csr)
    !
    if (sl.gt.0) then
       call primf(priml,ff)
       return
    endif
    !--------------------------------
    if (sr.lt.0) then
       call primf(primr,ff)
       return
    endif
    !--------------------------------
    call primf(priml,fl)
    call primf(primr,fr)
    call primu(priml,ul)
    call primu(primr,ur)
    !
    ff(:)=(sr*fl(:)-sl*fr(:)+sl*sr*(ur(:)-ul(:)))/(sr-sl)
    return
    !--------------------------------
    !
  end subroutine primfhll
  !
#endif
end subroutine hllfluxes
!======================================================================
