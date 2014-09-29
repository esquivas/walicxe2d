!=======================================================================
!   calculates HLLC fluxes from the primitive variables 
!   on all the domain
!   choice --> 1, uses u for the 1st half of timestep (first order)
!          --> 2, uses u and limiter for second order timestep
!=======================================================================
subroutine hllcfluxes(nb,choice,f,g)
#ifdef HLLC
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
           call primfhllc(priml,primr,ff)
           f(:,i,j)=ff(:)
           !   y direction 
           priml(:)=primit(nb,:,i ,j )
           primr(:)=primit(nb,:,i, jp)
           call swapy(priml,neq)
           call swapy(primr,neq)
           !
           call primfhllc(priml,primr,ff)
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
           call primfhllc(priml,primr,ff)
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
           call primfhllc(priml,primr,ff)
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
  subroutine primfhllc(priml,primr,ff)
    use parameters
    implicit none
    real, dimension(neq), intent(in)  :: priml, primr
    real, dimension(neq), intent(inout) :: ff
    real, dimension(neq)              :: uu, uuk
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
    slmul=sl-priml(2)
    srmur=sr-primr(2)
    rholul=priml(1)*priml(2)
    rhorur=primr(1)*primr(2)
    !
    sst = (srmur*rhorur-slmul*rholul-primr(4)+priml(4) )        &  
         / (srmur*primr(1)-slmul*priml(1) )
    !--------------------------------
    if (sst.ge.0.) then
       rhost=priml(1)*(slmul)/(sl-sst)
       ek= 0.5*priml(1)*(priml(2)**2+priml(3)**2)+cv*priml(4)
       !
       uuk(1)=rhost
       uuk(2)=rhost*sst
       uuk(3)=rhost*priml(3)
       uuk(4)=rhost*( ek/priml(1)+(sst-priml(2))*(sst+priml(4)/(priml(1)*slmul)) )
       !
#ifdef PASSIVES
       uuk(5:neqpas)=rhost*priml(5:neqpas)/priml(1)
#endif
       !
       call primf(priml,ff)
       call primu(priml,uu)
       ff(:)=ff(:) + sl*( uuk(:)-uu(:) )
       return
    endif
    !--------------------------------
    if (sst.le.0.) then
       rhost=primr(1)*(srmur)/(sr-sst)
       ek= 0.5*primr(1)*(primr(2)**2+primr(3)**2)+cv*primr(4)
       !
       uuk(1)=rhost
       uuk(2)=rhost*sst
       uuk(3)=rhost*primr(3)
       uuk(4)=rhost*( ek/primr(1)+(sst-primr(2))*(sst+primr(4)/(primr(4)*srmur)) )
       !
#ifdef PASSIVES
       uuk(5:neqpas)=rhost*primr(5:neqpas)/primr(1)
#endif
       !
       call primf(primr,ff)
       call primu(primr,uu)
       ff(:)=ff(:) + sr*( uuk(:)-uu(:) )
       return
    endif
    !--------------------------------
    !
    print'(a,4e10.3)', 'Error in hllc' ,priml(4),primr(4),priml(1),primr(1)
    stop
    !
  end subroutine primfhllc
  !
#endif
end subroutine hllcfluxes
!======================================================================
