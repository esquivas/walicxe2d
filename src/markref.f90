!====================================================================
!   marks all the leaf blocks for refinement or de-refinement
!====================================================================
subroutine markref(irefup,irefdown)
  use parameters
  use globals
  implicit none
  logical, intent(out),   dimension(nblocks*np) :: irefup,irefdown
  !
  integer :: i, nb
#ifdef MPIP
  integer :: source, dest1, dest2, nbup, nbdown, nbupp, nbdownp,err
  integer, dimension(nblocks) ::aux1,aux2, aux3
#endif
  !
  !--------------------------------------------------------------------
  !
#ifdef MPIP
  nbup=0
  nbdown=0
#endif
  !
  irefdown(:)=.false.
  irefup(:)=.false.
  do nb=nbmin,nbmax
     if( leaf(nb) .and. lp(nb,1).ge.1 ) then
        call criteria(nb,irefup,irefdown)
!***************************************
#ifdef MPIP
        !
        if (irefup(nb) ) then
           nbup=nbup+1
           aux1(nbup)=nb
           !print*,rank,nb,' up'
        end if
        if (irefdown(nb) ) then
           nbdown=nbdown+1
           aux2(nbdown)=nb
           !print*,rank,nb,' down'
        end if
#endif
!***************************************
     endif
  end do
  !
#ifdef MPIP
  !   take turns by processors to update flags
  do i=0,np-1
     nbupp=nbup
     nbdownp=nbdown
     call mpi_bcast(nbupp  ,1,mpi_integer,i,mpi_comm_world, err) 
     call mpi_bcast(nbdownp,1,mpi_integer,i,mpi_comm_world, err)
     !
     !  refining flag
     if (nbupp.ne.0) then
        aux3(1:nbupp)=aux1(1:nbupp)
        call mpi_bcast(aux3(1:nbupp),nbupp,mpi_integer,i   &
                                      ,mpi_comm_world, err)
        do nb=1,nbupp
           irefup(aux3(nb))=.true.
           !print*,aux3(nb),'up'
        end do
     end if
     !
     !   coarsening flag
     if (nbdownp.ne.0) then
        aux3(1:nbdownp)=aux2(1:nbdownp)
        call mpi_bcast(aux3(1:nbdownp),nbdownp,mpi_integer,i   &
             ,mpi_comm_world, err)
        do nb=1,nbdownp
           irefdown(aux3(nb))=.true.
           !print*,rank,aux3(nb),'down'
        end do
     end if
     !
     !print'(2i3,100l)',i,rank,irefdown(aux3(1:nbdownp))
     !print'(100i7)',i,rank,aux3(1:nbdownp)
     call mpi_barrier(mpi_comm_world,err)
  end do
#endif
  !

  !---------------------------------------------------------------------
  return
end subroutine markref
!=======================================================================
