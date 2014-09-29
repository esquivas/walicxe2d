!=======================================================================
!  refinement criterion to allow maximum difference of 1 level between
!  neighboring blocks when coarsening (i.e. only with irefdown)
!=======================================================================
subroutine critneighdown(irefdown)
  use parameters
  use globals
  implicit none
  logical, intent(inout), dimension(nblocks*np)   :: irefdown
  integer :: nb,n1
#ifdef MPIP
  integer ::  i, err
  integer, dimension(1:np) :: nbmaxp
  call mpi_allgather(nbmax, 1, mpi_integer, nbmaxp, 1, mpi_integer, &
       mpi_comm_world,err)
#endif 
  !
#ifdef MPIP
  do i=1,np
     do nb=(i-1)*nblocks+1,nbmaxp(i)
#else
     do nb=nbmin,nbmax
#endif
        if (irefdown(nb)) then
           !---------------------------------------------------
           !   Left
           n1=lp(nb,5)
           if (n1.ne.-1) then
              if ( .not.(leaf(n1)) ) irefdown(nb)=.false.
           endif
           !   Right
           n1=lp(nb,6)
           if (n1.ne.-1) then
              if ( .not.(leaf(n1)) ) irefdown(nb)=.false.
           endif
           !   Bottom
           n1=lp(nb,7)
           if (n1.ne.-1) then
              if ( .not.(leaf(n1)) ) irefdown(nb)=.false.
           endif
           !   Top
           n1=lp(nb,8)
           if (n1.ne.-1) then
              if ( .not.(leaf(n1)) ) irefdown(nb)=.false.
           endif
           !---------------------------------------------------
        end if
     end do
     !
#ifdef MPIP
  end do
  call  mpi_barrier(mpi_comm_world,err)
#endif
  !
end subroutine critneighdown
!=======================================================================
