!=======================================================================
!  refinement criterion to allow maximum difference of 1 level between
!  neighboring blocks when refining (i.e. only with irefup)
!=======================================================================
subroutine critneighup(irefup)
  use parameters
  use globals
  implicit none
  logical, intent(out), dimension(nblocks*np)   :: irefup
  !
  integer :: nl,nb,n1
#ifdef MPIP
  integer ::  i, err
  integer, dimension(1:np) :: nbmaxp
  call mpi_allgather(nbmax, 1, mpi_integer, nbmaxp, 1, mpi_integer, &
       mpi_comm_world,err)
#endif 
  !
  do nl=nlevs-1,1,-1
#ifdef MPIP
     do i=1,np
        do nb=(i-1)*nblocks+1,nbmaxp(i)
#else
     do nb=nbmin,nbmax
#endif
           if (irefup(nb) .and. lp(nb,1).eq.nl) then
              !---------------------------------------------------
              !   Left
              n1=lp(nb,5)
              if (n1.ne.-1) then
                 if( lp(n1,1).lt.lp(nb,1) ) irefup(n1)=.true.
              end if
              !---------------------------------------------------
              !   right
              n1=lp(nb,6)
              if (n1.ne.-1) then
                 if( lp(n1,1).lt.lp(nb,1) ) irefup(n1)=.true.
              end if
              !---------------------------------------------------
              !   bottom
              n1=lp(nb,7)
              if (n1.ne.-1) then
                 if( lp(n1,1).lt.lp(nb,1) ) irefup(n1)=.true.
              end if
              !---------------------------------------------------
              !   top
              n1=lp(nb,8)
              if (n1.ne.-1) then
                 if( lp(n1,1).lt.lp(nb,1) ) irefup(n1)=.true.
              end if
              !---------------------------------------------------

           end if
#ifdef MPIP
        end do
#endif
     end do
  end do
  !
  return
end subroutine critneighup
   !=======================================================================
