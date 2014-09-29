!=======================================================================
!   updates the mesh by refining and coarsening all leaf blocks
!=======================================================================
subroutine admesh
  use parameters
  use globals
  implicit none
  logical,dimension(nblocks*np) :: irefup, irefdown
  integer :: nb, nl,son1,son2,son3,son4,dad, i
  logical, dimension(3) :: brother
  integer, dimension(1:np) :: nbmaxp,nbmaxp1
#ifdef MPIP
  integer ::  err
  call mpi_allgather(nbmax, 1, mpi_integer, nbmaxp, 1, mpi_integer, &
       mpi_comm_world,err)
  nbmaxp1=nbmaxp
#else
  nbmaxp(1) =nbmax
  nbmaxp1(1)=nbmax
#endif
  !
  !   flag for refinement/coarsen
  call markref(irefup,irefdown)
  !   mark for refining on proximity of higher resolution grid
  call critneighup(irefup)
  !
  !   refine marked blocks
  !   (in order of levels to ensure neighbors are updated correctly)
  do nl=1,nlevs-1
     do i=1,np
        do nb=(i-1)*nblocks+1,nbmaxp1(i)
           !
           if (lp(nb,1).eq.nl .and. irefup(nb) ) then
              !+++++++++++++++++++++++++++++++++++++++++++++++++++
              if (lp(nb,3).eq.0) then
                 !  1st time a son is produced
                 son1=nbmaxp(i)+1
                 son2=son1+1
                 son3=son1+2
                 son4=son1+3
                 lp(nb,3)=son1
                 !------------------------------------------------
                 nbmaxp(i)=nbmaxp(i)+4
                 if (nbmax.ge.nbmaxtot) then
                    print*,'nblocks exceeded limit',nbmaxp(i),nbmaxtot
                    stop
                 end if
                 !------------------------------------------------
              else
                 son1=lp(nb,3)
                 son2=lp(son1,6)
                 son3=lp(son1,8)
                 son4=lp(son3,6)
              end if
              !+++++++++++++++++++++++++++++++++++++++++++++++++++
              !
              if (rank.eq.nb/nblocks) nbmax=nbmaxp(i)
              !
              call refine(nb,son1,son2,son3,son4)
              call updatelpup(nb,son1,son2,son3,son4,irefup,irefdown)
              !
           end if
           !
        end do
     end do
#ifdef MPIP
     call mpi_barrier(mpi_comm_world,err)
#endif
  end do
  !
  !  verify that the blocks to be coarsen will have max 1 lev diff
  call critneighdown(irefdown)
  !
  !   coarsen in order of levels
   do nl=nlevs,2,-1
      do i=1,np
         do nb=(i-1)*nblocks+1,nbmaxp1(i)
            !    
            if ( irefdown(nb) .and. (lp(nb,4) .eq.1)    &
                              .and. (lp(nb,1).eq. nl) ) then
               !+++++++++++++++++++++++++++++++++++++++++++++++++++
               dad = lp(nb,2)
               son1= nb
               son2= lp(nb,6)
               son3= lp(nb,8)
               son4= lp(lp(nb,6),8)
               brother(1) = irefdown(son2)
               brother(2) = irefdown(son3)
               brother(3) = irefdown(son4)
               !+++++++++++++++++++++++++++++++++++++++++++++++++++
               if ( all(brother(:)) ) then
                  !
                  call coarse(dad,son1,son2,son3,son4)
                  call updatelpdown(dad,son1,son2,son3,son4)
                  !
               end if
            end if
            !
            end do
      end do
      !
#ifdef MPIP
     call mpi_barrier(mpi_comm_world,err)
#endif
   end do
   !
#ifdef MPIP
     call mpi_barrier(mpi_comm_world,err)
#endif
   !
end subroutine admesh
!=======================================================================
