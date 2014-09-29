!======================================================================
!   locate neighbors to pass boundary conditions
!   it makes a list for the inner boundaries and a list with the outer
!   boundaries according to:
!
!   innerbounds(ninternal,3)  1: ID source,    2: dest ID, 3:direction
!
!   extbounds= (nextertot,5)  1: ID of source, 2: source rank
!                             3: ID of dest,   4: dest rank 
!                             5: direction
!
!   durection = -1 : left box boundary
!             = -2 : right box boundary
!             = -3 : bottom box boundary 
!             = -4 : top box boundary
!             =  1 : left   boundary at same   resolution
!             =  2 : left   boundary at higher resolution 1
!             =  3 : left   boundary at higher resolution 2 (only MPI)
!             =  4 : right  boundary at same   resolution
!             =  5 : right  boundary at higher resolution 1
!             =  6 : right  boundary at higher resolution 2 (only MPI)
!             =  7 : bottom boundary at same   resolution
!             =  8 : bottom boundary at higher resolution 1
!             =  9 : bottom boundary at higher resolution 2 (only MPI)
!             = 10 : top    boundary at same   resolution
!             = 11 : top    boundary at higher resolution 1
!             = 12 : top    boundary at higher resolution 2 (only MPI)
!======================================================================
subroutine locatebounds(time)
  use parameters
  use globals
  implicit none
  real, intent(in) :: time
  logical, dimension(nbmin:nbmaxtot,4)       :: bounds
  integer :: nb, n1, n2, dest1, dest2, source, i
#ifdef MPIP
  integer, dimension(0:np-1) :: recvcnts,displs
  logical, dimension(nbmin:nbmax,8) :: boundext
  integer, dimension(:,:), allocatable :: extboundsp
  integer :: nexternalp, err
  !
  source=rank
  !
  !   not sure if needed
  boundext(:,:)=.false.
  !
#endif
  !   reset boundary flag, lists and their size
  bounds(:,:)=.false.
  deallocate(innerbounds)
  allocate(innerbounds(nbleafs*4,3))
  ninternal=0
  !
#ifdef MPIP
  call mpi_barrier(mpi_comm_world,nb)
  !
  nexternalp=0
  allocate(extboundsp(nbleafs*4,5))
#endif
  !
  !   locate internal boundaries
  DO nb=nbmin,nbmax
     IF( leaf(nb) ) THEN
        !---------------- 
        !   left boundary
        !---------------- 
        if ( not(bounds(nb,1)) ) then
           !   if box boundary
           if (lp(nb,5).eq.-1) then
              bounds(nb,1)=.true.
              ninternal=ninternal+1
              innerbounds(ninternal,1)= nb
              innerbounds(ninternal,2)=  0
              innerbounds(ninternal,3)= -1
              !
              !   only if neighbor is at least at same level
           else if (lp(lp(nb,5),1).ge.lp(nb,1) ) then
              !
              n1=lp(nb,5)
              
              !
              ! if neighbor is at same resolution
              if ( leaf(n1) ) then
                 !
#ifdef MPIP
                 dest1=n1/nblocks
                 !   if boundary is in same processor
                 if (dest1.eq.rank) then
#endif
                    ninternal=ninternal+1
                    innerbounds(ninternal,1)= nb
                    innerbounds(ninternal,2)= n1
                    innerbounds(ninternal,3)=  1
                    bounds(nb,1)=.true.
                    bounds(n1,2)=.true.
#ifdef MPIP
                 else
                    !   if boundary in other processor
                    !   same resolution to the left
                    boundext(nb,1) = .true.
                    if (nb .le. n1) then
                       nexternalp=nexternalp+1
                       extboundsp(nexternalp,1)= nb
                       extboundsp(nexternalp,2)= rank
                       extboundsp(nexternalp,3)= n1
                       extboundsp(nexternalp,4)= dest1
                       extboundsp(nexternalp,5)= 1
                    end if
                 end if
#endif
                 !   if neighbor it at higher resolution
              else
                 bounds(nb,1)=.true.
#ifndef MPIP
                 !   if MPIP is not defined
                 !***************************************************
                 n1=lp(n1,3)+1
                 n2=n1+2
                 !
                 ninternal=ninternal+1
                 innerbounds(ninternal,1)= nb
                 innerbounds(ninternal,2)= n1
                 innerbounds(ninternal,3)=  2
                 bounds(n1,2)=.true.
                 bounds(n2,2)=.true.
                 !
#else
                 !   if MPIP is defined
                 !***************************************************
                 n1=lp(lp(n1,3),6)
                 n2=lp(n1      ,8)
                 !   look in which processor the neighbors are
                 dest1=n1/nblocks
                 dest2=n2/nblocks
                 !
                 bounds(nb,1)=.true.
                 !
                 if (dest1.ne.rank) then
                    nexternalp=nexternalp+1
                    extboundsp(nexternalp,1)= nb
                    extboundsp(nexternalp,2)= rank
                    extboundsp(nexternalp,3)= n1
                    extboundsp(nexternalp,4)= dest1
                    extboundsp(nexternalp,5)= 2
                    boundext(nb,1) = .true.
                 else
                    ninternal=ninternal+1
                    innerbounds(ninternal,1)= nb
                    innerbounds(ninternal,2)= n1
                    innerbounds(ninternal,3)=  2
                 end if
                 !
                 if (dest2.ne.rank) then
                    nexternalp=nexternalp+1
                    extboundsp(nexternalp,1)= nb
                    extboundsp(nexternalp,2)= rank
                    extboundsp(nexternalp,3)= n2
                    extboundsp(nexternalp,4)= dest2
                    extboundsp(nexternalp,5)= 3
                    boundext(nb,2) = .true.
                 else
                    ninternal=ninternal+1
                    innerbounds(ninternal,1)= nb
                    innerbounds(ninternal,2)= n2
                    innerbounds(ninternal,3)=  3
                 end if
                 !
#endif
              end if
           end if
        end if
        !---------------- 
        !   right boundary
        !---------------- 
        if ( not(bounds(nb,2)) ) then
           !   if box boundary
           if (lp(nb,6).eq.-1) then
              bounds(nb,2)=.true.
              ninternal=ninternal+1
              innerbounds(ninternal,1)= nb
              innerbounds(ninternal,2)=  0
              innerbounds(ninternal,3)= -2
              !
              !   only if neighbor is at least at same level
           else if (lp(lp(nb,6),1).ge.lp(nb,1) ) then
              !
              n1=lp(nb,6)
              !
              ! if neighbor is at same resolution
              if ( leaf(n1) ) then
                 !
#ifdef MPIP
                 dest1=n1/nblocks
                 !   if boundary is in same processor
                 if (dest1.eq.rank) then
#endif
                    ninternal=ninternal+1
                    innerbounds(ninternal,1)= nb
                    innerbounds(ninternal,2)= n1
                    innerbounds(ninternal,3)=  4
                    bounds(nb,2)=.true.
                    bounds(n1,1)=.true.
#ifdef MPIP
                 else
                    !   if boundary in other processor
                    !   same resolution to the right
                    boundext(nb,3) = .true.
                    if (nb .le. n1) then
                       nexternalp=nexternalp+1
                       extboundsp(nexternalp,1)= nb
                       extboundsp(nexternalp,2)= rank
                       extboundsp(nexternalp,3)= n1
                       extboundsp(nexternalp,4)= dest1
                       extboundsp(nexternalp,5)= 4
                    end if
                 end if
#endif
                 !   if neighbor it at higher resolution
              else
                 bounds(nb,2)=.true.
#ifndef MPIP
                 !   if MPIP is not defined
                 !***************************************************
                 n1=lp(n1,3)
                 n2=n1+2
                 !
                 ninternal=ninternal+1
                 innerbounds(ninternal,1)= nb
                 innerbounds(ninternal,2)= n1
                 innerbounds(ninternal,3)=  5
                 bounds(n1,1)=.true.
                 bounds(n2,1)=.true.
                 !
#else
                 !   if MPIP is defined
                 !***************************************************
                 n1=lp(n1,3)
                 n2=lp(n1,8)
                 !   look in which processor the neighbors are
                 dest1=n1/nblocks
                 dest2=n2/nblocks
                 !
                 bounds(nb,2)=.true.
                 !
                 if (dest1.ne.rank) then
                    nexternalp=nexternalp+1
                    extboundsp(nexternalp,1)= nb
                    extboundsp(nexternalp,2)= rank
                    extboundsp(nexternalp,3)= n1
                    extboundsp(nexternalp,4)= dest1
                    extboundsp(nexternalp,5)= 5
                    boundext(nb,3) = .true.
                 else
                    ninternal=ninternal+1
                    innerbounds(ninternal,1)= nb
                    innerbounds(ninternal,2)= n1
                    innerbounds(ninternal,3)=  5

                 end if
                 !
                 if (dest2.ne.rank) then
                    nexternalp=nexternalp+1
                    extboundsp(nexternalp,1)= nb
                    extboundsp(nexternalp,2)= rank
                    extboundsp(nexternalp,3)= n2
                    extboundsp(nexternalp,4)= dest2
                    extboundsp(nexternalp,5)= 6
                    boundext(nb,4) = .true.
                 else
                    ninternal=ninternal+1
                    innerbounds(ninternal,1)= nb
                    innerbounds(ninternal,2)= n2
                    innerbounds(ninternal,3)=  6
                 end if
                 !
#endif
              end if
           end if
        end if
        !----------------
        !   bottom boundary
        !---------------- 
        if ( not(bounds(nb,3)) ) then
           !   if box boundary
           if (lp(nb,7).eq.-1) then
              bounds(nb,3)=.true.
              ninternal=ninternal+1
              innerbounds(ninternal,1)= nb
              innerbounds(ninternal,2)=  0
              innerbounds(ninternal,3)= -3
              !
              !   only if neighbor is at least at same level
           else if (lp(lp(nb,7),1).ge.lp(nb,1) ) then
              !
              n1=lp(nb,7)
              !
              ! if neighbor is at same resolution
              if ( leaf(n1) ) then
                 !
#ifdef MPIP
                 dest1=n1/nblocks
                 !   if boundary is in same processor
                 if (dest1.eq.rank) then
#endif
                    ninternal=ninternal+1
                    innerbounds(ninternal,1)= nb
                    innerbounds(ninternal,2)= n1
                    innerbounds(ninternal,3)=  7
                    bounds(nb,3)=.true.
                    bounds(n1,4)=.true.
#ifdef MPIP
                 else
                    !   if boundary in other processor
                    !   same resolution to the left
                    boundext(nb,5) = .true.
                    if (nb .le. n1) then
                       nexternalp=nexternalp+1
                       extboundsp(nexternalp,1)= nb
                       extboundsp(nexternalp,2)= rank
                       extboundsp(nexternalp,3)= n1
                       extboundsp(nexternalp,4)= dest1
                       extboundsp(nexternalp,5)= 7
                    end if
                 end if
#endif
                 !   if neighbor it at higher resolution
              else
                 bounds(nb,3)=.true.
#ifndef MPIP
                 !   if MPIP is not defined
                 !***************************************************
                 n1=lp(n1,3)+2
                 n2=n1+1
                 !
                 ninternal=ninternal+1
                 innerbounds(ninternal,1)= nb
                 innerbounds(ninternal,2)= n1
                 innerbounds(ninternal,3)=  8
                 bounds(n1,4)=.true.
                 bounds(n2,4)=.true.
                 !
#else
                 !   if MPIP is defined
                 !***************************************************
                 n1=lp(lp(n1,3),8)
                 n2=lp(n1      ,6)
                 !   look in which processor the neighbors are
                 dest1=n1/nblocks
                 dest2=n2/nblocks
                 !
                 bounds(nb,3)=.true.
                 !
                 if (dest1.ne.rank) then
                    nexternalp=nexternalp+1
                    extboundsp(nexternalp,1)= nb
                    extboundsp(nexternalp,2)= rank
                    extboundsp(nexternalp,3)= n1
                    extboundsp(nexternalp,4)= dest1
                    extboundsp(nexternalp,5)= 8
                    boundext(nb,5) = .true.
                 else
                    ninternal=ninternal+1
                    innerbounds(ninternal,1)= nb
                    innerbounds(ninternal,2)= n1
                    innerbounds(ninternal,3)=  8
                 end if
                 !
                 if (dest2.ne.rank) then
                    nexternalp=nexternalp+1
                    extboundsp(nexternalp,1)= nb
                    extboundsp(nexternalp,2)= rank
                    extboundsp(nexternalp,3)= n2
                    extboundsp(nexternalp,4)= dest2
                    extboundsp(nexternalp,5)= 9
                    boundext(nb,6) = .true.
                 else
                    ninternal=ninternal+1
                    innerbounds(ninternal,1)= nb
                    innerbounds(ninternal,2)= n2
                    innerbounds(ninternal,3)=  9
                 end if
                 !
#endif
              end if
           end if
        end if
        !---------------- 
        !   top boundary
        !---------------- 
        if ( not(bounds(nb,4)) ) then
           !   if box boundary
           if (lp(nb,8).eq.-1) then
              bounds(nb,4)=.true.
              ninternal=ninternal+1
              innerbounds(ninternal,1)= nb
              innerbounds(ninternal,2)= 0
              innerbounds(ninternal,3)=-4
              !
              !   only if neighbor is at least at same level
           else if (lp(lp(nb,8),1).ge.lp(nb,1) ) then
              !
              n1=lp(nb,8)
              !
              ! if neighbor is at same resolution
              if ( leaf(n1) ) then
                 !
#ifdef MPIP
                 dest1=n1/nblocks
                 !   if boundary is in same processor
                 if (dest1.eq.rank) then
#endif
                    ninternal=ninternal+1
                    innerbounds(ninternal,1)= nb
                    innerbounds(ninternal,2)= n1
                    innerbounds(ninternal,3)= 10
                    bounds(nb,4)=.true.
                    bounds(n1,3)=.true.
#ifdef MPIP
                 else
                    !   if boundary in other processor
                    !   same resolution to the right
                    boundext(nb,7) = .true.
                    if (nb .le. n1) then
                       nexternalp=nexternalp+1
                       extboundsp(nexternalp,1)= nb
                       extboundsp(nexternalp,2)= rank
                       extboundsp(nexternalp,3)= n1
                       extboundsp(nexternalp,4)= dest1
                       extboundsp(nexternalp,5)= 10
                    end if
                 end if
#endif
                 !   if neighbor it at higher resolution
              else
                 bounds(nb,4)=.true.
#ifndef MPIP
                 !   if MPIP is not defined
                 !***************************************************
                 n1=lp(n1,3)
                 n2=n1+1
                 !
                 ninternal=ninternal+1
                 innerbounds(ninternal,1)= nb
                 innerbounds(ninternal,2)= n1
                 innerbounds(ninternal,3)= 11
                 bounds(n1,3)=.true.
                 bounds(n2,3)=.true.
                 !
#else
                 !   if MPIP is defined
                 !***************************************************
                 n1=lp(n1,3)
                 n2=lp(n1,6)
                 !   look in which processor the neighbors are
                 dest1=n1/nblocks
                 dest2=n2/nblocks
                 !
                 bounds(nb,4)=.true.
                 !
                 if (dest1.ne.rank) then
                    nexternalp=nexternalp+1
                    extboundsp(nexternalp,1)= nb
                    extboundsp(nexternalp,2)= rank
                    extboundsp(nexternalp,3)= n1
                    extboundsp(nexternalp,4)= dest1
                    extboundsp(nexternalp,5)= 11
                    boundext(nb,7) = .true.
                 else
                    ninternal=ninternal+1
                    innerbounds(ninternal,1)= nb
                    innerbounds(ninternal,2)= n1
                    innerbounds(ninternal,3)= 11
                 end if
                 !
                 if (dest2.ne.rank) then
                    nexternalp=nexternalp+1
                    extboundsp(nexternalp,1)= nb
                    extboundsp(nexternalp,2)= rank
                    extboundsp(nexternalp,3)= n2
                    extboundsp(nexternalp,4)= dest2
                    extboundsp(nexternalp,5)= 12
                    boundext(nb,8) = .true.
                 else
                    ninternal=ninternal+1
                    innerbounds(ninternal,1)= nb
                    innerbounds(ninternal,2)= n2
                    innerbounds(ninternal,3)= 12
                 end if
                 !
#endif
              end if
           end if
        end if
        !----------------
     END IF
  END DO
  !--------------------------------------------------------------------
  !
  !consolidate list of external boundaries
  !
#ifdef MPIP
  call mpi_allgather(nexternalp, 1, mpi_integer, recvcnts, 1, mpi_integer, &
       mpi_comm_world,err)
  !
  nexternal=sum(recvcnts(:))
  deallocate(extbounds)
  allocate(extbounds(nexternal,5))
  !
  displs(0)=0
  do i=1,np-1
     displs(i)=displs(i-1)+recvcnts(i-1)
  end do
  !
  do i=1,5
     call mpi_allgatherv(extboundsp(:,i) , nexternalp, mpi_integer, &
          extbounds (:,i)  , recvcnts    , displs,  &
          mpi_integer, mpi_comm_world,err)
  end do
  ! 
  deallocate(extboundsp)
  !
  !leafmpi(:)=.false.
  !do nb=1,nexternal
  !   if (rank.eq.extbounds(nb,2) ) leafmpi(extbounds(nb,1))=.true.
  !   if (rank.eq.extbounds(nb,4) ) leafmpi(extbounds(nb,3))=.true.
  !end do
  !do i=0,np-1
  !   if (rank.eq.i) then
  !      do nb=1,nexternal
  !         print'(5i5)',extbounds(nb,:)
  !      end do
  !   end if
  !   call mpi_barrier(mpi_comm_world,err)
  !end do
  !
#endif
  !--------------------------------------------------------------------
  !
  return
end subroutine locatebounds
!======================================================================
