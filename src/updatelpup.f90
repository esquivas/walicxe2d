!=======================================================================
!  update pointers after refining
!=======================================================================
subroutine updatelpup(nb, son1, son2, son3, son4, irefup, irefdown)
  use parameters
  use globals
  implicit none
  integer, intent(in) :: nb, son1, son2, son3, son4
  logical, intent(out), dimension(nblocks*np) :: irefup,irefdown
  integer :: nsr, nlev2
  !
  !  nlev2: level of refinement of the sons
  nlev2=lp(nb,1)+1
  !
#ifdef MPIP
  if (rank.eq.  nb/nblocks) nbleafs=nbleafs-1
  !
  if (rank.eq.son1/nblocks) nbleafs=nbleafs+1
  if (rank.eq.son2/nblocks) nbleafs=nbleafs+1
  if (rank.eq.son3/nblocks) nbleafs=nbleafs+1
  if (rank.eq.son4/nblocks) nbleafs=nbleafs+1
#endif
  !
  !---------------------------------------------------------------------
  !  son1
  lp(son1,1)=nlev2          ! level of refinement
  lp(son1,2)=nb             ! father
  !lp(son1,3)=0              ! son (only 1st=ID1)
  lp(son1,4)=1              ! my ID
  lp(son1,6)=son2           ! right
  lp(son1,8)=son3           ! up
  ! left  ***
  if (lp(nb,5).eq.-1.or.leaf(lp(nb,5))) then
     lp(son1,5)=lp(nb,5)
  else
#ifndef MPIP
     nsr=lp(lp(nb,5),3)+1
#else
     nsr=lp(lp(lp(nb,5),3),6)
#endif
     lp(son1,5)=nsr
     lp(nsr ,6)=son1
     !if (leafmpi(nb)) then
     !   if (rank.ne.nsr/nblocks) print*,'MPI LEAF LEFT',nb,nsr,son1
     !end if
  end if
  ! down  ****
  if (lp(nb,7).eq.-1.or.leaf(lp(nb,7))) then
     lp(son1,7)=lp(nb,7)
  else
#ifndef MPIP
     nsr=lp(lp(nb,7),3)+2
#else
     nsr=lp(lp(lp(nb,7),3),8)
#endif
     lp(son1,7)=nsr
     lp(nsr ,8)=son1
     !if (leafmpi(nb)) then
     !   if (rank.ne.nsr/nblocks) print*,'MPI LEAF DOWN',nb,nsr,son1
     !end if
  endif
  !---------------------------------------------------------------------
  !  son2
  lp(son2,1)=nlev2          ! level of refinement
  lp(son2,2)=nb             ! father
  !lp(son2,3)=0              ! son (only 1st=ID1)
  lp(son2,4)=2              ! my ID
  lp(son2,5)=son1           ! left
  lp(son2,8)=son4           ! up
  ! right ****
  if (lp(nb,6).eq.-1.or.leaf(lp(nb,6))) then
     lp(son2,6)=lp(nb,6)
  else
     nsr=lp(lp(nb,6),3)
     lp(son2,6)=nsr
     lp(nsr ,5)=son2
     !if (leafmpi(nb)) then
     !   if (rank.ne.nsr/nblocks)  print*,'MPI LEAF RIGHT',nb,nsr,son2
     !end if
  end if
  ! down  ****
  if (lp(nb,7).eq.-1.or.leaf(lp(nb,7))) then
     lp(son2,7)=lp(nb,7)
  else
#ifndef MPIP
     nsr=lp(lp(nb,7),3)+3
#else
     nsr=lp(lp(lp(lp(nb,7),3),8),6)
#endif
     lp(son2,7)=nsr
     lp(nsr ,8)=son2
     !if (leafmpi(nb)) then
     !   if (rank.ne.nsr/nblocks) print*,'MPI LEAF DOWN',nb,nsr,son2
     !end if
  end if
  !---------------------------------------------------------------------
  !  son3
  lp(son3,1)=nlev2          ! level of refinement
  lp(son3,2)=nb             ! father
  !lp(son3,3)=0              ! son (only 1st=ID1)
  lp(son3,4)=3              ! my ID
  lp(son3,6)=son4           ! right
  lp(son3,7)=son1           ! down
  ! left  ****
  if (lp(nb,5).eq.-1.or.leaf(lp(nb,5))) then
     lp(son3,5)=lp(nb,5)
  else
#ifndef MPIP
     nsr=lp(lp(nb,5),3)+3
#else
     nsr=lp(lp(lp(lp(nb,5),3),8),6)
#endif
     lp(son3,5)=nsr
     lp(nsr ,6)=son3
     !if (leafmpi(nb)) then
     !   if (rank.ne.nsr/nblocks) print*,'MPI LEAF LEFT',nb,nsr,son3
     !end if
  end if
  ! up    ****
  if (lp(nb,8).eq.-1.or.leaf(lp(nb,8))) then
     lp(son3,8)=lp(nb,8)
  else
     nsr=lp(lp(nb,8),3)
     lp(son3,8)=nsr
     lp(nsr ,7)=son3
     !if (leafmpi(nb)) then
     !   if (rank.ne.nsr/nblocks) print*,'MPI LEAF UP',nb,nsr,son3
     !end if
  end if
  !---------------------------------------------------------------------
  !  son4
  lp(son4,1 )=nlev2         ! level of refinement
  lp(son4,2)=nb             ! father
  !lp(son4,3)=0              ! son (only 1st=ID1)
  lp(son4,4)=4              ! my ID
  lp(son4,5)=son3           ! left
  lp(son4,7)=son2           ! down
  ! right ****
  if (lp(nb,6).eq.-1.or.leaf(lp(nb,6))) then
     lp(son4,6)=lp(nb,6)
  else
#ifndef MPIP
     nsr=lp(lp(nb,6),3)+2
#else
     nsr=lp(lp(lp(nb,6),3),8)
#endif
     lp(son4,6)=nsr
     lp(nsr ,5)=son4
     !if (leafmpi(nb)) then
     !   if (rank.ne.nsr/nblocks) print*,'MPI LEAF RIGHT',nb,nsr,son4
     !end if
  end if
  ! up    ****
  if (lp(nb,8).eq.-1.or.leaf(lp(nb,8))) then
     lp(son4,8)=lp(nb,8)
  else
#ifndef MPIP
     nsr=lp(lp(nb,8),3)+1
#else
     nsr=lp(lp(lp(nb,8),3),6)
#endif
     lp(son4,8)=nsr
     lp(nsr ,7)=son4
     !if (leafmpi(nb)) then
     !   if (rank.ne.nsr/nblocks) print*,'MPI LEAF UP',nb,nsr,son4
     !end if
  end if
  !---------------------------------------------------------------------
  ! NOTE that the neighbors marked **** are given on a coarser level
  ! if not available at same level
  !---------------------------------------------------------------------

    !  father no longer leaf
  leaf(nb)=.false.
  !  all sons are leafs now
  leaf(son1)=.true.
  leaf(son2)=.true.
  leaf(son3)=.true.
  leaf(son4)=.true.
  !
  !   de-mark for refining/coarsening
  irefup(nb)=.false.
  irefdown(nb)=.false.
  !
  irefdown(son1)=.false.
  irefdown(son2)=.false.
  irefdown(son3)=.false.
  irefdown(son4)=.false.
  !
end subroutine updatelpup
!=======================================================================
