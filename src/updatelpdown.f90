!=======================================================================
!  update pointers after coarsening
!=======================================================================
subroutine updatelpdown(dad, son1, son2, son3, son4)
  use parameters
  use globals
  implicit none
  integer, intent(in) :: dad, son1, son2, son3, son4
  integer :: nb1
  !
  !   update the father
  leaf(dad)=.true.
  !
  !   clear the sons:
  leaf(son1)=.false.
  leaf(son2)=.false.
  leaf(son3)=.false.
  leaf(son4)=.false.
  !
#ifdef MPIP
  if (rank.eq. dad/nblocks) nbleafs=nbleafs+1
  !
  if (rank.eq.son1/nblocks) nbleafs=nbleafs-1
  if (rank.eq.son2/nblocks) nbleafs=nbleafs-1
  if (rank.eq.son3/nblocks) nbleafs=nbleafs-1
  if (rank.eq.son4/nblocks) nbleafs=nbleafs-1
#endif
  !
  !  update neighbors
  !  left
  nb1=lp(dad,5)
  if (nb1.ne.-1 .and. not(leaf(nb1)) ) then
     nb1=lp(lp(nb1,3),6)
     lp(   nb1   ,6)=dad
     lp(lp(nb1,8),6)=dad
  end if
  !  right
  nb1=lp(dad,6)
  if (nb1.ne.-1 .and. not(leaf(nb1)) ) then
     nb1=lp(nb1,3)
     lp(   nb1   ,5)=dad
     lp(lp(nb1,8),5)=dad
  end if
  !  bottom
  nb1=lp(dad,7)
  if (nb1.ne.-1 .and. not(leaf(nb1)) ) then
     nb1=lp(lp(nb1,3),8)
     lp(   nb1   ,8)=dad
     lp(lp(nb1,6),8)=dad
  end if
  !  top
  nb1=lp(dad,8)
  if ( nb1.ne.-1 .and. not(leaf(nb1)) ) then
     nb1=lp(nb1,3)
     lp(   nb1   ,7)=dad
     lp(lp(nb1,6),7)=dad
  end if
  !
end subroutine updatelpdown
!=======================================================================
