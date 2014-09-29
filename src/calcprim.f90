!======================================================================
!   calculates the primitives on all the domain 
!   nb is the block ID
!======================================================================
subroutine calcprim(u1,prim1,neq,nxmin,nxmax,nymin,nymax)
  implicit none
  real,    intent(in),  dimension(neq,nxmin:nxmax,nymin:nymax) :: u1
  real,    intent(out), dimension(neq,nxmin:nxmax,nymin:nymax) :: prim1
  integer, intent(in) :: neq,nxmin,nxmax,nymin,nymax
  real    :: T
  integer :: i, j
  !
  do i=nxmin,nxmax
     do j=nymin,nymax
        !
        call uprim(prim1(:,i,j),u1(:,i,j),T)
        !
     end do
  end do
  !
end subroutine calcprim
!======================================================================

