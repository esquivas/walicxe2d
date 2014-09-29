!=======================================================================
!   This module contains all global variables
!=======================================================================
module globals
  implicit none
  !
  !   ID range of blocks on each proc
  integer :: nbmin, nbmaxtot
  !   number of active blocks and leafs   
  integer :: nbmax, nbleafs
  !    
  !   flow variables and fluxes
  real,    dimension(:,:,:,:), allocatable :: u, up, primit
  !
  !   pointers, refining, boundary  and coordinate info
  integer, dimension(:,:), allocatable :: lp
  logical, dimension( : ), allocatable :: leaf
  integer, dimension(:,:), allocatable :: icoord
  !
  !   spacing
  real,    dimension( : ), allocatable :: dx,dy
  !
  !
  !   arrays with boundaries info
  integer, dimension(:,:), allocatable :: innerbounds
  integer :: ninternal
#ifdef MPIP
  integer, dimension(:,:), allocatable :: extbounds
  integer :: nexternal
#endif
  !
  !   MPI rank and # of procs
  integer :: rank, np
  !
end module globals
!=======================================================================
