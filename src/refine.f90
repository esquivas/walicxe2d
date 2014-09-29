!=======================================================================
!  refine to next level
!=======================================================================
subroutine refine(nb,son1,son2,son3,son4)
  use parameters
  use globals
  implicit none
  integer, intent(in) :: nb,son1,son2,son3,son4
  integer :: i, j, nx0, ny0, nsr
  integer, parameter :: nx2=nx/2, ny2=ny/2
#ifdef MPIP
  integer :: err, status(MPI_STATUS_SIZE)
  integer :: source,dest1,dest2,dest3,dest4
  integer, dimension(2) :: nxy0
  integer, parameter :: scount=(neq*(nx2+2)*(ny2+2))
  real, dimension(:,:,:,:), allocatable :: recv
  source= nb/nblocks
  dest1 = son1/nblocks
  dest2 = son2/nblocks
  dest3 = son3/nblocks
  dest4 = son4/nblocks
#endif
  !
  !if(rank.eq.source) print'(a,i6,a,4i6)','refining   block  '       &
  !     , nb,' ->',son1,son2,son3,son4
  !
  !---------------------------------------------------------------------
  !
  !   update the u's
#ifndef MPIP
  do i=nxmin,nxmax
     do j=nymin,nymax
        u(son1,:,i,j)=u( nb,:, (i+1)/2    ,(j+1)/2     )
        u(son2,:,i,j)=u( nb,:, (i+1)/2+nx2,(j+1)/2     )
        u(son3,:,i,j)=u( nb,:, (i+1)/2    ,(j+1)/2+ny2 )
        u(son4,:,i,j)=u( nb,:, (i+1)/2+nx2,(j+1)/2+ny2 )
     end do
  end do
  !
  !   for the icoords
  nx0=icoord(nb,1)*2!(2**(nlevs2+1))
  ny0=icoord(nb,2)*2!(2**(nlevs2+1))
  !
#else
  !
  !   for the icoords
  if (rank.eq.source) then
     nx0=icoord(nb,1)*2!(2**(nlevs2+1))
     ny0=icoord(nb,2)*2!(2**(nlevs2+1))
     nxy0(1)=nx0
     nxy0(2)=ny0
  end if
  !
  !
  !   do blocks one by one
  !--------------------
  !   son1
  if (dest1.eq.source) then
     if (rank.eq.source) then
        do i=nxmin,nxmax
           do j=nymin,nymax
              u(son1,:,i,j)=u( nb,:, (i+1)/2     ,(j+1)/2     )
           end do
        end do
     end if
  else
     !   refine w/MPI
     !   send buffer
     if (rank.eq.source) then
        call mpi_send(   u(nb,:,0:nx2+1,0:ny2+1) ,scount,mpi_real_kind,& 
             dest1,0,mpi_comm_world, err)
        call mpi_send(nxy0, 2, mpi_integer,dest1,0,mpi_comm_world, err) 
     end if
     !   receive
     if (rank.eq.dest1) then
        !
        allocate(recv(1,neq,0:nx2+1,0:ny2+1))
        !
        call mpi_recv(recv(1 ,:,0:nx2+1,0:ny2+1), scount,mpi_real_kind, &
             source,0,mpi_comm_world, status, err )
        call mpi_recv(nxy0, 2 ,mpi_integer, &
             source,0,mpi_comm_world, status, err )
        nx0=nxy0(1)
        ny0=nxy0(2)
        !   refine
        do i=nxmin,nxmax
           do j=nymin,nymax
              u(son1,:,i,j)=recv( 1,:, (i+1)/2     ,(j+1)/2     )
           end do
        end do
        !
        deallocate(recv)
     end if
     !
     !print*,'son1: ', son1,' from ', source, 'on ',rank
     call mpi_barrier(mpi_comm_world,err)
  end if
  !--------------------
  !   son2
  if ( dest2.eq.source) then
     if (rank.eq.source) then
        do i=nxmin,nxmax
           do j=nymin,nymax
              u(son2,:,i,j)=u( nb,:, (i+1)/2+nx2,(j+1)/2     )
           end do
        end do
     end if
  else
     !   refine w/MPI
     !   send buffer
     if (rank.eq.source) then
        call mpi_send(   u(nb,:,nx2:nx+1,0:ny2+1) , scount, mpi_real_kind, & 
             dest2,0,mpi_comm_world, err)
        call mpi_send(nxy0, 2, mpi_integer,dest2,0,mpi_comm_world, err) 
     end if
     !   receive
     if (rank.eq.dest2) then
        allocate(recv(1,neq,nx2:nx+1,0:ny2+1))
        !
        call mpi_recv(recv, scount, mpi_real_kind, &
             source,0,mpi_comm_world, status, err )
        call mpi_recv(nxy0, 2 ,mpi_integer, &
             source,0,mpi_comm_world, status, err )
        nx0=nxy0(1)
        ny0=nxy0(2)
        !   refine
        do i=nxmin,nxmax
           do j=nymin,nymax
              u(son2,:,i,j)=recv( 1,:, (i+1)/2+nx2,(j+1)/2     )
           end do
        end do
        !
        deallocate(recv)
     end if
     !
     !print*,'son2: ', son2,' from ', source, 'on ',rank
     call mpi_barrier(mpi_comm_world,err)
  end if
  !--------------------
  !   son3
  if (dest3.eq.source) then
     if (rank.eq.source) then
        do i=nxmin,nxmax
           do j=nymin,nymax
              u(son3,:,i,j)=u( nb,:, (i+1)/2     ,(j+1)/2+ny2)
           end do
        end do
     end if
  else
     !   refine w/MPI
     !   send buffer
     if (rank.eq.source) then
        call mpi_send(   u(nb,:,0:nx2+1,ny2-1:ny) ,scount,mpi_real_kind, & 
             dest3,0,mpi_comm_world, err)
        call mpi_send(nxy0, 2, mpi_integer,dest3,0,mpi_comm_world, err) 
     end if
     !   receive
     if (rank.eq.dest3) then
        !
        allocate(recv(1,neq,0:nx2+1,ny2:ny+1))
        !
        call mpi_recv(recv(1 ,:,0:nx2+1,ny2:ny+1), scount,mpi_real_kind, &
             source,0,mpi_comm_world, status, err )
        call mpi_recv(nxy0, 2 ,mpi_integer, &
             source,0,mpi_comm_world, status, err )
        nx0=nxy0(1)
        ny0=nxy0(2)
        !   refine
        do i=nxmin,nxmax
           do j=nymin,nymax
              u(son3,:,i,j)=recv( 1,:, (i+1)/2     ,(j+1)/2+ny2)
           end do
        end do
        !
        deallocate(recv)
     end if
     !
     !print*,'son3: ', son3,' from ', source, 'on ',rank
     call mpi_barrier(mpi_comm_world,err)
  end if
  !--------------------
  !   son4
  if (dest4.eq.source) then
     if (rank.eq.source) then
        do i=nxmin,nxmax
           do j=nymin,nymax
              u(son4,:,i,j)=u( nb,:, (i+1)/2+nx2,(j+1)/2+ny2) 
           end do
        end do
     end if
  else
     !   refine w/MPI
     !   send buffer
     if (rank.eq.source) then
        call mpi_send(   u(nb,:,nx2:nx+1,ny2:ny+1) ,scount,mpi_real_kind, & 
             dest4,0,mpi_comm_world, err)
        call mpi_send(nxy0, 2, mpi_integer,dest4,0,mpi_comm_world, err) 
     end if
     !   receive
     if (rank.eq.dest4) then

        !
        allocate(recv(1,neq,nx2:nx+1,ny2:ny+1))
        !
        call mpi_recv(recv(1 ,:,nx2:nx+1,ny2:ny+1), scount,mpi_real_kind, &
             source,0,mpi_comm_world, status, err )
        call mpi_recv(nxy0, 2 ,mpi_integer, &
             source,0,mpi_comm_world, status, err )
        nx0=nxy0(1)
        ny0=nxy0(2)
        !   refine
        do i=nxmin,nxmax
           do j=nymin,nymax
              u(son4,:,i,j)=recv( 1,:, (i+1)/2+nx2,(j+1)/2+ny2) 
           end do
        end do
        !
        deallocate(recv)
     end if
     !
     !print*,'son4: ', son4,' from ', source, 'on ',rank
     call mpi_barrier(mpi_comm_world,err)
  end if
  !--------------------
  !
#endif
  !
  !   update the icoords
#ifdef MPIP
  if (rank.eq.dest1) then
#endif
     call calcprim(u(son1,:,:,:), primit(son1,:,:,:),   &
                   neq, nxmin, nxmax, nymin, nymax )
     icoord(son1,1)=nx0
     icoord(son1,2)=ny0
#ifdef MPIP
  end if
  if (rank.eq.dest2) then
#endif
     call calcprim(u(son2,:,:,:), primit(son2,:,:,:),   &
                   neq, nxmin, nxmax, nymin, nymax )
     icoord(son2,1)=nx0+nx
     icoord(son2,2)=ny0
#ifdef MPIP
  end if
  if (rank.eq.dest3) then
#endif
     call calcprim(u(son3,:,:,:), primit(son3,:,:,:),   &
                   neq, nxmin, nxmax, nymin, nymax )
     icoord(son3,1)=nx0
     icoord(son3,2)=ny0+ny
#ifdef MPIP
  end if
  if (rank.eq.dest4) then
#endif
     call calcprim(u(son4,:,:,:), primit(son4,:,:,:),   &
                   neq, nxmin, nxmax, nymin, nymax )
     icoord(son4,1)=nx0+nx
     icoord(son4,2)=ny0+ny
#ifdef MPIP
  end if
#endif
  !
  !--------------------------------------------------------------------
  !
end subroutine refine
!=======================================================================


