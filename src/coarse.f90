!=======================================================================
!  coarsen to parent level son1-4 --> dad
!  MPIP assumes that we have 2 ghost cells
!=======================================================================
subroutine coarse(dad,son1,son2,son3,son4)
  use parameters
  use globals
  implicit none
  integer, intent(in) :: dad,son1,son2,son3,son4
  integer :: i, j, i1, i2, j1, j2, ieq
  integer, parameter :: nx2=(nx/2), ny2=ny/2
  integer, parameter :: nxp1=nx+1,nyp1=ny+1,nx2p1=nx2+1,ny2p1=ny2+1
  integer            ::nx2i, ny2j
#ifdef MPIP
  integer:: source1,source2,source3,source4,dest
  real, dimension(1,neq,    0:nx2 ,    0:ny2 ) :: send1
  real, dimension(1,neq,nx2p1:nxp1,    0:ny2 ) :: send2
  real, dimension(1,neq,    0:nx2 ,ny2p1:nyp1) :: send3
  real, dimension(1,neq,nx2p1:nxp1,ny2p1:nyp1) :: send4
  integer, parameter :: scount=(nx2p1*ny2p1*neq)
  integer :: err, status(MPI_STATUS_SIZE)
  !integer, dimension(2) :: req
  !integer, dimension(MPI_STATUS_SIZE,2) :: status_array
#endif
  !
  !if (rank.eq.0) print'(a,2i6,a,4i6,3i6)','derefining blocks '       &
  !    ,rank,dad,' <-',son1,son2,son3,son4,nbmax,nbleafs
  
#ifndef MPIP
  !----------------------------------------------------------------
  !   the block itself
  do i=1,nx2
     do j  = 1,ny2
        i1 = 2*i-1
        i2 = i1+1
        j1 = 2*j-1
        j2 = j1+1
        nx2i = nx2+i
        ny2j = ny2+j
        
        do ieq=1,neq
           u(dad,ieq,i   ,j   ) = 0.25*sum( u(son1,ieq,i1:i2,j1:j2) )
           u(dad,ieq,nx2i,j   ) = 0.25*sum( u(son2,ieq,i1:i2,j1:j2) )
           u(dad,ieq,i   ,ny2j) = 0.25*sum( u(son3,ieq,i1:i2,j1:j2) )
           u(dad,ieq,nx2i,ny2j) = 0.25*sum( u(son4,ieq,i1:i2,j1:j2) )
        end do
        
     end do
  end do
  !
  !**************************************************************
  !   the ghost cells of the block
  !---------------------
  if (nghost.eq.1) then
     !   left
     do i=nxmin,0
        do j  = 1,ny2
           j1 = 2*j-1
           j2 = j1+1
           
           do ieq=1,neq
              u(dad,ieq,i ,j    ) = 0.5*sum( u(son1,ieq,i,j1:j2) )
              u(dad,ieq,i ,ny2+j) = 0.5*sum( u(son3,ieq,i,j1:j2) )
           end do
           
        end do
     end do
     !   right
     do i=nx+1,nxmax
        do j  = 1,ny2
           j1 = 2*j-1
           j2 = j1+1
           
           do ieq=1,neq
              u(dad,ieq,i ,j    ) = 0.5*sum( u(son2,ieq,i,j1:j2) )
              u(dad,ieq,i ,ny2+j) = 0.5*sum( u(son4,ieq,i,j1:j2) )
           end do
           
        end do
     end do
     !   bottom
     do i=1,nx2
        do j  = nymin,0
           i1 = 2*i-1
           i2 = i1+1
           
           do ieq=1,neq
              u(dad,ieq,i    ,j ) = 0.5*sum( u(son1,ieq,i1:i2,j) )
              u(dad,ieq,nx2+i,j ) = 0.5*sum( u(son2,ieq,i1:i2,j) )
           end do
           
        end do
     end do
     !   top
     do i=1,nx2
        do j  = ny+1,nymax
           i1 = 2*i-1
           i2 = i1+1
           
           do ieq=1,neq
              u(dad,ieq,i    ,j ) = 0.5*sum( u(son3,ieq,i1:i2,j) )
              u(dad,ieq,nx2+i,j ) = 0.5*sum( u(son4,ieq,i1:i2,j) )
           end do
           
        end do
     end do
     !---------------------
  else if (nghost.eq.2) then
     !   left
     do ieq=1,neq
        u(dad,ieq,0 ,0   ) =  0.25*sum( u(son1,ieq,-1:0,  -1:0    ) )
        u(dad,ieq,0 ,ny+1) =  0.25*sum( u(son3,ieq,-1:0,ny+1:nymax) )
     end do
     do i=nxmin,0
        do j  = 1,ny2
           j1 = 2*j-1
           j2 = j1+1
           
           do ieq=1,neq
              u(dad,ieq,i ,j    ) = 0.25*sum( u(son1,ieq,-1:0,j1:j2) )
              u(dad,ieq,i ,ny2+j) = 0.25*sum( u(son3,ieq,-1:0,j1:j2) )
           end do
           
        end do
     end do
     !   right
     do ieq=1,neq
        u(dad,ieq,nx+1 ,0   ) =  0.25*sum( u(son2,ieq,nx+1:nxmax,  -1:0    ) )
        u(dad,ieq,nx+1 ,ny+1) =  0.25*sum( u(son4,ieq,nx+1:nxmax,ny+1:nymax) ) 
     end do
     do i=nx+1,nxmax
        do j  = 1,ny2
           j1 = 2*j-1
           j2 = j1+1
           
           do ieq=1,neq
              u(dad,ieq,i ,j    ) = 0.25*sum( u(son2,ieq,nx+1:nxmax,j1:j2) )
              u(dad,ieq,i ,ny2+j) = 0.25*sum( u(son4,ieq,nx+1:nxmax,j1:j2) )
           end do
           
        end do
     end do
     !   bottom
     do ieq=1,neq
        u(dad,ieq, 0  ,0) = 0.25*sum( u(son1,ieq,  -1:0    ,-1:0) )
        u(dad,ieq,nx+1,0) = 0.25*sum( u(son2,ieq,nx+1:nxmax,-1:0) )
     end do
     do i=1,nx2
        do j  = nymin,0
           i1 = 2*i-1
           i2 = i1+1
           
           do ieq=1,neq
              u(dad,ieq,i    ,j ) = 0.25*sum( u(son1,ieq,i1:i2,-1:0) )
              u(dad,ieq,nx2+i,j ) = 0.25*sum( u(son2,ieq,i1:i2,-1:0) )
           end do
           
        end do
     end do
     !   top
     do ieq=1,neq
        u(dad,ieq, 0  ,ny+1) = 0.25*sum( u(son3,ieq,  -1:0    ,ny+1:nymax) )
        u(dad,ieq,nx+1,ny+1) = 0.25*sum( u(son4,ieq,nx+1:nxmax,ny+1:nymax) )
     end do
     do i=1,nx2
        do j  = ny+1,nymax
           i1 = 2*i-1
           i2 = i1+1
           
           do ieq=1,neq
              u(dad,ieq,i    ,j ) = 0.25*sum( u(son3,ieq,i1:i2,ny+1:nymax) )
              u(dad,ieq,nx2+i,j ) = 0.25*sum( u(son4,ieq,i1:i2,ny+1:nymax) )
           end do
           
        end do
     end do
     !---------------------
     call calcprim(u(dad,:,:,:), primit(dad,:,:,:),                   &
          neq, nxmin, nxmax, nymin, nymax )
  end if
  !----------------------------------------------------------------
#else
  !
  dest=dad/nblocks
  source1=son1/nblocks
  source2=son2/nblocks
  source3=son3/nblocks
  source4=son4/nblocks
  !----------------------------------------------------------------
  !   son1
  if (rank.eq.dest .or. rank.eq.source1) then
     if (source1.eq.dest) then
        !  both are in the same processor
        do i= 0,nx2
           do j =0,ny2
              do ieq=1,neq
                 u(dad,ieq,i,j) = &
                      0.25*sum( u(son1,ieq,2*i-1:2*i,2*j-1:2*j) )
              end do
           end do
        end do
     else
        !  different processor
        !print*,'coarsening w/MPI 1', dad,son1,rank
        if (rank.eq.source1) then
           !print*,'coarsening w/MPI 1', dad,son1
           do i= 0,nx2
              do j =0,ny2
                 do ieq=1,neq
                    send1(1,ieq,i,j) = &
                         0.25*sum( u(son1,ieq,2*i-1:2*i,2*j-1:2*j) )
                 end do
              end do
           end do
           call mpi_send(send1,scount,mpi_real_kind,dest,0,mpi_comm_world, err)
           !call mpi_Isend(send1,scount,mpi_real_kind,dest,0,mpi_comm_world,req(1), err)
           !call mpi_wait(req(1), status, err)
        end if
        if (rank.eq.dest) then
           call mpi_recv( u(dad,:,0:nx2,0:ny2)  &
                ,scount,mpi_real_kind,source1,0,mpi_comm_world, status, err)
           !call mpi_Irecv( u(dad,:,0:nx2,0:ny2)  &
           !     ,scount,mpi_real_kind,source1,0,mpi_comm_world, status, req(2), err)
           !call mpi_wait(req(2), status, err)
        end if
        !
        !call mpi_waitall( 2, req, status, err)
        !
     end if
  end if
  !print*,'done son1'
  !call mpi_barrier(mpi_comm_world,err)
  
  !   son2
  if (rank.eq.dest .or. rank.eq.source2) then
     if (source2.eq.dest) then
        !  both are in the same processor
        do i= 1,nx2+1
           do j =0,ny2
              do ieq=1,neq
                 u(dad,ieq,i+nx2,j) = &
                      0.25*sum( u(son2,ieq,2*i-1:2*i,2*j-1:2*j) )
              end do
           end do
        end do
     else
        !  different processor
        !print*,'coarsening w/MPI 2', dad,son2,rank
        if (rank.eq.source2) then
           !print*,'coarsening w/MPI 2', dad,son2
           do i= 1,nx2p1
              do j =0,ny2
                 do ieq=1,neq
                    send2(1,ieq,i+nx2,j) = &
                         0.25*sum( u(son2,ieq,2*i-1:2*i,2*j-1:2*j) )
                 end do
              end do
           end do
           call mpi_send(send2,scount,mpi_real_kind,dest,1,mpi_comm_world, err)

        end if
        if (rank.eq.dest) then
           call mpi_recv( u(dad,:,nx2p1:nxp1,0:ny2)  &
                ,scount,mpi_real_kind,source2,1,mpi_comm_world, status, err)

        end if
        !
     end if
  end if
  !print*,'done son2'
  !call mpi_barrier(mpi_comm_world,err)
  
  !   son3
  if (rank.eq.dest .or. rank.eq.source3) then
     if (source3.eq.dest) then
        !  both are in the same processor
        do i= 0,nx2
           do j =1,ny2+1
              do ieq=1,neq
                 u(dad,ieq,i,j+ny2) = &
                      0.25*sum( u(son3,ieq,2*i-1:2*i,2*j-1:2*j) )
              end do
           end do
        end do
     else
        !  different processor
        !print*,'coarsening w/MPI 3', dad,son3,rank
        if (rank.eq.source3) then
           !print*,'coarsening w/MPI 3', dad,son3
           do i= 0,nx2
              do j =1,ny2p1
                 do ieq=1,neq
                    send3(1,ieq,i,j+ny2) = &
                         0.25*sum( u(son3,ieq,2*i-1:2*i,2*j-1:2*j) )
                 end do
              end do
           end do
           call mpi_send(send3,scount,mpi_real_kind,dest,3,mpi_comm_world, err)
        end if
        if (rank.eq.dest) then
           call mpi_recv( u(dad,:,0:nx2,ny2p1:nyp1)  &
                ,scount,mpi_real_kind,source3,3,mpi_comm_world, status, err)
        end if
 
      end if
  end if
  !print*,'done son3'
  !call mpi_barrier(mpi_comm_world,err)

  !   son4
  if (rank.eq.dest .or. rank.eq.source4) then
     if (source4.eq.dest) then
        !  both are in the same processor
        do i= 1,nx2+1
           do j =1,ny2p1
              do ieq=1,neq
                 u(dad,ieq,i+nx2,j+ny2) = &
                      0.25*sum( u(son4,ieq,2*i-1:2*i,2*j-1:2*j) )
              end do
           end do
        end do
     else
        !  they are on a different processor
        !print*,'coarsening w/MPI 4', dad,son4,rank
        if (rank.eq.source4) then
           !print*,'coarsening w/MPI 4', dad,son4
           do i= 1,nx2+1
              do j =1,ny2+1
                 do ieq=1,neq
                    send4(1,ieq,i+nx2,j+ny2) = &
                         0.25*sum( u(son4,ieq,2*i-1:2*i,2*j-1:2*j) )
                 end do
              end do
           end do
           call mpi_send(send4,scount,mpi_real_kind,dest,4,mpi_comm_world, err)
        end if
        if (rank.eq.dest) then
           call mpi_recv( u(dad,:,nx2p1:nxp1,ny2p1:nyp1)  &
                ,scount,mpi_real_kind,source4,4,mpi_comm_world, status, err)
        end if

     end if
  end if
  !print*,'done son4'
  !call mpi_barrier(mpi_comm_world,err)

  !   fill the outer ghost cells
  if (rank.eq.dest) then
     do ieq=1,neq
        u(dad,ieq,nxmin,  :  ) = u(dad,ieq, 0  , :  )
        u(dad,ieq,nxmax,  :  ) = u(dad,ieq,nx+1, :  )
        u(dad,ieq,  :  ,nymin) = u(dad,ieq, :  , 0  )
        u(dad,ieq,  :  ,nymax) = u(dad,ieq, :  ,ny+1)
        !   the corners
        u(dad,ieq,nxmin,nymin) = u(dad,ieq,  0  ,nymin)
        u(dad,ieq,nxmax,nymin) = u(dad,ieq,nxmax, 0   )
        u(dad,ieq,nxmin,nymax) = u(dad,ieq,nxmin,ny+1 )
        u(dad,ieq,nxmax,nymax) = u(dad,ieq, nx+1,nymax)
     end do
     call calcprim(u(dad,:,:,:), primit(dad,:,:,:),                   &
          neq, nxmin, nxmax, nymin, nymax )
     !print*,rank,' just updated nb:',dad
  end if
  call mpi_barrier(mpi_comm_world,err)
  !----------------------------------------------------------------
#endif
  
end subroutine coarse
!=======================================================================
