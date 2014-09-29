!======================================================================
!   boundary conditions for the 1rst order half-timestep (applied to U)
!======================================================================
subroutine boundaryI(time)
  use parameters
  use globals
#ifdef OTHERB
  use user_mod
#endif
  implicit none
  real,    intent(in) :: time
  !
  integer :: nb, nbs, n1, n2, i, j, j1, j2, i1, i2, ny2j, nx2i, ieq
  integer, parameter :: nxp1=nx+1, nyp1=ny+1, nxm1=nx-1, nym1=ny-1
  integer, parameter :: nx2=nx/2,  ny2=ny/2
  integer, parameter :: nx2m1=nx2-1, ny2m1=ny2-1, nx2p1=nx2+1, ny2p1=ny2+1
#ifdef MPIP
  integer, parameter :: countx1=neq*(nx+2) , county1=neq*(ny+2)
  integer, parameter :: countx2=countx1/2  , county2=county1/2
  integer, parameter :: countx3=neq*(nx2+2), county3=neq*(ny2+2)
  integer ::  err, status(MPI_STATUS_SIZE)
  integer ::  source, dest1, dest2, nbd1, nbd2
  real, dimension(1,neq,1,    0:ny2p1) :: bx1l
  real, dimension(1,neq,1,  ny2:nyp1 ) :: bx2l
  real, dimension(1,neq,1,    0:ny2  ) :: bx1s
  real, dimension(1,neq,1,ny2p1:nyp1 ) :: bx2s
  real, dimension(1,neq,    0:nx2p1,1) :: by1l
  real, dimension(1,neq,  nx2:nxp1 ,1) :: by2l
  real, dimension(1,neq,    0:nx2  ,1) :: by1s
  real, dimension(1,neq,nx2p1:nxp1 ,1) :: by2s
#endif
  !
  !  internal boundaries
  do nb=1,ninternal
     select case (innerbounds(nb,3))
        !--------------------------------------------------------------
     case (-1)    ! LEFT BOX BC
        nbs=innerbounds(nb,1)
#ifdef REFXL
        !   reflecting BC'
        u(nbs,1       ,0,0:nyp1) =  u(nbs,1       ,1,0:nyp1) 
        u(nbs,2       ,0,0:nyp1) = -u(nbs,2       ,1,0:nyp1) 
        u(nbs,3:neqpas,0,0:nyp1) =  u(nbs,3:neqpas,1,0:nyp1) 

#endif
#ifdef OUTFXL
        !   outflow BC'
        u(nbs,:,0,0:nyp1) = u(nbs,:,1,0:nyp1) 
#endif
#ifdef INFXL
        !   inflow BC to be adapted to the particular problem'
#endif
        !--------------------------------------------------------------
     case (-2)    ! RIGHT BOX BC
        nbs=innerbounds(nb,1)
#ifdef REFXR
        !   reflecting BC'
        u(nbs,1       ,nxp1,0:nyp1) =  u(nbs,1       ,nx,0:nyp1) 
        u(nbs,2       ,nxp1,0:nyp1) = -u(nbs,2       ,nx,0:nyp1) 
        u(nbs,3:neqpas,nxp1,0:nyp1) =  u(nbs,3:neqpas,nx,0:nyp1) 
#endif
#ifdef OUTFXR
        !   outflow BC'
        u(nbs,:,nxp1,0:nyp1) = u(nbs,:,nx,0:nyp1) 
#endif
#ifdef INFXR
        !   inflow BC wind from the side
        dens=nenv
        u(nbs,1,nxp1,0:nyp1) = dens
        u(nbs,2,nxp1,0:nyp1) = -vism*dens
        u(nbs,3,nxp1,0:nyp1) = 0.
        u(nbs,4,nxp1,0:nyp1) = 0.5*dens*vism**2+cv*1.0001*dens*Tenv/Tempsc
        u(nbs,5,nxp1,0:nyp1) =0.9999*dens
        u(nbs,6,nxp1,0:nyp1) =-dens
#endif
        !--------------------------------------------------------------
     case (-3)    ! BOTTOM BOX BC
        nbs=innerbounds(nb,1)
#ifdef REFYB
        !   reflecting BC
        u(nbs,1:2     ,0:nxp1,0) =  u(nbs,1:2     ,0:nxp1,1) 
        u(nbs,3       ,0:nxp1,0) = -u(nbs,3       ,0:nxp1,1) 
        u(nbs,4:neqpas,0:nxp1,0) =  u(nbs,4:neqpas,0:nxp1,1)       
#endif
#ifdef OUTFYB
        !   outflow BC
        u(nb,:,0:nxp1,0) = u(nb,:,0:nxp1,1) 
#endif
#ifdef INFYB
        !   inflow BC to be adapted to the particular problem'
#endif
        !--------------------------------------------------------------
     case (-4)    ! TOP BOX BC
        nbs=innerbounds(nb,1)
#ifdef REFYT
        !   reflecting BC'
        u(nbs,1:2     ,0:nxp1,nyp1) =  u(nbs,1:2     ,0:nxp1,ny) 
        u(nbs,3       ,0:nxp1,nyp1) = -u(nbs,3       ,0:nxp1,ny) 
        u(nbs,4:neqpas,0:nxp1,nyp1) =  u(nbs,4:neqpas,0:nxp1,ny) 
#endif
#ifdef OUTFYT
        !   outflow BC'
        u(nbs,:,0:nxp1,nyp1) = u(nbs,:,0:nxp1,ny)
#endif
#ifdef INFYT
        !   inflow BC to be adapted to the particular problem'
#endif
        !--------------------------------------------------------------
     case ( 1)    ! LEFT @ same resolution
        nbs=innerbounds(nb,1)
        n1 =innerbounds(nb,2)
        ! own's left
        u(nbs,:,0   ,0:nyp1)=u(n1 ,:,nx,0:nyp1)
        ! neighbors's right
        u(n1 ,:,nxp1,0:nyp1)=u(nbs,:,1 ,0:nyp1)
        !--------------------------------------------------------------
     case ( 2)    ! LEFT @ higher resolution 1
        nbs=innerbounds(nb,1)
        n1 =innerbounds(nb,2)
#ifndef MPIP
        n2 =n1+2
#endif
        !own's left
        do j=1,ny2
           !
           j1 = 2*j-1
           j2 = j1+1
#ifndef MPIP
           ny2j=ny2+j
#endif
           do ieq=1,neq
              u(nbs,ieq,0,j   )= 0.25*sum(u(n1,ieq,nxm1:nx,j1:j2))
#ifndef MPIP
              u(nbs,ieq,0,ny2j)= 0.25*sum(u(n2,ieq,nxm1:nx,j1:j2))
#endif
           end do
        end do
        do ieq=1,neq
           u(nbs,ieq, 0, 0   )= 0.25*sum(u(n1,ieq,nxm1:nx,  -1:0    ))
#ifndef MPIP
           u(nbs,ieq, 0,nyp1 )= 0.25*sum(u(n2,ieq,nxm1:nx,nyp1:nymax))
#endif
        end do
        ! neighbor's right
        ! to finer level (0th order interpolation)
        do j=0,nyp1
           u(n1,:,nxp1,j)=u(nbs,:,1,(j+1)/2)
#ifndef MPIP
           u(n2,:,nxp1,j)=u(nbs,:,1,(j+1)/2+ny2)
#endif
        end do
        !--------------------------------------------------------------
#ifdef MPIP
     case ( 3)    ! LEFT @ higher resolution 2
        nbs=innerbounds(nb,1)
        n2 =innerbounds(nb,2)
        !
        !own's left
        do j=1,ny2
           !
           j1 = 2*j-1
           j2 = j1+1
           ny2j=ny2+j
           do ieq=1,neq
              u(nbs,ieq,0,ny2j)= 0.25*sum(u(n2,ieq,nxm1:nx,  j1:j2   ))
           end do
        end do
        do ieq=1,neq
           u(nbs,ieq, 0,nyp1 )= 0.25*sum(u(n2,ieq,nxm1:nx,nyp1:nymax))
        end do
        ! neighbor's right
        ! to finer level (0th order interpolation)
        do j=0,nyp1
           u(n2,:,nxp1,j)=u(nbs,:,1,(j+1)/2+ny2)
        end do
#endif
        !--------------------------------------------------------------
     case ( 4)    ! RIGHT @ same resolution
        nbs=innerbounds(nb,1)
        n1 =innerbounds(nb,2)
        ! own's right
        u(nbs,:,nxp1,0:nyp1)=u(n1 ,:,1 ,0:nyp1)
        ! neighbors 's left
        u(n1 ,:, 0  ,0:nyp1)=u(nbs,:,nx,0:nyp1)
        !--------------------------------------------------------------
     case ( 5)    ! RIGHT @ higher resolution 1
        nbs=innerbounds(nb,1)
        n1 =innerbounds(nb,2)
#ifndef MPIP
        n2 =n1+2
#endif
        !own's right
        do j=1,ny2
           !
           j1 = 2*j-1
           j2 = j1+1
#ifndef MPIP
           ny2j=ny2+j
#endif
           do ieq=1,neq
              u(nbs,ieq,nxp1,j   )= 0.25*sum(u(n1,ieq,1:2,j1:j2))
#ifndef MPIP              
              u(nbs,ieq,nxp1,ny2j)= 0.25*sum(u(n2,ieq,1:2,j1:j2))
#endif
           end do
        end do
        do ieq=1,neq
           u(nbs,ieq, nxp1, 0   )= 0.25*sum(u(n1,ieq,1:2,   -1:0   ))
#ifndef MPIP
           u(nbs,ieq, nxp1,nyp1 )= 0.25*sum(u(n2,ieq,1:2,nyp1:nymax))
#endif
        end do
        ! neighbor's left
        ! to finer level (0th order interpolation)
        do j=0,nyp1
           u(n1,:,0,j)=u(nbs,:,nx,(j+1)/2)
#ifndef MPIP
           u(n2,:,0,j)=u(nbs,:,nx,(j+1)/2+ny2)
#endif
        end do
        !--------------------------------------------------------------
#ifdef MPIP
     case ( 6)    ! RIGHT @ higher resolution 2
        nbs=innerbounds(nb,1)
        n2 =innerbounds(nb,2)
        !
        !own's right
        do j=1,ny2
           !
           j1 = 2*j-1
           j2 = j1+1
           ny2j=ny2+j
           do ieq=1,neq
              u(nbs,ieq,nxp1,ny2j )= 0.25*sum(u(n2,ieq,1:2,   j1:j2   ))
           end do
        end do
        do ieq=1,neq
           u(nbs,ieq, nxp1,nyp1   )= 0.25*sum(u(n2,ieq,1:2, nyp1:nymax))
        end do
        ! neighbor's left
        ! to finer level (0th order interpolation)
        do j=0,nyp1
           u(n2,:,0,j)=u(nbs,:,nx,(j+1)/2+ny2)
        end do
#endif
        !--------------------------------------------------------------
     case ( 7)    ! BOTTOM @ same resolution
        nbs=innerbounds(nb,1)
        n1 =innerbounds(nb,2)
        ! own's bottom
        u(nbs,:,0:nxp1,0   )=u(n1, :,0:nxp1,ny)
        ! neighbors's top
        u(n1 ,:,0:nxp1,nyp1)=u(nbs,:,0:nxp1,1 )
        !--------------------------------------------------------------
     case ( 8)    ! BOTTOM @ higher resolution 1
        nbs=innerbounds(nb,1)
        n1 =innerbounds(nb,2)
#ifndef MPIP
        n2 =n1+1
#endif
        !own's bottom
        do i=1,nx2
           !
           i1 = 2*i-1
           i2 = i1+1
#ifndef MPIP
           nx2i=nx2+i
#endif
           do ieq=1,neq
              u(nbs,ieq, i  ,0 )= 0.25*sum(u(n1,ieq,i1:i2,nym1:ny))
#ifndef MPIP
              u(nbs,ieq,nx2i,0 )= 0.25*sum(u(n2,ieq,i1:i2,nym1:ny))
#endif
           end do
        end do
        do ieq=1,neq
           u(nbs,ieq, 0  , 0)= 0.25*sum(u(n1,ieq,  -1:0    ,nym1:ny  ))
#ifndef MPIP
           u(nbs,ieq,nxp1, 0)= 0.25*sum(u(n2,ieq,nxp1:nxmax,nym1:ny  ))
#endif
        end do
        ! neighbor's top
        ! to finer level (0th order interpolation)
        do i=0,nxp1
           u(n1,:,i,nyp1)=u(nbs,:,(i+1)/2    ,1)
#ifndef MPIP
           u(n2,:,i,nyp1)=u(nbs,:,(i+1)/2+nx2,1)
#endif
        end do
        !--------------------------------------------------------------
#ifdef MPIP
     case ( 9)    ! BOTTOM @ higher resolution 2
        nbs=innerbounds(nb,1)
        n2 =innerbounds(nb,2)
        !
        !own's bottom
        do i=1,nx2
           !
           i1 = 2*i-1
           i2 = i1+1
           nx2i=nx2+i
           do ieq=1,neq
              u(nbs,ieq,nx2i,0 )= 0.25*sum(u(n2,ieq,   i1:i2   ,nym1:ny))
           end do
        end do
        do ieq=1,neq
           u(nbs,ieq,  nxp1 ,0 )= 0.25*sum(u(n2,ieq, nxp1:nxmax,nym1:ny))
        end do
        ! neighbor's top
        ! to finer level (0th order interpolation)
        do i=0,nxp1
           u(n2,:,i,nyp1)=u(nbs,:,(i+1)/2+nx2 ,1)
        end do
#endif
        !--------------------------------------------------------------
     case (10)    ! TOP @ same resolution
        nbs=innerbounds(nb,1)
        n1 =innerbounds(nb,2)
        ! own's top
        u(nbs,:,0:nxp1,nyp1)=u(n1, :,0:nxp1,1 )
        ! neighbors's bottom
        u(n1 ,:,0:nxp1, 0  )=u(nbs,:,0:nxp1,ny)
        !--------------------------------------------------------------
     case (11)    ! TOP @ higher resolution
        nbs=innerbounds(nb,1)
        n1 =innerbounds(nb,2)
#ifndef MPIP
        n2 =n1+1
#endif
        !own's top
        do i=1,nx2
           !
           i1 = 2*i-1
           i2 = i1+1
#ifndef MPIP
           nx2i=nx2+i
#endif
           do ieq=1,neq
              u(nbs,ieq, i  ,nyp1 )= 0.25*sum(u(n1,ieq,i1:i2,1:2))
#ifndef MPIP
              u(nbs,ieq,nx2i,nyp1 )= 0.25*sum(u(n2,ieq,i1:i2,1:2))
#endif
           end do
        end do
        do ieq=1,neq
           u(nbs,ieq, 0  , nyp1)= 0.25*sum(u(n1,ieq,  -1:0    ,1:2 ))
#ifndef MPIP
           u(nbs,ieq,nxp1, nyp1)= 0.25*sum(u(n2,ieq,nxp1:nxmax,1:2 ))
#endif
        end do
        ! neighbor's botom
        ! to finer level (0th order interpolation)
        do i=0,nxp1
           u(n1,:,i,0)=u(nbs,:,(i+1)/2    ,ny)
#ifndef MPIP
           u(n2,:,i,0)=u(nbs,:,(i+1)/2+nx2,ny)
#endif
        end do
        !--------------------------------------------------------------
#ifdef MPIP
     case (12)    ! TOP @ higher resolution
        nbs=innerbounds(nb,1)
        n2 =innerbounds(nb,2)
        !
        !own's bottom
        do i=1,nx2
           !
           i1 = 2*i-1
           i2 = i1+1
           nx2i=nx2+i
           do ieq=1,neq
              u(nbs,ieq,nx2i,nyp1 )= 0.25*sum(u(n2,ieq,  i1:i2   ,1:2))
           end do
        end do
        do ieq=1,neq
           u(nbs,ieq,   nxp1,nyp1 )= 0.25*sum(u(n2,ieq,nxp1:nxmax,1:2))
        end do
        ! neighbor's bottom
        ! to finer level (0th order interpolation)
        do i=0,nxp1
           u(n2,:,i,0)=u(nbs,:,(i+1)/2+nx2 ,ny)
        end do
#endif
        !--------------------------------------------------------------
     end select
 !-----------------------------------------------------------------
    !   inflow boundaries (wind)
    nbs=innerbounds(nb,1)
#ifdef OTHERB
    !
    call impose_user_bc(u(nbs,:,:,:),nbs,time)
    !
#endif
    !-----------------------------------------------------------------
    !
  end do
  !
#ifdef MPIP
  !
  !if (rank.eq.0) then
  !    do nb=1,nexternal
  !       print'(5i5)',extbounds(nb,:)
  !    end do
  !    print*,'==='
  !end if
  !   the external (MPI) boundaries
  !
   do nb=1,nexternal
     select case (extbounds(nb,5))
        !--------------------------------------------------------------
     case ( 1)    ! LEFT @ same resolution
        nbs    = extbounds(nb,1)
        source = extbounds(nb,2)
        nbd1   = extbounds(nb,3)
        dest1  = extbounds(nb,4)
        !
        if (rank.eq.source)                                      &
          call mpi_sendrecv(u(nbs,:,1,0:nyp1), county1, mpi_real_kind, dest1 ,0, &
                            u(nbs,:,0,0:nyp1), county1, mpi_real_kind, dest1, 0, &
                            mpi_comm_world, status , err)
        !
        if (rank.eq.dest1)                                       &
          call mpi_sendrecv(u(nbd1,:,nx  ,0:nyp1), county1, mpi_real_kind, source,0, &
                            u(nbd1,:,nxp1,0:nyp1), county1, mpi_real_kind, source,0, &
                            mpi_comm_world, status , err)
        !--------------------------------------------------------------
     case ( 2)    ! LEFT @ higher resolution 1
        nbs    = extbounds(nb,1)
        source = extbounds(nb,2)
        nbd1   = extbounds(nb,3)
        dest1  = extbounds(nb,4)
        !
        !   own's left
        if (rank.eq.source) then
           !   fill send buffer
           do ieq=1,neq 
              do j=0,ny2p1
                 bx1l(1,ieq,1,j ) = u(nbs,ieq,1,j )
              end do
           end do
           !   send/receive
           call mpi_sendrecv(bx1l, county3, mpi_real_kind, dest1 ,0, &
                             bx1s, county2, mpi_real_kind, dest1, 0, &
                             mpi_comm_world, status , err)
           !   fill u's with received buffer
           do ieq=1,neq
              do j=0,ny2
                 u(nbs,ieq,0,j )= bx1s(1,ieq,1,j )
              end do
           end do
        end if
        !
        !   neighbor's right
        if (rank.eq.dest1)  then
           !   fill send buffer
           do ieq=1,neq 
              bx1s(   1,ieq,1,0)=0.50*sum( u(nbd1,ieq,nxm1:nx,     0   ) )  
              do j=1,ny2
                 bx1s(1,ieq,1,j)=0.25*sum( u(nbd1,ieq,nxm1:nx,2*j-1:2*j) )
             end do
           end do
           !   send/receive
           call mpi_sendrecv(bx1s, county2, mpi_real_kind, source,0, &
                             bx1l, county3, mpi_real_kind, source,0, &
                             mpi_comm_world, status , err)
           !   fill u's with received buffer
           do ieq=1,neq
              do j=0,nyp1
                 u(nbd1,ieq,nxp1,j)= bx1l(1,ieq,1,(j+1)/2)
              end do
           end do
        end if
        !--------------------------------------------------------------
     case ( 3)    ! LEFT @ higher resolution 2
        nbs    = extbounds(nb,1)
        source = extbounds(nb,2)
        nbd2   = extbounds(nb,3)
        dest2  = extbounds(nb,4)
        !
        !   own's left
        if (rank.eq.source) then
           !   fill send buffer
           do ieq=1,neq 
              do j=ny2,nyp1
                 bx2l(1,ieq,1,j) = u(nbs,ieq,1,j)
              end do
           end do
           !   send/receive
           call mpi_sendrecv(bx2l, county3, mpi_real_kind, dest2 ,0, &
                             bx2s, county2, mpi_real_kind, dest2, 0, &
                             mpi_comm_world, status , err)
           !   fill u's with received buffer
           do ieq=1,neq
              do j=ny2p1,nyp1
                 u(nbs,ieq,0,j )=  bx2s(1,ieq,1,j )
              end do
           end do
        end if
        !
        !   neighbor's right
        if (rank.eq.dest2)  then
           !   fill send buffer
           do ieq=1,neq
              do j=1,ny2
                 bx2s(1,ieq,1,j+ny2)=0.25*sum( u(nbd2,ieq,nxm1:nx,2*j-1:2*j) )
              end do
              bx2s(   1,ieq,1,nyp1 )=0.50*sum( u(nbd2,ieq,nxm1:nx,nyp1     ) )
           end do
           !   send/receive
           call mpi_sendrecv(bx2s, county2, mpi_real_kind, source,0, &
                             bx2l, county3, mpi_real_kind, source,0, &
                             mpi_comm_world, status , err)
           !   fill u's with received buffer
           do ieq=1,neq
              do j=0,nyp1
                 u(nbd2,ieq,nxp1,j)= bx2l(1,ieq,1,(j+1)/2+ny2)
              end do
           end do
        end if
        !--------------------------------------------------------------
     case ( 4)    ! RIGHT @ same resolution
        nbs    = extbounds(nb,1)
        source = extbounds(nb,2)
        nbd1   = extbounds(nb,3)
        dest1  = extbounds(nb,4)
        !
        if (rank.eq.source)                                      &
          call mpi_sendrecv(u(nbs,:,nx  ,0:nyp1), county1, mpi_real_kind, dest1 ,0, &
                            u(nbs,:,nxp1,0:nyp1), county1, mpi_real_kind, dest1, 0, &
                            mpi_comm_world, status , err)
        !
        if (rank.eq.dest1)                                       &
          call mpi_sendrecv(u(nbd1,:,1 ,0:nyp1), county1, mpi_real_kind, source,0, &
                            u(nbd1,:,0 ,0:nyp1), county1, mpi_real_kind, source,0, &
                            mpi_comm_world, status , err)
        !--------------------------------------------------------------
     case ( 5)    ! RIGHT @ higher resolution 1
        nbs    = extbounds(nb,1)
        source = extbounds(nb,2)
        nbd1   = extbounds(nb,3)
        dest1  = extbounds(nb,4)
        !
        !   own's right
        if (rank.eq.source) then
           !   fill send buffer
           do ieq=1,neq 
              do j=0,ny2p1
                 bx1l(1,ieq,1,j ) = u(nbs,ieq,nx,j )
              end do
           end do
           !   send/receive
           call mpi_sendrecv(bx1l, county3, mpi_real_kind, dest1 ,0, &
                             bx1s, county2, mpi_real_kind, dest1, 0, &
                             mpi_comm_world, status , err)
           !   fill u's with received buffer
           do ieq=1,neq
              do j=0,ny2
                 u(nbs,ieq,nxp1,j )= bx1s(1,ieq,1,j )
              end do
           end do
        endif
        !
        !   neighbor's left
        if (rank.eq.dest1)  then
           !   fill send buffer
           do ieq=1,neq 
              bx1s(1,ieq,1,0)   =0.50*sum( u(nbd1,ieq,1:2,0        ) )  
              do j=1,ny2
                 bx1s(1,ieq,1,j)=0.25*sum( u(nbd1,ieq,1:2,2*j-1:2*j) )
              end do
           end do
           !   send/receive
           call mpi_sendrecv(bx1s, county2, mpi_real_kind, source,0, &
                             bx1l, county3, mpi_real_kind, source,0, &
                             mpi_comm_world, status , err)
           !   fill u's with received buffer
           do ieq=1,neq
              do j=0,nyp1
                 u(nbd1,ieq,0,j)=bx1l(1,ieq,1,(j+1)/2)
              end do
           end do
        end if
        !--------------------------------------------------------------
     case ( 6)    ! RIGHT @ higher resolution 2
        nbs    = extbounds(nb,1)
        source = extbounds(nb,2)
        nbd2   = extbounds(nb,3)
        dest2  = extbounds(nb,4)
        !
        !   own's right
        if (rank.eq.source) then
           !   fill send buffer
           do ieq=1,neq
              do j=ny2,nyp1
                 bx2l(1,ieq,1,j) = u(nbs,ieq,nx,j)
              end do
           end do
           !   send/receive
           call mpi_sendrecv(bx2l, county3, mpi_real_kind, dest2 ,0, &
                             bx2s, county2, mpi_real_kind, dest2, 0, &
                             mpi_comm_world, status , err)
           !   fill u's with received buffer
           do ieq=1,neq
              do j=ny2p1,nyp1
                 u(nbs,ieq,nxp1,j )= bx2s(1,ieq,1,j )
              end do
           end do
        endif
        !
        !   neighbor's left
        if (rank.eq.dest2)  then
           !   fill send buffer
           do ieq=1,neq
              do j=1,ny2
                 bx2s(1,ieq,1,j+ny2)=0.25*sum( u(nbd2,ieq,1:2,2*j-1:2*j) )
              end do
              bx2s(1,ieq,1,nyp1)    =0.50*sum( u(nbd2,ieq,1:2,nyp1     ) )
           end do
           !   send/receive
           call mpi_sendrecv(bx2s, county2, mpi_real_kind, source,0, &
                             bx2l, county3, mpi_real_kind, source,0, &
                             mpi_comm_world, status , err)
           !   fill u's with received buffer
           do ieq=1,neq
              do j=0,nyp1
                 u(nbd2,ieq,0,j)= bx2l(1,ieq,1,(j+1)/2+ny2)
              end do
           end do
        end if
        !--------------------------------------------------------------
     case ( 7)    ! BOTTOM @ same resolution
        nbs    = extbounds(nb,1)
        source = extbounds(nb,2)
        nbd1   = extbounds(nb,3)
        dest1  = extbounds(nb,4)
        !
        if (rank.eq.source)                                      &
          call mpi_sendrecv(u(nbs,:,0:nxp1,1 ), countx1, mpi_real_kind, dest1 ,0, &
                            u(nbs,:,0:nxp1,0 ), countx1, mpi_real_kind, dest1, 0, &
                            mpi_comm_world, status , err)
        !
        if (rank.eq.dest1)                                       &
          call mpi_sendrecv(u(nbd1,:,0:nxp1,ny  ), countx1, mpi_real_kind, source,0, &
                            u(nbd1,:,0:nxp1,nxp1), countx1, mpi_real_kind, source,0, &
                            mpi_comm_world, status , err)
        !--------------------------------------------------------------
     case ( 8)    ! BOTTOM @ higher resolution 1
        nbs    = extbounds(nb,1)
        source = extbounds(nb,2)
        nbd1   = extbounds(nb,3)
        dest1  = extbounds(nb,4)
        !
        !   own's bottom
        if (rank.eq.source) then
           !   fill send buffer
           do ieq=1,neq 
              do i=0,nx2p1
                 by1l(1,ieq,i ,1) = u(nbs,ieq,i ,1)
              end do
           end do
           !   send/receive
           call mpi_sendrecv(by1l, countx3, mpi_real_kind, dest1 ,0, &
                             by1s, countx2, mpi_real_kind, dest1, 0, &
                             mpi_comm_world, status , err)
           !   fill u's with received buffer
           do ieq=1,neq
              do i=0,nx2
                 u(nbs,ieq,i ,0)= by1s(1,ieq,i ,1)
              end do
           end do
        endif
        !
        !   neighbor's top
        if (rank.eq.dest1)  then
           !   fill send buffer
           do ieq=1,neq 
              by1s(1,ieq,0,1)   =0.50*sum( u(nbd1,ieq,     0   ,nym1:ny) )  
              do i=1,nx2
                 by1s(1,ieq,i,1)=0.25*sum( u(nbd1,ieq,2*i-1:2*i,nym1:ny) )
              end do
           end do
           !   send/receive
           call mpi_sendrecv(by1s, countx2, mpi_real_kind, source,0, &
                             by1l, countx3, mpi_real_kind, source,0, &
                             mpi_comm_world, status , err)
           !   fill u's with received buffer
           do ieq=1,neq
              do i=0,nxp1
                 u(nbd1,ieq,i,nyp1) = by1l(1,ieq,(i+1)/2,1)
              end do
           end do
        end if
        !--------------------------------------------------------------
     case ( 9)    ! BOTTOM @ higher resolution 2
        nbs    = extbounds(nb,1)
        source = extbounds(nb,2)
        nbd2   = extbounds(nb,3)
        dest2  = extbounds(nb,4)
        !
        !   own's bottom
        if (rank.eq.source) then
           !   fill send buffer
           do ieq=1,neq 
              do i=nx2,nxp1
                 by2l(1,ieq,i,1) = u(nbs,ieq,i,1)
              end do
           end do
           !   send/receive
           call mpi_sendrecv(by2l, countx3, mpi_real_kind, dest2 ,0, &
                             by2s, countx2, mpi_real_kind, dest2, 0, &
                             mpi_comm_world, status , err)
           !   fill u's with received buffer
           do ieq=1,neq
              do i=nx2p1,nxp1
                 u(nbs,ieq,i,0)=by2s(1,ieq,i,1)
              end do
           end do
        endif
        !
        !   neighbor's top
        if (rank.eq.dest2)  then
           !   fill send buffer
           do ieq=1,neq
              do i=1,nx2
                 by2s(1,ieq,i+nx2,1)=0.25*sum( u(nbd2,ieq,2*i-1:2*i,nym1:ny) )
              end do
              by2s(1,ieq,nxp1,1)    =0.50*sum( u(nbd2,ieq,  nxp1   ,nym1:ny) )
           end do
           !   send/receive
           call mpi_sendrecv(by2s, countx2, mpi_real_kind, source,0, &
                             by2l, countx3, mpi_real_kind, source,0, &
                             mpi_comm_world, status , err)
           !   fill u's with received buffer
           do ieq=1,neq
              do i=0,nxp1
                 u(nbd2,ieq,i,nyp1)= by2l(1,ieq,(i+1)/2+nx2,1)
              end do
           end do
        end if
        !--------------------------------------------------------------
     case (10)    ! TOP @ same resolution
        nbs    = extbounds(nb,1)
        source = extbounds(nb,2)
        nbd1   = extbounds(nb,3)
        dest1  = extbounds(nb,4)
        !
        if (rank.eq.source)                                      &
          call mpi_sendrecv(u(nbs,:,0:nxp1,ny  ), countx1, mpi_real_kind, dest1 ,0, &
                            u(nbs,:,0:nxp1,nxp1), countx1, mpi_real_kind, dest1, 0, &
                            mpi_comm_world, status , err)
        !
        if (rank.eq.dest1)                                       &
          call mpi_sendrecv(u(nbd1,:,0:nxp1,1), countx1, mpi_real_kind, source,0, &
                            u(nbd1,:,0:nxp1,0), countx1, mpi_real_kind, source,0, &
                            mpi_comm_world, status , err)
        !--------------------------------------------------------------
     case (11)    ! TOP @ higher resolution 1 
        nbs    = extbounds(nb,1)
        source = extbounds(nb,2)
        nbd1   = extbounds(nb,3)
        dest1  = extbounds(nb,4)
        !
        !   own's top
        if (rank.eq.source) then
           !   fill send buffer
           do ieq=1,neq 
              do i=0,nx2p1
                 by1l(1,ieq,i ,1) = u(nbs,ieq,i ,ny)
              end do
           end do
           !   send/receive
           call mpi_sendrecv(by1l, countx3, mpi_real_kind, dest1 ,0, &
                             by1s, countx2, mpi_real_kind, dest1, 0, &
                             mpi_comm_world, status , err)
           !   fill u's with received buffer
           do ieq=1,neq
              do i=0,nx2
                 u(nbs,ieq,i ,nyp1)= by1s(1,ieq,i ,1)
              end do
           end do
        endif
        !
        !   neighbor's bottom
        if (rank.eq.dest1)  then
           !   fill send buffer
           do ieq=1,neq 
              by1s(1,ieq,0,1)   =0.50*sum( u(nbd1,ieq,     0   ,1:2) )  
              do i=1,nx2
                 by1s(1,ieq,i,1)=0.25*sum( u(nbd1,ieq,2*i-1:2*i,1:2) )
              end do
           end do
           !   send/receive
           call mpi_sendrecv(by1s, countx2, mpi_real_kind, source,0, &
                             by1l, countx3, mpi_real_kind, source,0, &
                             mpi_comm_world, status , err)
           !   fill u's with received buffer
           do ieq=1,neq
              do i=0,nxp1
                 u(nbd1,ieq,i,0) =  by1l(1,ieq,(i+1)/2,1)
              end do
           end do
        end if
        !call mpi_finalize(err)
        !stop
        !--------------------------------------------------------------
     case (12)    ! TOP @ higher resolution 2
        nbs    = extbounds(nb,1)
        source = extbounds(nb,2)
        nbd2   = extbounds(nb,3)
        dest2  = extbounds(nb,4)
        !
        !   own's top
        if (rank.eq.source) then
           !   fill send buffer   
           do ieq=1,neq 
              do i=nx2,nxp1
                 by2l(1,ieq,i,1) = u(nbs,ieq,i,ny)
              end do
           end do
           !   send/receive
           call mpi_sendrecv(by2l, countx3, mpi_real_kind, dest2 ,0, &
                             by2s, countx2, mpi_real_kind, dest2, 0, &
                             mpi_comm_world, status , err)
           !   fill u's with received buffer
           do ieq=1,neq
              do i=nx2p1,nxp1
                 u(nbs,ieq,i,nyp1)=by2s(1,ieq,i,1)
              end do
           end do
        endif
        !
        !   neighbor's bottom
        if (rank.eq.dest2)  then
           !   fill send buffer
           do ieq=1,neq
              do i=1,nx2
                 by2s(1,ieq,i+nx2,1)=0.25*sum( u(nbd2,ieq,2*i-1:2*i,1:2) )
              end do
              by2s(1,ieq,nxp1,1)    =0.50*sum( u(nbd2,ieq,  nxp1   ,1:2) )
           end do
           !   send/receive
           call mpi_sendrecv(by2s, countx2, mpi_real_kind, source,0, &
                             by2l, countx3, mpi_real_kind, source,0, &
                             mpi_comm_world, status , err)
           !   fill u's with received buffer
           do ieq=1,neq
              do i=0,nxp1
                 u(nbd2,ieq,i,0)=by2l(1,ieq,(i+1)/2+nx2,1)
              end do
           end do
        end if
        !--------------------------------------------------------------
     end select
     !
  end do
  !
#endif
  !
  return
end subroutine boundaryI
!======================================================================
