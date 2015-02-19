!======================================================================
!   boundary conditions for the 2nd order half-timestep
!======================================================================
subroutine boundaryII(time)
  use parameters
  use globals
#ifdef OTHERB
  use user_mod
#endif
  implicit none
  real,    intent(in) :: time
  !
  integer :: nb, nbs, n1, n2, i, j, j1,j2, i1, i2, ny2j, nx2i, ieq
  integer, parameter :: nxp1=nx+1, nyp1=ny+1, nxm1=nx-1, nym1=ny-1
  integer, parameter :: nx2=nx/2,  ny2=ny/2
  integer, parameter :: nx2m1=nx2-1, ny2m1=ny2-1, nx2p1=nx2+1, ny2p1=ny2+1
  integer, parameter :: nx2p2=nx2+2, ny2p2=ny2+2
#ifdef MPIP
  integer, parameter :: countx1=2*neq*(nx+4) , county1=2*neq*(ny+4)
  integer, parameter :: countx2=2*neq*(nx2+1), countx3=neq*(nx2+2)
  integer, parameter :: county2=2*neq*(ny2+1), county3=neq*(ny2+2)
  integer :: err, status(MPI_STATUS_SIZE)
  integer :: source, dest1, dest2
  real, dimension(1,neq,     1,         0:ny2p1) :: bx1l
  real, dimension(1,neq,     1,       ny2:nyp1 ) :: bx2l
  real, dimension(1,neq,     2,         0:ny2  ) :: bx1s
  real, dimension(1,neq,     2,     ny2p1:nyp1 ) :: bx2s
  real, dimension(1,neq,    0:nx2p1,     1     ) :: by1l
  real, dimension(1,neq,  nx2:nxp1,      1     ) :: by2l
  real, dimension(1,neq,    0:nx2,       2     ) :: by1s
  real, dimension(1,neq,nx2p1:nxp1,      2     ) :: by2s
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
        j=2
        do i=-1,0
           up(nbs,1       ,i,:) =  up(nbs,1       ,j,:) 
           up(nbs,2       ,i,:) = -up(nbs,2       ,j,:) 
           up(nbs,3:neqpas,i,:) =  up(nbs,3:neqpas,j,:) 
           j=j-1
        end do
#endif
#ifdef OUTFXL
        !   outflow BC'
        j=2
        do i=-1,0
           up(nbs,:,i,:) = up(nbs,:,j,:)
           j=j-1
        end do
#endif
#ifdef INFXL
        !   inflow BC to be adapted to the particular problem'
#endif
        !--------------------------------------------------------------
     case (-2)    ! RIGHT BOX BC
        nbs=innerbounds(nb,1)
#ifdef REFXR
        !   reflecting BC'
        j=nx
        do i=nxp1,nxmax
           up(nbs,1       ,i,:) =  up(nbs,1       ,j,:) 
           up(nbs,2       ,i,:) = -up(nbs,2       ,j,:) 
           up(nbs,3:neqpas,i,:) =  up(nbs,3:neqpas,j,:) 
           j=j-1
        end do
#endif
#ifdef OUTFXR
        !   outflow BC'
        j=nx
        do i=nxp1,nxmax
           up(nbs,:,i,:) = up(nbs,:,j,:)
           j=j-1
        end do
#endif
#ifdef INFXR
        !   inflow BC wind from the side
        dens=nenv
        do i=nxp1,nxmax
           up(nbs,1,i,-1:nymax) = dens
           up(nbs,2,i,-1:nymax) = -vism*dens
           up(nbs,3,i,-1:nymax) = 0.
           up(nbs,4,i,-1:nymax) = 0.5*dens*vism**2+cv*1.0001*dens*Tenv/Tempsc
           up(nbs,5,i,-1:nymax) = 0.9999*dens
           up(nbs,6,i,-1:nymax) = -dens
        end do
#endif
        !--------------------------------------------------------------
     case (-3)    ! BOTTOM BOX BC
        nbs=innerbounds(nb,1)
#ifdef REFYB
        !   reflecting BC
        j=2
        do i=-1,0
           up(nbs,1:2     ,:,i) =  up(nbs,1:2     ,:,j) 
           up(nbs,3       ,:,i) = -up(nbs,3       ,:,j) 
           up(nbs,4:neqpas,:,i) =  up(nbs,4:neqpas,:,j)
           j=j-1
        end do
#endif
#ifdef OUTFYB
        !   outflow BC
        j=2
        do i=-1,0
           up(nbs,:,:,i) = up(nbs,:,:,j)
           j=j-1
        end do
#endif
#ifdef INFYB
        !   inflow BC to be adapted to the particular problem'
#endif
        !--------------------------------------------------------------
     case (-4)    ! TOP BOX BC
        nbs=innerbounds(nb,1)
#ifdef REFYT
        !   reflecting BC'
        j=ny
        do i=nyp1,nymax
           up(nbs,1:2     ,:,i) =  up(nbs,1:2     ,:,j) 
           up(nbs,3       ,:,i) = -up(nbs,3       ,:,j) 
           up(nbs,4:neqpas,:,i) =  up(nbs,4:neqpas,:,j)
           j=j-1
        end do
#endif
#ifdef OUTFYT
        !   outflow BC'
        j=ny
        do i=nyp1,nymax
           up(nbs,:,:,i) = up(nbs,:,:,j)
           j=j-1
        end do
#endif
#ifdef INFYT
        !   inflow BC to be adapted to the particular problem'
#endif
        !--------------------------------------------------------------
     case ( 1)    ! LEFT @ same resolution
        nbs=innerbounds(nb,1)
        n1 =innerbounds(nb,2)
        ! own's left
        up(nbs,:,   -1:0   ,:)=up(n1 ,:,nxm1:nx,:)
        ! neighbors's right
        up(n1 ,:,nxp1:nxmax,:)=up(nbs,:,   1:2 ,:)
        !--------------------------------------------------------------
     case ( 2)    ! LEFT @ higher resolution 1
        nbs=innerbounds(nb,1)
        n1 =innerbounds(nb,2)
#ifndef MPIP
        n2 =n1+2
#endif
        !own's left
        do j=1,ny2
#ifndef MPIP
           ny2j=ny2+j
#endif
           do i=-1,0   
              !
              i1 = 2*i-1+nx
              i2 = i1+1
              j1 = 2*j-1
              j2 = j1+1
              do ieq=1,neq
                 up(nbs,ieq,i,j   )= 0.25*sum(up(n1,ieq,i1:i2,j1:j2))
#ifndef MPIP
                 up(nbs,ieq,i,ny2j)= 0.25*sum(up(n2,ieq,i1:i2,j1:j2))
#endif
              end do
           end do
        end do
        do ieq=1,neq
           up(nbs,ieq, 0, 0 )= 0.25*sum(up(n1,ieq,nxm1:nx,   -1:0    ))
           up(nbs,ieq,-1, 0 )= 0.25*sum(up(n1,ieq,nx-3:nx-2 ,-1:0    ))
           up(nbs,ieq, 0,-1 )= up(nbs,ieq, 0, 0 )
           up(nbs,ieq,-1,-1 )= up(nbs,ieq,-1, 0 )
#ifndef MPIP
           up(nbs,ieq, 0,nyp1 )= 0.25*sum(up(n2,ieq,nxm1:nx,   nyp1:nymax))
           up(nbs,ieq,-1,nyp1 )= 0.25*sum(up(n2,ieq,nx-3:nx-2 ,nyp1:nymax))
           up(nbs,ieq, 0,nymax)= up(nbs,ieq, 0,nyp1)
           up(nbs,ieq,-1,nymax)= up(nbs,ieq,-1,nyp1)
#endif
        end do
        ! neighbor's right
        ! to finer level (0th order interpolation)
        do i=1,2
           do j=-1,nymax
              up(n1,:,i+nx,j)=up(nbs,:,(i+1)/2,(j+1)/2    )
#ifndef MPIP
              up(n2,:,i+nx,j)=up(nbs,:,(i+1)/2,(j+1)/2+ny2)
#endif
           end do
        end do
        !--------------------------------------------------------------
#ifdef MPIP
     case ( 3)    ! LEFT @ higher resolution 2
        nbs=innerbounds(nb,1)
        n2 =innerbounds(nb,2)
        !
        !own's left
        do i=-1,0
           do j=1,ny2
              !
              i1 = 2*i-1+nx
              i2 = i1+1
              j1 = 2*j-1
              j2 = j1+1
              ny2j=ny2+j
              do ieq=1,neq
                 up(nbs,ieq,i,ny2j)= 0.25*sum(up(n2,ieq,i1:i2,j1:j2))
              end do
           end do
        end do
        do ieq=1,neq
           up(nbs,ieq, 0,nyp1 )= 0.25*sum(up(n2,ieq,nxm1:nx  ,nyp1:nymax))
           up(nbs,ieq,-1,nyp1 )= 0.25*sum(up(n2,ieq,nx-3:nx-2,nyp1:nymax))
           up(nbs,ieq, 0,nymax)= up(nbs,ieq, 0,nyp1)
           up(nbs,ieq,-1,nymax)= up(nbs,ieq,-1,nyp1)
        end do
        ! neighbor's right
        ! to finer level (0th order interpolation)
        do i=1,2
           do j=-1,nymax
              up(n2,:,i+nx,j)=up(nbs,:,(i+1)/2,(j+1)/2+ny2)
           end do
        end do
#endif
        !--------------------------------------------------------------
     case ( 4)    ! RIGHT @ same resolution
        nbs=innerbounds(nb,1)
        n1 =innerbounds(nb,2)
        ! own's right
        up(nbs,:,nxp1:nxmax,:)=up(n1 ,:,   1:2 ,:)
        ! neighbors 's left
        up(n1 ,:,  -1:0    ,:)=up(nbs,:,nxm1:nx,:)
        !--------------------------------------------------------------
     case ( 5)    ! RIGHT @ higher resolution 1
        nbs=innerbounds(nb,1)
        n1 =innerbounds(nb,2)
#ifndef MPIP
        n2 =n1+2
#endif
        !own's right
        do i=nxp1,nxmax
           do j=1,ny2
              !
              i1 = 2*(i-nx)-1
              i2 = i1+1
              j1 = 2*j-1
              j2 = j1+1
#ifndef MPIP
              ny2j=ny2+j
#endif
              do ieq=1,neq
                 up(nbs,ieq,i,j   )= 0.25*sum(up(n1,ieq,i1:i2,j1:j2))
#ifndef MPIP              
                 up(nbs,ieq,i,ny2j)= 0.25*sum(up(n2,ieq,i1:i2,j1:j2))
#endif
              end do
           end do
        end do
        do ieq=1,neq
           up(nbs,ieq, nxp1, 0 )= 0.25*sum(up(n1,ieq,1:2,   -1:0   ))
           up(nbs,ieq,nxmax, 0 )= 0.25*sum(up(n1,ieq,3:4 ,  -1:0    ))
           up(nbs,ieq,nxp1 ,-1 )= up(nbs,ieq,nxp1 ,0 )  
           up(nbs,ieq,nxmax,-1 )= up(nbs,ieq,nxmax,0 )
#ifndef MPIP
           up(nbs,ieq, nxp1,nyp1 )= 0.25*sum(up(n2,ieq,1:2,nyp1:nymax))
           up(nbs,ieq,nxmax,nyp1 )= 0.25*sum(up(n2,ieq,3:4,nyp1:nymax))
           up(nbs,ieq,nxp1 ,nymax)= up(nbs,ieq,nxp1 ,nyp1)
           up(nbs,ieq,nxmax,nymax)= up(nbs,ieq,nxmax,nyp1)
#endif
        end do
        ! neighbor's left
        ! to finer level (0th order interpolation)
        do i=1,2
           do j=-1,nymax
              up(n1,:,i-2,j)=up(nbs,:,(i+1)/2+nxm1,(j+1)/2)
#ifndef MPIP
              up(n2,:,i-2,j)=up(nbs,:,(i+1)/2+nxm1,(j+1)/2+ny2)
#endif
           end do
        end do
        !--------------------------------------------------------------
#ifdef MPIP
     case ( 6)    ! RIGHT @ higher resolution 2
        nbs=innerbounds(nb,1)
        n2 =innerbounds(nb,2)
        !
        !own's right
        do i=nxp1,nxmax
           do j=1,ny2
              !
              i1 = 2*(i-nx)-1
              i2 = i1+1
              j1 = 2*j-1
              j2 = j1+1
              ny2j=ny2+j
              do ieq=1,neq
                 up(nbs,ieq,i,ny2j )= 0.25*sum(up(n2,ieq,i1:i2, j1:j2   ))
              end do
           end do
        end do
        do ieq=1,neq
           up(nbs,ieq,nxp1 ,nyp1 )= 0.25*sum(up(n2,ieq,1:2,nyp1:nymax))
           up(nbs,ieq,nxmax,nyp1 )= 0.25*sum(up(n2,ieq,3:4,nyp1:nymax))
           up(nbs,ieq,nxp1 ,nymax)= up(nbs,ieq,nxp1 ,nyp1)
           up(nbs,ieq,nxmax,nymax)= up(nbs,ieq,nxmax,nyp1)
        end do
        ! neighbor's left
        ! to finer level (0th order interpolation)
        do i=1,2
           do j=-1,nymax
              up(n2,:,i-2,j)=up(nbs,:,(i+1)/2+nxm1,(j+1)/2+ny2)
           end do
        end do
#endif
        !--------------------------------------------------------------
     case ( 7)    ! BOTTOM @ same resolution
        nbs=innerbounds(nb,1)
        n1 =innerbounds(nb,2)
        ! own's bottom
        up(nbs,:,:,  -1:0    )=up(n1, :,:,nym1:ny)
        ! neighbors's top
        up(n1 ,:,:,nyp1:nymax)=up(nbs,:,:,   1:2 )
        !--------------------------------------------------------------
     case ( 8)    ! BOTTOM @ higher resolution 1
        nbs=innerbounds(nb,1)
        n1 =innerbounds(nb,2)
#ifndef MPIP
        n2 =n1+1
#endif
        !own's bottom
        do i=1,nx2
#ifndef MPIP
           nx2i=nx2+i
#endif
           do j=-1,0
              !
              i1 = 2*i-1
              i2 = i1+1
              j1 = 2*j-1+ny
              j2 = j1+1
              do ieq=1,neq
                 up(nbs,ieq, i  ,j )= 0.25*sum(up(n1,ieq,i1:i2,j1:j2))
#ifndef MPIP
                 up(nbs,ieq,nx2i,j )= 0.25*sum(up(n2,ieq,i1:i2,j1:j2))
#endif
              end do
           end do
        end do
        do ieq=1,neq
           up(nbs,ieq, 0,   0 )= 0.25*sum(up(n1,ieq,  -1:0    ,nym1:ny  ))
           up(nbs,ieq, 0,  -1 )= 0.25*sum(up(n1,ieq,  -1:0    ,ny-3:ny-2))
           up(nbs,ieq,-1,   0 )= up(nbs,ieq, 0, 0 )
           up(nbs,ieq,-1,  -1 )= up(nbs,ieq,-1, 0 )
#ifndef MPIP
           up(nbs,ieq,nxp1,  0)= 0.25*sum(up(n2,ieq,nxp1:nxmax,nym1:ny  ))
           up(nbs,ieq,nxp1, -1)= 0.25*sum(up(n2,ieq,nxp1:nxmax,ny-3:ny-2))
           up(nbs,ieq,nxmax, 0)= up(nbs,ieq,nxp1, 0)
           up(nbs,ieq,nxmax,-1)= up(nbs,ieq,nxp1,-1)
#endif
        end do
        ! neighbor's top
        ! to finer level (0th order interpolation)
        do i=-1,nxmax
           do j=1,2
              up(n1,:,i,j+ny)=up(nbs,:,(i+1)/2    ,(j+1)/2)
#ifndef MPIP
              up(n2,:,i,j+ny)=up(nbs,:,(i+1)/2+nx2,(j+1)/2)
#endif
           end do
        end do
        !--------------------------------------------------------------
#ifdef MPIP
     case ( 9)    ! BOTTOM @ higher resolution 2
        nbs=innerbounds(nb,1)
        n2 =innerbounds(nb,2)
        !
        !own's bottom
        do i=1,nx2
           do j=-1,0
              !
              i1 = 2*i-1
              i2 = i1+1
              j1 = 2*j-1+ny
              j2 = j1+1
              nx2i=nx2+i
              do ieq=1,neq
                 up(nbs,ieq,nx2i,j )= 0.25*sum(up(n2,ieq,   i1:i2   ,j1:j2))
              end do
           end do
        end do
        do ieq=1,neq
           up(nbs,ieq,nxp1,  0)= 0.25*sum(up(n2,ieq,nxp1:nxmax,nym1:ny  ))
           up(nbs,ieq,nxp1, -1)= 0.25*sum(up(n2,ieq,nxp1:nxmax,ny-3:ny-2))
           up(nbs,ieq,nxmax, 0)= up(nbs,ieq,nxp1, 0)
           up(nbs,ieq,nxmax,-1)= up(nbs,ieq,nxp1,-1)
        end do
        ! neighbor's top
        ! to finer level (0th order interpolation)
        do i=-1,nxmax
           do j=1,2
              up(n2,:,i,j+ny)=up(nbs,:,(i+1)/2+nx2 ,(j+1)/2)
           end do
        end do
#endif
        !--------------------------------------------------------------
     case (10)    ! TOP @ same resolution
        nbs=innerbounds(nb,1)
        n1 =innerbounds(nb,2)
        ! own's top
        up(nbs,:,:,nyp1:nxmax)=up(n1, :,:,   1:2 )
        ! neighbors's bottom
        up(n1 ,:,:,  -1:0    )=up(nbs,:,:,nym1:ny)
        !--------------------------------------------------------------
     case (11)    ! TOP @ higher resolution 1
        nbs=innerbounds(nb,1)
        n1 =innerbounds(nb,2)
#ifndef MPIP
        n2 =n1+1
#endif
        !own's top
        do i=1,nx2
           do j=nyp1,nymax
              !
              i1 = 2*i-1
              i2 = i1+1
              j1 = 2*(j-ny)-1
              j2 = j1+1
#ifndef MPIP
              nx2i=nx2+i
#endif
              do ieq=1,neq
                 up(nbs,ieq, i  ,j )= 0.25*sum(up(n1,ieq,i1:i2,j1:j2))
#ifndef MPIP
                 up(nbs,ieq,nx2i,j )= 0.25*sum(up(n2,ieq,i1:i2,j1:j2))
#endif
              end do
           end do
        end do
        do ieq=1,neq
           up(nbs,ieq, 0,nyp1 )= 0.25*sum(up(n1,ieq,-1:0,1:2 ))
           up(nbs,ieq, 0,nymax)= 0.25*sum(up(n1,ieq,-1:0,3:4 ))
           up(nbs,ieq,-1,nyp1 )= up(nbs,ieq,0,nyp1 )
           up(nbs,ieq,-1,nymax)= up(nbs,ieq,0,nymax)
#ifndef MPIP
           up(nbs,ieq,nxp1 ,nyp1 )= 0.25*sum(up(n2,ieq,nxp1:nxmax,1:2 ))
           up(nbs,ieq,nxp1 ,nymax)= 0.25*sum(up(n2,ieq,nxp1:nxmax,3:4))
           up(nbs,ieq,nxmax,nyp1 )= up(nbs,ieq,nxp1 ,nyp1 )
           up(nbs,ieq,nxmax,nymax)= up(nbs,ieq,nxp1 ,nymax)
#endif
        end do
        ! neighbor's botom
        ! to finer level (0th order interpolation)
        do i=-1,nxmax
           do j=1,2
              up(n1,:,i,j-2)=up(nbs,:,(i+1)/2    ,(j+1)/2+nym1)
#ifndef MPIP
              up(n2,:,i,j-2)=up(nbs,:,(i+1)/2+nx2,(j+1)/2+nym1)
#endif
           end do
        end do
        !--------------------------------------------------------------
#ifdef MPIP
     case (12)    ! TOP @ higher resolution 2
        nbs=innerbounds(nb,1)
        n2 =innerbounds(nb,2)
        !
        !own's bottom
        do i=1,nx2
           do j=nyp1,nymax
              !
              i1 = 2*i-1
              i2 = i1+1
              j1 = 2*(j-ny)-1
              j2 = j1+1
              nx2i=nx2+i
              do ieq=1,neq
                 up(nbs,ieq,nx2i,j )= 0.25*sum(up(n2,ieq,  i1:i2   ,j1:j2))
              end do
           end do
        end do
        do ieq=1,neq
           up(nbs,ieq,nxp1, nyp1 )= 0.25*sum(up(n2,ieq,nxp1:nxmax,1:2))
           up(nbs,ieq,nxp1, nymax)= 0.25*sum(up(n2,ieq,nxp1:nxmax,3:4))
           up(nbs,ieq,nxmax,nyp1 )= up(nbs,ieq,nxp1 ,nyp1 )
           up(nbs,ieq,nxmax,nymax)= up(nbs,ieq,nxp1 ,nymax)
        end do
        ! neighbor's bottom
        ! to finer level (0th order interpolation)
        do i=-1,nxmax
           do j=1,2
              up(n2,:,i,j-2)=up(nbs,:,(i+1)/2+nx2 ,(j+1)/2+nym1)
           end do
        end do
#endif
        !--------------------------------------------------------------
     end select
     !-----------------------------------------------------------------
     !   inflow boundaries (wind)
     nbs=innerbounds(nb,1)
#ifdef OTHERB
    !
    call impose_user_bc(up(nbs,:,:,:),nbs,time)
    !
#endif
     !-----------------------------------------------------------------
     !
  end do
  !
#ifdef MPIP
  !
  !if (rank.eq.0) then
  !   do nb=1,nexternal
  !      print'(5i6)',extbounds(nb,:)
  !   end do
  !   print*,nx,ny
  !   print*,countx1,countx2,countx3 
  !   print*,county1,county2,county3 
  !   print*,'==='
  !end if
  !
  !   the external (MPI) boundaries
  !
  !print*,rank,nexternal,extbounds(nb,5)
  !call mpi_finalize(nb)
  !stop
  do nb=1,nexternal
     !print*,rank,nb,extbounds(nb,5)
     select case (extbounds(nb,5))
        !--------------------------------------------------------------
     case ( 1)    ! LEFT @ same resolution
        nbs    = extbounds(nb,1)
        source = extbounds(nb,2)
        n1     = extbounds(nb,3)
        dest1  = extbounds(nb,4)
        !
        if (rank.eq.source)                                            &
          call mpi_sendrecv(up(nbs,:, 1:2,-1:nymax), county1, mpi_real_kind, &
                                                     dest1 ,0,          &
                            up(nbs,:,-1:0,-1:nymax), county1, mpi_real_kind, &
                                                     dest1, 0,          &
                            mpi_comm_world, status , err)
        !
        if (rank.eq.dest1)                                             &
          call mpi_sendrecv(up(n1,:,nxm1:nx   ,-1:nymax), county1, mpi_real_kind, &
                                                          source,0,          &
                            up(n1,:,nxp1:nxmax,-1:nymax), county1, mpi_real_kind, &
                                                          source,0,          &
                            mpi_comm_world, status , err)
        !--------------------------------------------------------------
     case ( 2)    ! LEFT @ higher resolution 1
        nbs    = extbounds(nb,1)
        source = extbounds(nb,2)
        n1     = extbounds(nb,3)
        dest1  = extbounds(nb,4)
        !
        !   own's left
        if (rank.eq.source) then
           !   fill send buffer
           do ieq=1,neq 
              do j=0,ny2p1
                 bx1l(1,ieq,1,j ) = up(nbs,ieq,1,j )
              end do
           end do
           !   send/receive
           call mpi_sendrecv(bx1l, county3, mpi_real_kind, dest1 ,0, &
                             bx1s, county2, mpi_real_kind, dest1, 0, &
                             mpi_comm_world, status , err)
           !   fill u's with received buffer
           do i=-1,0
              i1=i+2
              do j=0,ny2
                 do ieq=1,neq
                    up(nbs,ieq,i,j )= bx1s(1,ieq,i1,j )
                 end do
              end do
           end do
           do ieq=1,neq
              up(nbs,ieq, 0,-1 )= up(nbs,ieq, 0, 0 )
              up(nbs,ieq,-1,-1 )= up(nbs,ieq,-1, 0 )
           end do
        end if
        !
        !   neighbor's right
        if (rank.eq.dest1)  then
           !   fill send buffer
           do j=0,ny2
              j1=2*j-1
              j2=j1+1
              do i=-1,0
                 i1=2*i-1+nx
                 i2 = i1+1
                 do ieq=1,neq
                    bx1s(1,ieq,i+2,j)=0.25*sum( up(n1,ieq,i1:i2,j1:j2) )
                 end do
              end do
           end do
           !   send/receive
           call mpi_sendrecv(bx1s, county2, mpi_real_kind, source,0, &
                             bx1l, county3, mpi_real_kind, source,0, &
                             mpi_comm_world, status , err)
           !   fill u's with received buffer
           do i=nxp1,nxmax
              do j=-1,nymax
                 do ieq=1,neq
                    up(n1,ieq,i,j)= bx1l(1,ieq,1,(j+1)/2)
                 end do
              end do
           end do
        end if
        !--------------------------------------------------------------
     case ( 3)    ! LEFT @ higher resolution 2
        nbs    = extbounds(nb,1)
        source = extbounds(nb,2)
        n2   = extbounds(nb,3)
        dest2  = extbounds(nb,4)
        !
        !   own's left
        if (rank.eq.source) then
           !   fill send buffer
           do ieq=1,neq 
              do j=ny2,nyp1
                 bx2l(1,ieq,1,j) = up(nbs,ieq,1,j)
              end do
           end do
           !   send/receive
           call mpi_sendrecv(bx2l, county3, mpi_real_kind, dest2 ,0, &
                             bx2s, county2, mpi_real_kind, dest2, 0, &
                             mpi_comm_world, status , err)
           !   fill u's with received buffer
           do i=-1,0
              i1=i+2
              do j=ny2p1,nyp1
                 do ieq=1,neq
                    up(nbs,ieq,i,j )=  bx2s(1,ieq,i1,j )
                 end do
              end do
           end do
          do ieq=1,neq
            up(nbs,ieq, 0,nymax )= up(nbs,ieq, 0,nyp1 )  
            up(nbs,ieq,-1,nymax )= up(nbs,ieq,-1,nyp1 )
          end do
        end if
        !
        !   neighbor's right
        if (rank.eq.dest2)  then
           !   fill send buffer
           do j=1,ny2p1
              j1=2*j-1
              j2=j1+1
              ny2j=ny2+j
              do i=-1,0
                 i1=2*i-1+nx
                 i2 = i1+1
                 do ieq=1,neq
                    bx2s(1,ieq,i+2,ny2j)=0.25*sum( up(n2,ieq,i1:i2,j1:j2) )
                 end do
              end do
           end do
           !   send/receive
           call mpi_sendrecv(bx2s, county2, mpi_real_kind, source,0, &
                             bx2l, county3, mpi_real_kind, source,0, &
                             mpi_comm_world, status , err)
           !   fill u's with received buffer
           do i=nxp1,nxmax
              do j=-1,nymax
                 do ieq=1,neq
                    up(n2,ieq,i,j)= bx2l(1,ieq,1,(j+1)/2+ny2)
                 end do
              end do
           end do
        end if
        !--------------------------------------------------------------
     case ( 4)    ! RIGHT @ same resolution
        nbs    = extbounds(nb,1)
        source = extbounds(nb,2)
        n1   = extbounds(nb,3)
        dest1  = extbounds(nb,4)
        !
        if (rank.eq.source)                                      &
          call mpi_sendrecv(up(nbs,:,nxm1:nx,   -1:nymax), county1, mpi_real_kind, &
                                                           dest1 ,0,          &
                            up(nbs,:,nxp1:nxmax,-1:nymax), county1, mpi_real_kind, &
                                                           dest1, 0,          &
                            mpi_comm_world, status , err)
        !
        if (rank.eq.dest1)                                       &
          call mpi_sendrecv(up(n1,:, 1:2 ,-1:nymax), county1, mpi_real_kind, & 
                                                              source,0, &
                            up(n1,:,-1:0 ,-1:nymax), county1, mpi_real_kind, &
                                                              source,0, &
                            mpi_comm_world, status , err)
        !--------------------------------------------------------------
     case ( 5)    ! RIGHT @ higher resolution 1
        nbs    = extbounds(nb,1)
        source = extbounds(nb,2)
        n1   = extbounds(nb,3)
        dest1  = extbounds(nb,4)
        !
        !   own's right
        if (rank.eq.source) then
           !   fill send buffer
           do ieq=1,neq 
              do j=0,ny2p1
                 bx1l(1,ieq,1,j ) = up(nbs,ieq,nx,j )
              end do
           end do
           !   send/receive
           call mpi_sendrecv(bx1l, county3, mpi_real_kind, dest1 ,0, &
                             bx1s, county2, mpi_real_kind, dest1, 0, &
                             mpi_comm_world, status , err)
           !   fill u's with received buffer
           do i=nxp1,nxmax
              i1=i-nx
              do j=0,ny2
                 do ieq=1,neq
                    up(nbs,ieq,i,j )= bx1s(1,ieq,i1,j )
                 end do
              end do
           end do
           do ieq=1,neq
              up(nbs,ieq,nxp1 ,-1 )= up(nbs,ieq, nxp1,  0 )
              up(nbs,ieq,nxmax,-1 )= up(nbs,ieq, nxmax, 0 )
           end do
        endif
        !
        !   neighbor's left
        if (rank.eq.dest1)  then
           !   fill send buffer
           do j=0,ny2
              j1=2*j-1
              j2=j1+1
              do i=1,2
                 i1=2*i-1
                 i2 = i1+1
                 do ieq=1,neq
                    bx1s(1,ieq,i,j)=0.25*sum( up(n1,ieq,i1:i2,j1:j2) )
                 end do
              end do
           end do
           !   send/receive
           call mpi_sendrecv(bx1s, county2, mpi_real_kind, source,0, &
                             bx1l, county3, mpi_real_kind, source,0, &
                             mpi_comm_world, status , err)
           !   fill u's with received buffer
           do i=-1,0
              do j=-1,nymax
                 do ieq=1,neq
                    up(n1,ieq,i,j)=bx1l(1,ieq,1,(j+1)/2)
                 end do
              end do
           end do
        end if
        !--------------------------------------------------------------
     case ( 6)    ! RIGHT @ higher resolution 2
        nbs    = extbounds(nb,1)
        source = extbounds(nb,2)
        n2   = extbounds(nb,3)
        dest2  = extbounds(nb,4)
        !
        !   own's right
        if (rank.eq.source) then
           !   fill send buffer
           do ieq=1,neq
              do j=ny2,nyp1
                 bx2l(1,ieq,1,j) = up(nbs,ieq,nx,j)
              end do
           end do
           !   send/receive
           call mpi_sendrecv(bx2l, county3, mpi_real_kind, dest2 ,0, &
                             bx2s, county2, mpi_real_kind, dest2, 0, &
                             mpi_comm_world, status , err)
           !   fill u's with received buffer
           do i=nxp1,nxmax
              i1=i-nx
              do j=ny2p1,nyp1
                 do ieq=1,neq
                    up(nbs,ieq,i,j )= bx2s(1,ieq,i1,j )
                 end do
              end do
           end do
           do ieq=1,neq
              up(nbs,ieq,nxp1 ,nymax )= up(nbs,ieq, nxp1,  nyp1 )
              up(nbs,ieq,nxmax,nymax )= up(nbs,ieq, nxmax, nyp1 )
           end do
        endif
        !
        !   neighbor's left
        if (rank.eq.dest2)  then
           !   fill send buffer
           do j=1,ny2p1
              j1=2*j-1
              j2=j1+1
              ny2j=ny2+j
              do i=1,2
                 i1=2*i-1
                 i2 = i1+1
                 do ieq=1,neq
                    bx2s(1,ieq,i,ny2j)=0.25*sum( up(n2,ieq,i1:i2,j1:j2) )
                 end do
              end do
           end do
           !   send/receive
           call mpi_sendrecv(bx2s, county2, mpi_real_kind, source,0, &
                             bx2l, county3, mpi_real_kind, source,0, &
                             mpi_comm_world, status , err)
           !   fill u's with received buffer
           do i=-1,0
              do j=-1,nymax
                 do ieq=1,neq
                    up(n2,ieq,i,j)= bx2l(1,ieq,1,(j+1)/2+ny2)
                 end do
              end do
           end do
        end if
         !--------------------------------------------------------------
     case ( 7)    ! BOTTOM @ same resolution
        nbs    = extbounds(nb,1)
        source = extbounds(nb,2)
        n1     = extbounds(nb,3)
        dest1  = extbounds(nb,4)
        !
        if (rank.eq.source)                                      &
          call mpi_sendrecv(up(nbs,:,-1:nxmax, 1:2 ), countx1, mpi_real_kind, &
                                                      dest1 ,0,          &
                            up(nbs,:,-1:nxmax,-1:0 ), countx1, mpi_real_kind, &
                                                      dest1, 0,          &
                            mpi_comm_world, status , err)
        !
        if (rank.eq.dest1)                                       &
          call mpi_sendrecv(up(n1,:,-1:nxmax,nym1:ny   ), countx1, mpi_real_kind, &
                                                          source,0, &
                            up(n1,:,-1:nxmax,nyp1:nymax), countx1, mpi_real_kind, &
                                                          source,0, &
                            mpi_comm_world, status , err)
        !--------------------------------------------------------------
     case ( 8)    ! BOTTOM @ higher resolution 1
        nbs    = extbounds(nb,1)
        source = extbounds(nb,2)
        n1     = extbounds(nb,3)
        dest1  = extbounds(nb,4)
        !
        !   own's bottom
        if (rank.eq.source) then
           !   fill send buffer
           do ieq=1,neq 
              do i=0,nx2p1
                 by1l(1,ieq,i ,1) = up(nbs,ieq,i ,1)
              end do
           end do
           !   send/receive
           call mpi_sendrecv(by1l, countx3, mpi_real_kind, dest1 ,0, &
                             by1s, countx2, mpi_real_kind, dest1, 0, &
                             mpi_comm_world, status , err)
           !   fill u's with received buffer
           do j=-1,0
              j1=j+2
              do i=0,nx2
                 do ieq=1,neq
                    up(nbs,ieq,i ,j)= by1s(1,ieq,i ,j1)
                 end do
              end do
           end do
           do ieq=1,neq
              up(nbs,ieq,-1,  0)= up(nbs,ieq, 0, 0 )
              up(nbs,ieq,-1, -1)= up(nbs,ieq, 0,-1 )
           end do
        endif
        !
        !   neighbor's top
        if (rank.eq.dest1)  then
           !   fill send buffer
           do j=-1,0
              j1=2*j-1+ny
              j2 = j1+1
              do i=0,nx2
                 i1=2*i-1
                 i2=i1+1
                 do ieq=1,neq 
                    by1s(1,ieq,i,j+2)=0.25*sum( up(n1,ieq,i1:i2, j1:j2) )
                 end do
              end do
           end do
           !   send/receive
           call mpi_sendrecv(by1s, countx2, mpi_real_kind, source,0, &
                             by1l, countx3, mpi_real_kind, source,0, &
                             mpi_comm_world, status , err)
           !   fill u's with received buffer
           do j=nyp1,nymax
              do i=-1,nxmax
                 do ieq=1,neq
                    up(n1,ieq,i,j) = by1l(1,ieq,(i+1)/2,1)
                 end do
              end do
           end do
        end if
        !--------------------------------------------------------------
     case ( 9)    ! BOTTOM @ higher resolution 2
        nbs    = extbounds(nb,1)
        source = extbounds(nb,2)
        n2   = extbounds(nb,3)
        dest2  = extbounds(nb,4)
        !
        !   own's bottom
        if (rank.eq.source) then
           !   fill send buffer
           do ieq=1,neq 
              do i=nx2,nxp1
                 by2l(1,ieq,i,1) = up(nbs,ieq,i,1)
              end do
           end do
           !   send/receive
           call mpi_sendrecv(by2l, countx3, mpi_real_kind, dest2 ,0, &
                             by2s, countx2, mpi_real_kind, dest2, 0, &
                             mpi_comm_world, status , err)
           !   fill u's with received buffer
           do j=-1,0
              j1=j+2
              do i=nx2p1,nxp1
                 do ieq=1,neq
                    up(nbs,ieq,i,j)=by2s(1,ieq,i,j1)
                 end do
              end do
           end do
           do ieq=1,neq
              up(nbs,ieq,-1,nyp1  )= up(nbs,ieq,0,nyp1  )  
              up(nbs,ieq,-1,nymax )= up(nbs,ieq,0,nymax )
           end do
        endif
        !
        !   neighbor's top
        if (rank.eq.dest2)  then
           !   fill send buffer
           do j=-1,0
              j1=2*j-1+ny
              j2 = j1+1
              do i=1,nx2p1
                 i1=2*i-1
                 i2=i1+1
                 nx2i=nx2+i
                 do ieq=1,neq
                    by2s(1,ieq,nx2i,j+2)=0.25*sum( up(n2,ieq,i1:i2,j1:j2) )
                 end do
              end do
           end do
           !   send/receive
           call mpi_sendrecv(by2s, countx2, mpi_real_kind, source,0, &
                             by2l, countx3, mpi_real_kind, source,0, &
                             mpi_comm_world, status , err)
           !   fill u's with received buffer
           do j=nyp1,nymax
              do i=-1,nxmax
                 do ieq=1,neq
                    up(n2,ieq,i,j)= by2l(1,ieq,(i+1)/2+nx2,1)
                 end do
              end do
           end do
        end if
        !--------------------------------------------------------------
     case (10)    ! TOP @ same resolution
        nbs    = extbounds(nb,1)
        source = extbounds(nb,2)
        n1   = extbounds(nb,3)
        dest1  = extbounds(nb,4)
        !
        if (rank.eq.source)                                      &
          call mpi_sendrecv(up(nbs,:,-1:nxmax,nym1:ny   ), countx1, mpi_real_kind, &
                                                                    dest1 ,0, &
                            up(nbs,:,-1:nxmax,nyp1:nymax), countx1, mpi_real_kind, &
                                                                    dest1, 0, &
                            mpi_comm_world, status , err)
        !
        if (rank.eq.dest1)                                       &
          call mpi_sendrecv(up(n1,:,-1:nxmax, 1:2), countx1, mpi_real_kind, &
                                                    source,0, &
                            up(n1,:,-1:nxmax,-1:0), countx1, mpi_real_kind, &
                                                    source,0, &
                            mpi_comm_world, status, err)
        !--------------------------------------------------------------
     case (11)    ! TOP @ higher resolution 1 
        nbs    = extbounds(nb,1)
        source = extbounds(nb,2)
        n1     = extbounds(nb,3)
        dest1  = extbounds(nb,4)
        !
        !   own's top
        if (rank.eq.source) then
           !   fill send buffer
           do ieq=1,neq 
              do i=0,nx2p1
                 by1l(1,ieq,i ,1) = up(nbs,ieq,i ,ny)
              end do
           end do
           !   send/receive
           call mpi_sendrecv(by1l, countx3, mpi_real_kind, dest1 ,0, &
                             by1s, countx2, mpi_real_kind, dest1, 0, &
                             mpi_comm_world, status , err)
           !   fill u's with received buffer
           do j=nyp1,nymax
              j1=j-ny
              do i=0,nx2
                 do ieq=1,neq
                    up(nbs,ieq,i ,j )= by1s(1,ieq,i ,j1)
                 end do
              end do
           end do
           do ieq=1,neq
              up(nbs,ieq,-1,nyp1  )= up(nbs,ieq, 0, nyp1  )
              up(nbs,ieq,-1,nymax )= up(nbs,ieq, 0, nymax )
           end do
        endif
        !
        !   neighbor's bottom
        if (rank.eq.dest1)  then
           !   fill send buffer
           do j=1,2
              j1=2*j-1
              j2 = j1+1
              do i=0,nx2
                 i1=2*i-1
                 i2=i1+1
                 do ieq=1,neq 
                    by1s(1,ieq,i,j)=0.25*sum( up(n1,ieq,i1:i2,j1:j2) )
                 end do
              end do
           end do
           !   send/receive
           call mpi_sendrecv(by1s, countx2, mpi_real_kind, source,0, &
                             by1l, countx3, mpi_real_kind, source,0, &
                             mpi_comm_world, status , err)
           !   fill u's with received buffer
           do i=-1,nxmax
              do j=-1,0
                 do ieq=1,neq
                    up(n1,ieq,i,j) =  by1l(1,ieq,(i+1)/2,1)
                 end do
              end do
           end do
        end if
        !--------------------------------------------------------------
     case (12)    ! TOP @ higher resolution 2
        nbs    = extbounds(nb,1)
        source = extbounds(nb,2)
        n2   = extbounds(nb,3)
        dest2  = extbounds(nb,4)
        !
        !   own's top
        if (rank.eq.source) then
           !   fill send buffer   
           do ieq=1,neq 
              do i=nx2,nxp1
                 by2l(1,ieq,i,1) = up(nbs,ieq,i,ny)
              end do
           end do
           !   send/receive
           call mpi_sendrecv(by2l, countx3, mpi_real_kind, dest2 ,0, &
                             by2s, countx2, mpi_real_kind, dest2, 0, &
                             mpi_comm_world, status , err)
           !   fill u's with received buffer
           do j=nyp1,nymax
              j1=j-ny
              do i=nx2p1,nxp1
                 do ieq=1,neq
                    up(nbs,ieq,i,j)=by2s(1,ieq,i,j1)
                 end do
              end do
           end do
           do ieq=1,neq
              up(nbs,ieq,nxmax,nyp1  )= up(nbs,ieq, nxp1, nyp1  )
              up(nbs,ieq,nxmax,nymax )= up(nbs,ieq, nxp1, nymax )
           end do
        endif
        !
        !   neighbor's bottom
        if (rank.eq.dest2)  then
           !   fill send buffer
           do j=1,2
              j1=2*j-1
              j2 = j1+1
              do i=1,nx2p1
                 i1=2*i-1
                 i2=i1+1
                 nx2i=nx2+i
                 do ieq=1,neq
                    by2s(1,ieq,nx2i,j)=0.25*sum( up(n2,ieq,i1:i2,j1:j2) )
                 end do
              end do
           end do
           !   send/receive
           call mpi_sendrecv(by2s, countx2, mpi_real_kind, source,0, &
                             by2l, countx3, mpi_real_kind, source,0, &
                             mpi_comm_world, status , err)
           !   fill u's with received buffer
           do j=-1,0
              do i=-1,nxmax
                 do ieq=1,neq
                    up(n2,ieq,i,j)=by2l(1,ieq,(i+1)/2+nx2,1)
                 end do
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
end subroutine boundaryII
!======================================================================
