!========================================================================
!   redistributes the load among all available processors, according to a 
!   Hilbert ordering
!========================================================================
subroutine loadbalance
#ifdef MPIP
  use parameters
  use globals
  implicit none
  integer :: err, status(MPI_STATUS_SIZE)
  integer :: nb,i, j, leafspp,nbtot, ilevs,ip, nbleafstot
  integer :: mynb,nb2,source,dest,count
  integer, dimension(0:np-1) :: recvcnts,displs
  integer, dimension(:),   allocatable :: xh,yh,indexH,idnew,idold
  integer, dimension(:),   allocatable :: hilbert,hilbertp
  integer, dimension(:,:), allocatable :: icoordp,lpp
  logical, dimension(:),   allocatable :: leafp
  !
  allocate(hilbert(nbmin:nbmax))
  allocate(xh(nbmin:nbmax))
  allocate(yh(nbmin:nbmax))
  !
  !  identify  coordinate of each block
  do nb=nbmin,nbmax
     xh(nb)=icoord(nb,1)*2**(nlevs-lp(nb,1))/nx
     yh(nb)=icoord(nb,2)*2**(nlevs-lp(nb,1))/ny
  end do
  !
  ! number of active blocks in each processor
  mynb=nbmax-nbmin+1
  !
  !   assign a Hilbert number to each block
  call  hilbert2d( xh,yh,hilbert,nlevs+1,mynb)
  !
  !   add the leafs and total blocks of all processors
  call mpi_allreduce(nbleafs, nbleafstot, 1, mpi_integer, mpi_sum,     &
       mpi_comm_world,err)
  leafspp=nbleafstot/np
  !
  call mpi_allreduce(mynb, nbtot, 1, mpi_integer, mpi_sum,             &
       mpi_comm_world,err)
  !
  !   consolidate list of Hilbert numbers
  call mpi_allgather(mynb, 1, mpi_integer, recvcnts, 1, mpi_integer, &
       mpi_comm_world,err)
  displs(0)=0
  do i=1,np-1
     displs(i)=displs(i-1)+recvcnts(i-1)
  end do
  !
  allocate(indexH(nbtot))
  allocate(idnew(nbtot))
  allocate(idold(nbtot))
  allocate(hilbertp(nbtot))
  !
  call mpi_gatherv(hilbert ,mynb, mpi_integer, hilbertp,recvcnts,    &
       displs ,mpi_integer, master, mpi_comm_world,err)
  !
  !print'(7i8)',rank,mynb,nbmin,nbmax,nbleafs,nbtot,nbleafstot!,recvcnts,displs
  !
  !   clear up's to use them as temporary buffer
  up(:,:,:,:)=0
  allocate(    lpp(nblocks*np,8))
  allocate(icoordp(nbmin:nbmaxtot,2))
  allocate(  leafp(nblocks*np  ))
  icoordp(:,:)=0
  lpp(:,:)=0
  leafp(:)=.false.
  !
  !----------------------------------------
  !   use index table to assign new block #
  !   and update pointers
  IF (rank.eq.master)  THEN
     !
     !   create an index table  based on the Hilbert number
     call indexx(nbtot,hilbertp,indexH)
     !
     j=1
     do ip=0,np-1
        do nb=(ip)*nblocks+1,recvcnts(ip)+(ip)*nblocks
           idold(J)=nb
           j=j+1
        end do
     end do
     !
     i=1
     j=1
     nb2=0
     idnew(1)=1
     !
     do nb=1,nbtot
        !
        idnew(indexH(nb))=(i-1)*nblocks+nb-nb2
        if ( leaf( idold(indexH(nb)) ) ) j=j+1
        if (j.gt.(leafspp).and.i.lt.np) then
           i=i+1
           j=1
           nb2=nb
        end if
        if (idnew(indexH(nb))/i .eq. nblocks) then
           print'(a,3i3)','Error: reached maximum number of blocks for rank ', i
           stop
        end if
        !
     end do
     !
     do nb=1,nbtot
        !
        !   update the leaf tag 
        leafp(idnew(nb)) = leaf(idold(nb))
        !
        !   update level and ID pointers
        lpp(idnew(nb),1) = lp(idold(nb),1)
        lpp(idnew(nb),4) = lp(idold(nb),4)
        !
        !   father and sons
        if(lp(idold(nb),1).eq.1) then
           !   root blocks
           !   father
           lpp(idnew(nb),2)=0
           !   son
           lpp(idnew(nb),3)=idnew(locate(idold,nbtot,lp(idold(nb),3)))
        else
           !   the rest
           !   father
           lpp(idnew(nb),2)= idnew(locate(idold,nbtot,lp(idold(nb),2)))
           !   son
           lpp(idnew(nb),3)= idnew(locate(idold,nbtot,lp(idold(nb),3)))
            !
        end if
        !
        !   neighbors
        !   left
        if(lp(idold(nb),5).ne.-1) then
           lpp(idnew(nb),5)= idnew(locate(idold,nbtot,lp(idold(nb),5)))
        else
           lpp(idnew(nb),5)= -1
        end if
        !   right
        if(lp(idold(nb),6).ne.-1) then
           lpp(idnew(nb),6)= idnew(locate(idold,nbtot,lp(idold(nb),6)))
        else
           lpp(idnew(nb),6)= -1
        end if
        !   bottom
        if(lp(idold(nb),7).ne.-1) then
           lpp(idnew(nb),7)= idnew(locate(idold,nbtot,lp(idold(nb),7)))
        else
           lpp(idnew(nb),7)= -1
        end if
        !   top
        if(lp(idold(nb),8).ne.-1) then
           lpp(idnew(nb),8)= idnew(locate(idold,nbtot,lp(idold(nb),8)))
        else
           lpp(idnew(nb),8)= -1
        end if
        !
        !if (leafp(idnew(nb))) then
           !print'(2i7,a,i7,a,8i7,a,3i7)',nb,idold(nb),'---> *',idnew(nb) &
           !     ,':',lpp(idnew(nb),:)!,';',lp(idold(nb),3),locate(idold,nbtot,lp(idold(nb),3))
        !else
        !   print'(2i7,a,i7,a,8i7,a,3i7)',nb,idold(nb),'--->  ',idnew(nb) &
        !        ,':',lpp(idnew(nb),:)!,';',lp(idold(nb),3),locate(idold,nbtot,lp(idold(nb),3))
        !end if
        !
     end do
     !
  END IF
  !
  call mpi_bcast(idnew,nbtot       ,mpi_integer,0,mpi_comm_world, err)
  call mpi_bcast(idold,nbtot       ,mpi_integer,0,mpi_comm_world, err)
  call mpi_bcast(lpp,np*nblocks*8  ,mpi_integer,0,mpi_comm_world, err)
  call mpi_bcast(leafp,np*nblocks  ,mpi_logical,0,mpi_comm_world, err)
  lp=lpp
  leaf=leafp
  !
  !   update nbmax
  nbmax=nbmin-1
  nbleafs=0
  do nb=1,nbtot
     if(idnew(nb)/nblocks.eq.rank) nbmax=nbmax+1
  end do
  !---------------------------------------------------------------------
  !     redistribute the workload
  do nb=1,nbtot
     !
     source=idold(nb)/nblocks
     dest=idnew(nb)/nblocks
     if(dest.ge.np) dest=np-1
     !
     if (source.eq.dest.and.rank.eq.source)  then
        up(idnew(nb),:,:,: ) =      u(idold(nb),:,:,:)
        icoordp(idnew(nb),:     ) = icoord(idold(nb),:    )
     else
        !   u -> up   --------------------
        if(rank.eq.source) then
           count=neq*(nxmax-nymin+1)*(nxmax-nymin+1)
           call mpi_send(  u(idold(nb)     ,:,:,:), count, mpi_real_kind    &
                , dest, 0, mpi_comm_world,err)
        end if
        if(rank.eq.dest  ) then
           count=neq*(nxmax-nymin+1)*(nxmax-nymin+1)
           call mpi_recv( up(idnew(nb),:,:,:), count, mpi_real_kind   &
                , source,0, mpi_comm_world, status, err)
        end if
        !   icoord -> icoordp   ----------
        if(rank.eq.source) then
           call mpi_send(  icoord(idold(nb)     ,:), 2 , mpi_integer &
                , dest, 0, mpi_comm_world,err)
        end if
        if(rank.eq.dest  ) then
           call mpi_recv( icoordp(idnew(nb),:), 2 , mpi_integer &
                , source,0, mpi_comm_world, status, err)
        end if
        !---------------------------------
     end if
     !
     call mpi_barrier(mpi_comm_world,err)
     !
  end do
  !---------------------------------------------------------------------
  !   update arrays with the temporary buffers
  u=up
  icoord=icoordp
  !
  !   update the leafs
  nbleafs=0
  do nb=nbmin,nbmax
     if( leaf(nb) ) then 
        nbleafs=nbleafs+1
        call calcprim(u(nb,:,:,:), primit(nb,:,:,:),   &
                      neq, nxmin, nxmax, nymin, nymax ) 
     end if
  end do
  !
  !print'(i7,a,i7,a,i7,a,2i7)',rank,' nbleafs:', nbleafs, ' leafspp:'   &
  !     , leafspp,' bounds:',nbmin,nbmax
  !
  call mpi_barrier(mpi_comm_world,err)
  !
  deallocate(xh,yh,hilbert,hilbertp,icoordp,leafp,lpp,idnew,idold,indexH)
  !
contains
  !
!========================================================================
!   This routine assigns a Hilbert number to a 2D map
!========================================================================
subroutine hilbert2d(x,y,order,bit_length,npoint)
  implicit none
  integer, INTENT(IN)                     ::bit_length,npoint
  integer, INTENT(IN) ,dimension(1:npoint)::x,y
  integer, INTENT(OUT),dimension(1:npoint)::order
  !
  logical,dimension(0:2*bit_length-1)::i_bit_mask 
  logical,dimension(0:1*bit_length-1)::x_bit_mask,y_bit_mask
  integer,dimension(0:3,0:1,0:3)::state_diagram
  integer::i,ip,cstate,nstate,b0,b1,sdigit,hdigit

  if(bit_length>bit_size(bit_length))then
     write(*,*)'Maximum bit length=',bit_size(bit_length)
     write(*,*)'stop in hilbert2d'
     !call clean_stop
  endif
 
  state_diagram = RESHAPE( (/ 1, 0, 2, 0, &
                            & 0, 1, 3, 2, &
                            & 0, 3, 1, 1, &
                            & 0, 3, 1, 2, &
                            & 2, 2, 0, 3, &
                            & 2, 1, 3, 0, &
                            & 3, 1, 3, 2, &
                            & 2, 3, 1, 0  /), &
                            & (/ 4, 2, 4 /) )

  do ip=1,npoint

     ! convert to binary
     do i=0,bit_length-1
        x_bit_mask(i)=btest(x(ip),i)
        y_bit_mask(i)=btest(y(ip),i)
     enddo
     
     ! interleave bits
     do i=0,bit_length-1
        i_bit_mask(2*i+1)=x_bit_mask(i)
        i_bit_mask(2*i  )=y_bit_mask(i)
     end do
     
     ! build Hilbert ordering using state diagram
     cstate=0
     do i=bit_length-1,0,-1
        b1=0 ; if(i_bit_mask(2*i+1))b1=1
        b0=0 ; if(i_bit_mask(2*i)  )b0=1
        sdigit=b1*2+b0
        nstate=state_diagram(sdigit,0,cstate)
        hdigit=state_diagram(sdigit,1,cstate)
        i_bit_mask(2*i+1)=btest(hdigit,1)
        i_bit_mask(2*i  )=btest(hdigit,0)
        cstate=nstate
     enddo

     ! save Hilbert key as double precision real
     order(ip)=0.
     do i=0,2*bit_length-1
        b0=0 ; if(i_bit_mask(i))b0=1
        order(ip)=order(ip)+b0*2**i!+dble(b0)*dble(2)**i
     end do

  end do

end subroutine hilbert2d
!=======================================================================

!=======================================================================
!   this routine indexes an array arr(1:n), i.e. outputs
!   indx(1:n) such that arr(indx(j)) is in ascending order for
!   j=1, 2, 3 ...n.  arr and n are not changed
!   Taken from NR
!=======================================================================
subroutine indexx(n,arr,indx)
  implicit none
  integer, intent(in)   :: n
  integer, intent(in),  dimension(n) :: arr
  integer, intent(out), dimension(n) :: indx
  integer :: i,indxt,ir,itemp,j,jstack,k,l,a
  integer, parameter :: m=5, nstack=50
  integer, dimension(nstack) :: istack
  !
  do j=1,n
     indx(j)=j
  end do
  jstack=0
  l=1
  ir=n
  do
     if(ir-l.lt.m)then
        do j=l+1,ir
           indxt=indx(j)
           a=arr(indxt)
           do i=j-1,l,-1
              if(arr(indx(i)).le.a) exit
              indx(i+1)=indx(i)
           end do
           indx(i+1)=indxt
        end do
        if(jstack.eq.0) RETURN
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
     else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l)).gt.arr(indx(ir))) then
           itemp=indx(l)
           indx(l)=indx(ir)
           indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(ir))) then
           itemp=indx(l+1)
           indx(l+1)=indx(ir)
           indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(l+1))) then
           itemp=indx(l)
           indx(l)=indx(l+1)
           indx(l+1)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l+1)
        a=arr(indxt)
        do
           do
              i=i+1
              if(arr(indx(i)).ge.a) exit
           end do
           do
              j=j-1
              if(arr(indx(j)).le.a) exit
           end do
           if(j.lt.i) exit   
           itemp=indx(i)
           indx(i)=indx(j)
           indx(j)=itemp
        end do
        indx(l+1)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.nstack) pause 'NSTACK to small to indexx'
        if(ir-i+1.ge.j-l)then
           istack(jstack)=ir
           istack(jstack-1)=i
           ir=j-1
        else
           istack(jstack)=j-1
           istack(jstack-1)=l
           l=i
        end if
     end if
  end do
  !
end subroutine indexx
!=======================================================================

!=======================================================================
!   given an array xx(1:n), and a given value of x, returns a value j
!   such that x is between xx(j) and xx(j+1). xx(1:n) must be monotonic,
!   eithet increasing or decreasing. j=0 or j=n is returned to indicate
!   that x is out of range.
!   Taken from NR (modified to take arrays f integers, and the size n
!   explicitly.
!=======================================================================
function locate(xx,n,x)
  implicit none
  integer :: locate
  integer, intent(in), dimension(n) :: xx
  integer, intent(in):: x, n
  integer            :: jl,jm,ju
  logical            :: ascnd
  !
  ascnd = (xx(n) >= xx(1))
  jl=0
  ju=n+1
  !  
  do
     if(ju-jl <= 1) exit
     jm=(ju+jl)/2
     if (ascnd .eqv. (x >= xx(jm) ) ) then
        jl=jm
     else
        ju=jm
     end if
  end do
  if (x == xx(1)) then
     locate=1
  else if (x==xx(n)) then
     locate=n!-1
  else
     locate=jl
  end if
  !
end function locate
  !
#endif  
end subroutine loadbalance
!========================================================================
