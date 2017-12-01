!=======================================================================
!   Thuis subroutine generates the initial grid
!=======================================================================
subroutine basegrid(itprint)
  use parameters
  use globals
  use user_mod
  implicit none
  integer, intent(inout)  :: itprint
  !
  logical,dimension(nblocks*np) :: irefup,irefdown
  integer :: err,unit, nb,nl,nb2, nl2,nbr,i,j,k, rankorig
  integer :: nbminp,nbmaxp,nbleafsp,unitin,nbstep,nblocksorig
  character (len=128) file1,file2
  !
  if(iwarm.eq.0) then
    !------------1st level-----------------------
    ! sibling numbering follows the convention:
    !                       Y
    !       ---------       ^
    !       |ID3|ID4|       |
    !       ---------       |
    !       |ID1|ID2|       +---> X
    !       ---------
    !   initialized root blocks
    lp(1,1)=1           ! level of refinement
    lp(1,2)=0           ! father
    lp(1,3)=0           ! son (only 1st=ID1)
    lp(1,4)=1           ! my sibling ID
    lp(1,5)=-1          ! neighbor left
    lp(1,6)= 2          ! neighbor right
    lp(1,7)=-1          ! neighbor down
    lp(1,8)=-1          ! neighbor up
    !                     -1 = grid boundary
    icoord(1,1)=0
    icoord(1,2)=0
    leaf(1)=.true.
    call initial_conditions(u(1,:,:,:), 1, 0.) ! nb=1, time=0.
    call calcprim(u(1,:,:,:), primit(1,:,:,:),                        &
                  neq, nxmin, nxmax, nymin, nymax)
    irefup(1)=.true.
    !---
    lp(2,1)=1           ! level of refinement
    lp(2,2)=0           ! father
    lp(2,3)=0           ! son (only 1st=ID1)
    lp(2,4)=1           ! my sibling ID
    lp(2,5)= 1          ! neighbor left
    lp(2,6)= 3          ! neighbor right
    lp(2,7)=-1          ! neighbor down
    lp(2,8)=-1          ! neighbor up
    !                     -1 = grid boundary
    icoord(2,1)=nx
    icoord(2,2)=0
    leaf(2)=.true.
    call initial_conditions(u(2,:,:,:), 2, 0.) 
    call calcprim(u(2,:,:,:), primit(2,:,:,:),                        &
                 neq, nxmin, nxmax, nymin, nymax )
    irefup(2)=.true.
    !---
    lp(3,1)=1           ! level of refinement
    lp(3,2)=0           ! father
    lp(3,3)=0           ! son (only 1st=ID1)
    lp(3,4)=1           ! my sibling ID
    lp(3,5)= 2          ! neighbor left
    lp(3,6)= 4          ! neighbor right
    lp(3,7)=-1          ! neighbor down
    lp(3,8)=-1          ! neighbor up
    !                     -1 = grid boundary
    icoord(3,1)=2*nx
    icoord(3,2)=0
    leaf(3)=.true.
    call initial_conditions(u(3,:,:,:), 3, 0.) 
    call calcprim(u(3,:,:,:), primit(3,:,:,:),                        &
                 neq, nxmin, nxmax, nymin, nymax )
    irefup(3)=.true.
    !---
    lp(4,1)=1           ! level of refinement
    lp(4,2)=0           ! father
    lp(4,3)=0           ! son (only 1st=ID1)
    lp(4,4)=1           ! my sibling ID
    lp(4,5)= 3          ! neighbor left
    lp(4,6)=-1          ! neighbor right
    lp(4,7)=-1          ! neighbor down
    lp(4,8)=-1          ! neighbor up
    !                     -1 = grid boundary
    icoord(4,1)=3*nx
    icoord(4,2)=0
    leaf(4)=.true.
    call initial_conditions(u(4,:,:,:), 4, 0.) 
    call calcprim(u(4,:,:,:), primit(4,:,:,:),                        &
                 neq, nxmin, nxmax, nymin, nymax )
    irefup(4)=.true.
    !---
    !
    !   The resulting number of blocks and leafs
    nbmax  = 4
    nbleafs= 4
    !
    do nbr=1,nbroot
      lp(nbr,3)=nbmax+1
      call refine(    nbr,nbmax+1,nbmax+2,nbmax+3,nbmax+4)
      call updatelpup(nbr,nbmax+1,nbmax+2,nbmax+3,nbmax+4, &
                       irefup,irefdown)
      nbmax=nbmax+4
     end do
     !
     !   apply initial conditions to level 2 sons
     do nb=5,20
        call initial_conditions(u(nb,:,:,:), nb, 0.)
        call calcprim(u(nb,:,:,:), primit(nb,:,:,:),                   &
                      neq, nxmin, nxmax, nymin, nymax ) 
     end do
     !
     print*,rank,'---------'
     !--------------------------------------------
     do nl=2,nlevs-1
        !   mark for refinement
        do nb=3,nbmax
           irefup(nb)=.false.
           if( leaf(nb) .and. lp(nb,1).eq.nl ) then
              !  Force refinement in the lower left corner 
              !if ( icoord(nb,1) <= 150 )  then  !   100 pixels in X
                 if ( icoord(nb,2) <= 100   ) then  !   100 pixels in y
                    irefup(nb)=.true.
                    irefdown(nb)=.false.
                 endif
              !endif
           end if
        end do
        !
        !   mark for refining on proximity of higher resolution grid
        call critneighup0
        !..............................................................
        !
        !   refine marked blocks and apply initial conditions
        do nb=nbmin,nbmax
           !if (lp(nb,1).lt.4) then
           !   irefup(nb)=.true.
           !   irefdown(nb)=.false.
           !end if
           !
           if (irefup(nb)) then
              !print*,'********************************************',nb
              !   to do it in order of levels
              do nl2=2,nlevs-1
                 if (lp(nb,1).eq.nl2) then
                    !print'(a,i6,2i6,a,2i5)','refining block', nb           &
                    !     ,icoord(nb,:),':',nbmax+1,nbmax+4
                    lp(nb,3)=nbmax+1
                    call refine(nb,nbmax+1,nbmax+2,nbmax+3,nbmax+4)
                    call updatelpup(nb,nbmax+1,nbmax+2,nbmax+3,nbmax+4,&
                                    irefup,irefdown)
                    nbmax=nbmax+4
                    !   apply init conditions
                    do nb2=lp(nb,3),lp(nb,3)+3
                      call initial_conditions(u(nb2,:,:,:), nb2, 0.)
                      call calcprim(u(nb2,:,:,:), primit(nb2,:,:,:),   &
                                      neq, nxmin, nxmax, nymin, nymax ) 
                    end do
                 endif
              end do
              !
           end if
        end do
     end do
     !
  else
     !   warm start
#ifdef MPIP
     write(file2,'(a,i3.3,a,i3.3,a)')                                 &
           trim(outputpath)//'BIN/blocks'  ,rank,'.',itprint,'.bin'
     unit=10+rank
#else
     write(file1,'(a,i3.3,a)') trim(outputpath)//'BIN/pointers',itprint,'.bin'
     write(file2,'(a,i3.3,a)') trim(outputpath)//'BIN/blocks'  ,itprint,'.bin'
     unit=10
#endif        
     !
     !   If MPI Take turns to read all the pointer files
#ifdef MPIP
     do i=0,np-1
        if (rank.eq.i) then
           do j=0,np-1
              write(file1,'(a,i3.3,a,i3.3,a)')                        &
                    trim(outputpath)//'BIN/pointers',j,'.',itprint,'.bin'
#endif
              open(unit=unitin,file=file1,status='unknown',form='unformatted', &
                   convert='LITTLE_ENDIAN',recordtype='stream',access='sequential')       
              read(unitin) nbminp,nbmaxp,nbleafsp,nblocksorig
              nbstep=(nblocks)*j+1-nbminp
              !
              read(unitin) lp(nbminp+nbstep:nbmaxp+nbstep,:)
              !
              if (nblocks.ne.nblocksorig) then
                 !   fix the pointers if the number of blocks
                 !   has changed from the last run
                 do nb=nbminp+nbstep,nbmaxp+nbstep 
                    if (maxval(lp(nb,:)).gt.nblocksorig) then
                       !print'(8i8)',lp(nb,:)
                       !print'(8i8)',lp(nb,:)/nblocksorig
                       do k=2,8
                          rankorig=(lp(nb,k)/nblocksorig)
                          if (rankorig .gt. 0) then
                             
                             lp(nb,k)= lp(nb,k) -rankorig*nblocksorig+ &
                                       (nblocks)*rankorig
                          end if
                       end do
                       !print'(8i8)',lp(nb,:)
                       !print*,'--'
                    end if
                 end do
                 !if (j==0) stop
              endif
              read(unitin) leaf(nbminp+nbstep:nbmaxp+nbstep)
              !
#ifdef MPIP
              if (rank.eq.j) then
#endif
                 nbleafs = nbleafsp
                 nbmax=nbmaxp+nbstep
                 !print'(5i8)',rank,nbmin,nbmax,nbminp,nbmaxp
                 read(unitin) icoord(nbmin:nbmax,:)
#ifdef MPIP
              end if
#endif
              close(unitin)
#ifdef MPIP
           end do
        end if
        call mpi_barrier(mpi_comm_world,err)
     end do
#endif
     !
     !   points (only leaf blocks)
     open(unit=unitin,file=file2,status='unknown',form='unformatted',  &
          convert='LITTLE_ENDIAN',recordtype='stream',access='sequential')
     do nb=nbmin,nbmax
        if(leaf(nb)) then
           read(unitin) u(nb,:,:,:)
           call calcprim(u(nb,:,:,:), primit(nb,:,:,:),                &
                         neq, nxmin, nxmax, nymin, nymax )

        end if
     end do
     close(unitin)
     !
     print*,rank,' just read ',trim(file1),' & ' ,trim(file2)
     itprint=itprint+1
     !
  endif
  !
contains
  !=======================================================================
  !  refinement criterion to allow maximum difference of 1 level between
  !  neighboring blocks when refining (i.e. only with irefup)
  !=======================================================================
  subroutine critneighup0
    use parameters
    use globals
    implicit none
    integer :: nb,n1
    !
    do nb=nbmin,nbmax
       if (irefup(nb)) then
          !---------------------------------------------------
          !   Left
          n1=lp(nb,5)
          if (n1.ne.-1) then
             if ( lp(n1,1).lt.lp(nb,1) ) then
                irefup(n1)=.true.
                if (lp(lp(n1,5),1).lt.lp(n1,1) ) irefup(lp(n1,5))=.true.
                if (lp(lp(n1,6),1).lt.lp(n1,1) ) irefup(lp(n1,6))=.true.
                if (lp(lp(n1,7),1).lt.lp(n1,1) ) irefup(lp(n1,7))=.true.
                if (lp(lp(n1,8),1).lt.lp(n1,1) ) irefup(lp(n1,8))=.true.   
             end if
          end if
          !   Right
          n1=lp(nb,6)
          if (n1.ne.-1) then
             if ( lp(n1,1).lt.lp(nb,1) ) then
                irefup(n1)=.true.
                if (lp(lp(n1,5),1).lt.lp(n1,1) ) irefup(lp(n1,5))=.true.
                if (lp(lp(n1,6),1).lt.lp(n1,1) ) irefup(lp(n1,6))=.true.
                if (lp(lp(n1,7),1).lt.lp(n1,1) ) irefup(lp(n1,7))=.true.
                if (lp(lp(n1,8),1).lt.lp(n1,1) ) irefup(lp(n1,8))=.true.   
             end if
          end if
          !   Bottom
          n1=lp(nb,7)
          if (n1.ne.-1) then
             if ( lp(n1,1).lt.lp(nb,1) ) then
                irefup(n1)=.true.
                if (lp(lp(n1,5),1).lt.lp(n1,1) ) irefup(lp(n1,5))=.true.
                if (lp(lp(n1,6),1).lt.lp(n1,1) ) irefup(lp(n1,6))=.true.
                if (lp(lp(n1,7),1).lt.lp(n1,1) ) irefup(lp(n1,7))=.true.
                if (lp(lp(n1,8),1).lt.lp(n1,1) ) irefup(lp(n1,8))=.true.
             end if
          end if
          !   Top
          n1=lp(nb,8)
          if (n1.ne.-1) then
             if ( lp(n1,1).lt.lp(nb,1) ) then
                irefup(n1)=.true.
                if (lp(lp(n1,5),1).lt.lp(n1,1) ) irefup(lp(n1,5))=.true.
                if (lp(lp(n1,6),1).lt.lp(n1,1) ) irefup(lp(n1,6))=.true.
                if (lp(lp(n1,7),1).lt.lp(n1,1) ) irefup(lp(n1,7))=.true.
                if (lp(lp(n1,8),1).lt.lp(n1,1) ) irefup(lp(n1,8))=.true.   
             end if
          end if
          !---------------------------------------------------
       end if
    end do
    !
  end subroutine critneighup0
  !====================================================================
  !
end subroutine basegrid
!======================================================================

