!=======================================================================
!   writes the output to a file
!   The file format(s) is set in the makefile
!=======================================================================
subroutine output(itprint)
  use parameters
  use globals
  implicit none
  integer, intent(in) :: itprint
#ifdef MPIP
  integer :: err, status(MPI_STATUS_SIZE)
#endif
  character (len=128) file1, file2
  character (len=50) cbuffer
  character(1), parameter  :: lf = char(10) 
  real, dimension(neq) :: uu, prim
  real  :: T
  integer :: nb, npoints, i, j, k, n, unitout,ncells
  integer, dimension(:), allocatable :: nbl
  !
  allocate(nbl(nbleafs))
  !
  i=1
  do nb=nbmin,nbmax
     if(leaf(nb)) then
        nbl(i)=nb
        i=i+1
     end if
  end do
  !
  print*,'there are',nbleafs,' leafs, on rank ',rank, nbmax
  !print*,'with ID:',nbl
  !
  !
  !-----------------------------------------------------------------
  !   output *.bin (points y pointers that are used for warm start)
  !-----------------------------------------------------------------
  !
#ifdef OUTBIN
#ifdef MPIP
    write(file1,'(a,i3.3,a,i3.3,a)') trim(outputpath)//'BIN/pointers',rank,'.',itprint,'.bin'
    write(file2,'(a,i3.3,a,i3.3,a)') trim(outputpath)//'BIN/blocks'  ,rank,'.',itprint,'.bin'
  unitout=rank+10
#else
  write(file1,'(a,i3.3,a)') trim(outputpath)//'BIN/pointers',itprint,'.bin'
  write(file2,'(a,i3.3,a)') trim(outputpath)//'BIN/blocks'  ,itprint,'.bin'
  unitout=10
#endif
  !   mesh info (nbmin,nbmax,icoord,leaf & pointers)
  open(unit=unitout,file=file1,status='unknown',form='unformatted', &
       convert='LITTLE_ENDIAN',recordtype='stream',access='sequential')
  write(unitout) nbmin,nbmax,nbleafs,nblocks
  write(unitout) lp(nbmin:nbmax,:)
  write(unitout) leaf(nbmin:nbmax)
  write(unitout) icoord(nbmin:nbmax,:)
  close(unitout)
  !
  !   points (only leaf blocks)
  open(unit=unitout,file=file2,status='unknown',form='unformatted', &
       convert='LITTLE_ENDIAN',recordtype='stream',access='sequential')
  do nb=nbmin,nbmax
     if(leaf(nb)) then
        write(unitout) u(nb,:,:,:)
     end if
  end do
  close(unitout)
#endif
  !
  !-----------------------------------------------------------------
  !   output *.vtk
  !-----------------------------------------------------------------
  !
#ifdef OUTVTK
#ifdef MPIP
  !*****************************************************************
  !   write to .visit file to include several subsets
  if (rank.eq.0) then
     if (itprint.eq.0)then
        open(unit=7,file=trim(outputpath)//'.visit',status='unknown',form='formatted')
        write(7,'(a,i)'),'!NBLOCKS ',np
        do i=0,np-1
           write(7,'(a,i3.3,a,i3.3,a)') './VTK/0',i,'.',itprint,'.vtk'
        end do
        close(7)
     else
        open(unit=7,file=trim(outputpath)//'.visit',status='unknown',form='formatted'    &
             ,position='append')
        do i=0,np-1
           write(7,'(a,i3.3,a,i3.3,a)') './VTK/0',i,'.',itprint,'.vtk'
        end do
        close(7)
     end if
  end if
  !*****************************************************************
  write(file1,'(a,i3.3,a,i3.3,a)') trim(outputpath)//'VTK/0',rank,'.',itprint,'.vtk'
  unitout=rank+10
#else
  write(file1,'(a,i3.3,a)') trim(outputpath)//'VTK/out',itprint,'.vtk'
  unitout=10
#endif
  open(unit=unitout,file=file1,status='unknown',   &
           form='unformatted',recordtype='STREAM', &
           action='write',convert='LITTLE_ENDIAN',    &
           access='sequential' )
  !
  !   write the header
  npoints=nbleafs*(nx+1)*(ny+1)
  ncells= nbleafs*nx*ny
  !
  write(unitout) '# vtk DataFile Version 4.2 ',lf
  write(unitout) 'output from walicxe',lf
  write(unitout) 'BINARY',lf
  write(unitout) 'DATASET UNSTRUCTURED_GRID',lf
  write(cbuffer,'("POINTS", i, " float")'),npoints
  write(unitout) trim(cbuffer),lf
  !
  do nb=1,nbleafs
     do i=0,nx
        do j=0,ny
           write(unitout)                                             &
                ( float(icoord(nbl(nb),1)+i) )*dx(lp(nbl(nb),1))*rsc, &
                ( float(icoord(nbl(nb),2)+j) )*dy(lp(nbl(nb),1))*rsc,0.
        end do
     end do
  end do
  write(unitout) lf
  write(cbuffer,'("CELLS ",2(i,1x))'), NCELLS, ncells*5
  write(unitout) trim(cbuffer),lf
  !
  k=0
  do nb=1,nbleafs
     do i=0,nx
        do j=0,ny
           if(i.lt.nx.and.j.lt.ny) write(unitout) 4,k,k+1,k+ny+1,k+ny+2
           k=k+1
        end do
     end do
  end do
  write(unitout) lf
  write(cbuffer,fmt='(a,i)') 'CELL_TYPES', ncells
  write(unitout) trim(cbuffer),lf
  do nb=1,nbleafs
     do i=1,nx
        do j=1,ny
           write(unitout) 8
        end do
     end do
  end do
  write(unitout) lf
  write(cbuffer,'(a,i)')'CELL_DATA', ncells
  write(unitout) trim(cbuffer),lf
  !
  !   Density
  write(cbuffer,'(a)')  'SCALARS rho float 1'
  write(unitout) trim(cbuffer),lf
  write(cbuffer,'(a)')  'LOOKUP_TABLE default'
  write(unitout) trim(cbuffer),lf
  do nb=1,nbleafs
     do i=1,nx
        do j=1,ny
           write(unitout) u(nbl(nb),1,i,j)*rhosc
        end do
     end do
  end do
  write(unitout) lf
  !
  !   Gas Pressure
  write(cbuffer,'(a)')  'SCALARS P float 1'
  write(unitout) trim(cbuffer),lf
  write(cbuffer,'(a)')  'LOOKUP_TABLE default'
  write(unitout) trim(cbuffer),lf
  do nb=1,nbleafs
     do i=1,nx
        do j=1,ny
           write(unitout) primit(nbl(nb),4,i,j)
        end do
     end do
  end do
  write(unitout) lf
  !
  !   Density of neutrals
  write(cbuffer,'(a)')  'SCALARS rho_HI float 1'
  write(unitout) trim(cbuffer),lf
  write(cbuffer,'(a)')  'LOOKUP_TABLE default'
  write(unitout) trim(cbuffer),lf
  do nb=1,nbleafs
     do i=1,nx
        do j=1,ny
           write(unitout) primit(nbl(nb),5,i,j)
        end do
     end do
  end do
  write(unitout) lf
  !
  !   Temperature
  write(cbuffer,'(a)')  'SCALARS T float 1'
  write(unitout) trim(cbuffer),lf
  write(cbuffer,'(a)')  'LOOKUP_TABLE default'
  write(unitout) trim(cbuffer),lf
  do nb=1,nbleafs
     do i=1,nx
        do j=1,ny
           uu(:)=u(nbl(nb),:,i,j)
           call uprim (prim,uu,T)
           write(unitout) T
        end do
     end do
  end do
  write(unitout) lf
  !
  !   Block ID
  write(cbuffer,'(a)')  'SCALARS blockID float 1'
  write(unitout) trim(cbuffer),lf
  write(cbuffer,'(a)')  'LOOKUP_TABLE default'
  write(unitout) trim(cbuffer),lf
  do nb=1,nbleafs
     do i=1,nx
        do j=1,ny
           write(unitout) nbl(nb)
        end do
     end do
  end do
  write(unitout) lf
  !  
  !   Passive scalar
  write(cbuffer,'(a)')  'SCALARS Passive float 1'
  write(unitout) trim(cbuffer),lf
  write(cbuffer,'(a)')  'LOOKUP_TABLE default'
  write(unitout) trim(cbuffer),lf
  do nb=1,nbleafs
     do i=1,nx
        do j=1,ny
           write(unitout) primit(nbl(nb),6,i,j)
        end do
     end do
  end do
  write(unitout) lf
  !
  !  Velocity
  write(cbuffer,'(a)')  'VECTORS vel float'
  write(unitout) trim(cbuffer),lf
  do nb=1,nbleafs
     do i=1,nx
        do j=1,ny
           write(unitout)  primit(nbl(nb),2,i,j)*sqrt(vsc2)  &
                            ,primit(nbl(nb),3,i,j)*sqrt(vsc2)  &
                            ,0.
        end do
     end do
  end do
  write(unitout) lf
  !
  close(unitout)
  !
  deallocate(nbl)
#endif
!
#ifdef MPIP
  call mpi_barrier(mpi_comm_world,err)
#endif
  !
end subroutine output
!========================================================================
