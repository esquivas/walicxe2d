!=======================================================================
!  initializes mpi, grid spacing and time integration
!=======================================================================
subroutine initmain(time, tprint, itprint)
  use parameters
  use globals
  use user_mod
#ifdef COOLINGDMC
  use dmc_module
#endif
  implicit none
  real,    intent(out) :: time, tprint
  integer, intent(out) :: itprint
  integer :: ilevs, err
#ifdef COOLINGDMC
  integer :: i
  real (kind=8) :: a, b
#endif
  !
  !--------------------------------------------------------------------
  !initializes mpi
#ifdef MPIP
  call mpi_init (err)
  call mpi_comm_rank (mpi_comm_world,rank,err)
  call mpi_comm_size (mpi_comm_world,np ,err)
#else
  rank=0
  np=1
#endif
  if(rank.eq.master) then
     print'(a)',"***********************************************"
     print'(a)',"                 __________                   "
     print'(a)',"___      _______ ___  /__(_)_________  ______ "
     print'(a)',"__ | /| / /  __ `/_  /__  /_  ___/_  |/_/  _ \"
     print'(a)',"__ |/ |/ // /_/ /_  / _  / / /__ __>  < /  __/"
     print'(a)',"____/|__/ \__,_/ /_/  /_/  \___/ /_/|_| \___/ "
     print'(a)',"                                              "
#ifdef MPIP
     print'(a,i3,a)',"*      running with mpi on", np , " processors      *"
     print'(a)',"***********************************************"
#else
     print'(a)',"***********************************************"
     print'(a)',"*       running on a single processor         *"
     print'(a)',"***********************************************"
#endif
  end if
#ifdef MPIP
  call mpi_barrier(mpi_comm_world, err)
  print '(a,i3,a)', 'processor ', rank,' ready'
#endif
  !
  nbmin   =nblocks*rank+1
  nbmaxtot=nblocks*(rank+1)
  nbmax=nbmin
  nbleafs=0
  !--------------------------------------------------------------------
  !   grid spacing at each level
  allocate(dx(nlevs),dy(nlevs))
  do ilevs=1,nlevs
     dx(ilevs) =xmax/nx/( 2.**(ilevs-1) )/float(nbroot)
     dy(ilevs) =ymax/ny/( 2.**(ilevs-1) )
  end do
  !
  !   initialize time integration
  if (iwarm.eq.0) then
     itprint=0
     time=0.
     tprint=0.
  else
     itprint=itprint0
     time=itprint*dtprint
     tprint=time+dtprint
  end if
  !
  !   allocate memory for big arrays
  allocate (       u(nbmin:nbmaxtot,neq,nxmin:nxmax,nymin:nymax))
  allocate (      up(nbmin:nbmaxtot,neq,nxmin:nxmax,nymin:nymax))
  allocate (  primit(nbmin:nbmaxtot,neq,nxmin:nxmax,nymin:nymax))
  allocate (      lp(nblocks*np, 8 ) )
  allocate (    leaf(nblocks*np    ) )
  allocate (  icoord(nbmin:nbmaxtot, 2 ) )
  !
#ifdef COOLINGDMC
  if(rank.eq.master) then
     open(unit=10,file='../src/DMClib/coolingDMC.tab',status='old')
     print'(a)'
     print'(a)', 'reading DMC cooling table'
     do i=1,41
        read(10,*) a, b
        cooltab(1,i)=10.d0**(a)
        cooltab(2,i)=10.d0**(-b)
     end do
     close(unit=10)
  endif
#ifdef MPIP
  call mpi_bcast(cooltab,82,mpi_double_precision,0,mpi_comm_world,err)
#endif
#endif
  !
  !  User input initialization, it is called always,
  !  it has to be there, even empty
  call init_user_mod()
  !
  !--------------------------------------------------------------------
  !   write report of compilation parameters
#ifdef MPIP
  call mpi_barrier(mpi_comm_world, err)
  if(rank.eq.master) then
     print'(a)',''
#endif
#ifdef DOUBLEP
     print'(a)', 'Double precision used (reals are 8 bytes long)'
     print'(a)', ''
#else
     print'(a)', 'Single precision used (reals are 4 bytes long)'
     print'(a)', ''
#endif
#ifdef HLLC
     print'(a)', 'The Riemann solver is HLLC'
     print'(a)', ''
#endif
#ifdef HLL
     print'(a)', 'The Riemann solver is HLL'
     print'(a)', ''
#endif
#ifdef HLL_HLLC
     print'(a)', 'The Riemann solver is HLL-HLLC (hybrid)'
     print'(a)', ''
#endif
#ifdef ADIABATIC
     print'(a)', 'The code is in ADIABATIC mode'
     print'(a)', ''
#endif
#ifdef COOLINGH
     print'(a)', 'Radiative cooling ON (w/parametrized cooling curve)'
     print'(a)', ''
#endif
#ifdef COOLINGBBG
     print'(a)', 'Radiative cooling ON (w/ Benjamin Benson & Cox 2003 prescription)'
     print'(a)', ''
#endif
#ifdef COOLINGDMC
     print'(a)', 'Radiative cooling ON (w/ Dalgarno & Mc Cray, coronal eq.)'
     print'(a)', ''
#endif
     print'(a)', '-----  OUTPUT -----------------------'
     print'(a)', 'path: '//trim(outputpath)
     print'(a)', 'in the following format(s):'
#ifdef OUTBIN
     print'(a)', '*.bin (raw unformatted)'
#endif
#ifdef OUTVTK
     print'(a)', '*.vtk (binary VTK)'
#endif
     print'(a)', ''
     print'(a)', '----- BOUNDARY CONDITIONS -----------'
#ifdef PERIODX
     print'(a)', 'LEFT & RIGHT: PERIODIC'
#endif
#ifdef PERIODY
     print'(a)', 'BOTTOM & TOP: PERIODIC'
#endif
#ifdef REFXL
     print'(a)', 'LEFT:   REFLECTIVE (CLOSED)'
#endif
#ifdef OUTFXL
     print'(a)', 'LEFT:   OUTFLOW    (OPEN)'
#endif
#ifdef INFXL
     print'(a)', 'LEFT:   INFLOW     (user defined)'
#endif
#ifdef REFXR
     print'(a)', 'RIGHT:  REFLECTIVE (CLOSED)'
#endif
#ifdef OUTFXR
     print'(a)', 'RIGHT:  OUTFLOW    (OPEN)'
#endif
#ifdef INFXR
     print'(a)', 'RIGHT:  INFLOW     (user defined)'
#endif
#ifdef REFYB
     print'(a)', 'BOTTOM: REFLECTIVE (CLOSED)'
#endif
#ifdef OUTFYB
     print'(a)', 'BOTTOM: OUTFLOW    (OPEN)'
#endif
#ifdef INFYB
     print'(a)', 'BOTTOM: INFLOW     (user defined)'
#endif
#ifdef REFYT
     print'(a)', 'TOP:    REFLECTIVE (CLOSED)'
#endif
#ifdef OUTFYT
     print'(a)', 'TOP:    OUTFLOW    (OPEN)'
#endif
#ifdef INFYT
     print'(a)', 'TOP:    INFLOW     (user defined)'
#endif
#ifdef OTHERB
     print'(a)', 'Other boundaries enabled (otherbounds.f90)'
#endif
     print'(a)', ''
#if LIMITER==-1
     print'(a)', 'No average in the limiter (reduces to 1st order)'
#endif
#if LIMITER==0
     print'(a)', 'No limiter'
#endif
#if LIMITER==1
     print'(a)', 'MINMOD limiter -most diffusive-'
#endif
#if LIMITER==2
     print'(a)', 'Falle Limiter (Van Leer)'
#endif
#if LIMITER==3
     print'(a)', 'Van Albada Limiter'
#endif
#if LIMITER==4
     print'(a)', 'UMIST limiter -least diffusive-'
#endif
#if LIMITER==5
     print'(a)', 'Woodward Limiter (MC-limiter; monotonized central difference)'
#endif
#if LIMITER==6
     print'(a)', 'SUPERBEE limiter (tends to flatten circular waves)'
#endif
     print'(a)', ''
     print'(a)','***********************************************'
#ifdef MPIP
  end if
  call mpi_barrier(mpi_comm_world, err)
#endif
  !--------------------------------------------------------------------
  !
end subroutine initmain
!=======================================================================
