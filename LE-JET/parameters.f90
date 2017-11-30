!=======================================================================
!   This module contains all variables and constants that are global
!=======================================================================
module parameters
  implicit none 
#ifdef MPIP
  include "mpif.h"
#endif
  !   output path
  character (len=128),parameter ::                                     &
       outputpath='./Cool'
  integer, parameter :: neq=5,ndim=2   ! number of eqs (+scal) and dimensions
  integer, parameter :: npas=1         ! number of passive scalars
  integer, parameter :: nghost=2       ! number of ghost cells
  integer, parameter :: nblocks=10000  ! max number of blocks allowed p/proc
  integer, parameter :: nx=32  ,ny=32  ! block size
  integer, parameter :: nlevs=7        ! number of levels
  integer, parameter :: nbroot=4       ! number of root blocks
  !   this way maximum resolution is 2**(nlevs-1)*(nx,ny)
  !-------------------------------------------------------------------
  !   mpi parameters
#ifdef MPIP
  ! number of time iterations before load balancing
  integer, parameter :: nload=10
#endif
  !
  !  some constants
  real, parameter :: pc=3.0857E18, AU=1.496e13
  real, parameter :: clight=2.99E10
  real, parameter :: Msun=1.99E33, Rsun=6.955e10
  real, parameter :: amh=1.66e-24,Kb=1.38e-16,Rg=8.3145e7
  real, parameter :: yr=3.156E7
  real, parameter :: pi=acos(-1.)
  real, parameter :: amol= 1.0
  !
  !--------------------------------------------------------------------
  !   the following options could go on an runtime input file
  !--------------------------------------------------------------------
  real, parameter :: cv=1.5, gamma=(cv+1.)/cv 
  !
  !  box size in normalized units
  real, parameter :: xmax=1.0, ymax=0.25
  !
  !  box size in physical (cgs) units
  real, parameter :: xphys=6.E17 !cm 
  !--------------------------------------------------------------------
  !
  !   scaling factors to physical (cgs) units
  real, parameter :: T0=1.e4                ! reference temperature (for cs)
  real, parameter :: rsc=xphys/xmax         ! distance scaling
  real, parameter :: rhosc= amh*amol          ! mass density scaling
  real, parameter :: Tempsc=T0*gamma        ! Temperature scaling
  real, parameter :: vsc2 = gamma*Rg*T0/amol  ! Velocity scaling squared
  real, parameter :: vsc = sqrt(vsc2)       ! Velocity scaling
  real, parameter :: Psc = rhosc*vsc2       ! Pressure scaling
  real, parameter :: tsc =rsc/vsc    ! time scaling
  real, parameter :: bsc = sqrt(4.0*pi*Psc) !magnetic field scaling
  !
  !--------------------------------------------------------------------
  !
  !   maximum time of integration and output interval
  real, parameter :: tmax   = 1000.*yr/tsc
  real, parameter :: dtprint=  5.*yr/tsc
  !   Courant number
  real, parameter :: cfl=0.5
  !   Viscosity
  real, parameter :: eta=0.001
  ! 
  !   for an iwarm (0->from t=0,1->starts from t=itprint*dtprint
  integer, parameter :: iwarm=0
  integer, parameter :: itprint0=100
  !--------------------------------------------------------------------
  !   some derived parameters (no need of user's input below this line)
  !--------------------------------------------------------------------
  integer, parameter :: neqpas=4+npas
  integer, parameter :: master=0
  !
  !   array bounds
  integer, parameter :: nxmin=1-nghost,    nxmax=nx+nghost,            &
                        nymin=1-nghost,    nymax=ny+nghost
#ifndef MPIP
  ! number of time iterations before load balancing
  integer, parameter :: nload=1
#endif
  !
  !   set floating point precision (kind) for MPI messages
#ifdef MPIP
#ifdef DOUBLEP
  integer, parameter :: mpi_real_kind=mpi_real8
#else
  integer, parameter :: mpi_real_kind=mpi_real4
#endif
#endif
  !--------------------------------------------------------------------
  !
end module parameters
!=======================================================================
