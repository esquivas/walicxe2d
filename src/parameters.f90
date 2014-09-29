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
       outputpath='/datos/esquivel/W2D-test/'
  integer, parameter :: neq=6,ndim=2   ! number of eqs (+scal) and dimensions
  integer, parameter :: npas=2         ! number of passive scalars
  integer, parameter :: nghost=2       ! number of ghost cells
  integer, parameter :: nblocks=10000  ! max number of blocks allowed p/proc
  integer, parameter :: nx=16  ,ny=16  ! block size
  integer, parameter :: nlevs=7        ! number of levels
  integer, parameter :: nbroot=4       ! number of root blocks
  !   this way maximum resolution is 2**(nlevs-1)*nx*ny
  !-------------------------------------------------------------------
  !   mpi parameters
#ifdef MPIP
  ! number of time iterations before load balancing
  integer, parameter :: nload=50
#endif
  !
  !  some constants
  real, parameter :: pc=3.0857E18, AU=1.496e13
  real, parameter :: clight=2.99E10
  real, parameter :: Msun=1.99E33, Rsun=6.955e10
  real, parameter :: amh=1.66e-24,Kb=1.38e-16,Rg=8.3145e7
  real, parameter :: yr=3.156E7
  real, parameter :: pi=acos(-1.)
  real, parameter :: amol= 0.6

  !--------------------------------------------------------------------
  !   the following options could go on an runtime input file
  !--------------------------------------------------------------------
  real, parameter :: cv=1.5, gamma=(cv+1.)/cv 
  !
  !  box size in normalized units
  real, parameter :: xmax=1.0, ymax=0.25
  !
  !  box size in physical (cgs) units
  real, parameter :: xphys=3.E17 !cm 
  !--------------------------------------------------------------------
  !
  real, parameter :: T0= 1.e4  ! a reference Temperature for scalings
  !   scaling factors to physical (cgs) units
  real, parameter :: rsc=xphys/xmax, rhosc=amh*amol
  real, parameter :: Tempsc=T0*gamma
  real, parameter :: Psc = gamma*amh*Rg*T0
  real, parameter :: vsc2 =gamma*T0*Rg/amol
  real, parameter :: tsc = rsc/sqrt(vsc2)
  !--------------------------------------------------------------------
  !
  !   maximum time of integration and output interval
  real, parameter :: tmax   = 355.*yr/tsc
  real, parameter :: dtprint=  10.*yr/tsc
  !   Courant number
  real, parameter :: cfl=0.7
  !   Viscosity
  real, parameter :: eta=0.001
  ! 
  !   for an iwarm (0->from t=0,1->starts from t=itprint*dtprint
  integer, parameter :: iwarm=0
  integer, parameter :: itprint0=29
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
