!=======================================================================
!  User imput module
!  This is an attempt to have all input needed from the user in a single 
!  file
!  This module should load additional modules (i.e. star, jet, sn), to
!  impose initial and boundary conditions (such as sources)
!=======================================================================
module user_mod
  !  here we should load additional modules
  use Jet
  !use SN
  implicit none
  !  Variables global to the module can be here
  ! real :: variablename

contains

!--------------------------------------------------------------------
!   Here the parameters of the User Defined module are initialized, 
!   and scaled to code units
!   This routine is alway called initmain.f90
!   It has to be present, even if empty
!--------------------------------------------------------------------
subroutine init_user_mod()
  implicit none
    
  !  initialize modules loaded by user if needed
  call init_jet()

end subroutine init_user_mod
!--------------------------------------------------------------------


!--------------------------------------------------------------------
!   Here the domain (in one block) is initialized at t=0 
!   In general takes U, nb and time as argument, if time is not
!   required just ignore it...
!--------------------------------------------------------------------
subroutine initial_conditions(u,nb,time)
  use parameters, only : neq, nxmin, nxmax, nymin, nymax   &
                       , Tempsc, pc, cv
  implicit none
  real,    intent(out) :: u(neq,nxmin:nxmax,nymin:nymax)
  integer, intent(in)  :: nb
  real,    intent(in)  :: time
  !   Surrounding ISM
  real, parameter :: nenv=15.0, Tenv=1.e3
  integer ::  i, j

  !  fill ISM
  do i=nxmin,nxmax
     do j=nymin,nymax
      u(1,i,j)= nenv
      u(2,i,j)= 0.
      u(3,i,j)= 0.
      u(4,i,j)= cv*1.0001*nenv*Tenv/Tempsc
      u(5,i,j)= nenv*0.9999  !  neutral fraction
      u(6,i,j)=-nenv
   end do
  end do

  !  Impose Jet, overwties ISM
  call impose_jet(u, nb, time)
  !call impose_SN(u,nb, 0.,0.)
    
end subroutine initial_conditions
!--------------------------------------------------------------------


!--------------------------------------------------------------------
!   User Defined Boundary conditions (block with ID nb)
!--------------------------------------------------------------------
#ifdef OTHERB
subroutine impose_user_bc(u,nb,time)
  use parameters, only : neq, nxmin, nxmax, nymin, nymax 
  implicit none
  real,    intent(out) :: u(neq,nxmin:nxmax,nymin:nymax)
  integer, intent(in)  :: nb
  real,    intent(in)  :: time
  
  call impose_jet(u, nb, time)
  
  end subroutine impose_user_bc
#endif
!--------------------------------------------------------------------
end module user_mod
