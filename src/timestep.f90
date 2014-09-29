!========================================================================
!   calculates the timestep as permited by the CFL criterion
!   (courant number=cfl defined in parameters.f90)
!   For coronal cooling it restricts the timestep with the cooling time
!   as well (this should be updated to incude all cooling
!========================================================================
subroutine timestep(dt)
  use parameters
  use globals
  use dmc_module
  implicit none
  real,    intent(out) :: dt
  real                 :: t, cs, dtp, dtcool
  integer              :: i, j, nb, err,lev
#ifdef COOLINGDMC
  real (kind=8) :: Aloss
#endif
  !
  dtp=1.e30
  dtcool=1e30
  do nb=nbmin,nbmax
     !
     if ( leaf(nb) ) then
        !
        lev=lp(nb,1)
        !
        do i=1,nx
           do j=1,ny
              call sound(primit(nb,4,i,j),primit(nb,1,i,j),cs)
              dtp=min(dtp,dx(lev)/(abs(primit(nb,2,i,j))+cs))  
              dtp=min(dtp,dy(lev)/(abs(primit(nb,3,i,j))+cs))
              !
              !   restrict the timestep w/t_cool
#ifdef COOLINGDMC
              T=(primit(nb,4,i,j)/primit(nb,1,i,j))*Tempsc
              if(T.gt.1e4) then 
                 Aloss=cooldmc(T)*(primit(nb,1,i,j)**2)
                 !
                 dtcool=min(dtcool, cv*primit(nb,4,i,j)/aloss )
              endif
#endif
              !
           end do
        end do
     end if
     !
  end do
  !
  !
  dtp=min(cfl*dtp,0.4*dtcool*Psc/tsc)
#ifdef MPIP
  call mpi_allreduce(dtp, dt, 1, mpi_real_kind, mpi_min, mpi_comm_world,err)
#else
  dt=dtp
#endif
  !
  !if(rank.eq.1) print*,dt,dtcool,dtp
  !call mpi_finalize(err)
  !stop
  !
end subroutine timestep
!========================================================================

