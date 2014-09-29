!=======================================================================
!   Integrates- Euler equations in two dimensions 
!   with adaptive mesh refinement
!   the integration method is set in the makefile file
!
!   the flow variables:
!   ieq=1 : rho  (total)
!       2 : rho u
!       3 : rho v
!       4 : Internal energy (thermal+kinetic)
!
!       additional variables advected into the flow:
!   with parametrized (H) cooling
!       5 : rho (neutrals)
!       6:  passive scalar * rho
!   with BBC cooling (additional scalars should be added at the end)
!       5:  n_HI
!       6:  n_HII
!       7:  n_HeI
!       8:  n_HeII
!       9:  n_HeIII
!      10:  rho*zbar
!      11:  ne   
!=======================================================================
program walicxe
  use parameters
  use globals
  implicit none
  integer :: err
  integer :: itprint,iload
  real    :: time, dt, tprint
  !--------------------------------------------------------------------
  !
  !   initializes mpi, grid spacing and other variables
  call initmain(time, tprint, itprint)
  !
  !--------------------------------------------------------------------
  !   set initial conditions
  if (iwarm.eq.0) then
     !   if iwarm eq 0 only master processor, load balance
     !   will distrubute the work
     if(rank.eq.0) then 
        call basegrid(itprint)
     else
        nbmax=nbmin-1
        nbleafs=0
     end if
  else
     call basegrid(itprint)
  end if
  !
  allocate(innerbounds(nbleafs*4,3))
#ifdef MPIP
  allocate(extbounds(nbleafs*4,5))
#endif
  !   balance the load if mpi is enabled
#ifdef MPIP
  call loadbalance
#endif
  call admesh
  call admesh
  call admesh
  call admesh
  if (rank.eq.0) print*,'---'
  !
  !
  !******************************************************************
  !   time integration (main loop)
  do while (time.le.tmax)
     !
     do iload=1,nload
        !------------------------------
        !   output at intervals tprint
        if(time.ge.tprint) then
           call output(itprint)
           if (rank.eq.0) then 
              print'(a,e10.2,a,e10.2,a,f5.2,a)',                        &
                   't:',time,'  dt:',dt,'  ',(time/tmax)*100,'%'
              print'(a,i4)','------------- wrote output ------------- :'&
                   , itprint
           end if
           tprint=tprint+dtprint
           itprint=itprint+1
        end if
        !------------------------------
        !   computes the timestep  
        call timestep(dt)
        time=time + dt
        if (rank.eq.1) print'(i5,a,es12.3,a,es12.3,a,es12.3,a)'             &
             ,iload,' time=',time*tsc/3.156e7,' yr  dt=', dt*tsc/3.156e7,   &
                    ' yr tprint=',tprint*tsc/3.156e7,' yr'
        !
        !
        !   advances the solution
        call tstep(time,dt)
        !
        !  update mesh
        call admesh
        !
      end do
     !
#ifdef MPIP
     call loadbalance
#endif
     if (rank.eq.0) print*,"---"
     !
  end do
  !******************************************************************
  !   finishes
  deallocate(u,up,primit,lp,leaf,icoord,dx, dy)
  if (rank.eq.0) print'(a)',"--- My work here is done, have a nice day ---"
#ifdef MPIP
  call mpi_finalize(err)
#endif
  stop
end program  walicxe
!====================================================================
