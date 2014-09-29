!======================================================================
!   integration from t to t+dt with approximate Riemnann solver
!======================================================================
subroutine tstep(time,dt)
  use parameters
  use globals
  implicit none
  real, intent(in) :: time, dt
  real,    dimension(neq,nxmin:nxmax,nymin:nymax) :: f, g
  real :: dtm
  integer :: nb, lev
  !
  !   locate the boundaries
  call locatebounds(time)
  !
  !--------------------------------------------------------------------
  !  1st half timestep
  !--------------------------------------------------------------------
  !
  dtm=dt/2.
  !
  DO nb=nbmin,nbmax
     IF( leaf(nb) ) THEN
        !
        !   calculate the fluxes using the primitives
        !   (piecewise constant)
        !   ---> f_i^n=Riemann[p_i^n,p_{i+1}^n]
#ifdef HLL
        call hllfluxes(nb,1,f,g)
#endif
#ifdef HLLC
        call hllcfluxes(nb,1,f,g)
#endif
#ifdef HLL_HLLC
        call hll_hllcfluxes(nb,1,f,g)
#endif
        !
        !   upwind timestep (block nb)
        !   up=u^{n+1/2}=u_i^n-(dt/2dx)[f_i^n - f_{i-1}^n]
        call step(nb,dtm,f,g)
        !
        !   add viscosity
        !call viscosity(nb,nbmin,nbmaxtot,u,up)
        !
     END IF
  END DO
  !
  !--------------------------------------------------------------------
  !  2nd half timestep
  !--------------------------------------------------------------------
  !
  !   boundary conditions on up
    call boundaryII(time)!,ninternal,nexternal,innerbounds,extbounds)
  !
  DO nb=nbmin,nbmax
     IF( leaf(nb) ) THEN
        !
        !    update the primitives with up primit=primit^{n+1/2}
        call calcprim(up(nb,:,:,:), primit(nb,:,:,:),              &
             neq, nxmin, nxmax, nymin, nymax ) 
        !
        !
        !   calculate the fluxes using the primitives
        !   with linear reconstruction (piecewise linear)
        !   p_{R/L}=p_{i(+/-)1/2}^{n+1/2}
        !   ---> f_{i+1/2}^(n+1/2)=Riemann[p_{R,i}^{n+1/2},p_{L+1}^{n+1/2}]
#ifdef HLL
        call hllfluxes(nb,2,f,g)
#endif
#ifdef HLLC
        call hllcfluxes(nb,2,f,g)
#endif
#ifdef HLL_HLLC
        call hll_hllcfluxes(nb,2,f,g)
#endif
        !
        !   upwind timestep (block nb)
        !   up=u^{n+1}=u_i^n-(dt/dx)[f_{i+1/2}^{n+1/2}-f_{i-1/2}^{n+1/2}]
        call step(nb,dt,f,g)
        !
        !   add viscosity
        if (lp(nb,1).ge.(nlevs-1)) call viscosity(nb,nbmin,nbmaxtot,u,up)
        !
        !   Replace the u's with the ups
        u(nb,:,:,:)=up(nb,:,:,:)
        !
     END IF
  END DO
  !--------------------------------------------------------------------
  !
  !   boundary contiditions on u
  call boundaryI(time)
  !
  !   if cooling's off just update the primitives
#ifdef ADIABATIC
  DO nb=nbmin,nbmax
     IF( leaf(nb) ) THEN
        call calcprim(u(nb,:,:,:), primit(nb,:,:,:),              &
             neq, nxmin, nxmax, nymin, nymax ) 
     END IF
  END DO
#endif
  !
  !****************************************************
  !   apply cooling
  !****************************************************
#ifdef COOLINGH
  !   add cooling to the conserved variables
  call coolingh(dt*tsc)
  !  update the primitives with u
  DO nb=nbmin,nbmax
     IF( leaf(nb) ) THEN
        call calcprim(u(nb,:,:,:), primit(nb,:,:,:),              &
                      neq, nxmin, nxmax, nymin, nymax ) 
     END IF
  END DO
#endif
  !****************************************************
#ifdef COOLINGDMC
  !   the primitives are updated in the cooling routine
  call coolingdmc(dt*tsc)
#endif
  !****************************************************
  !   BBC cooling not implemented in the 2D version yet
  !****************************************************
end subroutine tstep
!======================================================================
