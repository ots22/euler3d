module m_flux
  use m_eos, only: cons_to_prim
  use m_state
  
contains

  function flux(cs,ps,dirn)
    real cs(nq), ps(nq), flux(nq)
    integer dirn
    real v(3), p
    real rho, mom(3), rhoE, rho_lambda
    real rho_fl, mom_fl(3), rhoE_fl, rho_lambda_fl

    v = prim_get_v(ps)
    p = prim_get_p(ps)
    
    rho = cons_get_rho(cs)
    mom = cons_get_mom(cs)
    rhoE = cons_get_rhoE(cs)
    rho_lambda = cons_get_rho_lambda(cs)

    rho_fl = mom(dirn)
    mom_fl = v(dirn) * mom
    mom_fl(dirn) = mom_fl(dirn) + p
    rhoE_fl = v(dirn) * (rhoE + p)
    rho_lambda_fl = v(dirn) * rho_lambda

    flux = [rho_fl, mom_fl, rhoE_fl, rho_lambda_fl]
  end function flux

  function lax_friedrichs_flux(L, R, fl_L, fl_R, dx_dt) result(fl)
    real, dimension(nq) :: L, R, fl_L, fl_R, fl
    real dx_dt
    fl = 0.5 * (fl_L + fl_R + dx_dt * (L - R))
  end function lax_friedrichs_flux

  function richtmeyer_flux(L, R, fl_L, fl_R, dx_dt, dirn) result(fl)
    real, dimension(nq) :: L, R, fl_L, fl_R, fl, ri_state, ri_pstate
    real dx_dt
    integer dirn
    ri_state = 0.5 * (L + R + (1.0/dx_dt) * (fl_L - fl_R))
    ri_pstate = cons_to_prim(ri_state)
    fl = flux(ri_state, ri_pstate, dirn)
  end function richtmeyer_flux

  function force_flux(L, R, fl_L, fl_R, dx_dt, dirn) result(fl)
    real, dimension(nq) :: L, R, fl_L, fl_R, fl
    real dx_dt
    integer dirn
    fl = 0.5 * (lax_friedrichs_flux(L,R,fl_L,fl_R,dx_dt) &
         &    + richtmeyer_flux(L,R,fl_L,fl_R,dx_dt,dirn))
  end function force_flux
end module m_flux
