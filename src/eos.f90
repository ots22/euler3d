module m_eos
  real, parameter :: gamma=1.4
contains
  function cons_to_prim(cs) result(ps)
    use m_state
    real rho, v(3), p, lambda
    real, intent(in) :: cs(nq)
    real ps(nq)
    rho = cons_get_rho(cs)
    v = cons_get_mom(cs)/rho
    p = (gamma-1) * (cons_get_rhoE(cs) - 0.5*rho*dot_product(v,v))
    lambda = cons_get_rho_lambda(cs)/rho
    ps = [rho, v, p, lambda]
  end function cons_to_prim

  function prim_to_cons(ps) result(cs)
    use m_state
    real rho, v(3), mom(3), rho_lambda, rhoE
    real, intent(in) :: ps(nq)
    real cs(nq)
    rho = prim_get_rho(ps)
    v = prim_get_v(ps)
    mom = rho * v
    rhoE = 0.5 * rho * dot_product(v,v) + prim_get_p(ps)/(gamma-1)
    rho_lambda = rho * prim_get_lambda(ps)
    cs = [rho, mom, rhoE, rho_lambda]
  end function prim_to_cons
end module m_eos
