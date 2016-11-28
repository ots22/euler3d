module m_state
  integer, parameter :: nq = 6
  integer, parameter :: prim_rho=1, prim_v(3)=[2,3,4], prim_p=5, prim_lambda=6
  integer, parameter :: cons_rho=1, cons_mom(3)=[2,3,4], cons_rhoE=5, cons_rho_lambda=6
contains
  pure function prim_get_rho(u) result(rho)
    real, intent(in) :: u(nq)
    real rho
    rho = u(prim_rho)
  end function prim_get_rho

  pure function prim_get_v(u) result (v)
    real, intent(in) :: u(nq)
    real v(3)
    v = u(prim_v)
  end function prim_get_v

  pure function prim_get_p(u) result (p)
    real, intent(in) :: u(nq)
    real p
    p = u(prim_p)
  end function prim_get_p

  pure function prim_get_lambda(u) result (lambda)
    real, intent(in) :: u(nq)
    real lambda
    lambda = u(prim_lambda)
  end function prim_get_lambda

  pure function cons_get_rho(u) result(rho)
    real, intent(in) :: u(nq)
    real rho
    rho = u(cons_rho)
  end function cons_get_rho

  pure function cons_get_mom(u) result (mom)
    real, intent(in) :: u(nq)
    real mom(3)
    mom = u(cons_mom)    
  end function cons_get_mom
  
  pure function cons_get_rhoE(u) result(rhoE)
    real, intent(in) :: u(nq)
    real rhoE
    rhoE = u(cons_rhoE)
  end function cons_get_rhoE

  pure function cons_get_rho_lambda(u) result(rho_lambda)
    real, intent(in) :: u(nq)
    real rho_lambda
    rho_lambda = u(cons_rho_lambda)
  end function cons_get_rho_lambda
end module m_state
