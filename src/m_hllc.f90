module m_hllc
  use m_state
  use m_eos
  private
  public :: wave_speeds, hllc_flux
contains
  ! Only the fastest and slowest speeds are given (and needed)
  subroutine wave_speeds(speeds,L,R,primL,primR,dirn)
    real, intent(out) :: speeds(2)
    real, intent(in), dimension(nq) :: L, R, primL, primR
    integer, intent(in) :: dirn
    real, parameter :: min_density = 1e-8, min_pressure = 1e-8
    logical vacuum_left, vacuum_right

    real :: rhoL, pL, vL, rhoEL
    real :: rhoR, pR, vR, rhoER

    ! acoustic wave speeds
    real aL, aR
    ! enthalpies
    real HL, HR
    ! Roe-averaged states
    real v_Roe, H_Roe, a_Roe

    rhoL = L(cons_rho); rhoR = R(cons_rho)
    pL = primL(prim_p); pR = primR(prim_p)
    vL = primL(prim_v(dirn)); vR = primR(prim_v(dirn))
    rhoEL = L(cons_rhoE); rhoER = R(cons_rhoE)

    vacuum_left = rhoL < min_density .or. pL < min_pressure
    vacuum_right = rhoR < min_density .or. pR < min_pressure

    if (vacuum_left.and.vacuum_right) then
       speeds = 0
       return
    end if

    if (vacuum_left) then
       aR = sqrt(gamma*pR/rhoR)
       speeds(1) = vR - 2*aR/(gamma-1)
       speeds(2) = vR + aR
    else if (vacuum_right) then
       aL = sqrt(gamma*pL/rhoL)
       speeds(1) = vL - aL
       speeds(2) = vL + 2*aL/(gamma-1)
    else
       aL = sqrt(gamma*pL/rhoL)
       HL = (rhoEL + pL)/rhoL
       aR = sqrt(gamma*pR/rhoR)
       HR = (rhoER + pR)/rhoR

       v_Roe = (sqrt(rhoL)*vL + sqrt(rhoR)*vR)/(sqrt(rhoL) + sqrt(rhoR))
       H_Roe = (sqrt(rhoL)*HL + sqrt(rhoR)*HR)/(sqrt(rhoL) + sqrt(rhoR))
       a_Roe = sqrt((gamma-1)*(H_Roe - 0.5*v_Roe*v_Roe))

       speeds(1) = min(v_Roe-a_Roe, vL-aL)
       speeds(2) = max(v_Roe+a_Roe, vR+aR)
!       speeds = [v_Roe-a_Roe, v_Roe+a_Roe]

    end if
  end subroutine wave_speeds

  subroutine hllc_flux(F,S,L,R,primL,primR,dirn)
    use m_flux, only: flux
    real, intent(out) :: F(nq), S(2) ! flux and wave speed estimates
    real, dimension(nq), intent(in) :: L,R,primL,primR
    integer, intent(in) :: dirn
    real Sstar ! wave speed of star region
    real, dimension(nq) :: FL,FR, FLstar, FRstar
    real, dimension(nq) :: starL, starR

    real  :: rhoL, rhoR
    real  :: pL, pR
    real, dimension(3) :: vL, vR
    real, dimension(3) :: momL, momR
    real  :: rhoEL, rhoER
    real  :: rho_lambdaL, rho_lambdaR

    real Lfactor, Rfactor

    rhoL = L(cons_rho); rhoR = R(cons_rho)
    pL = primL(prim_p); pR = primR(prim_p)
    vL = primL(prim_v); vR = primR(prim_v)
    momL = L(cons_mom); momR = R(cons_mom)
    rhoEL = L(cons_rhoE); rhoER = R(cons_rhoE)
    rho_lambdaL = L(cons_rho_lambda); rho_lambdaR = R(cons_rho_lambda)

    call wave_speeds(S,L,R,primL,primR,dirn)

    if (S(1)>=0) then
       F = flux(L,primL,dirn)
    else if (0>=S(2)) then
       F = flux(R,primR,dirn)
    else
       Sstar = pR - pL + rhoL*vL(dirn)*(S(1) - vL(dirn)) - rhoR*vR(dirn)*(S(2) - vR(dirn))
       Sstar = Sstar / (rhoL*(S(1) - vL(dirn)) - rhoR*(S(2) - vR(dirn)))
       F = HUGE(F)
       if (S(1)<0.and.0<=Sstar) then
          Lfactor = ((S(1) - vL(dirn))/(S(1) - Sstar))
          starL(cons_rho)        = rhoL
          starL(cons_mom)        = momL
          starL(cons_mom(dirn))  = rhoL * Sstar
          starL(cons_rhoE)       = rhoEL + (rhoL*Sstar - rhoL*vL(dirn)) * (Sstar + pL/(S(1)-vL(dirn)))
          starL(cons_rho_lambda) = rho_lambdaL

          starL = starL * Lfactor
          FL = flux(L,primL,dirn)
          FLstar = FL + S(1) * (starL - L)
          F = FLstar
       else if (Sstar<=0.and.0<S(2)) then
          Rfactor = ((S(2) - vR(dirn))/(S(2) - Sstar))
          starR(cons_rho)       = rhoR
          starR(cons_mom)       = momR
          starR(cons_mom(dirn)) = rhoR * Sstar
          starR(cons_rhoE)      = rhoER + (rhoR*Sstar - rhoR*vR(dirn)) * (Sstar + pR/(S(2)-vR(dirn)))
          starR(cons_rho_lambda) = rho_lambdaR
          starR = starR * Rfactor
          FR = flux(R,primR,dirn)
          FRstar = FR + S(2) * (starR - R)
          F = FRstar
       end if
    end if
  end subroutine hllc_flux

end module m_hllc
