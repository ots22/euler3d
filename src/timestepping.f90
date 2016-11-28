module m_timestepping
contains
  subroutine advance_solution_1d(U, max_wave_speed, dx_dt, dirn)
    use m_reconstruct
    use m_domain, only: nx, nb
    use m_state, only: nq
    use m_eos
    use m_flux
    use m_hllc
    real, intent(inout) :: U(:,-nb+1:)
    real, intent(inout) :: max_wave_speed
    real, intent(in) :: dx_dt
    integer, intent(in) :: dirn

    ! the numerical fluxes, same shape as U, (but not all used)
    real, dimension(nq, lbound(U,2):ubound(U,2)) :: F, UL, UR
    real, dimension(nq) :: Li, Ri, primLi, primRi, fl_Li, fl_Ri
    real :: S(2) ! fastest left and right wave speed estimates

    integer ix

    ! plain FORCE uses the trivial reconstruction:
!    UL = U(:,:,1)
!    UR = U(:,:,1)

    call linear_reconstruct(U(:,:),UL,UR)

    do ix=0,nx(dirn)
!       Li = UR(:,ix)   ! right face reconstruction of the left cell
!       Ri = UL(:,ix+1) ! left face reconstruction of the right cell

! reuse these variables - change this
       Li = UL(:,ix)
       Ri = UR(:,ix)

       primLi = cons_to_prim(Li)
       primRi = cons_to_prim(Ri)

       fl_Li = flux(Li, primLi, dirn)
       fl_Ri = flux(Ri, primRi, dirn)

!       F(:,ix) = force_flux(Li, Ri, fl_Li, fl_Ri, dx_dt, dirn)

       ! half-dt step, with centred fluxes
       UL(:,ix) = UL(:,ix) + 0.5 * (fl_Li - fl_Ri) / dx_dt
       UR(:,ix) = UR(:,ix) + 0.5 * (fl_Li - fl_Ri) / dx_dt
    end do

!   Repeat with the half-dt reconstructed states
    do ix=0,nx(dirn)
       Li = UR(:,ix)   ! right face reconstruction of the left cell
       Ri = UL(:,ix+1) ! left face reconstruction of the right cell
       primLi = cons_to_prim(Li)
       primRi = cons_to_prim(Ri)

       fl_Li = flux(Li, primLi, dirn)
       fl_Ri = flux(Ri, primRi, dirn)

       call hllc_flux(F(:,ix),S,Li,Ri,primLi,primRi,dirn)

       max_wave_speed = maxval([max_wave_speed, abs(S)])
    end do

!   apply the computed flux
    do ix=1,nx(dirn)
       U(:,ix) = U(:,ix) + (F(:,ix-1) - F(:,ix)) / dx_dt
    end do
  end subroutine advance_solution_1d
end module m_timestepping
