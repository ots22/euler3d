module m_reconstruct
  private
  public linear_reconstruct
contains
  elemental function phi_vanleer(d1, d2) result(phi)
    real phi, r
    real, intent(in) :: d1, d2
    r = d1/merge(d2, 1.0, abs(d2).gt.1E-9)
    !phi = (r + abs(r)) / (1 + abs(r))
    if (r.le.0) then
       phi = 0
    else
       phi = min(2*r/(1+r), 2/(1+r))
    endif
  end function phi_vanleer

  subroutine linear_reconstruct(U, UL, UR)
    use m_state, only: cons_rho
    real, intent(in) :: U(:,0:)
    real, intent(out) :: UL(:,0:), UR(:,0:)
    integer nq, nx
    integer iq, ix
    real slopes(size(U,1),size(U,2)-2)
    nq=size(U,1); nx=size(U,2)-2
    slopes(:,1:nx) = 0.5 * (U(:,2:nx+1) - U(:,0:nx-1))
    
    ! limit each component on itself
    do ix=1,nx; do iq=1,nq       
       slopes(iq,ix) = slopes(iq,ix) * phi_vanleer(U(iq,ix)-U(iq,ix-1), U(iq,ix+1)-U(iq,ix))
    end do; end do

    UL = U
    UR = U
    UL(:,1:nx) = U(:,1:nx) - 0.5 * slopes
    UR(:,1:nx) = U(:,1:nx) + 0.5 * slopes
  end subroutine linear_reconstruct

end module m_reconstruct
