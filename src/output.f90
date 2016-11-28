module m_output
  private
  public visit_output

contains
  
  subroutine visit_output(unit,it,t,U)
    use m_domain
    use m_state
    use m_eos

    integer, intent(in) :: unit, it
    real, intent(in) :: t, U(:,-nb+1:,-nb+1:,-nb+1:)

    real rho, p
    integer ix,iy,iz

    write (unit,'(A)') '# vtk DataFile Version 3.0'
    write (unit,'(A)') 'vtk output'
    write (unit,'(A)') 'ASCII'
    write (unit,'(A)') 'DATASET RECTILINEAR_GRID'
    write (unit,'(A)') 'FIELD FieldData 2'
    write (unit,'(A)') 'TIME 1 1 double'
    write (unit,'(E16.7E3)') t
    write (unit,'(A)') 'CYCLE 1 1 int'
    write (unit,'(I3)') it
    write (unit,'(A,I7,I7,I7)') 'DIMENSIONS ', nx(1)+1, nx(2)+1, nx(3)+1

    write (unit,'(A,I7,1X,A)') 'X_COORDINATES', nx(1)+1, 'FLOAT'
    write (unit,'(E16.7E3)') REAL(0)
    do ix=1,nx(1)
       write (unit,'(E16.7E3)') ix*dx
    end do

    write (unit,'(A,I7,1X,A)') 'Y_COORDINATES', nx(2)+1, 'FLOAT'
    write (unit,'(E16.7E3)') REAL(0)
    do ix=1,nx(2)
       write (unit,'(E16.7E3)') ix*dx
    end do

    write (unit,'(A,I7,1X,A)') 'Z_COORDINATES', nx(3)+1, 'FLOAT'
    write (unit,'(E16.7E3)') REAL(0)
    do ix=1,nx(3)
       write (unit,'(E16.7E3)') ix*dx
    end do

    write (unit,'(A,I15)') 'CELL_DATA', product(nx)

    call write_scalars_header('rho')
    do iz=1,nx(3); do iy=1,nx(2); do ix=1,nx(1)
       rho = cons_get_rho(U(:,ix,iy,iz))
       write (unit,'(E16.7E3)') rho
    end do; end do; end do;

    call write_scalars_header('pressure')
    do iz=1,nx(3); do iy=1,nx(2); do ix=1,nx(1)
       p = prim_get_p(cons_to_prim(U(:,ix,iy,iz)))
       write (unit,'(E16.7E3)') p
    end do; end do; end do;   

    write (unit,'(A)') 'VECTORS momentum FLOAT'
    write (unit,'(3E16.7E3)') U(cons_mom,1:nx(1),1:nx(2),1:nx(3))

    call write_scalars_header('rhoE')
    write (unit,'(E16.7E3)') U(cons_rhoE,1:nx(1),1:nx(2),1:nx(3))

    call write_scalars_header('rho_lambda')
    write (unit,'(E16.7E3)') U(cons_rho_lambda,1:nx(1),1:nx(2),1:nx(3))
       
  contains
    subroutine write_scalars_header(name)
      character(*), intent(in) :: name
      write (unit,'(A,I7)') 'SCALARS ' // trim(name) // ' FLOAT', 1
      write (unit,'(A)') 'LOOKUP_TABLE default'
    end subroutine write_scalars_header
  end subroutine visit_output

end module m_output
