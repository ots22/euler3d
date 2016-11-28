module m_config
  real cfl
  integer outstep
  real tmax
  integer, parameter :: MAX_FILENAME_LEN=1000
  character(len=MAX_FILENAME_LEN) output_filename_stem
  namelist/CONFIG/cfl,outstep,tmax,output_filename_stem

contains
  subroutine init_config(unit)
    integer, intent(in) :: unit
    cfl = 0.8
    outstep = 1
    tmax = 0
    output_filename_stem = 'output_'

    read (unit, nml=CONFIG)
  end subroutine init_config
end module m_config

module m_simulation_data
  ! simulation time
  real t
  ! timestep number
  integer it
  ! size of current timestep
  real dt

  ! two arrays of conserved solution data
  ! dimensions are as follows:
  !      1   - nq of conserved solution component
  !      2:4 - nx(:)+nb, three spatial dimensions + ghost cells
  !      5   - 2 copies of solution array
  real, allocatable, dimension(:,:,:,:) :: U
contains
  subroutine init_simulation_data
    use m_domain, only: nx,nb
    use m_state, only: nq
    it = 0
    t = 0
    allocate(U(nq, -nb+1:nx(1)+nb, -nb+1:nx(2)+nb, -nb+1:nx(3)+nb))
  end subroutine init_simulation_data
end module m_simulation_data

program euler3d
  use m_domain
  use m_simulation_data
  use m_config
  use m_timestepping
! exit the main timestepping loop?
  logical stopflag
! wall-clock start and current time
  real clock_start, time
! fastest wave speed in the entire domain
  real max_wave_speed

  real dx_dt
  integer dirn
  logical dirns(3)

  call init_domain(5)
  call init_config(5)
  call init_simulation_data

  call apply_IC
  call apply_BC(do_x=.true., do_y=.true., do_z=.true.)
  call output

  max_wave_speed = get_initial_max_wave_speed()

  stopflag = .false.
  call clock(clock_start)
  do
     dt = (cfl/max_wave_speed) * dx
     if (it.lt.5) dt = dt/5
     dx_dt = dx/dt

     if (t + dt > tmax) then
        dt = tmax - t
        if (dt < 1e-5) exit
        stopflag = .true.
     end if

     ! max_wave_speed is updated in place by the calls to 'sweep'
     max_wave_speed = 0
     call sweep(dirn=1)
     call apply_BC(do_x=.false., do_y=.true., do_z=.true.)

     call sweep(dirn=2)
     call apply_BC(do_x=.true., do_y=.false., do_z=.true.)

     call sweep(dirn=3)
     call apply_BC(do_x=.true., do_y=.true., do_z=.true.)

     t = t + dt
     it = it + 1

     if (stopflag) exit
     if (mod(it,outstep).eq.0) call output

     call clock(time)
     write (0,*) it, t, dt, time-clock_start
  end do
  call output

contains
  subroutine apply_ic
    use m_eos, only: prim_to_cons
    integer ix,iy,iz
    real r(2), r2
    do iz=-nb+1,nx(3)+nb; do iy=-nb+1,nx(2)+nb; do ix=-nb+1,nx(1)+nb
       U(:,ix,iy,iz) = prim_to_cons([1.0, 0.75, 0.0, 0.0, 1.0, 0.0])

       r = real((/ ix,iy /) - real(nx(1:2))/2) /nx(1:2)
       r2 = dot_product(r,r)
       if (ix.ge.0.3*nx(1)) then
          U(:,ix,iy,iz) = prim_to_cons([0.125, 0.0, 0.0, 0.0, 0.1, 1.0])
       end if

       ! U(:,ix,iy,iz,1) = prim_to_cons([1.0, 0.0, 1.0, 0.0, 1.0, 0.0])
       ! U(:,ix,iy,iz,2) = prim_to_cons([1.0, 0.0, 1.0, 0.0, 1.0, 0.0])

       ! if (iy.ge.nx(2)/3.and.iy.le.2*nx(2)/3) then
       !    U(:,ix,iy,iz,1) = prim_to_cons([1.0, 0.0, 1.0, 0.0, 1.0, 1.0])
       !    U(:,ix,iy,iz,2) = prim_to_cons([1.0, 0.0, 1.0, 0.0, 1.0, 1.0])
       ! endif

    end do; end do; end do
  end subroutine apply_ic

  ! apply boundary conditions on the specified directions
  subroutine apply_bc(do_x,do_y,do_z)
    logical, intent(in) :: do_x,do_y,do_z
! periodic boundaries
    ! if (do_x) U(:,-nb+1:0,:,:,:) = U(:,nx(1)+1-nb:nx(1),:,:,:)
    ! if (do_y) U(:,:,-nb+1:0,:,:) = U(:,:,nx(2)+1-nb:nx(2),:,:)
    ! if (do_z) U(:,:,:,-nb+1:0,:) = U(:,:,:,nx(3)+1-nb:nx(3),:)

    ! if (do_x) U(:,nx(1)+1:nx(1)+nb,:,:,:) = U(:,1:nb,:,:,:)
    ! if (do_y) U(:,:,nx(2)+1:nx(2)+nb,:,:) = U(:,:,1:nb,:,:)
    ! if (do_z) U(:,:,:,nx(3)+1:nx(3)+nb,:) = U(:,:,:,1:nb,:)

! transmissive boundaries
    if (do_x) U(:,-nb+1:0,:,:) = U(:,1:nb,:,:)
    if (do_y) U(:,:,-nb+1:0,:) = U(:,:,1:nb,:)
    if (do_z) U(:,:,:,-nb+1:0) = U(:,:,:,1:nb)

    if (do_x) U(:,nx(1)+1:nx(1)+nb,:,:) = U(:,nx(1)-nb+1:nx(1),:,:)
    if (do_y) U(:,:,nx(2)+1:nx(2)+nb,:) = U(:,:,nx(2)-nb+1:nx(2),:)
    if (do_z) U(:,:,:,nx(3)+1:nx(3)+nb) = U(:,:,:,nx(3)-nb+1:nx(3))
  end subroutine apply_bc

  subroutine output
    use m_output
    integer, parameter :: visit_output_unit=16
    character(len=MAX_FILENAME_LEN) fname
    write (fname,'(A,I0,A)') trim(output_filename_stem), it, '.vtk'
    open(unit=visit_output_unit,file=fname)
    call visit_output(visit_output_unit,it,t,U(:,:,:,:))
    close(16)
  end subroutine output

  function get_initial_max_wave_speed() result(max_wave_speed)
    real max_wave_speed
    max_wave_speed = HUGE(max_wave_speed)
  end function get_initial_max_wave_speed

  ! performs a flux-sweep solution update along the specified
  ! coordinate direction
  subroutine sweep(dirn)
    integer dirn
    integer ix,iy,iz
    select case (dirn)
    case (1)
       do iy=1,nx(2); do iz=1,nx(3)
         call advance_solution_1d(U(:,:,iy,iz), max_wave_speed, dx_dt, dirn)
       end do; end do
    case (2)
       do ix=1,nx(1); do iz=1,nx(3)
          call advance_solution_1d(U(:,ix,:,iz), max_wave_speed, dx_dt, dirn)
       end do; end do
    case (3)
       do ix=1,nx(1); do iy=1,nx(2)
          call advance_solution_1d(U(:,ix,iy,:), max_wave_speed, dx_dt, dirn)
       end do; end do
    end select
  end subroutine sweep

  subroutine clock(t)
!$  use omp_lib
    real, intent(out) :: t
    call cpu_time(t)
!   use openmp clock if it is enabled:
!$  t = omp_get_wtime()
  end subroutine clock
end program euler3d
