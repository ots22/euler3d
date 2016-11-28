module m_domain
  ! length of the domain (in first dimension)
  real domain_length
  ! number of cells in each spatial dimension
  integer nx(3)
  ! number of boundary cells
  integer, parameter :: nb=4
  ! grid resolution
  real dx
  
  namelist/DOMAIN/domain_length,nx
contains
  subroutine init_domain(unit)
    use m_error
    integer, intent(in) :: unit
    nx=0; domain_length=0
    read(unit, nml=DOMAIN)
    call assert(all(nx.ge.0),"number of grid cells (nx) must be positive")
    call assert(domain_length.gt.0, "domain_length must be positive")
    dx = domain_length/nx(1)
  end subroutine init_domain
end module m_domain
