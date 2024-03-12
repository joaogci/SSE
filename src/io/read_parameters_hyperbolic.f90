subroutine read_parameters_hyperbolic(p,q,nl,beta,epsilon, &
  therm_cycles,n_bins,mc_sweeps,S,J_perp,J_par) &
  bind(c,name='read_parameters_hyperbolic')

  use,intrinsic :: iso_c_binding,only:c_char,c_int,c_double,c_long
  implicit none

  integer(kind=c_int),    intent(inout) :: p
  integer(kind=c_int),    intent(inout) :: q
  integer(kind=c_int), intent(inout) :: nl

  real(kind=c_double),    intent(inout) :: beta
  real(kind=c_double),    intent(inout) :: epsilon
  integer(kind=c_long),   intent(inout) :: therm_cycles
  integer(kind=c_int),    intent(inout) :: n_bins
  integer(kind=c_long),   intent(inout) :: mc_sweeps

  real(kind=c_double),    intent(inout) :: S
  real(kind=c_double),    intent(inout) :: J_perp
  real(kind=c_double),    intent(inout) :: J_par

  namelist /Lattice/ p,q,nl
  namelist /Simulation/ beta,epsilon,therm_cycles,n_bins,mc_sweeps
  namelist /Hamiltonian/ S,J_perp,J_par

  open(unit=100,file='parameters',status='old')
  read(unit=100,nml=Lattice)
  read(unit=100,nml=Simulation)
  read(unit=100,nml=Hamiltonian)
  close(unit=100)
  
endsubroutine read_parameters_hyperbolic

